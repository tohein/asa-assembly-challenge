import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashMap;

/**
 * Start de novo assembly using de Bruijn graphs.
 */
public class LaunchAssembly {

    public static void main(String args[]) {
        boolean verbose = true;

        if (args.length == 0) {
            System.out.println("Assembly Challenge - de novo assembly using de Bruijn graphs in Java");
            System.out.println();
            System.out.println("Usage: java LaunchAssembly [reads.fasta] [output.fasta] [size] [k-mer]");
            System.out.println();
            System.out.println(" [reads.fasta] - input fasta file");
            System.out.println(" [output.fasta] - output file (fasta)");
            System.out.println(" [size] - k-mer size (odd integer)");
            System.out.println(" [k-mer] - (optional) print incoming and outgoing k-mers of [k-mer]");
            System.out.println();
            return;
        }
        if (args.length < 3) {
            System.err.println("Insufficient number of arguments.");
            System.out.println("Usage: java LaunchAssembly [reads.fasta] [output.fasta] [size] [k-mer]");
            return;
        }

        // parse k-mer size
        int k;
        try {
            k = Integer.parseInt(args[2]);
        } catch (NumberFormatException e) {
            System.err.println("Invalid k-mer size.");
            return;
        }
        if (k % 2 == 0 || k < 3) {
            System.err.println("k should be odd integer >= 3.");
            return;
        }

        // read file names
        String inputFile = args[0];
        String outputFile = args[1];

        // read in input fasta file
        DNAString[] inputs = null;
        System.out.print("Reading input file ... ");
        try {
            inputs = DNAStringUtils.readFasta(inputFile);
        } catch (FileNotFoundException e1) {
            System.err.println("Could not find input file.");
        } catch (IOException e2) {
            System.err.println("Failed to read input file.");
        }
        System.out.println("done.");

        if (args.length == 4) {
            // find k-mer neighbors
            DNAString kmer = new DNAString(args[3]);
            if (kmer.length() != k) {
                System.err.println("Invalid k-mer. k-mer has to be of specified length.");
                return;
            }
            // solve task (a) manually (de Bruijn graph will be built on error corrected reads and their rc)
            LinkedHashMap<DNAString, Integer> counts = DNAStringUtils.kmerCounts(inputs, k, true);
            System.out.print("Incoming k-mers - ");
            boolean hasEdge = false;
            for (DNAString s : kmer.subSequence(0, k - 1).allVariations(-1)) {
                if (counts.containsKey(s)) {
                    hasEdge = true;
                    System.out.print(s.toString() + " ");
                }
            }
            if (!hasEdge) System.out.println("None");
            else System.out.println();
            System.out.print("Outgoing k-mers - ");
            hasEdge = false;
            for (DNAString s : kmer.subSequence(1, k).allVariations(k)) {
                if (counts.containsKey(s)) {
                    hasEdge = true;
                    System.out.print(s.toString() + " ");
                }
            }
            if (!hasEdge) System.out.println("None");
            System.out.println();
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ReadCorrector rcor = new ReadCorrector();
        rcor.setReads(inputs, k);
        DBGraph G = null;
        if (inputs.length < 2) {
            // error correction will only be used if there are at least two reads
            System.out.println("(single read -> no error correction, cutoff = 0)");
            System.out.println("BUILD DE BRUIJN GRAPH");
            G = new DBGraph(inputs, rcor.getCounts(), k, verbose);
            G.simplify(false, 0, verbose);
        } else {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            System.out.println("\nPHASE I : AUTO CUTOFF DETECTION (this might take a while)\n");
            // build de Bruijn graph without read error correction and coverage cutoff
            G = new DBGraph(inputs, rcor.getCounts(), k, verbose);
            G.simplify(true, 0, verbose);
            System.out.println();
            // analyze k-mer distribution to find cutoff
            CovCutoffFinder cutoffFinder = new CovCutoffFinder(G.getAdjustedCovData());
            int cutoff = 2;
            try {
                cutoff = cutoffFinder.findCutoff(verbose);
            } catch (Exception e) {
                e.printStackTrace();
                System.out.println("\nFailed to determine cutoff. Use default cutoff: 2");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            System.out.println("\nPHASE II : READ ERROR CORRECTION\n");

            // correct reads using spectral alignment and the more conservative bestReplacement algorithm
            if (cutoff > 0) {
                rcor.setCutoff(cutoff + 1);
                rcor.computeSAC(verbose);
                rcor.setCutoff(cutoff);
                rcor.computeSAC(verbose);
                rcor.computeBestReplacement(verbose);
            } else {
                System.out.println("(no read error correction)");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            System.out.println("\nPHASE III : BUILD FINAL DE BRUIJN GRAPH\n");
            if (cutoff > 0) {
                // rebuild de Bruijn graph on corrected reads and use coverage cutoff to obtain final contigs
                G = new DBGraph(inputs, rcor.getCounts(), k, verbose);
                G.simplify(true, cutoff, verbose);
            } else {
                System.out.println("(no re-build necessary)");
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // display summary and save output to fasta

        System.out.println("\nSUMMARY\n");
        DNAString[] contigs = G.getSequences(true, 0);
        int max = G.getMaxSequenceLength();
        System.out.println("Number of contigs: " + G.getSize() / 2);
        System.out.println("Max contig length: " + max);
        System.out.println("Avg contig length: " + G.getAvgSequenceLength());
        System.out.println();

        System.out.print("Saving to " + outputFile + " ... ");
        try {
            DNAStringUtils.writeFasta(outputFile, contigs);
        } catch (IOException e) {
            System.out.println();
            System.err.println("Could not save to output file.");
            return;
        }
        System.out.println("done.");
    }
}
