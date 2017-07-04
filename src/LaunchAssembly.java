import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashMap;

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
        int k = 0;
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
            // solve task (a) manually (de Bruijn graph will be built on error corrected reads)
            LinkedHashMap<DNAString, Integer> counts = DNAStringUtils.kmerCounts(inputs, k, false);
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

        // TODO find cutoff
        int cutoff = 2;
        // error correction will only be used if there are at least two reads
        boolean correctErrors = true;
        if (inputs.length < 2) {
            correctErrors = false;
            cutoff = 0;
        }

        System.out.println();
        System.out.println("Read error correction: ");
        ReadCorrector rcor = new ReadCorrector();
        rcor.setReads(inputs, k);

        if (correctErrors) {
            rcor.setCutoff(cutoff + 1);
            rcor.computeSAC(true);
            rcor.setCutoff(cutoff);
            rcor.computeSAC(true);
            rcor.computeBestReplacement(true);
        }

        System.out.println();
        DBGraph G = new DBGraph(inputs, rcor.getCounts(), k, true);
        G.simplify(correctErrors, cutoff, true);
        DNAString[] contigs = G.getSequences(true, 0);

        System.out.println();
        int max = G.getMaxSequenceLength();
        System.out.println("Max contig length: " + max);
        System.out.println("Avg contig length: " + G.getAvgSequenceLength());
        System.out.println();

        // save to output fasta
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
