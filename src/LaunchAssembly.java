import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashSet;

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
			System.out.println(" [size] - k-mer size");
			System.out.println(" [k-mer] - (optional) print incoming and outgoing k-mers of [k-mer]");
			System.out.println();
			return;
		}
		if (args.length < 3) {
			System.err.println("Insufficient number of arguments.");
			System.out.println("Usage: java LaunchAssembly [reads.fasta] [output.fasta] [size] [k-mer]");
			return;
		}
		
		// read command line arguments
		int k = 0;
		try {
			k = Integer.parseInt(args[2]);			
		} catch (NumberFormatException e) {
			System.err.println("Invalid k-mer size.");
			return;
		}
		String inputFile = args[0];
		String outputFile = args[1];
		
		String[] inputs = null;
		try {
			inputs = SeqUtils.readFasta(inputFile);
		} catch (FileNotFoundException e1) {
			System.err.println("Could not find input file.");
		} catch (IOException e2) {
			System.err.println("Failed to read input file.");
		}
		
		// solve task (a) manually (de Bruijn graph will be built on error corrected reads)
		if (args.length == 4) {
			String kmer = args[3];
			if (kmer.length() != k) {
				System.err.println("Invalid k-mer. k-mer has to be of specified length.");
				return;
			}			
			LinkedHashSet<String> kmers = SeqUtils.allKmers(inputs, k);
			System.out.print("Incoming k-mers - ");
			boolean hasEdge = false;
			for (String s : SeqUtils.extendFront(kmer.substring(0, k-1))) {
				if (kmers.contains(s)) {
					hasEdge = true;
					System.out.print(s + " ");
				}
			}
			if (!hasEdge) System.out.println("None");
			else System.out.println();
			System.out.print("Outgoing k-mers - ");
			hasEdge = false;
			for (String s : SeqUtils.extendBack(kmer.substring(1, k))) {
				if (kmers.contains(s)) {
					hasEdge = true;
					System.out.print(s + " ");
				}
			}			
			if (!hasEdge) System.out.println("None");
			System.out.println();
		}
		
		int cutoff = 10;
		// correct reads
		inputs = SeqUtils.spectralAlignment(inputs, cutoff, k, verbose);
		// compute contigs
		System.out.println();
		DBGraph G = new DBGraph(inputs, k, true);		
		String[] contigs = G.findContigs(true, true, true, cutoff);		
		System.out.println();
		System.out.println("Max contig length: " + G.getMaxContig());
		System.out.println();
		System.out.print("Saving to " + outputFile + " ... ");
		try {
			SeqUtils.writeFasta(outputFile, contigs);
		} catch (IOException e) {
			System.out.println();
			System.err.println("Could not save to output file.");
			return;
		}	
		System.out.println("done.");
	}
}
