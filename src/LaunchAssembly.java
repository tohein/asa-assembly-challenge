import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashSet;

public class LaunchAssembly {
	
	public static void main(String args[]) {
		//boolean verbose = true;
		
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


		int[] cutoff1 = {10};
		int[] cutoff2 = {1,2,3,4,5,7,10,15};
		int maxContigLength = 0;
		int topc1 = -1;
		int topc2 = -1;

		// save k-mer counts
		try {
			String countfile = "/home/tohei/Downloads/counts.txt";
			System.out.println(countfile);
			SeqUtils.saveKmerCounts(countfile, SeqUtils.kmerCounts(inputs, 21));
		} catch (IOException e) {
			e.printStackTrace();
		}

		// correct reads
		inputs = SeqUtils.spectralAlignment(inputs, 30, k, true);

		// find best parameter settings
		System.out.println("Testing different paramters ... ");
		for (int i = 0; i < cutoff1.length; i++) {
			for (int j = 0; j < cutoff2.length; j++) {
				System.out.println("Cutoffs: " + cutoff1[i] + " (rmv tips), " + cutoff2[j] + " (rmv low cov).");
				DBGraph G = new DBGraph(inputs, k, false);
				String[] contigs = G.findContigs(true, false, true, cutoff1[i], cutoff2[j]);
				System.out.println();
				int max = G.getMaxContig();
				System.out.println("Max contig length: " + max);
				System.out.println();
				if (max > maxContigLength) {
					maxContigLength = max;
					topc1 = i;
					topc2 = j;
				}
			}
		}
		// repeat for best parameters
		System.out.println("Top parameters: " + cutoff1[topc1] + " (rmv tips), " + cutoff2[topc2] + " (rmv low cov).");
		DBGraph G = new DBGraph(inputs, k, true);
		String[] contigs = G.findContigs(true, true, true, cutoff1[topc1], cutoff2[topc2]);
		System.out.println();
		int max = G.getMaxContig();
		System.out.println("Max contig length: " + max);
		System.out.println();

		// save
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
