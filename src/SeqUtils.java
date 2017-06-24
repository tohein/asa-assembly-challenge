import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;


/**
 * Various utility functions for reading and processing
 * sequence reads.
 * @author tohei
 * 
 */
public final class SeqUtils {
	private static final HashMap<Character, Character> COMPLEMENTS;
	static {
        COMPLEMENTS = new HashMap<Character, Character>();
        COMPLEMENTS.put('A', 'T');
        COMPLEMENTS.put('T', 'A');
        COMPLEMENTS.put('G', 'C');
        COMPLEMENTS.put('C', 'G');
    }
	
	private SeqUtils() {
        throw new AssertionError();
    }
	
	/**
	 * Compute the set of k-mers which are a substring of at least one 
	 * String in reads.               
	 * <p>
	 * @param reads array of reads.      
	 * @param k k-mer length.
	 * @return LinkedHashSet containing all k-mers.
	 */
	public static LinkedHashSet<String> allKmers(String[] reads, int k) {
		LinkedHashSet<String> kmers = new LinkedHashSet<String>(reads.length);
		for (String s : reads) {
			for (int i = 0; i < s.length() - k + 1; i++) {
				kmers.add(s.substring(i, i + k));
			}
		}
		return kmers;
	}

	/**
	 * Computes for every k-mer m the  number of reads which
	 * contain m. 
	 * <p>
	 * @param reads array of reads.      
	 * @param k k-mer length.
	 * @return LinkedHashSet mapping k-mers to counts.
	 */
	public static LinkedHashMap<String, Integer> kmerCounts(String[] reads, int k) {
		return kmerCounts(reads, k, false);
	}

    /**
     * Computes for every k-mer m the  number of reads which
     * contain m.
     * <p>
     * @param reads array of reads.
     * @param k k-mer length.
     * @param countRC count reverse complements of reads.
     * @return LinkedHashSet mapping k-mers to counts.
     */
    public static LinkedHashMap<String, Integer> kmerCounts(String[] reads, int k, boolean countRC) {
        LinkedHashMap<String, Integer> kmerCounts = new LinkedHashMap<String, Integer>(reads.length);
        for (String s : reads) {
            // updated contains all k-mers found in read s
            LinkedHashSet<String> updated = new LinkedHashSet<String>();
            for (int i = 0; i < s.length() - k + 1; i++) {
                String kmer = s.substring(i, i + k);
                if (!updated.contains(kmer)) {
                    int count = 1;
                    if (kmerCounts.containsKey(kmer)) {
                        count += kmerCounts.get(kmer);
                    }
                    kmerCounts.put(kmer, count);
                }
                updated.add(kmer);
            }

            if (countRC) {
                // updated contains all k-mers found in reverse complement of s
                String sRC = reverseComplement(s);
                updated = new LinkedHashSet<String>();
                for (int i = 0; i < sRC.length() - k + 1; i++) {
                    String kmer = sRC.substring(i, i + k);
                    if (!updated.contains(kmer)) {
                        int count = 1;
                        if (kmerCounts.containsKey(kmer)) {
                            count += kmerCounts.get(kmer);
                        }
                        kmerCounts.put(kmer, count);
                    }
                    updated.add(kmer);
                }
            }
        }
        return kmerCounts;
    }
	
	/**
	 * Compute the reverse complement of a sequence.
	 * <p>
	 * @param s String comprised of the letters 'A', 'T', 'G', 'C'.
	 * @return Reverse complement of input string.
	 */
	public static String reverseComplement(String s) {
		String revComp = "";
		for (int i = 0; i < s.length(); i++) {
			revComp = COMPLEMENTS.get(s.charAt(i)) + revComp;
		}
		return revComp;
	}

    /**
     * Compute the complement of a character.
     * <p>
     * @param c One of 'A', 'T', 'G', 'C'.
     * @return Complement of input string.
     */
    public static char complement(char c) {
        return COMPLEMENTS.get(c);
    }
	
	/**
	 * Compute four extensions of input string s by adding one 
	 * additional base ('A', 'T', 'G', 'C') to the end of s.
	 * <p>
	 * @param s input sequence.
	 * @return String array of length 4 containing extensions.
	 */
	public static String[] extendBack(String s) {
		String[] extB = new String[COMPLEMENTS.size()];
		int k = 0;
		for (char c : COMPLEMENTS.keySet()) {
			extB[k] = s + c;
			k++;
		}
		return extB;
	}
	
	/**
	 * Compute four extensions of input string s by adding one 
	 * additional base ('A', 'T', 'G', 'C') to the beginning of s.
	 * <p>
	 * @param s input sequence.
	 * @return String array of length 4 containing extensions.
	 */
	public static String[] extendFront(String s) {
		String[] extF = new String[COMPLEMENTS.size()];
		int k = 0;
		for (char c : COMPLEMENTS.keySet()) {
			extF[k] = c + s;
			k++;
		}
		return extF;
	}
		
	/**
	 * Reads in a FASTA file and returns reads in a string array.
	 * @param fileName String containing full path to FASTA file.
	 * @return String array containing reads.
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static String[] readFasta(String fileName) throws FileNotFoundException, IOException {
		String line;
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		ArrayList<String> sequences = new ArrayList<String>();
		String seq = "";
		while((line = in.readLine()) != null) {
		   if(line.startsWith(">")) {
			   if (seq.length() > 0) {
				   sequences.add(seq);
			   } 
			   seq = "";			   
		   } else {
			   seq = seq + line;
		   }
		}
		if (seq.length() > 0) {
			sequences.add(seq);
		} 
		in.close();
		return sequences.toArray(new String[sequences.size()]);
	}
	
	/**
	 * Writes string array to FASTA file.
	 * @param fileName path to output file.
	 * @param reads String array containing reads.
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static void writeFasta(String fileName, String[] reads) throws IOException {
		final int MAX_LENGTH = 80;
		
		PrintWriter writer = new PrintWriter(fileName, "UTF-8");
		for (int i = 0; i < reads.length; i++) {
			writer.println(">CONTIG_" + i);
			String seq = reads[i];
			// each line of a sequence should have fewer than 80 characters
			for (int j = 0; j < seq.length(); j += MAX_LENGTH) {
				writer.println(seq.substring(j, Math.min(seq.length(), j + MAX_LENGTH)));
		    }
		}
		writer.close();
	}

    /**
     * Saves k-mer counts map to tab separated file.
     * @param fileName path to output file.
     * @param counts Map as generated by kmerCounts.
     * @throws IOException
     */
	public static void saveKmerCounts(String fileName, LinkedHashMap<String, Integer> counts) throws IOException {
        PrintWriter writer = new PrintWriter(fileName, "UTF-8");
        for (Map.Entry<String, Integer> entry : counts.entrySet()) {
            writer.println(entry.getKey() + "\t" + entry.getValue());
        }
        writer.close();
    }
	
	/**
	 * Function for error correction of reads using spectral alignment.
	 * <p>
	 * @param reads array of reads. 
	 * @param cutoff integer threshold for coverage cutoff. Every k-mer with
	 * lower coverage will be replaced by its best Hamming neighbor.
	 * @param k k-mer length.
	 * @param verbose be verbose.
	 * @return String array of error corrected reads.
	 */
	public static String[] spectralAlignment(String[]reads, int cutoff, int k, boolean verbose) {
		long startTime = System.currentTimeMillis();
		if (verbose) System.out.print("Computing spectral alignment (cutoff = " + cutoff + ") ... ");

		LinkedHashMap<String, Integer> kmerCounts = kmerCounts(reads, k, true);
		int numOfReplacements = 0;
		for (int i = 0; i < reads.length; i++) {
			String read = reads[i];
			for (int j = 0; j < read.length() - k + 1; j++) {
				String v = read.substring(j, j + k);
				if (kmerCounts.get(v) < cutoff) {
					numOfReplacements ++;
					//System.out.println("Low cov string " + v);
					reads[i] = read.substring(0, j) + bestHammingNeighbor(v, cutoff, kmerCounts) + read.substring(j + k);
				}
			}
		}
		if (verbose) {
			long endTime = System.currentTimeMillis();
			System.out.println("done (" + (endTime - startTime)/1000 + " seconds).");
			System.out.println("(Total of " + numOfReplacements + " replacements)");
		}
		return reads;
		
	}

	/**
	 * Find best Hamming neighbor of v.
	 * @param v input sequence.
	 * @param cutoff initial count threshold. Only consider candidates with higher coverage.
	 * @param kmerCounts LinkedHashMap mapping k-mers to counts.
	 * @return best Hamming neighbor of v.
	 */
	public static String bestHammingNeighbor(String v, int cutoff, LinkedHashMap<String, Integer> kmerCounts) {
		String bestSeq = v;
		int bestCount = cutoff - 1;
		for (int i = 0; i < v.length(); i++) {
			for (char c : COMPLEMENTS.keySet()) {
				String candidate = v.substring(0, i) + c + v.substring(i + 1);
				//System.out.println("Candidate: "+candidate);
				if (kmerCounts.containsKey(candidate) && (kmerCounts.get(candidate) > bestCount)) {
					bestSeq = candidate;
					bestCount = kmerCounts.get(candidate);
				}
			}
		}
		return bestSeq;
	}
}
