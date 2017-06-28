import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Class for error correction of DNAString reads.
 * <p>
 * All error correction methods are variations of spectral alignment.
 * <p>
 * (1) Majority Alignment Correction (MAC) - correction by majority vote
 * (2) Consensus Alignment Correction (CAC) - like MAC, but relying on unanimous consensus
 * (3) Spectral Alignment Correction (SAC) - basic spectral alignment
 * <p>
 * Created by tohei on 6/26/17.
 */
public class ReadCorrector {

    /**
     * k-mer length.
     */
    private int k;

    /**
     * Coverage cutoff to distinguish between solid (coverage above cutoff) and weak (coverage below cutoff) k-mers.
     */
    private int cutoff;

    /**
     * Reads to correct.
     */
    private DNAString[] reads;

    /**
     * Coverage counts of all k-mers in reads.
     */
    private LinkedHashMap<DNAString, Integer> kmerCounts;

    /**
     * Limit number of replacements (majority vote)
     */
    private final boolean LIMIT_REPLACEMENTS = false;
    /**
     * Maximum number of replacements per read (majority vote).
     */
    private final int MAX_REPLACEMENTS = 4;

    /**
     * Get counts of k-mers and their reverse complements.
     *
     * @return LinkedHashSet mapping k-mers to counts.
     */
    public Map<DNAString, Integer> getCounts() {
        return Collections.unmodifiableMap(kmerCounts);
    }

    /**
     * Associate this ReadCorrector with the given reads.
     *
     * @param reads reads to work on.
     * @param k     k-mer length.
     */
    public void setReads(DNAString[] reads, int k) {
        this.k = k;
        this.reads = reads;
        kmerCounts = DNAStringUtils.kmerCounts(reads, k);
    }

    /**
     * Setter for coverage cutoff.
     *
     * @param cutoff integer coverage cutoff.
     */
    public void setCutoff(int cutoff) {
        this.cutoff = cutoff;
    }

    /**
     * Convenience class for saving replacement votes of k-mers in a read.
     */
    private class ReplacementList {
        /**
         * Array of maps relating bases to votes.
         */
        private LinkedHashMap<Byte, Integer>[] candidateReplacements;

        /**
         * Old read.
         */
        private DNAString old;

        /**
         * Create new ReplacementList.
         *
         * @param old read to correct.
         */
        public ReplacementList(DNAString old) {
            this.old = old;
            this.candidateReplacements = new LinkedHashMap[old.length()];
            for (int i = 0; i < old.length(); i++) {
                this.candidateReplacements[i] = new LinkedHashMap<Byte, Integer>();
            }
        }

        /**
         * Get the length of the old read.
         *
         * @return
         */
        public int getListLength() {
            return old.length();
        }

        /**
         * Add a new candidate to the list (if does not exist) or update the count of an existing
         * candidate (if it does).
         *
         * @param pos position to update.
         * @param c   base voted for.
         */
        public void addCandidate(int pos, byte c) {
            addCandidate(pos, c, 1);
        }

        /**
         * Add a new candidate to the list (if does not exist) or update the count of an existing
         * candidate (if it does).
         *
         * @param pos    position to update.
         * @param c      base voted for.
         * @param weight weight of the vote.
         */
        public void addCandidate(int pos, byte c, int weight) {
            Integer count = candidateReplacements[pos].get(c);
            if (count == null) {
                candidateReplacements[pos].put(c, weight);
            } else {
                candidateReplacements[pos].put(c, count + weight);
            }
        }

        /**
         * Compute the consensus replacement of all bases in the read.
         *
         * @return corrected read (DNAString).
         */
        public DNAString getConsensusReplacement() {
            byte[] seq = new byte[getListLength()];
            for (int i = 0; i < getListLength(); i++) {
                seq[i] = getConsensusCandidate(i);
            }
            return new DNAString(seq);
        }

        /**
         * The consensus replacement at pos is the majority candidate (if the vote is unanimous) or
         * the old base (if it is not)
         *
         * @param pos position to evaluate.
         * @return consensus base.
         */
        private byte getConsensusCandidate(int pos) {
            LinkedHashMap<Byte, Integer> posCandidates = candidateReplacements[pos];

            if (posCandidates.keySet().size() == 1) {
                return posCandidates.keySet().iterator().next();
            } else {
                return this.old.byteAt(pos);
            }
        }

        /**
         * Compute the majority replacement of all bases in the read.
         *
         * @return corrected read (DNAString).
         */
        public DNAString getMajorityReplacement() {
            int replLen = 0;
            if (LIMIT_REPLACEMENTS) {
                replLen = MAX_REPLACEMENTS;
            } else {
                replLen = getListLength();
            }
            int[] topCandidatePos = new int[replLen];
            // default initialization is all zeros
            int[] topCandidateScores = new int[replLen];

            for (int i = 0; i < getListLength(); i++) {
                byte c = getMajorityCandidate(i);
                int score = candidateReplacements[i].get(c);
                for (int j = 0; j < topCandidatePos.length; j++) {
                    if (score > topCandidateScores[j]) {
                        topCandidatePos[j] = i;
                        topCandidateScores[j] = score;
                        break;
                    }
                }
            }
            byte[] repl = old.toByteArray();
            for (int i = 0; i < topCandidatePos.length; i++) {
                int pos = topCandidatePos[i];
                repl[pos] = getMajorityCandidate(pos);
            }
            return new DNAString(repl);
        }

        /**
         * Find the majority candidate at a given position.
         *
         * @param pos position to evaluate.
         * @return majority base.
         */
        private byte getMajorityCandidate(int pos) {
            LinkedHashMap<Byte, Integer> posCandidates = candidateReplacements[pos];
            int maxCount = -1;
            byte topRepl = 0;
            for (Byte c : posCandidates.keySet()) {
                if (posCandidates.get(c) > maxCount) {
                    topRepl = c;
                }
            }
            return topRepl;
        }

        /**
         * Construct a string representation of the votes (for debugging purposes).
         *
         * @return a String representation of this ReplacementList.
         */
        @Override
        public String toString() {
            String s = "  ";
            for (int i = 0; i < getListLength(); i++) {
                s = s + i + " ";
            }
            final char[] alphabet = {'A', 'C', 'G', 'T'};
            for (char c : alphabet) {
                s = s + "\n";
                s = s + c + " ";
                for (int i = 0; i < getListLength(); i++) {
                    Integer count = candidateReplacements[i].get(c);
                    if (count == null) {
                        s = s + "0 ";
                    } else {
                        s = s + count.toString() + " ";
                    }
                }
            }
            return s;
        }
    }

    /**
     * (M)ajority (A)lignment (C)orrection.
     * <p>
     * Every k-mer in a read 'votes' for a replacement of itself by its best Hamming neighbor. The vote is weighted by
     * the coverage of that k-mer. Every base is replaced by the majority candidate of the respective best Hamming
     * neighbors of all overlapping k-mers.
     * k-mers.
     *
     * @param verbose be verbose.
     * @return the number of corrected reads.
     */
    public int computeMAC(boolean verbose) {
        long startTime = System.currentTimeMillis();
        if (verbose) System.out.print("Computing majority alignment (cutoff = " + cutoff + ") ... ");
        int numOfRepl = 0;

        for (int i = 0; i < reads.length; i++) {
            DNAString read = reads[i];
            ReplacementList cand = new ReplacementList(read);
            for (int j = 0; j < read.length() - k + 1; j++) {
                DNAString v = read.subSequence(j, j + k);
                if (kmerCounts.get(v) < cutoff) {
                    v = bestHammingNeighbor(v);
                }
                int weight = kmerCounts.get(v);
                for (int l = 0; l < k; l++) {
                    cand.addCandidate(j + l, v.byteAt(l), weight);
                }
            }
            reads[i] = cand.getMajorityReplacement();
            if (verbose) {
                /*System.out.println("\n" + read + ": ");
                System.out.println(cand.toString());*/
                if (!read.equals(reads[i])) {
                    numOfRepl++;
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime) / 1000 + " seconds).");
            System.out.println(numOfRepl + " of " + reads.length + " reads were corrected.");
        }
        kmerCounts = DNAStringUtils.kmerCounts(reads, k);
        return numOfRepl;
    }

    /**
     * (C)onsensus (A)lignment (C)orrection.
     * <p>
     * Every k-mer in a read 'votes' for a replacement of itself by its best Hamming neighbor. A base is replaced, if
     * all overlapping k-mers agree on that base.
     * k-mers.
     *
     * @param verbose be verbose.
     * @return the number of corrected reads.
     */
    public int computeCAC(boolean verbose) {
        long startTime = System.currentTimeMillis();
        if (verbose) System.out.print("Computing majority alignment (cutoff = " + cutoff + ") ... ");
        int numOfRepl = 0;

        for (int i = 0; i < reads.length; i++) {
            DNAString read = reads[i];
            ReplacementList cand = new ReplacementList(read);
            for (int j = 0; j < read.length() - k + 1; j++) {
                DNAString v = read.subSequence(j, j + k);
                if (kmerCounts.get(v) < cutoff) {
                    v = bestHammingNeighbor(v);
                }
                int weight = kmerCounts.get(v);
                for (int l = 0; l < k; l++) {
                    cand.addCandidate(j + l, v.byteAt(l));
                }
            }
            reads[i] = cand.getConsensusReplacement();
            if (verbose) {
                /*System.out.println("\n" + read + ": ");
                System.out.println(cand.toString());*/
                if (!read.equals(reads[i])) {
                    numOfRepl++;
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime) / 1000 + " seconds).");
            System.out.println(numOfRepl + " of " + reads.length + " reads were corrected.");
        }
        kmerCounts = DNAStringUtils.kmerCounts(reads, k);
        return numOfRepl;
    }

    /**
     * (S)ectral (A)lignment (C)orrection.
     * <p>
     * Low coverage kmers will be replaced by their best Hamming neighbor.
     *
     * @param verbose be verbose.
     * @return number of replacements made.
     */
    public int computeSAC(boolean verbose) {
        long startTime = System.currentTimeMillis();
        if (verbose) System.out.print("Computing spectral alignment (cutoff = " + cutoff + ") ... ");

        int numOfRepl = 0;
        for (int i = 0; i < reads.length; i++) {
            DNAString read = reads[i];
            for (int j = 0; j < read.length() - k + 1; j++) {
                DNAString v = read.subSequence(j, j + k);
                if (kmerCounts.get(v) < cutoff) {
                    numOfRepl++;
                    reads[i] = read.subSequence(0, j).concat(bestHammingNeighbor(v)).concat(read.subSequence(j + k, read.length()));
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime) / 1000 + " seconds).");
            System.out.println("(Total of " + numOfRepl + " replacements)");
        }
        kmerCounts = DNAStringUtils.kmerCounts(reads, k);
        return numOfRepl;
    }

    /**
     * Find best Hamming neighbor of v.
     * <p>
     * Searches for the solid k-mer with the highest count which only differs in one base from v. If no such
     * k-mer exists, v will be returned.
     *
     * @param seq input sequence.
     * @return best Hamming neighbor of v.
     */
    private DNAString bestHammingNeighbor(DNAString seq) {
        DNAString bestSeq = seq;
        int bestCount = cutoff - 1;
        for (int i = 0; i < seq.length(); i++) {
            for (DNAString candidate : seq.allVariations(i)) {
                if (kmerCounts.containsKey(candidate) && (kmerCounts.get(candidate) > bestCount)) {
                    bestSeq = candidate;
                    bestCount = kmerCounts.get(candidate);
                }
            }
        }
        return bestSeq;
    }

    public static void main(String[] args) {
        String[] reads = {"ACTA", "TAGT", "ACTA", "ACTT"};
        DNAString[] reads2 = new DNAString[reads.length];
        for (int i = 0; i < reads.length; i++) {
            reads2[i] = new DNAString(reads[i]);
        }
        System.out.println("Old reads:");
        for (DNAString read : reads2) {
            System.out.print(read.toString() + " ");
        }
        ReadCorrector r = new ReadCorrector();
        r.setReads(reads2, 3);
        r.setCutoff(2);
        r.computeSAC(true);
        for (DNAString read : reads2) {
            System.out.print(read.toString() + " ");
        }
    }
}
