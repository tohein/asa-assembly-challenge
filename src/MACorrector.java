import java.io.IOException;
import java.util.LinkedHashMap;

/**
 * (M)ajority (A)lignment Corrector.
 * Created by tohei on 6/23/17.
 */
public class MACorrector {
    private int k;
    private int cutoff;
    private String[] reads;
    public LinkedHashMap<String, Integer> kmerCounts;

    public MACorrector(int k) {
        this.k = k;
    }
    public MACorrector(){}

    public void setReads(String[] reads, int k) {
        this.k = k;
        this.setReads(reads);
    }

    public void setReads(String[] reads) {
        this.reads = reads;
        kmerCounts = SeqUtils.kmerCounts(reads, k, false);
    }

    public void setCutoff(int cutoff) {
        this.cutoff = cutoff;
    }

    private class ReplacementList {
        private LinkedHashMap<Character, Integer>[] candidateReplacements;
        private int length;

        public ReplacementList(int length) {
            this.length = length;
            this.candidateReplacements = new LinkedHashMap[length];
            for (int i = 0; i < length; i++) {
                this.candidateReplacements[i] = new LinkedHashMap<Character, Integer>();
            }
        }

        public void addCandidate(int pos, char c) {
            Integer count = candidateReplacements[pos].get(c);
            if (count == null) {
                candidateReplacements[pos].put(c, 0);
            } else {
                candidateReplacements[pos].put(c, count + 1);
            }
        }

        public String getMajorityReplacement() {
            String s = "";
            for (int i = 0; i < this.length; i++) {
                s = s + getMajorityCandidate(i);
            }
            return s;
        }

        private char getMajorityCandidate(int pos) {
            LinkedHashMap<Character, Integer> posCandidates = candidateReplacements[pos];
            int maxCount = -1;
            char topRepl = 0;
            for (Character c : posCandidates.keySet()) {
                if (posCandidates.get(c) > maxCount) {
                    topRepl = c;
                }
            }
            return topRepl;
        }
    }

    public String[] correctReads(boolean verbose) {
        long startTime = System.currentTimeMillis();
        if (verbose) System.out.print("Computing majority alignment (cutoff = " + cutoff + ") ... ");
        int numOfRepl = 0;

        for (int i = 0; i < reads.length; i++) {
            String read = reads[i];
            ReplacementList cand = new ReplacementList(read.length());
            for (int j = 0; j < read.length() - k + 1; j++) {
                String v = read.substring(j, j + k);
                if (kmerCounts.get(v) < cutoff) {
                    v = bestHammingNeighbor(v);
                }
                for (int l = 0; l < k; l++) {
                    cand.addCandidate(j + l, v.charAt(l));
                }
            }
            reads[i] = cand.getMajorityReplacement();
            if (verbose) {
                if (!read.equals(reads[i])) {
                    numOfRepl ++;
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime)/1000 + " seconds).");
            System.out.println(numOfRepl + " of " + reads.length + " reads were corrected.");
        }
        return reads;
    }

    /**
     * Find best Hamming neighbor of v.
     * @param v input sequence.
     * @return best Hamming neighbor of v.
     */
    private String bestHammingNeighbor(String v) {
        final char[] alphabet = {'A', 'C', 'G', 'T'};
        String bestSeq = v;
        int bestCount = cutoff - 1;
        for (int i = 0; i < v.length(); i++) {
            for (char c : alphabet) {
                String candidate = v.substring(0, i) + c + v.substring(i + 1);
                if (kmerCounts.containsKey(candidate) && (kmerCounts.get(candidate) > bestCount)) {
                    bestSeq = candidate;
                    bestCount = kmerCounts.get(candidate);
                }
            }
        }
        return bestSeq;
    }

    public static void main(String[] args) {
        String[] reads = null;
        try {
            reads = SeqUtils.readFasta("/home/tohei/Data/ASAData/reads_complex.fasta");
        } catch (IOException e) {
            e.printStackTrace();
        }
        MACorrector mac  = new MACorrector();
        System.out.print("Counting k-mers ... ");
        mac.setReads(reads, 21);
        System.out.println("done.");
        mac.setCutoff(10);
        mac.correctReads(true);
    }

}
