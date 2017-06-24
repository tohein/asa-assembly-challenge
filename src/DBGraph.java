import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.NoSuchElementException;


/**
 * Simple implementation of a de Brujin graph for
 * sequence assembly.
 * <p>
 *
 * @author tohei
 */
public class DBGraph {
    /**
     * Stores all blocks in the graph.
     */
    private LinkedHashMap<String, Block> blocks;
    /**
     * k-mer size
     */
    private int k;

    // TODO figure out a way to store short sequences

    /**
     * Blocks are nodes in the de Bruijn graph.
     * Every block corresponds to a k-mer or a simplified node
     * and stores a pointer to its twin block (reverse complement).
     *
     * @author tohei
     */
    private class Block {
        private LinkedList<Block> outgoing;
        private Block twin;
        private String seq;
        private int count;
        private int mult;

        public Block(String seq, int count, int mult) {
            outgoing = new LinkedList<Block>();
            this.twin = null;
            this.seq = seq;
            this.count = count;
            this.mult = mult;
        }

        public Block(String seq) {
            this(seq, 0, 0);
        }

        public int getIndegree() {
            int indegree = 0;
            if (this.twin != null) {
                indegree = this.twin.outgoing.size();
            }
            return indegree;
        }

        public int getOutdegree() {
            return this.outgoing.size();
        }

        public boolean addOutgoing(Block n) {
            if (!outgoing.contains(n)) {
                outgoing.add(n);
                return true;
            }
            return false;
        }

        public boolean rmvOutgoing(Block n) {
            return outgoing.remove(n);
        }

        public boolean addIncoming(Block n) {
            if (!twin.outgoing.contains(n)) {
                twin.outgoing.add(n);
                return true;
            }
            return false;
        }

        public boolean rmvIncoming(Block n) {
            return this.twin.outgoing.remove(n);

        }

        public String getSeq() {
            return seq;
        }

        public String getRCSeq() {
            return twin.seq;
        }

        public String getShortSeq() {
            return seq.substring(k - 1);
        }

        public float getCov() {
            return count / (float) mult;
        }

        @Override
        public String toString() {
            String twinSeq = "";
            if (twin == null) {
                twinSeq = "No twin";
            } else {
                twinSeq = twin.seq;
            }
            String s = "[Block: " + seq + "]";
            s = s + "\n[Twin:  " + twinSeq + "]";
            s = s + "\nOut:[";
            for (Block n : outgoing) {
                s = s + n.seq + " ";
            }
            s = s + "]";
            return s;
        }

    }

    public int getSize() {
        return blocks.size();
    }

    public int getK() {
        return k;
    }

    public DBGraph(String[] inputs, int k) {
        this(inputs, k, false);
    }

    public DBGraph(String[] inputs, int k, boolean verbose) {

        long startTime = System.currentTimeMillis();
        if (verbose) {
            System.out.println("Building de Bruijn graph from " + k
                    + "-mers of " + inputs.length + " reads ...");
        }
        LinkedHashMap<String, Integer> kmers = SeqUtils.kmerCounts(inputs, k);
        blocks = new LinkedHashMap<String, Block>(kmers.size());
        this.k = k;

        if (verbose) {
            System.out.println(" adding nodes ... ");
        }
        // add blocks
        for (String kmer : kmers.keySet()) {
            String kmerRC = SeqUtils.reverseComplement(kmer);
            // TODO count kmers and their complements separately
            int count = kmers.get(kmer);
            if (kmers.containsKey(kmerRC)) {
                count += kmers.get(kmerRC);
            }
            Block n1;
            Block n2;
            if (!blocks.containsKey(kmer)) {
                n1 = new Block(kmer, count, 1);
                n2 = new Block(kmerRC, count, 1);
                n1.twin = n2;
                n2.twin = n1;

                blocks.put(kmer, n1);
                blocks.put(kmerRC, n2);
            }
        }
        if (verbose) {
            System.out.println(" adding edges ... ");
        }
        // add edges
        for (Map.Entry<String, Block> entry : blocks.entrySet()) {
            String seq = entry.getKey();
            for (String s : SeqUtils.extendBack(seq.substring(1, k))) {
                Block n = blocks.get(s);
                // do not add edges to twin node
                if (n != null && n != entry.getValue().twin) {
                    entry.getValue().addOutgoing(n);
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime) / 1000 +
                    " seconds). " + getSize() + " nodes.\n");
        }
    }

    public boolean mergeBubbles(boolean verbose) {
        return false;
    }

    public boolean removeLowCoverageBlocks(boolean verbose, int cutoff) {
        if (verbose) {
            System.out.print("Removing low coverage blocks. ");
        }
        LinkedList<Block> blocksToRemove = new LinkedList<DBGraph.Block>();
        for (Map.Entry<String, Block> entry : blocks.entrySet()) {
            Block n = entry.getValue();
            if (n.getCov() < cutoff) {
                blocksToRemove.add(n);
            }
        }
        if (blocksToRemove.size() > 0) {
            System.out.print("Graph size - before: " + this.getSize());
            for (Block block : blocksToRemove) {
                this.rmvBlock(block);
            }
            System.out.println(", after: " + this.getSize());
            return true;
        } else {
            System.out.println("No changes");
            return false;
        }
    }

    public boolean removeTips(boolean verbose, int cutoff) {
        if (verbose) {
            System.out.print("Removing Tips. ");
        }
        LinkedList<Block> blocksToRemove = new LinkedList<DBGraph.Block>();
        for (Map.Entry<String, Block> entry : blocks.entrySet()) {
            Block n = entry.getValue();
            if ((n.getSeq().length() < 2 * k) && (n.getCov() < cutoff)) {
                if (n.getIndegree() == 0 || n.getOutdegree() == 0) {
                    blocksToRemove.add(n);
                }
            }
        }
        if (blocksToRemove.size() > 0) {
            System.out.print("Graph size - before: " + this.getSize());
            for (Block block : blocksToRemove) {
                this.rmvBlock(block);
            }
            System.out.println(", after: " + this.getSize());
            return true;
        } else {
            System.out.println("No changes");
            return false;
        }
    }

    private void rmvBlock(Block n) {
        for (Block block : n.outgoing) {
            block.twin.rmvOutgoing(n.twin);
        }
        for (Block block : n.twin.outgoing) {
            block.twin.rmvOutgoing(n);
        }
        blocks.remove(n.getSeq());
        blocks.remove(n.getRCSeq());
    }

    public boolean collapse(boolean verbose) {
        long startTime = System.currentTimeMillis();
        int blocksPriorCollapse = this.getSize();

        LinkedList<Block> startingBlocks = new LinkedList<Block>();
        for (Map.Entry<String, Block> entry : blocks.entrySet()) {
            Block n = entry.getValue();
            if (n.getOutdegree() == 1) {
                if (n.getIndegree() != 1) {
                    // multiple or no ancestors
                    startingBlocks.add(n);
                } else {
                    // single ancestor
                    Block ancestor = (n.twin.outgoing.iterator().next()).twin;
                    if (ancestor.getOutdegree() > 1) {
                        startingBlocks.add(n);
                    }
                }
            }
        }
        if (verbose) {
            System.out.println("Collapsing graph. Number of starting blocks: " + startingBlocks.size());
        }
        boolean collapsed = false;
        for (Block start : startingBlocks) {
            if (blocks.containsKey(start.getSeq())) {
                boolean startBlockCollapse = collapseLinearStretch(start);
                collapsed = collapsed || startBlockCollapse;
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print("Collapsed de Bruijn graph (" + (endTime - startTime) / 1000 + " seconds). ");
            if (collapsed) {
                System.out.println("Graph size - before: " + blocksPriorCollapse + ", after: " + getSize());
            } else {
                System.out.println("No changes.");
            }
        }
        return collapsed;
    }

    private boolean collapseLinearStretch(Block start) {
        boolean collapsed = false;
        LinkedList<Block> linearStretch = new LinkedList<Block>();
        linearStretch.add(start);

        // all start blocks have exactly one outgoing edge
        Block n = start.outgoing.iterator().next();
        while ((n.getIndegree() == 1) && !(linearStretch.contains(n) || linearStretch.contains(n.twin))) {
            linearStretch.add(n);
            if (n.getOutdegree() == 1) {
                n = n.outgoing.iterator().next();
            } else {
                break;
            }

        }
        if (linearStretch.size() > 1) {
            collapsed = true;
            collapseBlocks(linearStretch);
        }
        return collapsed;
    }

    private void collapseBlocks(LinkedList<Block> linearStretch) {
        Block first = linearStretch.getFirst();
        Block last = linearStretch.getLast();

        String shortSeq = "";
        int count = 0;
        int mult = 0;
        for (Block block : linearStretch) {
            count = count + block.count;
            mult = mult + block.mult;
            if (block == first) continue;
            shortSeq = shortSeq + block.getShortSeq();
        }
        String shortSeqRC = SeqUtils.reverseComplement(shortSeq);
        String newSeq = first.getSeq() + shortSeq;
        String newSeqRC = shortSeqRC + first.getRCSeq();
        Block newBlock = new Block(newSeq, count, mult);
        Block newTwin = new Block(newSeqRC, count, mult);
        newBlock.twin = newTwin;
        newTwin.twin = newBlock;

        newBlock.outgoing = last.outgoing;
        newTwin.outgoing = first.twin.outgoing;

        // update incoming to n
        for (Block block : newTwin.outgoing) {
            block.twin.rmvOutgoing(first);
            block.twin.addOutgoing(newBlock);
        }
        // update outgoing from m
        for (Block block : newBlock.outgoing) {
            block.twin.rmvOutgoing(last.twin);
            block.twin.addOutgoing(newTwin);
        }

        // remove old blocks and add new blocks
        for (Block block : linearStretch) {
            blocks.remove(block.getSeq());
            blocks.remove(block.getRCSeq());
        }
        blocks.put(newSeq, newBlock);
        blocks.put(newSeqRC, newTwin);
    }

    public String[] getIncoming(String s) {
        Block m = blocks.get(s);
        if (m == null) {
            throw new NoSuchElementException("String not found.");
        } else {
            m = m.twin;
            String[] out = new String[m.getOutdegree()];
            int j = 0;
            for (Block n : m.outgoing) {
                out[j] = n.getRCSeq();
                j++;
            }
            return out;
        }
    }

    public String[] getOutgoing(String s) {
        Block m = blocks.get(s);
        if (m == null) {
            throw new NoSuchElementException("String not found.");
        } else {
            String[] out = new String[m.getOutdegree()];
            int j = 0;
            for (Block n : m.outgoing) {
                out[j] = n.getSeq();
                j++;
            }
            return out;
        }
    }

    public boolean contains(String s) {
        return blocks.containsKey(s);
    }

    public void simplify(boolean verbose, boolean correctErrors, int cutoffTips, int cutoffCov) {
        if (verbose) {
            String s = null;
            if (correctErrors) s = "with";
            else s = "without";
            System.out.println("Computing contigs " + s + " error correction:");
        }
        int iter = 1;
        boolean modified = true;

        while (modified) {
            if (verbose) {
                System.out.println("\t------------- Iteration " + iter + " -------------");
            }
            modified = this.collapse(verbose);
            if (correctErrors) {
                modified = this.removeTips(verbose, cutoffTips) || modified;
                modified = this.removeLowCoverageBlocks(verbose, cutoffCov) || modified;
            }
            iter++;
        }
    }

    public String[] getContigs(boolean includeRC) {
        LinkedHashSet<String> visited = new LinkedHashSet<String>(blocks.size());
        String[] contigs = null;
        if (includeRC) {
            contigs = new String[blocks.size()];
        } else {
            contigs = new String[blocks.size() / 2];
        }
        int i = 0;
        for (Map.Entry<String, Block> entry : blocks.entrySet()) {
            Block block = entry.getValue();
            if (!visited.contains(block.getSeq())) {
                contigs[i] = block.getSeq();
                i++;
                if (includeRC) {
                    contigs[i] = block.getRCSeq();
                    i++;
                }
                visited.add(block.getSeq());
                visited.add(block.getRCSeq());
            }
        }
        return contigs;
    }

    public int getMaxContigLength() {
        int max = 0;
        for (String seq : blocks.keySet()) {
            if (seq.length() > max) {
                max = seq.length();
            }
        }
        return max;
    }

    public float getAvgContigLength() {
        int sum = 0;
        for (String seq : blocks.keySet()) {
                sum += seq.length();
        }
        return (float)sum/blocks.size();
    }

    public int getMaxCount() {
        int max = 0;
        for (String seq : blocks.keySet()) {
            if (blocks.get(seq).count > max) {
                max = seq.length();
            }
        }
        return max;
    }

    @Override
    public String toString() {
        String s = "";
        for (Map.Entry<String, Block> entry : blocks.entrySet()) {
            s = s + entry.getValue().toString() + "\n";
        }
        return s;
    }
}
