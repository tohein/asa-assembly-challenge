import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Simple implementation of a de Brujin graph for
 * sequence assembly.
 * <p>
 *
 * @author tohei
 */
public class DBGraph {

    //TODO choose different data structure to avoid problems with duplicate nodes (node.seq = node.twin.seq) created by collapse
    /**
     * Stores all nodes in the graph.
     */
    private LinkedHashMap<DNAString, DBGNode> nodesMap;

    /**
     * k-mer size
     */
    private int k;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Nodes in the de Bruijn graph.
     * <p>
     * Every node corresponds to a k-mer or a simplified linear stretch.
     * and stores a pointer to its twin block (reverse complement). The duality of nodes edges
     * leads to a bipartite graph.
     *
     * @author tohei
     */
    private class DBGNode implements Comparable<DBGNode> {
        /**
         * outgoing edges (incoming edges are parallel to the outgoing edges of the twin node)
         */
        private LinkedList<DBGEdge> outgoing;
        /**
         * twin node
         */
        private DBGNode twin;
        /**
         * node sequence
         */
        private DNAString seq;
        /**
         * coverage
         */
        private int count;

        /**
         * used by search algorithm for collapsing bubbles
         */
        private float dist;

        /**
         * Create new DBGNode from a sequence and its read count.
         *
         * @param seq   sequence to store in the new node.
         * @param count number of reads containing any of the k-mers represented in this node.
         */
        public DBGNode(DNAString seq, int count) {
            outgoing = new LinkedList<>();
            this.twin = null;
            this.seq = seq;
            this.count = count;
            this.dist = 0;
        }

        /**
         * Get the in-degree of this node (which is equal to the out-degree of its twin).
         *
         * @return the number of incoming edges.
         */
        private int getIndegree() {
            int indegree = 0;
            if (this.twin != null) {
                indegree = this.twin.outgoing.size();
            }
            return indegree;
        }

        /**
         * Get the out-degree of this node.
         *
         * @return number of outgoing edges.
         */
        private int getOutdegree() {
            return this.outgoing.size();
        }

        /**
         * Compute the size of this DBGNode, i.e. the number of k-mers represented by this node.
         *
         * @return the number of k-mers represented by this node.
         */
        private int getSize() {
            return seq.length() - k + 1;
        }

        /**
         * Computes the coverage of this node.
         * <p>
         * If this node only represents a single k-mer, the coverage is equal to count. Otherwise
         * it is the average of the counts of all k-mers represented in this node.
         *
         * @return the coverage of the contained sequence.
         */
        private float getCov() {
            return count / (float) getSize();
        }

        /**
         * Adds a new edge to the specified node or updates an existing one.
         * <p>
         * As all edges are mirrored by a twin edge, this will also create the reverse edge from n.twin
         * to this nodes's twin. If the edges already exists, this method will either increase or
         * overwrite their multiplicity depending on the value specified in update.
         *
         * @param targetNode   DBGNode to create an edge to.
         * @param multiplicity multiplicity of edge to create.
         * @param update       in case the edge already exists, add to edge multiplicity (true) or overwrite (false)
         */
        private void addEdgeTo(DBGNode targetNode, int multiplicity, boolean update) {
            boolean edgeExists = false;
            for (DBGEdge edge : outgoing) {
                if (edge.target == targetNode) {
                    if (update) {
                        edge.multiplicity += multiplicity;
                        edge.twin.multiplicity += multiplicity;
                    } else {
                        edge.multiplicity = multiplicity;
                        edge.twin.multiplicity = multiplicity;
                    }
                    edgeExists = true;
                    break;
                }
            }
            if (!edgeExists) {
                DBGEdge forward = new DBGEdge();
                DBGEdge backward = new DBGEdge();
                forward.twin = backward;
                forward.target = targetNode;
                forward.multiplicity = multiplicity;
                backward.twin = forward;
                backward.target = this.twin;
                backward.multiplicity = multiplicity;
                this.outgoing.add(forward);
                targetNode.twin.outgoing.add(backward);
            }
        }

        /**
         * Removes an edge and its twin reverse-edge.
         *
         * @param targetNode target node.
         * @return true if the edge existed, false otherwise.
         */
        private boolean rmvEdgeTo(DBGNode targetNode) {
            DBGEdge forward = null;
            for (DBGEdge edge : outgoing) {
                if (edge.target == targetNode) {
                    forward = edge;
                    break;
                }
            }
            if (forward == null) {
                return false;
            } else {
                outgoing.remove(forward);
                targetNode.twin.outgoing.remove(forward.twin);
                return true;
            }
        }

        /**
         * Removes all incoming and outgoing edges from this node and its twin.
         */
        private void isolate() {
            while (outgoing.size() > 0) {
                // remove edge
                DBGEdge forward = outgoing.remove();
                // remove reverse edge
                forward.target.twin.outgoing.remove(forward.twin);
            }
            while (twin.outgoing.size() > 0) {
                // remove edge
                DBGEdge forward = twin.outgoing.remove();
                // remove reverse edge
                forward.target.twin.outgoing.remove(forward.twin);
            }
        }

        //TODO find better way of computing extended sequences
        private DNAString getShortSeq() {
            return seq.subSequence(k - 1, seq.length());
        }


        /**
         * Generates a String representation of this node.
         *
         * @return String representation of this node.
         */
        @Override
        public String toString() {
            String twinSeq;
            if (twin == null) {
                twinSeq = "No twin";
            } else {
                twinSeq = twin.seq.toString();
            }
            String s = "[Node: " + seq + "]";
            s = s + "\n[Twin:  " + twinSeq + "]";
            s = s + "\nOut:[";
            for (DBGEdge edge : outgoing) {
                s = s + edge.target.seq.toString() + " ";
            }
            s = s + "]";
            return s;
        }

        @Override
        public int compareTo(DBGNode n) {
            return Float.compare(this.dist, n.dist);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Edges in the de Bruijn graph.
     * <p>
     * DBGEdges are directed edges which connect overlapping sequences. Every edge stores a pointer to its twin
     * edge which connects the target twin node to the source twin. The duality of nodes edges
     * leads to a bipartite graph.
     *
     * @author tohei
     */
    private class DBGEdge {
        /**
         * Target node of the directed DBGEdge.
         */
        private DBGNode target;
        /**
         * Number of reads using this edge.
         */
        private int multiplicity;
        /**
         * Twin edge.
         */
        private DBGEdge twin;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Create new de Bruijn Graph (DBGraph) object from an array of DNAString reads for a given k-mer length k.
     * <p>
     * The graph is created in 3 steps:
     * (1) unique k-mers are counted in the set of reads and turned into new DBGNodes
     * (2) nodes are connected based on k-1 overlap of their sequences
     * (3) node multiplicities are computed by iterating over the read array
     *
     * @param reads   array of DNAString reads.
     * @param kmers   LinkedHashSet mapping k-mers to counts.
     * @param k       k-mer size for creating this new graph.
     * @param verbose be verbose.
     */
    public DBGraph(DNAString[] reads, Map<DNAString, Integer> kmers, int k, boolean verbose) {

        // track time needed to create graph
        long startTime = System.currentTimeMillis();
        if (verbose) {
            System.out.println("Building de Bruijn graph from " + k
                    + "-mers of " + reads.length + " reads ...");
        }

        // initialize instance variables
        this.k = k;
        this.nodesMap = new LinkedHashMap<>(kmers.size());

        if (verbose) {
            System.out.println(" adding nodes ... ");
        }

        // add nodesMap
        for (DNAString kmer : kmers.keySet()) {
            // (if kmer is in nodesMap, so is its complement)
            if (!nodesMap.containsKey(kmer)) {
                DNAString kmerRC = kmer.reverseComplement();
                int count = kmers.get(kmer);

                // create new node and its twin
                DBGNode newNode = new DBGNode(kmer, count);
                DBGNode newTwin = new DBGNode(kmerRC, count);
                newNode.twin = newTwin;
                newTwin.twin = newNode;
                nodesMap.put(kmer, newNode);
                nodesMap.put(kmerRC, newTwin);
            }
        }

        if (verbose) {
            System.out.println(" adding edges ... ");
        }

        // add edges
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode node = nodesMap.get(nodeSeq);
            for (DNAString s : nodeSeq.subSequence(1, k).allVariations(k)) {
                DBGNode n = nodesMap.get(s);
                // do not add edges to twin node or itself
                if (n != null && n != node && n != node.twin) {
                    node.addEdgeTo(n, 1, false);
                }
            }
        }

        if (verbose) {
            System.out.println(" computing edge multiplicities ... ");
        }

        // follow read paths in the graph and update edge multiplicities
        for (DNAString s : reads) {
            DNAString endSeq = s.subSequence(0, k);
            for (int i = 0; i < s.length() - k; i++) {
                DNAString startSeq = endSeq;
                endSeq = s.subSequence(i + 1, i + 1 + k);
                DBGNode start = nodesMap.get(startSeq);
                DBGNode end = nodesMap.get(endSeq);
                if (start != end && start != end.twin) {
                    start.addEdgeTo(end, 1, true);
                }
            }
        }

        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime) / 1000 +
                    " seconds). " + getSize() + " nodes.\n");
        }
    }

    /**
     * Get number of nodes in the graph.
     *
     * @return number of nodes in the graph.
     */
    public int getSize() {
        return nodesMap.size();
    }

    /**
     * Get k-mer length.
     *
     * @return k-mer length used to create the graph.
     */
    public int getK() {
        return k;
    }

    /**
     * Remove a node and its twin from the graph.
     * <p>
     * This method first isolates the block (removes all incoming and outgoing edges). Then the node and its twin
     * are removed from nodesMap.
     *
     * @param n DBGNode to remove.
     */
    private void removeBlock(DBGNode n) {
        n.isolate();
        nodesMap.remove(n.seq);
        nodesMap.remove(n.twin.seq);
    }

    /**
     * Check if the graph contains a node corresponding to the given sequence.
     *
     * @param s DNAString to look for.
     * @return true if the graph contains a node with the argument sequence.
     */
    public boolean contains(DNAString s) {
        return nodesMap.containsKey(s);
    }

    /**
     * Get all node sequences from the graph.
     *
     * @param includeRC         include the reverse complement of every sequence.
     * @param minSequenceLength minimum sequence length.
     * @return DNAString array containing all node sequences.
     */
    public DNAString[] getSequences(boolean includeRC, int minSequenceLength) {
        LinkedHashSet<DNAString> visited = new LinkedHashSet<>(nodesMap.size());
        LinkedList<DNAString> contigs = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (!visited.contains(n.seq) && n.seq.length() > minSequenceLength) {
                contigs.add(n.seq);
                if (includeRC) {
                    contigs.add(n.twin.seq);
                }
                visited.add(n.seq);
                visited.add(n.twin.seq);
            }
        }
        return contigs.toArray(new DNAString[contigs.size()]);
    }

    /**
     * Get all node sequences from the graph.
     *
     * @return DNAString array containing all node sequences.
     */
    public DNAString[] getSequences() {
        return getSequences(true, 0);
    }

    /**
     * Find the length of the largest sequence.
     *
     * @return maximum sequence length.
     */
    public int getMaxSequenceLength() {
        int max = 0;
        for (DNAString nodeSeq : nodesMap.keySet()) {
            if (nodeSeq.length() > max) {
                max = nodeSeq.length();
            }
        }
        return max;
    }

    /**
     * Compute the average length of the sequences in this graph.
     *
     * @return average sequence length
     */
    public float getAvgSequenceLength() {
        int sum = 0;
        for (DNAString nodeSeq : nodesMap.keySet()) {
            sum += nodeSeq.length();
        }
        return (float) sum / nodesMap.size();
    }

    /**
     * Compute the maximum coverage of any node in the graph.
     *
     * @return maximum coverage of all nodes in the graph.
     */
    public float getMaxCov() {
        float max = 0;
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.getCov() > max) {
                max = n.getCov();
            }
        }
        return max;
    }

    /**
     * Remove bubbles from the graph.
     *
     * @param verbose be verbose
     * @return true, if bubbles were resolved, false otherwise.
     */
    public boolean removeBubbles(boolean verbose) {
        long startTime = System.currentTimeMillis();
        int blocksPriorCollapse = this.getSize();

        if (verbose) {
            System.out.print("Removing bubbles. ");
        }
        LinkedList<DBGNode> startingNodes = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.getOutdegree() > 1) {
                startingNodes.add(n);
            }
        }
        if (verbose) {
            System.out.println("Removing bubbles. Number of starting nodes: " + startingNodes.size());
        }

        boolean bubbleCollapsed = false;
        startingNodeLoop:
        for (DBGNode start : startingNodes) {
            // check that start was not yet modified
            if (!nodesMap.containsKey(start.seq) || start.getOutdegree() < 2) {
                continue startingNodeLoop;
            }
            for (DBGEdge edge1 : start.outgoing) {
                for (DBGEdge edge2 : start.outgoing) {
                    if (edge1 != edge2) {
                        DBGNode neighbor1 = edge1.target, neighbor2 = edge2.target;
                        if (neighbor1.getOutdegree() == 1 && neighbor2.getOutdegree() == 1) {
                            if (neighbor1.outgoing.getFirst().target == neighbor2.outgoing.getFirst().target) {
                                if (removeSimpleBubble(start, neighbor1, neighbor2)) {
                                    bubbleCollapsed = true;
                                    continue startingNodeLoop;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print(" Removed bubbles (" + (endTime - startTime) / 1000 + " seconds). ");
            if (bubbleCollapsed) {
                System.out.println("Graph size - before: " + blocksPriorCollapse + ", after: " + getSize());
            } else {
                System.out.println("No changes.");
            }
        }
        return bubbleCollapsed;
    }

    /**
     * Removes a basic bubble consisting of four nodes.
     * <p>
     *
     * @param top top path node.
     * @param bot bottom path node.
     * @return true if bubble was collapsed.
     */
    private boolean removeSimpleBubble(DBGNode start, DBGNode top, DBGNode bot) {
        DNAString topSeq = top.seq;
        DNAString botSeq = bot.seq;
        int topEdgeScore = 0;
        int botEdgeScore = 0;
        for (DBGEdge edge : start.outgoing) {
            if (edge.target == top) topEdgeScore = edge.multiplicity;
            if (edge.target == bot) botEdgeScore = edge.multiplicity;
        }
        DBGNode nodeToRemove = top.getSize() / topEdgeScore < bot.getSize() / botEdgeScore ? bot : top;
        if (top.getSize() > bot.getSize()) { // top longer
            topSeq = top.seq.subSequence(0, bot.seq.length());
            botSeq = bot.seq;
        } else if (bot.getSize() > top.getSize()) { // bottom longer
            botSeq = bot.seq.subSequence(0, top.seq.length());
            topSeq = top.seq;
        }
        int maxDist = (int) (0.1 * topSeq.length());
        if (DNAStringUtils.LevDistance(topSeq, botSeq) < maxDist) {
            removeBlock(nodeToRemove);
            return true;
        }
        return false;
    }

    /**
     * Remove low coverage nodes and their twins (chimeric connections).
     *
     * @param cutoff  coverage cutoff.
     * @param verbose be verbose.
     * @return true if nodes were removed, false otherwise.
     */
    public boolean removeLowCoverageBlocks(int cutoff, boolean verbose) {
        if (verbose) {
            System.out.print("Removing low coverage nodes. ");
        }
        LinkedList<DBGNode> nodesToRemove = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.getCov() < cutoff) {
                nodesToRemove.add(n);
            }
        }
        if (nodesToRemove.size() > 0) {
            System.out.print("Graph size - before: " + this.getSize());
            for (DBGNode node : nodesToRemove) {
                if (nodesMap.containsKey(node.seq)) {
                    removeBlock(node);
                }
            }
            System.out.println(", after: " + this.getSize());
            return true;
        } else {
            System.out.println("No changes");
            return false;
        }
    }

    /**
     * Remove tips from the graph.
     * <p>
     * A node-twin pair with only one outgoing edge is considered an erroneous tip if its sequence is shorter
     * than 2k and there exists an edge going from the junction to the graph with higher multiplicity than the
     * edge going from the junction to the tip.
     *
     * @param verbose be verbose.
     * @return true if nodes were removed, false otherwise.
     */
    public boolean removeTips(boolean verbose) {
        if (verbose) {
            System.out.print("Removing Tips. ");
        }
        LinkedList<DBGNode> blocksToRemove = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.seq.length() < 2 * k) { // probably one case is sufficient
                if (n.getIndegree() == 0 && n.getOutdegree() == 1) {
                    DBGEdge pathFromTip = n.outgoing.getFirst();
                    for (DBGEdge edge : pathFromTip.target.outgoing) {
                        if (edge.multiplicity > pathFromTip.multiplicity) {
                            blocksToRemove.add(n);
                        }
                    }
                } else if (n.getOutdegree() == 0 && n.getIndegree() == 1) {
                    DBGEdge pathFromTip = n.twin.outgoing.getFirst();
                    for (DBGEdge edge : pathFromTip.target.outgoing) {
                        if (edge.multiplicity > pathFromTip.multiplicity) {
                            blocksToRemove.add(n);
                        }
                    }
                }
            }
        }
        if (blocksToRemove.size() > 0) {
            System.out.print("Graph size - before: " + this.getSize());
            for (DBGNode node : blocksToRemove) {
                this.removeBlock(node);
            }
            System.out.println(", after: " + this.getSize());
            return true;
        } else {
            System.out.println("No changes");
            return false;
        }
    }

    /**
     * Collapse all linear stretches in the graph.
     *
     * @param verbose be verbose
     * @return true if a linear stretch was collapsed.
     */
    public boolean collapse(boolean verbose) {
        long startTime = System.currentTimeMillis();
        int blocksPriorCollapse = this.getSize();

        LinkedList<DBGNode> startingBlocks = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.getOutdegree() == 1) {
                if (n.getIndegree() != 1) {
                    // multiple or no ancestors
                    startingBlocks.add(n);
                } else {
                    // single ancestor
                    DBGEdge backward = n.twin.outgoing.getFirst();
                    DBGNode ancestor = backward.target.twin;
                    // check if ancestor is a branching node
                    if (ancestor.getOutdegree() > 1) {
                        startingBlocks.add(n);
                    }
                }
            }
        }
        if (verbose) {
            System.out.println("Collapsing graph. Number of starting nodes: " + startingBlocks.size());
        }
        boolean collapsed = false;
        for (DBGNode start : startingBlocks) {
            if (nodesMap.containsKey(start.seq)) {
                collapsed = collapseLinearStretch(start) || collapsed;
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print(" Collapsed de Bruijn graph (" + (endTime - startTime) / 1000 + " seconds). ");
            if (collapsed) {
                System.out.println("Graph size - before: " + blocksPriorCollapse + ", after: " + getSize());
            } else {
                System.out.println("No changes.");
            }
        }
        return collapsed;
    }


    /**
     * Check if a given node is the beginning of a linear stretch and collapse if possible.
     *
     * @param start first node of a potential linear stretch with exactly one outgoing edge.
     * @return true if start was the beginning of a linear stretch.
     */
    private boolean collapseLinearStretch(DBGNode start) {
        boolean collapsed = false;
        LinkedList<DBGNode> linearStretch = new LinkedList<>();
        linearStretch.add(start);

        // all start nodesMap have exactly one outgoing edge
        DBGNode n = start.outgoing.getFirst().target;
        while ((n.getIndegree() == 1) && !(linearStretch.contains(n) || linearStretch.contains(n.twin))) {
            linearStretch.add(n);
            if (n.getOutdegree() == 1) {
                n = n.outgoing.getFirst().target;
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

    /**
     * Collapse list of linear nodes into single supernode.
     *
     * @param linearStretch LinkedList of DBGNodes to collapse.
     */
    private void collapseBlocks(LinkedList<DBGNode> linearStretch) {
        DBGNode first = linearStretch.getFirst();
        DBGNode last = linearStretch.getLast();

        // TODO: use byte array
        DNAString shortSeq = new DNAString();
        int count = 0;
        for (DBGNode block : linearStretch) {
            count = count + block.count;
            if (block == first) continue;
            shortSeq = shortSeq.concat(block.getShortSeq());
        }
        DNAString shortSeqRC = shortSeq.reverseComplement();
        DNAString newSeq = first.seq.concat(shortSeq);
        DNAString newSeqRC = shortSeqRC.concat(first.twin.seq);

        DBGNode newNode = new DBGNode(newSeq, count);
        DBGNode newTwin = new DBGNode(newSeqRC, count);
        newNode.twin = newTwin;
        newTwin.twin = newNode;

        // connect new node to the rest of the graph
        // (reverse edges are added automatically)
        for (DBGEdge edge : first.twin.outgoing) {
            newTwin.addEdgeTo(edge.target, edge.multiplicity, false);
        }
        for (DBGEdge edge : last.outgoing) {
            newNode.addEdgeTo(edge.target, edge.multiplicity, false);
        }

        // detach first and last node from the rest of the graph
        first.isolate();
        last.isolate();

        // remove old nodes and add new nodes
        for (DBGNode node : linearStretch) {
            nodesMap.remove(node.seq);
            nodesMap.remove(node.twin.seq);
        }

        // should be impossible as every node seq starts with unique k-mer
        if (nodesMap.containsKey(newSeq)) {
            System.err.println("ERROR: OVERWRITING EXISTING NODE");
        }
        if (newSeq.equals(newSeqRC)) {
            System.err.println("ERROR: NODE EQUALS TWIN");
        }

        nodesMap.put(newSeq, newNode);
        nodesMap.put(newSeqRC, newTwin);
    }

    /**
     * Simplifies the graph by repeatedly applying removeTips, removeLowCoverageBlocks and collapse.
     *
     * @param correctErrors apply removeTips and removeLowCoverageBlocks and collapse (true) or just collapse (false)
     * @param cutoffCov     cutoff threshold for removeLowCoverageBlocks
     * @param verbose       be verbose
     */
    public void simplify(boolean correctErrors, int cutoffCov, boolean verbose) {
        if (verbose) {
            String s;
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
                modified = this.removeTips(verbose) || modified;
                modified = this.removeBubbles(verbose) || modified;
                modified = this.removeLowCoverageBlocks(cutoffCov, verbose) || modified;
            }
            iter++;
        }
    }

    /**
     * Generate String representation of the graph (for debugging purposes).     *
     *
     * @return String representation of this DBGraph.
     */
    @Override
    public String toString() {
        String s = "";
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            s = s + n.toString() + "\n";
        }
        return s;
    }


    public static void main(String[] args) {
        String inputFile = "/home/tohei/Data/ASAData/reads_complex.fasta";
        String outputFile = "/home/tohei/Data/ASAData/out_new.fasta";

        System.out.println("Reading files ... ");
        DNAString[] inputs = null;
        try {
            inputs = DNAStringUtils.readFasta(inputFile);
        } catch (FileNotFoundException e1) {
            System.err.println("Could not find input file.");
        } catch (IOException e2) {
            System.err.println("Failed to read input file.");
        }
        int k = 21;

        System.out.println("Error correction ....");
        ReadCorrector rcor = new ReadCorrector();
        rcor.setReads(inputs, k);

        rcor.setCutoff(3);
        rcor.computeSAC(true);
        rcor.setCutoff(2);
        rcor.computeSAC(true);

        DBGraph G = new DBGraph(inputs, rcor.getCounts(), k, true);

        G.simplify(true, 2, true);
        DNAString[] contigs = G.getSequences(false, k);

        System.out.println();
        int max = G.getMaxSequenceLength();
        System.out.println("Max contig length: " + max);
        System.out.println("Avg contig length: " + G.getAvgSequenceLength());
        System.out.println();

        // save
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
