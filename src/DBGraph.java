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

    /**
     * Maximum percentage of sequence length that is allowed to differ for
     * bubble removal (in terms of edit distance)
     */
    private static final double MAX_BUBBLE_DIFF = 0.1;

    /**
     * Maximum length of bubble path.
     */
    private static final int MIN_BUBBLE_LEN = 150;

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
     * and stores a pointer to its twin block (reverse complement).
     *
     * @author tohei
     */
    private class DBGNode {
        /**
         * outgoing edges (incoming edges are parallel or equal to the outgoing edges of the twin node)
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
        }

        /**
         * Get the in-degree of this node (which is equal to the out-degree of its twin).
         *
         * @return the number of incoming edges.
         */
        private int getIndegree() {
            int indegree = 0;
            if (twin != null) {
                indegree = twin.outgoing.size();
            }
            return indegree;
        }

        /**
         * Checks whether or not this node is a tip.
         * <p>
         * A tip is a linear node with in-degree or out-degree 1 and a sequence shorter than 2 * k.
         *
         * @return true if node is a tip.
         */
        private boolean isTip() {
            if (seq.length() > 2 * k) return false;
            else {
                if (getOutdegree() == 1 && getIndegree() == 0) return true;
                if (getIndegree() == 1 && getOutdegree() == 0) return true;
                return false;
            }
        }

        /**
         * Get the out-degree of this node.
         *
         * @return number of outgoing edges.
         */
        private int getOutdegree() {
            return outgoing.size();
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
         * If this node only represents a single k-mer, the coverage is equal to count (the number of reads
         * containing this node's sequence). Otherwise it is the average of the counts of all k-mers
         * represented in this node.
         *
         * @return the coverage of the contained sequence.
         */
        private float getCov() {
            return count / (float) getSize();
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
                if (forward.twin != null) targetNode.twin.outgoing.remove(forward.twin);
                return true;
            }
        }

        /**
         * Adds a new edge to the specified node or updates an existing one.
         * <p>
         * If addReverseEdge is true, this will also create the reverse edge from n.twin to this nodes's twin.
         * If the edges already exists, this method will either increase or overwrite their multiplicity
         * depending on the value specified in update. This method will not check if the edges to create are
         * valid connections in the de Bruijn graph!
         *
         * @param targetNode     DBGNode to create an edge to.
         * @param multiplicity   multiplicity of edge to create.
         * @param update         in case the edge already exists, add to edge multiplicity (true) or overwrite (false).
         * @param addReverseEdge add reverse edge (true) or not (false).
         */
        private void addEdgeTo(DBGNode targetNode, int multiplicity, boolean update, boolean addReverseEdge) {
            boolean edgeExists = false;
            for (DBGEdge edge : outgoing) {
                if (edge.target == targetNode) {
                    if (update) {
                        edge.multiplicity += multiplicity;
                        if (edge.twin != null) edge.twin.multiplicity += multiplicity;
                    } else {
                        edge.multiplicity = multiplicity;
                        if (edge.twin != null) edge.twin.multiplicity = multiplicity;
                    }
                    edgeExists = true;
                    break;
                }
            }
            if (!edgeExists) {
                DBGEdge forward = new DBGEdge();
                forward.target = targetNode;
                forward.multiplicity = multiplicity;
                DBGEdge backward = null;
                if (addReverseEdge) {
                    // add reverse edge
                    backward = new DBGEdge();
                    backward.twin = forward;
                    backward.target = this.twin;
                    backward.multiplicity = multiplicity;
                    targetNode.twin.outgoing.add(backward);
                }
                forward.twin = backward;
                this.outgoing.add(forward);
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
                if (forward.twin != null) forward.target.twin.outgoing.remove(forward.twin);
            }
            while (twin.outgoing.size() > 0) {
                // remove edge
                DBGEdge forward = twin.outgoing.remove();
                // remove reverse edge
                if (forward.twin != null) forward.target.twin.outgoing.remove(forward.twin);
            }
        }

        /**
         * Returns suffix of this node's sequence starting at position k-1.
         *
         * @return new DNAString containing the suffix at position k-1 of this node's sequence.
         */
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
            String s = "[Node: " + seq + "]" + " Cov: " + getCov();
            s = s + "\n Out:[";
            for (DBGEdge edge : outgoing) {
                s = s + edge.target.seq.toString() + "[" + edge.multiplicity + "]" + " ";
            }
            s = s + "]";
            return s;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Edges in the de Bruijn graph.
     * <p>
     * DBGEdges are directed edges used to connect nodes with overlapping sequences. Every edge
     * stores a pointer to its twin edge which connects the target twin node to the source twin.
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
                    + "-mers of " + reads.length + " reads:");
        }

        // initialize instance variables
        this.k = k;
        this.nodesMap = new LinkedHashMap<>(kmers.size());

        if (verbose) {
            System.out.println(" adding nodes ... ");
        }

        // add nodes
        for (DNAString kmer : kmers.keySet()) {
            // (if kmer is in nodesMap, so will be its complement)
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
                if (n != null) {
                    if (n == node.twin) {
                        node.addEdgeTo(n, 1, false, false);
                    } else {
                        node.addEdgeTo(n, 1, false, true);
                    }
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
                if (end == start.twin) {
                    start.addEdgeTo(end, 1, true, false);
                } else {
                    start.addEdgeTo(end, 1, true, true);
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
     * Produces a hash map of all the DNAStrings (nodes) in this graph and their coverage.
     *
     * @return LinkedHashMap mapping DNAstrings to coverage.
     */
    public LinkedHashMap<DNAString, Integer> getCovData() {
        LinkedHashMap<DNAString, Integer> counts = new LinkedHashMap<>(getSize() / 2);
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            counts.put(nodeSeq, (int) Math.ceil(n.getCov()));
        }
        return counts;
    }

    /**
     * Produces a SimpleVec vector holding all coverage values of nodes in this graph adjusted
     * by node sizes.
     * <p>
     * For every DBGNode n in this graph, this methods adds n.getSize() entries containing its average
     * coverage to the output array.
     *
     * @return node coverages amplified by node size.
     */
    public SimpleVec getAdjustedCovData() {
        LinkedList<Double> counts = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            for (int i = 0; i < n.getSize(); i++) {
                counts.add((double) n.getCov());
            }
        }
        double[] countsArray = new double[counts.size()];
        int i = 0;
        for (Double d : counts) {
            countsArray[i] = d;
            i++;
        }
        return new SimpleVec(countsArray);
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
     * Determine the maximum coverage of any node in the graph.
     *
     * @return maximum coverage of any node in the graph.
     */
    public float getMaxCov() {
        float max = 0;
        for (DBGNode node : nodesMap.values()) {
            if (node.getCov() > max) {
                max = node.getCov();
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
     * Remove bubbles from the graph.
     * <p>
     * This method removes bubbles which consist of a short non-branching path from a source node s to a target node t
     * and a second potentially nonlinear path from s to t. The algorithm creates a set of linear candidate paths and
     * calls findAndCutBubble on each of them.
     *
     * @param verbose be verbose
     * @return true, if bubbles were resolved, false otherwise.
     */
    public boolean removeBubbles(boolean verbose) {
        if (verbose) System.out.print("Removing bubbles ");
        long startTime = System.currentTimeMillis();
        int blocksPriorCollapse = getSize();

        LinkedList<DBGNode> simplePaths = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            // find linear paths between two junction nodes
            if (n.getIndegree() == 1 && n.getOutdegree() == 1) {
                if (n.outgoing.getFirst().target.getIndegree() > 1) {
                    if (n.twin.outgoing.getFirst().target.twin.getOutdegree() > 1) {
                        if (n.getSize() < MIN_BUBBLE_LEN) simplePaths.add(n);
                    }
                }
            }
        }

        boolean bubbleCollapsed = false;
        for (DBGNode path : simplePaths) {
            if (nodesMap.containsKey(path.seq)) {
                bubbleCollapsed = findAndCutBubble(path) || bubbleCollapsed;
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print("[" + (endTime - startTime) / 1000 + " seconds] ");
            if (bubbleCollapsed) {
                System.out.println("Graph size - before: " + blocksPriorCollapse + ", after: " + getSize());
            } else {
                System.out.println("No changes.");
            }
        }
        return bubbleCollapsed;
    }

    /**
     * Determine if the given linear node is part of a bubble and delete the path if possible.
     * <p>
     * This method searches for an alternative path connecting the ancestor of the given node to its
     * outgoing neighbor using Dijkstra's algorithm. The length of a directed edge is the number of
     * k-mers represented by the target node. If an alternative path is found, the given path node
     * is deleted if it has lower coverage and if the edit distance between the two paths is less
     * than 10% of the length of the shorter one.
     *
     * @param path node with one incoming and one outgoing edge.
     * @return true, if bubble was cut and false otherwise.
     */
    private boolean findAndCutBubble(DBGNode path) {
        DBGNode start = path.twin.outgoing.getFirst().target.twin;
        DBGNode goal = path.outgoing.getFirst().target;

        // find path from start to goal using Dijkstra's shortest path algorithm
        // TODO use memory efficient data structures
        LinkedList<DBGNode> toVisit = new LinkedList<>();
        LinkedHashMap<DBGNode, Integer> distances = new LinkedHashMap<>();
        LinkedHashMap<DBGNode, DBGNode> previous = new LinkedHashMap<>();

        for (DBGEdge edge : start.outgoing) {
            if (edge.target != path) {
                toVisit.add(edge.target);
                distances.put(edge.target, edge.target.getSize());
                previous.put(edge.target, start);
            }
        }
        distances.put(start, 0);
        previous.put(start, null);

        int minDistToVisit = 0;
        while ((toVisit.size() > 0) && (minDistToVisit < path.getSize() + 2)) {
            // find node with smallest distance
            int minIdx = 0;
            for (int i = 0; i < toVisit.size(); i++) {
                if (distances.get(toVisit.get(i)) < distances.get(toVisit.get(minIdx))) {
                    minIdx = i;
                }
            }
            DBGNode n = toVisit.remove(minIdx);
            minDistToVisit = distances.get(n);

            for (DBGEdge edge : n.outgoing) {
                if (!previous.containsKey(edge.target)) {
                    toVisit.add(edge.target);
                    distances.put(edge.target, distances.get(n) + edge.target.getSize());
                    previous.put(edge.target, n);
                } else {
                    if (distances.get(edge.target) > distances.get(n) + edge.target.getSize()) {
                        distances.put(edge.target, distances.get(n) + edge.target.getSize());
                        previous.put(edge.target, n);
                    }
                }
                // if alternative path is found, backtrack sequence
                if (edge.target == goal) {
                    DBGNode altPathWalker = n;
                    DNAString altSequence;
                    if (previous.get(altPathWalker) == start) {
                        altSequence = altPathWalker.seq;
                    } else {
                        altSequence = altPathWalker.getShortSeq();
                    }
                    float avgCoverage = altPathWalker.getCov();
                    int numOfNodes = 1;
                    while (previous.get(altPathWalker) != start) {
                        altPathWalker = previous.get(altPathWalker);
                        numOfNodes++;
                        avgCoverage += altPathWalker.getCov();
                        if (previous.get(altPathWalker) == start) {
                            altSequence = altPathWalker.seq.concat(altSequence);
                        } else {
                            altSequence = altPathWalker.getShortSeq().concat(altSequence);
                        }
                    }
                    avgCoverage = avgCoverage / numOfNodes;
                    // if sequences are similar enough, delete the simple path.
                    int maxDist = (int) (MAX_BUBBLE_DIFF * Math.min(path.seq.length(), altSequence.length()));
                    if (DNAStringUtils.LevDistance(path.seq, altSequence) < maxDist) {
                        if (avgCoverage > path.getCov()) {
                            removeBlock(path);
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    /**
     * Remove low coverage nodes and their twins.
     *
     * @param cutoff  coverage cutoff.
     * @param verbose be verbose.
     * @return true if nodes were removed, false otherwise.
     */
    public boolean removeLowCoverageBlocks(int cutoff, boolean verbose) {
        if (verbose) System.out.print("Apply cov cutoff ");
        long startTime = System.currentTimeMillis();

        LinkedList<DBGNode> nodesToRemove = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.getCov() < cutoff) {
                nodesToRemove.add(n);
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print("[" + (endTime - startTime) / 1000 + " seconds] ");
        }
        if (nodesToRemove.size() > 0) {
            if (verbose) System.out.print("Graph size - before: " + getSize());
            for (DBGNode node : nodesToRemove) {
                if (nodesMap.containsKey(node.seq)) {
                    removeBlock(node);
                }
            }
            if (verbose) System.out.println(", after: " + getSize());
            return true;
        } else {
            if (verbose) System.out.println("No changes");
            return false;
        }
    }

    /**
     * Remove tips from the graph.
     * <p>
     * A node-twin pair with only one outgoing edge is considered an erroneous tip if its sequence is shorter
     * than 2k and there exists an edge going from a junction to the graph with higher multiplicity than the
     * edge going from the junction to the tip.
     *
     * @param verbose be verbose.
     * @return true if nodes were removed, false otherwise.
     */
    public boolean removeTips(boolean verbose) {
        if (verbose) System.out.print("Removing tips    ");
        long startTime = System.currentTimeMillis();

        LinkedList<DBGNode> blocksToRemove = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (n.isTip()) {
                DBGEdge pathFromTip;
                DBGNode junction;
                if (n.getIndegree() == 1) {
                    pathFromTip = n.twin.outgoing.getFirst();
                    junction = pathFromTip.target.twin;
                } else {
                    pathFromTip = n.outgoing.getFirst();
                    junction = pathFromTip.target;
                }
                for (DBGEdge edge : junction.outgoing) {
                    if (edge.multiplicity > pathFromTip.multiplicity) {
                        blocksToRemove.add(n);
                        break;
                    }
                }
            }
        }
        int blocksPriorTipRemoval = getSize();
        boolean removedTips = false;
        for (DBGNode node : blocksToRemove) {
            // check that node still exists and is a tip
            if (nodesMap.containsKey(node.seq) && node.isTip()) {
                removedTips = true;
                removeBlock(node);
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print("[" + (endTime - startTime) / 1000 + " seconds] ");
        }
        if (removedTips) {
            if (verbose) System.out.println("Graph size - before: " + blocksPriorTipRemoval + ", after: " + getSize());
            return true;
        } else {
            if (verbose) System.out.println("No changes.");
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
        if (verbose) System.out.print("Collapsing graph ");
        long startTime = System.currentTimeMillis();
        int blocksPriorCollapse = getSize();

        // find potential first nodes in a linear sequence
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
                    if (ancestor.getOutdegree() > 1 || ancestor == n.twin) {
                        startingBlocks.add(n);
                    }
                }
            }
        }

        boolean collapsed = false;
        for (DBGNode start : startingBlocks) {
            if (nodesMap.containsKey(start.seq)) {
                LinkedList<DBGNode> linearStretch = expandLinearStretch(start);
                if (linearStretch.size() > 1) {
                    collapseBlocks(linearStretch);
                    collapsed = true;
                }
            }
        }
        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.print("[" + (endTime - startTime) / 1000 + " seconds] ");
            if (collapsed) {
                System.out.println("Graph size - before: " + blocksPriorCollapse + ", after: " + getSize());
            } else {
                System.out.println("No changes.");
            }
        }
        return collapsed;
    }


    /**
     * Expand linear stretch starting with the given node as far as possible.
     *
     * @param start first node of a potential linear stretch with exactly one outgoing edge.
     * @return LinkedList containing the found linear stretch.
     */
    private LinkedList<DBGNode> expandLinearStretch(DBGNode start) {
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
        return linearStretch;
    }

    /**
     * Collapse list of linear nodes into single super-node.
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
            if (edge.target == first) {
                newTwin.addEdgeTo(newNode, edge.multiplicity, false, false);
            } else {
                newTwin.addEdgeTo(edge.target, edge.multiplicity, false, true);
            }
        }
        for (DBGEdge edge : last.outgoing) {
            if (edge.target == last.twin) {
                newNode.addEdgeTo(newTwin, edge.multiplicity, false, false);
            } else {
                newNode.addEdgeTo(edge.target, edge.multiplicity, false, true);
            }
        }

        // detach first and last node from the rest of the graph
        first.isolate();
        last.isolate();

        // remove old nodes
        for (DBGNode node : linearStretch) {
            nodesMap.remove(node.seq);
            nodesMap.remove(node.twin.seq);
        }

        // should be impossible as every node seq starts with unique k-mer
        if (nodesMap.containsKey(newSeq)) {
            System.err.println("WARNING: overwriting existing node");
        }
        // should be impossible when starting with odd k-mer size
        if (newSeq.equals(newSeqRC)) {
            System.err.println("WARNING: node equals its own twin");
        }

        // add new nodes
        nodesMap.put(newSeq, newNode);
        nodesMap.put(newSeqRC, newTwin);
    }

    /**
     * Simplifies the graph by repeatedly applying removeTips, removeBubbles, removeLowCoverageBlocks and collapse.
     *
     * @param correctErrors use error correction (true) or just collapse (false).
     * @param cutoffCov     cutoff threshold for removeLowCoverageBlocks.
     * @param verbose       be verbose.
     */
    public void simplify(boolean correctErrors, int cutoffCov, boolean verbose) {
        if (verbose) {
            String s;
            if (correctErrors) s = "with";
            else s = "without";
            System.out.println("Simplifying graph (" + s + " error correction):");
        }
        int iter = 1;
        boolean modified = true;

        while (modified) {
            if (verbose) {
                System.out.println("\t------------ Iteration " + iter + " ------------");
            }
            modified = collapse(verbose);
            if (correctErrors) {
                modified = removeTips(verbose) || modified;
                modified = removeBubbles(verbose) || modified;
            }
            iter++;
        }
        if (verbose) {
            System.out.println("\t-------------------------------------");
        }
        if (correctErrors) {
            removeLowCoverageBlocks(cutoffCov, verbose);
            collapse(verbose);
        }
        if (verbose) {
            System.out.println("\nDone. Number of twin-nodes: " + getSize() / 2);
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
        DNAString[] sequences = getSequences(false, 0);
        for (DNAString nodeSeq : sequences) {
            DBGNode n = nodesMap.get(nodeSeq);
            s = s + "-----------------------------------\n";
            s = s + n.toString() + "\n";
            s = s + n.twin.toString() + "\n";
        }
        s = s + "-----------------------------------\n";
        return s;
    }
}
