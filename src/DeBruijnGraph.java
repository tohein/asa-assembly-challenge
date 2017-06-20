import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.NoSuchElementException;


public class DeBruijnGraph {
	private LinkedHashMap<String, Node> nodes;
	private int k; 
	// TODO add reads array and store positions in nodes
	
	private class Node {
		private LinkedList<Node> outgoing;
		private Node twin;
		private String seq;
		
		public Node(String seq) {
			outgoing = new LinkedList<Node>();
			this.twin = null;
			this.seq = seq;
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
		
		public boolean addOutgoing(Node n) {
			if (!outgoing.contains(n)) {
				outgoing.add(n);
				return true;
			}
			return false;
		}
		
		public boolean rmvOutgoing(Node n) {
			return outgoing.remove(n);
		}
		
		public boolean addIncoming(Node n) {
			if (!twin.outgoing.contains(n)) {
				twin.outgoing.add(n);
				return true;
			}
			return false;
		}
		
		public boolean rmvIncoming(Node n) {
			return this.twin.outgoing.remove(n);

		}
		
		public String getSeq() {
			return seq;
		}
		
		public String getRCSeq() {
			return twin.seq;
		}
		
		@Override
		public String toString() {
			String twinSeq = "";
			if (twin == null) {
				twinSeq = "No twin";
			} else {
				twinSeq = twin.seq;
			}
			String s = "[Node: " + seq + "]";
			s = s + "\n[Twin: " + twinSeq + "]";
			s = s + "\nOut:[";
			for (Node n : outgoing) {
				s = s + n.seq + " ";
			}
			s = s + "]";
			return s;
		}
	}
	
	public int getSize() {
		return nodes.size();
	}
	
	public int getK() {
		return k;
	}
	
//////////////////////////////////////////////////////////////////////////

	public DeBruijnGraph(String[] inputs, int k) {
		this(inputs, k, false);
	}
	
	public DeBruijnGraph(String[] inputs, int k, boolean verbose) {
		// get unique k-mers
		// TODO remove reverse complements 
		long startTime = System.currentTimeMillis();
		if (verbose) {
			System.out.print("Building de Bruijn graph from " + k + "-mers of " + inputs.length + " reads ...");			
		}
		LinkedHashSet<String> kmers = SeqUtils.allKmers(inputs, k);
		nodes = new LinkedHashMap<String, Node>(kmers.size());
		this.k = k;

		if (verbose) {
			System.out.print(" adding nodes ...");			
		}
		// add nodes and rc twins
		for (String kmer : kmers) {
			String kmerRC = SeqUtils.reverseComplement(kmer);
			Node n1; Node n2;
			if (!nodes.containsKey(kmer)) {				
				n1 = new Node(kmer);
				n2 = new Node(kmerRC);
				n1.twin = n2;
				n2.twin = n1;
				
				nodes.put(kmer, n1);
				nodes.put(kmerRC, n2);
			}
		}		
		if (verbose) {
			System.out.print(" adding edges ...");			
		}
		// add edges
		for (Map.Entry<String, Node> entry : nodes.entrySet()) {
			String seq = entry.getKey();
			for (String s : SeqUtils.extendBack(seq.substring(1, k))) {
				Node n = nodes.get(s);
				if (n != null && n != entry.getValue().twin) {
					entry.getValue().addOutgoing(n);
				}
			}
		}
		if (verbose) {
			long endTime = System.currentTimeMillis();
			System.out.print(" done (" + (endTime - startTime)/1000 + 
					" seconds).\nTotal of " + getSize() + " nodes.\n");			
		}
	}

//////////////////////////////////////////////////////////////////////////
	
	public boolean collapse(boolean verbose) {
		long startTime = System.currentTimeMillis();
		int nodesPriorCollapse = this.getSize();
		
		LinkedList<Node> startingNodes = new LinkedList<Node>();
		for (Map.Entry<String, Node> entry : nodes.entrySet()) {
			Node n = entry.getValue();
			if (n.getOutdegree() == 1) {
				if (n.getIndegree() != 1) { 
					// multiple or no ancestors
					startingNodes.add(n);
				} else { 
					// single ancestor
					Node ancestor = (n.twin.outgoing.iterator().next()).twin;
					if (ancestor.outgoing.size() > 1) {
						startingNodes.add(n);
					}
				}
			}
		}
		if (verbose) {
			System.out.println("Collapsing graph. Number of starting nodes: "+startingNodes.size());
		}
		boolean collapsed = false;
		for (Node start : startingNodes) {
			if (nodes.containsKey(start.seq)) {		
				boolean startNodeCollapse = collapseLinearStretches(start);	
				collapsed = collapsed || startNodeCollapse;					
			}
		}
		if (verbose) {
			long endTime = System.currentTimeMillis();
			System.out.println("Collapsed de Bruijn graph in " + (endTime - startTime)/1000 + 
					" seconds. Graph size - before: "+nodesPriorCollapse+", after: "+getSize());
		}
		return collapsed;
	}
	
	private boolean collapseLinearStretches(Node start) {
		boolean collapsed = false;
		LinkedList<Node> linearStretch = new LinkedList<Node>();
		linearStretch.add(start);
		
		// all start nodes have exactly one outgoing edge
		Node n = start.outgoing.iterator().next();
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
			collapseNodes(linearStretch);
		}
		return collapsed;
	}	
	
	private void collapseNodes(LinkedList<Node> linearStretch) {
		Node n = linearStretch.removeFirst();
		Node m = linearStretch.getLast();
		
		String newSeq = n.seq;
		for (Node node : linearStretch) {
			newSeq = newSeq +  node.seq.substring(k-1);
		}
		String newTwinSeq = SeqUtils.reverseComplement(newSeq);
		Node newNode = new Node(newSeq);
		Node newTwin = new Node(newTwinSeq);
		newNode.twin = newTwin;
		newTwin.twin = newNode;
		// test for loop from n.twin to n
		//boolean ntwinloop = n.twin.rmvOutgoing(n);
		// test for loop from m to m.twin
		//boolean mtwinloop = m.rmvOutgoing(m.twin);
		
		newNode.outgoing = m.outgoing;
		newTwin.outgoing = n.twin.outgoing;

		// update incoming to n
		for (Node node : newTwin.outgoing) {
			node.twin.rmvOutgoing(n);
			node.twin.addOutgoing(newNode);
		}		
		// update outgoing from m
		for (Node node : newNode.outgoing) {
			node.twin.rmvOutgoing(m.twin);
			node.twin.addOutgoing(newTwin);
		}
		
/*		if (ntwinloop) {
			newNode.addIncoming(newTwin);
		}
		if (mtwinloop) {
			newNode.addOutgoing(newTwin);
		}*/
		
		// remove old nodes and add new nodes
		nodes.remove(n.seq);
		nodes.remove(n.twin.seq);
		for (Node node : linearStretch) {
			nodes.remove(node.seq);
			nodes.remove(node.twin.seq);			
		}
		nodes.put(newSeq, newNode);
		nodes.put(newTwinSeq, newTwin);
	}
	
//////////////////////////////////////////////////////////////////////////

	public String[] getIncoming(String s) {
		Node m = nodes.get(s);		
		if (m == null) {
			throw new NoSuchElementException("String not found.");
		} else {
			m = m.twin;
			String[] out = new String[m.getOutdegree()];
			int j = 0;
			for (Node n : m.outgoing) {
				out[j] = n.getRCSeq();
				j++;
			}
			return out;
		}
	}
	
	public String[] getOutgoing(String s) {
		Node m = nodes.get(s);
		if (m == null) {
			throw new NoSuchElementException("String not found.");
		} else {
			String[] out = new String[m.getOutdegree()];
			int j = 0;
			for (Node n : m.outgoing) {
				out[j] = n.getSeq();
				j++;
			}
			return out;
		}
	}
	
	public boolean contains(String s) {
		return nodes.containsKey(s);
	}
	
	public String[] getContigs(boolean includeRC) {
		LinkedHashSet<String> visited = new LinkedHashSet<String>(nodes.size());
		String[] contigs = null;
		if (includeRC) {
			contigs = new String[nodes.size()];
		} else {
			contigs = new String[nodes.size()/2];
		}
		int i = 0;
		for (Map.Entry<String, Node> entry : nodes.entrySet()) {				
			Node node = entry.getValue();
			if (!visited.contains(node.getSeq())) {
				contigs[i] = node.getSeq();
				i++;
				if (includeRC) {
					contigs[i] = node.twin.getSeq();
					i++;
				}
				visited.add(node.getSeq());
				visited.add(node.twin.getSeq());				
				}
		}
		return contigs;
	}
	
	public int getMaxContig() {
		int max = 0;
		for (String seq: nodes.keySet()) {
			if (seq.length() > max) {
				max = seq.length();
			}
		}
		return max;
	}
	
	@Override
	public String toString() {
		String s = "";
		for (Map.Entry<String, Node> entry : nodes.entrySet()) {
			s = s + entry.getValue().toString() + "\n";
		}
		return s;
	}
}
