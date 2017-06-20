import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.NoSuchElementException;


/**
 * Simple implementation of a de Brujin graph for 
 * sequence assembly.
 * <p>
 * @author tohei
 * 
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
	 * @author tohei
	 *
	 */
	private class Block {
		private LinkedList<Block> outgoing;
		private Block twin;
		private String seq;
		private int cov;
		private int mult;
		
		public Block(String seq, int cov, int mult) {
			outgoing = new LinkedList<Block>();
			this.twin = null;
			this.seq = seq;
			this.cov = cov;
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
			return seq.substring(k-1);
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
			System.out.print("Building de Bruijn graph from " + k 
					+ "-mers of " + inputs.length + " reads ...");			
		}
		LinkedHashSet<String> kmers = SeqUtils.allKmers(inputs, k);
		blocks = new LinkedHashMap<String, Block>(kmers.size());
		this.k = k;

		if (verbose) {
			System.out.print(" adding blocks ...");			
		}
		// add blocks 
		for (String kmer : kmers) {
			String kmerRC = SeqUtils.reverseComplement(kmer);
			Block n1; Block n2;
			if (!blocks.containsKey(kmer)) {				
				n1 = new Block(kmer);
				n2 = new Block(kmerRC);
				n1.twin = n2;
				n2.twin = n1;
				
				blocks.put(kmer, n1);
				blocks.put(kmerRC, n2);
			}
		}		
		if (verbose) {
			System.out.print(" adding edges ...");			
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
			System.out.print(" done (" + (endTime - startTime)/1000 + 
					" seconds).\nTotal of " + getSize() + " blocks.\n");			
		}
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
			System.out.println("Collapsing graph. Number of starting blocks: "+startingBlocks.size());
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
			System.out.println("Collapsed de Bruijn graph in " + (endTime - startTime)/1000 + 
					" seconds. Graph size - before: "+blocksPriorCollapse+", after: "+getSize());
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
		for (Block block : linearStretch) {
			if (block == first) continue;
			shortSeq = shortSeq +  block.getShortSeq();
		}
		String shortSeqRC = SeqUtils.reverseComplement(shortSeq);
		String newSeq = first.getSeq() + shortSeq;
		String newSeqRC = shortSeqRC + first.getRCSeq();
		Block newBlock = new Block(newSeq);
		Block newTwin = new Block(newSeqRC);
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
	
	public String[] getContigs(boolean includeRC) {
		LinkedHashSet<String> visited = new LinkedHashSet<String>(blocks.size());
		String[] contigs = null;
		if (includeRC) {
			contigs = new String[blocks.size()];
		} else {
			contigs = new String[blocks.size()/2];
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
	
	public int getMaxContig() {
		int max = 0;
		for (String seq: blocks.keySet()) {
			if (seq.length() > max) {
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
