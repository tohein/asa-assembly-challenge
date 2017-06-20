import java.util.LinkedHashSet;


public class LaunchAssembly {
	
	public static void main(String args[]) {
		
		// (a)
		String[] inputs = {"AACGTACGTAGGTACGTCG"};
		DBGraph G = new DBGraph(inputs, 5);
		System.out.println("(a)");
		LinkedHashSet<String> kmers = SeqUtils.allKmers(inputs, 5);
		String[] inc = G.getIncoming("CGTAC");
		String[] out = G.getOutgoing("CGTAC");
		System.out.print("Incoming k-mers - ");
		boolean existingEdge = false;
		for (String s : inc) {
			if (kmers.contains(s)) {
				System.out.print(s + " ");
				existingEdge = true;
			}
		}
		if (!existingEdge) {
			System.out.print("None");
		}
		existingEdge = false;
		System.out.print("\nOutgoing k-mers - ");
		for (String s : out) {
			if (kmers.contains(s)) {
				System.out.print(s + " ");
				existingEdge = true;
			}
		}
		if (!existingEdge) {
			System.out.print("None");
		}
		System.out.println("\n");
		
		// (b)
		System.out.println("(b)");
		int k = 21;
		try {
			inputs = SeqUtils.readFasta("/home/tohei/Downloads/reads_simple.fasta");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		G = new DBGraph(inputs, k, true);
		G.collapse(true);
		System.out.println("Collapsed de Bruijn graph: ");
		System.out.println(G.toString());				
		
		// (c)
		System.out.println("(c)");
		try {
			inputs = SeqUtils.readFasta("/home/tohei/Downloads/reads_simple_error.fasta");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int cutoff = 10;		
		System.out.print("Computing spectral alignment (cutoff = " + cutoff + ") ... ");
		inputs = SeqUtils.spectralAlignment(inputs, cutoff, k);
		G = new DBGraph(inputs, k, true);
		G.collapse(true);
		System.out.println("Collapsed de Bruijn graph: ");
		System.out.println(G.toString());		
		
		// (d)
		System.out.println("(d)");
		long startTime = System.currentTimeMillis();
		try {			
			inputs = SeqUtils.readFasta("/home/tohei/Downloads/reads_complex.fasta");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long endTime = System.currentTimeMillis();
		System.out.println("Read file in " + (endTime - startTime)/1000 + " seconds.");

		//inputs = Arrays.copyOfRange(inputs, 0, inputs.length); 
		System.out.print("Computing spectral alignment (cutoff = " + cutoff + ") ... ");
		startTime = System.currentTimeMillis();
		inputs = SeqUtils.spectralAlignment(inputs, cutoff, k);
		endTime = System.currentTimeMillis();
		System.out.println("done (" + (endTime - startTime)/1000 + " seconds).");
		System.out.println("Number of reads: "+inputs.length);
		
		G = new DBGraph(inputs, k, true);	
		//System.out.println(G.toString());
		
		G.collapse(true);			
		boolean collapsed = G.collapse(true);			
		System.out.println(collapsed);		
		//System.out.println(G.toString());	
		
	
		System.out.println("Max contig length: " + G.getMaxContig());
/*		String[] contigs = G.getContigs(true);
		try {
			SeqUtils.writeFasta("/home/tohei/Downloads/contigs.fasta", contigs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String[] contigs2 = null;
		try {			
			//inputs = SeqUtils.readFasta("/home/tohei/Downloads/reads_simple_error.fasta");
			contigs2 = SeqUtils.readFasta("/home/tohei/Downloads/reads_complex.unitigs.fa");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		System.out.println("CONTIGS: " + contigs.length);
		System.out.println("CONTIGS2: " + contigs2.length);
		
		LinkedHashSet<String> cont2 = new LinkedHashSet<String>(contigs.length);
		for (String s : contigs2) {
			cont2.add(s);
		}
		boolean[] sim = new boolean[contigs.length/2];
		for (int i = 0; i < contigs.length-1; i=i+2) {
			String c1 = contigs[i];
			String c2 = contigs[i+1];
			boolean contained1 = cont2.contains(c1);
			boolean contained2 = cont2.contains(c2);
			if (!(contained1 || contained2)) { 
				System.out.println(contigs[i]);;			
			}
		}
		
		System.out.println("\n_______________________\n");
		
		for (int i = 0; i < contigs2.length; i++) {
			String c1 = contigs2[i];			
			if (!G.contains(c1)) System.out.println(c1);
		}*/
	}
}
