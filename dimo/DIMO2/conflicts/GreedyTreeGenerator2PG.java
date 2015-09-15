
package conflicts;

public class GreedyTreeGenerator2PG {

	private static void generateTree2PG(String filename, int method, String outputfilename, String  outputfilename2) {
		FileProblem problem = new FileProblem(filename);
		FilePopulation pop = new FilePopulation(problem);
		DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(pop, 1);
		dmga.computeTree2PG(method, outputfilename, outputfilename2);
	}
	
	
	/**
	 * computes a tree-like visualization of the hierarchical clustering based
	 * greedy algorithm for both delta-MOSS and kEMOSS. 
	 * 
	 * usage: "GreedyTreeGenerator2PG filename -o outputfilename -o2 outputfilenamefor2PG"
	 * 
	 * the optional argument '-o outputfilename' indicates that all output is
	 *   written to the file 'outputfilename'
         *                       '-o2 outputfilenamefor2PG' indicates that output for 2PG framework is
         *   written to the file 'outputfilenamefor2PG'
	 */
	public static void main(String[] args) {
		if (args == null || !(args.length == 1 || args.length == 5)) {
			System.out.println("Wrong usage.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   GreedyTreeGenerator2PG filename -o outputfilename -o2 outputfilename2PG");
			System.out.println();
			System.out.println("Adding '-o outputfilename' as second argument will result");
			System.out.println("   in writing all output to the file 'outputfilename'");
			System.out.println("   instead of writing to standard output.");
			System.out.println("Adding '-o2 outputfilenamefor2PG' as last argument will result");
			System.out.println("   in writing 2PG output to the file 'outputfilenamefor2PG' ");

		} else {
			String filename = args[0];			
			String outputfilename = ""; // standard: output written to stdout
                        String outputfilename2PG = ""; // standard: output written to 2PG framework

			if (args.length == 5) {
				if (args[1].compareTo("-o") == 0) {
					outputfilename = args[2];
				} else {
					System.out.println("Error: second argument is not '-o'.");
					System.exit(-1);
                                }
				if (args[3].compareTo("-o2") == 0) {
					outputfilename2PG = args[4];
				} else {
					System.out.println("Error: last argument is not '-o2'.");
					System.exit(-1);
				}

			}
			generateTree2PG(filename, 0, outputfilename, outputfilename2PG);
		}
	}


}
