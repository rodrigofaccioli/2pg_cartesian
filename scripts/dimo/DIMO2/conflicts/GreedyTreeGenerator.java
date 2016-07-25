
package conflicts;

public class GreedyTreeGenerator {

	private static void generateTree(String filename, int method, String outputfilename) {
		FileProblem problem = new FileProblem(filename);
		FilePopulation pop = new FilePopulation(problem);
		DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(pop, 1);
		dmga.computeTree(method, outputfilename);
	}
	
	
	/**
	 * computes a tree-like visualization of the hierarchical clustering based
	 * greedy algorithm for both delta-MOSS and kEMOSS. 
	 * 
	 * usage: "GreedyTreeGenerator filename [-o outputfilename]"
	 * 
	 * the optional argument '-o outputfilename' indicates that all output is
	 *   written to the file 'outputfilename'
	 */
	public static void main(String[] args) {
		if (args == null || !(args.length == 1 || args.length == 3)) {
			System.out.println("Wrong usage.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   GreedyTreeGenerator filename [-o outputfilename]");
			System.out.println();
			System.out.println("Adding '-o outputfilename' as last argument will result");
			System.out.println("   in writing all output to the file 'outputfilename'");
			System.out.println("   instead of writing to standard output.");
		} else {
			String filename = args[0];
			
			String outputfilename = ""; // standard: output written to stdout
			if (args.length == 3) {
				if (args[args.length - 2].compareTo("-o") == 0) {
					outputfilename = args[args.length - 1];
				} else {
					System.out.println("Error: second last argument is not '-o'.");
					System.exit(-1);
				}
			}
			generateTree(filename, 0, outputfilename);
		}
	}


}
