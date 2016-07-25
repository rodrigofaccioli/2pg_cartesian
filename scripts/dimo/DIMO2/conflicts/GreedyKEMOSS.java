
package conflicts;

import java.util.Vector;

import conflicts.sets.ObjectiveSet;

public class GreedyKEMOSS {

	private FileProblem problem;
	private Population pop;
	private Controller con;
	
	public void start(String filename, int k, String outputfilename) {		
		init(filename);
		Vector<String> lines = new Vector<String>(); // to store output
		
		/* new Controller(this.problem, searchSpace={0,1}^n, this.population,
		                   weak dominance relation, epsilon=0): */
		this.con = new Controller(this.problem, 1, this.pop, 1, 0);
		
		ObjectiveSet output = con.greedyAlgorithmForGivenK(k);
		boolean[] elements = output.getElements();
		for (int i=0; i<elements.length; i++) {
			if (elements[i]) {
				lines.add(Integer.toString(i));
			}
		}
		
		/* finally output everything */
		Output.print(lines, outputfilename);
	}

	
	private void init(String filename){
		this.problem = new FileProblem(filename);
		this.pop = new FilePopulation(problem);
	}

	/**
	 * Performs the greedy algorithm for kEMOSS and writes the objectives contained
	 *    in the computed set line-by-line either to a specified output file or to stdout.
	 *
	 * @param args
	 * 			args[0]: name of file, with information about the individuals
	 * 						data format: "id objectivevalue1 objectivevalue2 ..."
	 * 			args[1]: given k 
	 *          an additional optional argument '-o outputfilename' will indicate
	 *             that all output is written to the file 'outputfilename'
	 * 
	 */
	public static void main(String[] args) {
		if (args == null || args.length <2 || args.length > 4) {
			System.out.println("Wrong usage.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   GreedyKEMOSS filename k [-o outputfilename]");
			System.out.println();
			System.out.println("Adding '-o outputfilename' as last argument will result");
			System.out.println("   in writing all output to the file 'outputfilename'");
			System.out.println("   instead of writing to standard output.");
		} else {
			GreedyKEMOSS gkemoss = new GreedyKEMOSS();
			String outputfilename = ""; // standard: output written to stdout
			if (args.length == 4) {
				if (args[args.length - 2].compareTo("-o") == 0) {
					outputfilename = args[args.length - 1];
				} else {
					System.out.println("Error: second last argument is not '-o'.");
					System.exit(-1);
				}
			} else if (args.length == 3) {
				System.out.println("Error: 3 arguments not possible");
				System.exit(-1);
			}
			gkemoss.start(args[0], (new Integer(args[1])).intValue(), outputfilename);
		} 
	}

}
