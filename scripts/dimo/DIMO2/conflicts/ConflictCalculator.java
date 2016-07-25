package conflicts;

import java.util.Vector;

import conflicts.sets.IntSet;
import conflicts.sets.ObjectiveSet;

public class ConflictCalculator {

	private FileProblem problem;
	private Population pop;
	private Controller con;

	public void start(String filename) {		
		init(filename);
		
		Vector<String> lines = new Vector<String>(); // to store output
		
		this.con = new Controller(this.problem, 1, this.pop, 1, 0);
		java.util.Calendar calendar = new java.util.GregorianCalendar();
		long milliseconds = calendar.getTimeInMillis();
		IntSet output = con.exactAlgorithm();
		lines.add("One possible minimum objective set: " + output.toString());
		int size = output.size();
		lines.add(size + " objectives needed for " + filename + " and weak dominance");
		java.util.Calendar calendar2 = new java.util.GregorianCalendar();
		long millis = calendar2.getTimeInMillis();
		lines.add("Elapsed time during computation: " + (millis-milliseconds) + " milliseconds");
		
		/* finally output everything to stdout */
		Output.print(lines, "");
	}

	public void start(String filename, int type, double value, boolean greedy, String outputfilename) {		
		init(filename);
		// Controller:
		this.con = new Controller(this.problem, 1, this.pop, 1, 0);
		
		Vector<String> lines = new Vector<String>(); // to store output
		
		java.util.Calendar calendar = new java.util.GregorianCalendar();
		long milliseconds = calendar.getTimeInMillis();
		if (type == 1) {
			if (greedy) {
				ObjectiveSet output = con.greedyAlgorithmForGivenDelta(value);
				lines.add("delta-MOSS problem with given delta= " + value);
				lines.add("greedy algorithm");
				lines.add("------------------------------------------------------");
				lines.add("One possible minimal objective set: " + output.toString());
				lines.add(output.size() + " objectives needed for " + filename + " with failure " + output.getDelta());
			} else {
				ObjectiveSet output = con.exactAlgorithmForGivenDelta(value);
				lines.add("delta-MOSS problem with given delta= " + value);
				lines.add("exact algorithm");
				lines.add("------------------------------------------------------");
				lines.add("One possible minimal objective set: " + output.toString());
				lines.add(output.size() + " objectives needed for " + filename + " with failure " + output.getDelta());
			}
		} else {
			int intvalue = (new Double(value)).intValue();
			if (greedy) {
				ObjectiveSet output = con.greedyAlgorithmForGivenK(intvalue);
				lines.add("delta-MOSS problem with given k= " + intvalue);
				lines.add("greedy algorithm");
				lines.add("------------------------------------------------------");
				lines.add("One possible objective set with minimal delta: " + output.toString());
				lines.add(output.size() + " objectives needed for " + filename + " with failure " + output.getDelta());
			} else {				
				ObjectiveSet output = con.exactAlgorithmForGivenK(intvalue);
				lines.add("delta-MOSS problem with given k= " + intvalue);
				lines.add("exact algorithm");
				lines.add("------------------------------------------------------");
				lines.add("One possible objective set with minimal delta: " + output.toString());
				lines.add(output.size() + " objectives needed for " + filename + " with failure " + output.getDelta());
			}
		}
		java.util.Calendar calendar2 = new java.util.GregorianCalendar();
		long millis = calendar2.getTimeInMillis();
		lines.add("Elapsed time during computation: " + (millis-milliseconds) + " milliseconds");
		
		/* finally output everything */
		Output.print(lines, outputfilename);
	}

	
	private void init(String filename){
		this.problem = new FileProblem(filename);
		this.pop = new FilePopulation(problem);
	}


	/**
	 * Computes the size of the minimum non-redundant set for a given
	 * set of individuals in a file named data.
	 *
	 * @param args
	 * 			args[0]: name of file, with information about the individuals
	 * 						data format: "id objectivevalue1 objectivevalue2 ..."
	 * 			args[1]/
	 * 			args[2]: type/value for the two different delta-MOSS problems:
	 * 						a) the case of given delta value (type = 1, value = delta) and
	 *						b) the case of given k (type = 2, value = k).
	 * 			if only one argument is given, the exact algorithm for MOSS is performed
	 * 			args[3]: if == 'g' then the greedy algorithm is performed
	 * 					 if != 'g' or omitted, the exact algorithm is used
	 *          moreover, an additional optional argument '-o outputfilename' will indicate
	 *          that all output is written to the file 'outputfilename'
	 * 
	 */
	public static void main(String[] args) {
		ConflictCalculator cc = new ConflictCalculator();
		
		if (args == null || args.length == 2 || args.length > 7) {
			System.out.println("Computes the size of a minimum non-redundant set for a given");
			System.out.println("  set of individuals in a file named data.");
			System.out.println();
			System.out.println("usage of ConflictCalculator: ");
			System.out.println("   ConflictCalculator filename type value [-o outputfilename]");
			System.out.println("   or");
			System.out.println("   ConflictCalculator filename type value g [-o outputfilename]");
			System.out.println();
			System.out.println("   where type = 1:  delta-MOSS");
			System.out.println("   and   type = 2:  k-EMOSS");
			System.out.println("   and value=delta, resp. value=k ");
			System.out.println("   If the single character g is used as optional forth");
			System.out.println("     argument, the greedy heuristic is used.");
			System.out.println("   If only one argument is given (filename), the exact");
			System.out.println("     algorithm for MOSS is performed");
			System.out.println("   Adding '-o outputfilename' as last argument will result in");
			System.out.println("     writing all output to the file 'outputfilename' instead of");
			System.out.println("     writing to standard output.");
		} else {
			String filename = args[0];
			String outputfilename = ""; // standard: output written to stdout
			if (args.length == 1) {
				cc.start(filename);
			} else if (args.length == 3) {
				int type = new Integer(args[1]).intValue();
				double value = new Double(args[2]).doubleValue();
				cc.start(filename, type, value, false, outputfilename);
			} else if (args.length == 4) {
				int type = new Integer(args[1]).intValue();
				double value = new Double(args[2]).doubleValue();
				char c = args[3].charAt(0); 
				if (c == 'g') {
					cc.start(filename, type, value, true, outputfilename);
				} else {
					cc.start(filename, type, value, false, outputfilename);
				}
			} else if (args.length == 5) {
				if (args[3].compareTo("-o") == 0) {
					int type = new Integer(args[1]).intValue();
					double value = new Double(args[2]).doubleValue();
					outputfilename = args[4];
					cc.start(filename, type, value, false, outputfilename);
				} else {
					System.out.println("Error: forth argument is not '-o'.");
				}
			} else if (args.length ==6) {
				if (args[4].compareTo("-o") == 0) {
					int type = new Integer(args[1]).intValue();
					double value = new Double(args[2]).doubleValue();
					char c = args[3].charAt(0); 
					outputfilename = args[5];
					if (c == 'g') {
						cc.start(filename, type, value, true, outputfilename);
					} else {
						cc.start(filename, type, value, false, outputfilename);
					}
				} else {
					System.out.println("Error: fifth argument is not '-o'.");
				}
			}
		}
	}

}
