
package conflicts;

import java.util.Vector;

public class FrontCalculator {
	
	double[][] points;  // objective values of points ([individual id][objective number]) 
	int m=0;            // number of individuals
	int k=0;            // number of objectives + 1
	
	public FrontCalculator(double[][] points) {
		this.m = points.length;
		this.k = points[0].length;
		this.points = new double[m][k];
		for (int i=0; i<points.length; i++) {
			for (int j=0; j<points[0].length; j++) {
				this.points[i][j] = points[i][j];
			}
		}
	}
	
	/**
	 * Prints all dominated (dominatedPoints==true) or non-dominated
	 * (dominatedPoints=false) points with objective values in points[][]
	 * with respect to the weak dominance relation and minimization. 
	 * 
	 * pre: points[][] should not contain indifferent solution pairs!
	 */
	public void printDominatedNonDominated(boolean dominatedPoints, String outputfilename){
		Vector<String> lines = new Vector<String>(); // to store output
		
		for (int i=0; i<this.m; i++) {		 
			boolean nondominated = true;
			for (int j=0; j<this.m; j++) {
				if (i != j) {
					boolean dominatedByJ = true;
					/* starting from column 1 because column 0 contains the
					 * id's of the individuals: */
					for (int k=1; k<this.k; k++) {
						if (this.points[i][k] < this.points[j][k]) {
							dominatedByJ = false;
						}
					}
					if (dominatedByJ) {
						nondominated = false;
					} 
				}
			}
			if (dominatedPoints) {
				if (!nondominated) {
					lines.add(print(i));
				}
			}
			else {
				if (nondominated) {
					lines.add(print(i));
				}
			}
				
		}
		
		/* finally output everything */
		Output.print(lines, outputfilename);
	}
	
	private String print(int i) {
		String str = new Double(points[i][0]).toString();				
		String cs = ".0";
		String cs2 = "";
		str = str.replace(cs, cs2);
		String returnString = str + " ";
		for (int s=1; s<this.k; s++) {
			returnString = returnString + this.points[i][s] + " ";
		}
		
		return returnString;
	}

	/**
	 * @param args args[0]: filename
	 * 			   args[1]: if == 0 or omitted, the nondominated points are printed
	 * 					    if == 1, the dominated points are printed
	 * 	           an additional optional argument '-o outputfilename' will indicate
	 *             that all output is written to the file 'outputfilename'
	 */
	public static void main(String[] args) {
		String outputfilename = ""; // standard: output written to stdout
		if (args == null || args.length < 1 || args.length > 4) {
			System.out.println("Wrong usage.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   FrontCalculator filename dominated");
			System.out.println("   or");
			System.out.println("   FrontCalculator filename dominated -o outputfilename");
			System.out.println();
			System.out.println("where");
			System.out.println("   the nondominated points in filename are printed if");
			System.out.println("      dominated = 0");
			System.out.println("   and the dominated points in filename are printed if");
			System.out.println("      dominated = 1");
			System.out.println();
			System.out.println("Adding '-o outputfilename' as last argument will result");
			System.out.println("   in writing all output to the file 'outputfilename'");
			System.out.println("   instead of writing to standard output.");
		} else if (args.length == 4) {
			if (args[args.length - 2].compareTo("-o") == 0) {
				outputfilename = args[args.length - 1];
				String filename = args[0];
				FileProblem problem = new FileProblem(filename);
				FrontCalculator fc = new FrontCalculator(problem.getPoints());
				int dominated = new Integer(args[1]).intValue();
				if (dominated == 0) {
					fc.printDominatedNonDominated(false, outputfilename);
				} else {
					fc.printDominatedNonDominated(true, outputfilename);
				}
			} else {
				System.out.println("Error: second last argument is not '-o'.");
			}
		} else if (args.length == 3) {
			if (args[args.length - 2].compareTo("-o") == 0) {
				outputfilename = args[args.length - 1];
				String filename = args[0];
				FileProblem problem = new FileProblem(filename);
				FrontCalculator fc = new FrontCalculator(problem.getPoints());
				fc.printDominatedNonDominated(false, outputfilename);
			} else {
				System.out.println("Error: second last argument is not '-o'.");
			}
		} else if (args.length == 2) {
			String filename = args[0];
			FileProblem problem = new FileProblem(filename);
			FrontCalculator fc = new FrontCalculator(problem.getPoints());
			int dominated = new Integer(args[1]).intValue();
			if (dominated == 0) {
				fc.printDominatedNonDominated(false, outputfilename);
			} else {
				fc.printDominatedNonDominated(true, outputfilename);
			}
		} else if (args.length == 1){
			String filename = args[0];
			FileProblem problem = new FileProblem(filename);
			FrontCalculator fc = new FrontCalculator(problem.getPoints());
			fc.printDominatedNonDominated(false, outputfilename);
		}
	}

}
