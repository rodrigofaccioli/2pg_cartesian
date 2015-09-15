
package conflicts;

import java.util.Vector;

/**
 * Computes the maximal spread of a given population:
 * 
 *    max_spread = max_i [ max_p (f_i(p)) - min_p (f_i(p)) ]
 * 
 * 
 * Attention: only works for population of at least 2 individuals
 * and at least 2 objectives
 * 
 * Usage:
 * 
 *    MaxSpreadCalculator filename [-o outputfilename]
 *    
 * the optional argument '-o outputfilename' indicates that all output is
 * written to the file 'outputfilename'
 * 
 */
public class MaxSpreadCalculator {

	public static void main(String[] args) {
		if (args == null || !(args.length == 1 || args.length == 3)) {
			System.out.println("Wrong usage.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   MaxSpreadCalculator filename [-o outputfilename]");
			System.out.println();
			System.out.println("Adding '-o outputfilename' as last argument will result");
			System.out.println("   in writing all output to the file 'outputfilename'");
			System.out.println("   instead of writing to standard output.");
		} else {
			String outputfilename = ""; // standard: output written to stdout
			if (args.length == 3) {
				if (args[args.length - 2].compareTo("-o") == 0) {
					outputfilename = args[args.length - 1];
				} else {
					System.out.println("Error: second last argument is not '-o'.");
					System.exit(-1);
				}
			}			
			String filename = args[0];
			FileProblem problem = new FileProblem(filename);
			FilePopulation pop = new FilePopulation(problem);
			DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(pop, 1);
			double[][] m = dmga.values;
			double[] max = new double[m[0].length];
			double[] min = new double[m[0].length];
			for (int i=0; i<max.length; i++) {
				max[i] = Double.NEGATIVE_INFINITY;
				min[i] = Double.POSITIVE_INFINITY;
			}
			for (int i=0; i<m.length; i++) {
				for (int j=0; j<m[0].length; j++) {
					max[j] = Math.max(max[j], m[i][j]);
					min[j] = Math.min(min[j], m[i][j]);
				}
			}
			double maximum = Double.NEGATIVE_INFINITY;
			double minimum = Double.POSITIVE_INFINITY;
			double[] diff = new double[max.length];
			for (int i=0; i<max.length; i++) {
				diff[i] = Math.abs(max[i] - min[i]);
				if (max[i] > maximum) {
					maximum = max[i];
				}
				if (min[i] < minimum) {
					minimum = min[i];
				}
			}
			double diffmax = Double.NEGATIVE_INFINITY;
			for (int i=0; i<diff.length; i++) {
				diffmax = Math.max(diffmax, diff[i]);
			}
			
			/* output results */
			Vector<String> toPrint = new Vector<String>(); // to store output
			toPrint.add("The minimal and maximal values in " + filename + " are");
			toPrint.add("Min= " + minimum);
			toPrint.add("Max= " + maximum);
			toPrint.add("and the maximal spread of the population is " + diffmax);
			Output.print(toPrint, outputfilename);
		}
	}

}
