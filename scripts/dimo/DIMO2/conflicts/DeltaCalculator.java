package conflicts;

import java.util.Vector;

import conflicts.sets.ObjectiveSet;

public class DeltaCalculator {

	private static void calculateDelta(String filename, int[] set1, int[] set2, String outputfilename) {
		FileProblem problem = new FileProblem(filename);
		FilePopulation pop = new FilePopulation(problem);
		DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(pop, 1);
		ObjectiveSet os1 = new ObjectiveSet(set1, problem.os_dim);
		ObjectiveSet os2 = new ObjectiveSet(set2, problem.os_dim);
	
		Vector<String> toWrite = new Vector<String>();
		toWrite.add("delta equals " + dmga.getDeltaMinFor(os1, os2));
		Output.print(toWrite, outputfilename);
	}
	
	
	/**
	 * computes the delta value for solutions in filename according to
	 * the sets set1 and set2 which are given by the arguments
	 * arg[1]-arg[s-1] and arg[s+1]-arg[arg.length-1]
	 * where s is a separator argument including the character 's'.
	 * 
	 * usage: "DeltaCalculator filename 1 3 5 8 s 2 4 "
	 *        or
	 *        "DeltaCalculator filename 1 3 5 8 s 2 4 -o outputfilename"
	 * 
	 * the additional last argument '-o outputfilename' will indicate
	 * that all output is written to the file 'outputfilename'
	 * 
	 */
	public static void main(String[] args) {
		if (args == null || args.length < 4) {
			System.out.println("Computes the delta error for all solutions in filename");
			System.out.println("  according to the sets set1 and set2 which are given");
			System.out.println("  by the arguments arg[1]-arg[s-1] and arg[s+1]-arg[arg.length-1]");
			System.out.println("  where s is a separator argument including the character s.");
			System.out.println();
			System.out.println("example usage of DeltaCalculator for calculating the delta");
			System.out.println("  error between the sets {f_1, f_3, f_5} and {f_3, f_4, f_7}:");
			System.out.println();
			System.out.println("  DeltaCalculator filename 1 3 5 s 3 4 7");	
			System.out.println();
			System.out.println();
			System.out.println("  Adding '-o outputfilename' as last argument will result in");
			System.out.println("    writing all output to the file 'outputfilename' instead of");
			System.out.println("    writing to standard output.");
		} else {
			String filename = args[0];
			int i=0;
			for (i=1;i<args.length; i++) {
				if (args[i].charAt(0) == 's') {
					break;
				}
			}
			int[] set1 = new int[i-1];
			for (int j=0; j<set1.length; j++) {
				set1[j] = new Integer(args[j+1]).intValue();
			}
			int[] set2 = new int[args.length-i-1];
			
			String outputfilename = "";
			if (args[args.length - 2].compareTo("-o") == 0) {
				outputfilename = args[args.length - 1];
				for (int j=0; j<set2.length-2; j++) {
					set2[j] = new Integer(args[i+j+1]).intValue();
				}
			} else {
				for (int j=0; j<set2.length; j++) {
					set2[j] = new Integer(args[i+j+1]).intValue();
				}
			}

			calculateDelta(filename, set1, set2, outputfilename);
		}
	}

}
