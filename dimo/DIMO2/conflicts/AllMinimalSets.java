package conflicts;

import java.util.Iterator;
import java.util.Vector;
import conflicts.sets.ObjectiveSet;
import conflicts.sets.SetOfObjectiveSets;

public class AllMinimalSets {

	private FileProblem problem;
	private Population pop;
	private Controller con;

	public void start(String filename, int type, double value, String outputfilename) {		
		init(filename);
		// Controller:
		this.con = new Controller(this.problem, 1, this.pop, 1, 0);
		
		Vector<String> lines = new Vector<String>(); // to store output
		java.util.Calendar calendar = new java.util.GregorianCalendar();
		long milliseconds = calendar.getTimeInMillis();
		SetOfObjectiveSets output;
		if (type == 1) {
				output = con.allMinimalSetsDelta(value);
				lines.add("delta-MOSS problem with given delta= " + value);
				lines.add("exact algorithm");
				lines.add("------------------------------------------------------");
				
		} else {
			int intvalue = (new Double(value)).intValue();
				output = con.allMinimalSetsK(intvalue);
				lines.add("k-EMOSS problem with given k= " + intvalue);
				lines.add("exact algorithm");
				lines.add("------------------------------------------------------");
				
		}
		Iterator<ObjectiveSet> iter = output.getElements().iterator();
		int i=1;
		while (iter.hasNext()) {
			ObjectiveSet os = iter.next();
			lines.add(i + " " + os.size() + "  " + os.toString());
			i++;
		}
		lines.add("------------------------------------------------------");
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
	 * 						(data format: "id objectivevalue1 objectivevalue2 ...")
	 * 			args[1]/
	 * 			args[2]: type/value for the two different delta-MOSS problems:
	 * 						a) the case of given delta value (type = 1, value = delta) and
	 *						b) the case of given k (type = 2, value = k).
	 *          optional:
	 *          with '-o' the last argument gives a file name into which the result is
	 *          written instead of writing to standard output.
	 * 			
	 * 
	 */
	public static void main(String[] args) {
		AllMinimalSets ams = new AllMinimalSets();
		if (args == null || !(args.length == 3 || args.length == 5)) {
			System.out.println("computes all minimal objective sets for a Pareto set");
			System.out.println("  approximation given in file filename");
			System.out.println();
			System.out.println("usage:");
			System.out.println("   AllMinimalSets filename type value");
			System.out.println("   or");
			System.out.println("   AllMinimalSets filename type value -o outputfilename");
			System.out.println();
			System.out.println("   where type = 1:  delta-MOSS");
			System.out.println("   and   type = 2:  k-EMOSS");
			System.out.println("   and value=delta, resp. value=k ");
			System.out.println("   if '-o outputfilename' option is used, the result will be");
			System.out.println("   written to outputfilename and not to standard out");
			
		} else {
			String filename = args[0];
			int type = new Integer(args[1]).intValue();
			double value = new Double(args[2]).doubleValue();
			String outputfilename = "";
			if (args.length == 5) {
				if (args[3].compareTo("-o") == 0) {
					outputfilename = args[4];
				} else {
					System.out.println("Error: forth argument is not '-o'.");
				}
			}
			
			ams.start(filename, type, value, outputfilename);
		}
	}

}
