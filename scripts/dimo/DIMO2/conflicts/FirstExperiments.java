
package conflicts;

import java.util.Vector;

import cern.jet.random.engine.MersenneTwister;
import edu.cornell.lassp.houle.RngPack.RandomElement;

public class FirstExperiments {

	RandomElement re;
	
	void firstExperiments(int n, int k_min, int k_max, int k_step, int iterations, int problem, int algo, double relation, String outputfilename) {
		RandomElement newseed = new MersenneTwister((int)System.currentTimeMillis());
		this.firstExperiments(n, k_min, k_max, k_step, iterations, problem, algo, relation, newseed.choose(Integer.MIN_VALUE, Integer.MAX_VALUE), outputfilename);	
	}
	
	void firstExperiments(int n, int k_min, int k_max, int k_step, int iterations, int problem, int algo, double relation, int seed, String outputfilename) {
		double average = 0;
		this.re = new MersenneTwister(seed);
		int dominanceRelation = 1;
		
		Vector<String> lines = new Vector<String>(); // to store output
		
		switch (problem) {
			case 0: lines.add("RandomProblem with objective values in [0,1]");
			        break;
			case 1: lines.add("RandomProblem with objective values in {0,1}");
			        break;
		}
		switch (algo) {
			case 0: lines.add("exact algorithm");
					break;
			case 1: lines.add("greedy algorithm");
		        	break;
		}
		if (relation == -1) {
			lines.add("Dominance relation: weak dominance");
			dominanceRelation = Controller.WEAK_DOMINANCE;
		} else {
			if (relation >= 0) {
				lines.add("Dominance relation: " + relation + "-dominance");
				dominanceRelation = Controller.EPSILON_DOMINANCE;
			} else {
				lines.add("ERROR OCCURED: epsilon < 0");
			}
		}
		if (n<0) {
			lines.add("Number of Individuals is always " + Math.pow(2,-n) + " and the population equals the whole decision space");
		} else {
			lines.add("Number of Individuals is always " + n);
		}
		lines.add("k_min: " + k_min);
		lines.add("k_max: " + k_max);
		lines.add("k_step: " + k_step);
		lines.add("Random seed for the MersenneTwister: " + seed);
		lines.add("Each experiment is iterated " + iterations + " times with the same parameters but different random seeds");
		lines.add("______________________________________________________");
		lines.add("------------------------------------------------------");
		for (int k=k_min; k< k_max; k+= k_step) {
			Problem_BS prob;
			Population pop;
			Controller con;
			lines.add("number of objectives: " + k);
			average = 0;
				java.util.Calendar calendar = new java.util.GregorianCalendar();
				long milliseconds = calendar.getTimeInMillis();
			for (int i=0; i<iterations; i++) {
				switch (problem) {
					case 0: prob = new RandomProblem_BS(k, re);
							break;
					case 1: prob = new RandomProblem_BS_discrete(k, re);
							break;
					default: prob = new RandomProblem_BS(k, re);
							pop = new Population_BS(n, prob);
							lines.add("Default...");
				}
				if (n<0) {
					pop = new Population_BS(-n, prob);
				} else {
					pop = new Population_BS(n, (((int)(Math.log(n)/Math.log(2)))+1), prob);
				}							
				con = new Controller(prob, 1, pop, dominanceRelation, relation);

				if (algo == 0) {					
					average = average + con.exactAlgorithm().size();
				} else if (algo == 1) {
					average = average + con.greedyAlgorithm().size();
				}

			}
			
			average = average / iterations;
			lines.add("average number of objectives needed: " + average);
			java.util.Calendar calendar2 = new java.util.GregorianCalendar();
			long millis = calendar2.getTimeInMillis();
			lines.add("Elapsed time during computation: " + (millis-milliseconds) + " milliseconds");
			lines.add("--------------------------------------------------------");
		}
		
		/* finally output everything */
		Output.print(lines, outputfilename);
	}
	
	/**
	 * Performs the algorithms of bz2007d for the results on the random problem.
	 * 
	 * usage `FirstExperiments n k_min k_max k_step iterations problem algo relation seed [-o outputfilename]`
	 * or    `FirstExperiments n k_min k_max k_step iterations problem algo relation [-o outputfilename]`
	 * 
	 * where
	 * 
	 *  - n: decision space dimension
	 *    if n<0 then the population is the whole search space of dimension |n|
	 *    if n>0 then the population size is n
	 *  - k: objective space dimension from k_min up to k_max in steps of k_step
	 *  - each parameter set is iterated for iterations times
	 *  - problem=0: RandomProblem with values in [0,1]
	 *    problem=1: RandomProblem with values in {0,1}
	 *  - algo=0: exact algorithm
	 *    algo=1: greedy algorithm
	 *  - relation == -1: weak dominance relation
	 *    relation == epsilon >= 0: epsilon-dominance 
	 *  - seed: random seed for the Mersenne Twister
	 *  - the optional last argument '-o outputfilename' results in  writing all
	 *    output to the file 'outputfilename' instead of writing to standard output.
	 * 
	 */ 	
	public static void main(String[] args) {
		String outputfilename = ""; // standard: output written to stdout
		if (args == null || args.length < 8 || args.length > 11) {
			System.out.println("Wrong syntax.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   FirstExperiments n k_min k_max k_step iterations problem algo relation seed [-o outputfilename]");
			System.out.println("   or");
			System.out.println("   FirstExperiments n k_min k_max k_step iterations problem algo relation [-o outputfilename]");
			System.out.println();
			System.out.println("   if n<0 then the population is the whole search space of dimension |n|");
			System.out.println("   if n>0 then the population size is n");
			System.out.println("   k: objective space dimension from k_min up to k_max in steps of k_step");
			System.out.println("   each parameter set is iterated for iterations times");
			System.out.println("   problem==0: RandomProblem with objective values in [0,1]");
			System.out.println("   problem==1: RandomProblem with objective values in {0,1}");
			System.out.println("   algo==0: exact algorithm");
			System.out.println("   algo==1: greedy algorithm");	
			System.out.println("   relation is a double value, iff =-1 the dominance relation is the weak dominance relation");
			System.out.println("      if relation = epsilon > 0 then the relation is the epsilon-dominance relation");
			System.out.println("   seed: random seed for the Mersenne Twister");
			System.out.println("   Adding '-o outputfilename' as last argument will result in writing all output to the file");
			System.out.println("   'outputfilename' instead of writing to standard output.");
			
		} else if (args.length == 11){
			if (args[args.length - 2].compareTo("-o") == 0) {
				outputfilename = args[args.length - 1];
				FirstExperiments fi = new FirstExperiments();
				fi.firstExperiments((new Integer(args[0]).intValue()), (new Integer(args[1]).intValue()),(new Integer(args[2]).intValue()), (new Integer(args[3]).intValue()),(new Integer(args[4]).intValue()), (new Integer(args[5]).intValue()), (new Integer(args[6]).intValue()), (new Double(args[7]).doubleValue()), (new Integer(args[8]).intValue()), outputfilename);
			} else {
				System.out.println("Error: second last argument is not '-o'.");
			}
		} else if (args.length == 10) {
			if (args[args.length - 2].compareTo("-o") == 0) {
				outputfilename = args[args.length - 1];
				FirstExperiments fi = new FirstExperiments();
				fi.firstExperiments((new Integer(args[0]).intValue()), (new Integer(args[1]).intValue()),(new Integer(args[2]).intValue()), (new Integer(args[3]).intValue()),(new Integer(args[4]).intValue()), (new Integer(args[5]).intValue()), (new Integer(args[6]).intValue()), (new Double(args[7]).doubleValue()), outputfilename);
			} else {
				System.out.println("Error: second last argument is not '-o'.");
			}
		} else if (args.length == 9) {
			FirstExperiments fi = new FirstExperiments();
			fi.firstExperiments((new Integer(args[0]).intValue()), (new Integer(args[1]).intValue()),(new Integer(args[2]).intValue()), (new Integer(args[3]).intValue()),(new Integer(args[4]).intValue()), (new Integer(args[5]).intValue()), (new Integer(args[6]).intValue()), (new Double(args[7]).doubleValue()), (new Integer(args[8]).intValue()), outputfilename);
		} else if (args.length == 8) {
			FirstExperiments fi = new FirstExperiments();
			fi.firstExperiments((new Integer(args[0]).intValue()), (new Integer(args[1]).intValue()),(new Integer(args[2]).intValue()), (new Integer(args[3]).intValue()),(new Integer(args[4]).intValue()), (new Integer(args[5]).intValue()), (new Integer(args[6]).intValue()), (new Double(args[7]).doubleValue()), outputfilename);
		}
	}

}
