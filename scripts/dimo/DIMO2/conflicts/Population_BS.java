
package conflicts;

import edu.cornell.lassp.houle.RngPack.RandomElement;

public class Population_BS extends Population {
	 	 
	 /**
	 * @param dim the dimension of the decision space
	 * @param problem the problem which contains the fitness functions
	 */
	public Population_BS(int dim, Problem_BS problem) {
		super();
		
		this.ds_dim = dim;
		this.mu = (new Double(Math.pow(2,ds_dim))).intValue();
		this.problem = problem;

		/* Constructing all possible decision vectors */
		boolean[] vector = new boolean[dim];
		Individual_BS individuum;
		for (int i=0;i<mu;i++) {
			individuum = new Individual_BS(i, ds_dim, problem.getObjectiveSpaceDimension(), vector, problem.getFitness(vector));
			ind.add(i, individuum);
			vector = increase(vector);
		}
	}
	
	  
	 /**
	  * manipulates the boolean vector a. 
	  * Pre: a is the bit representation of the number n
	  * Post: a ist the bit representation of number n+1 
	  * 
	  * @param a a boolean vector
	  */
	 private boolean[] increase(boolean[] a) {
		 boolean[] b = new boolean[a.length];
		 for (int i=0;i<a.length;i++) {
			 b[i] = a[i];			 
		 }
		 int i = b.length-1;
		 while (i>=0) {
			 if (b[i]) {
				 b[i] = false;
			 } else {
				 b[i] = true;
				 break;
			 }				 
			 i--;
		 }
		 return b;		 
	 }
	
	/**
	 * @param mu number of individuals in the population
	 * @param dim dimension of the decision space
	 * @param problem the problem which contains the fitness functions 
	 */
	public Population_BS(int mu, int dim, Problem_BS problem) {
		super();
		
		this.ds_dim = dim;
		this.mu = mu;
		this.problem = problem;
		
		Individual_BS individuum;
		
		for (int i=0;i<mu;i++) {
			boolean[] vector = new boolean[dim];
			for (int j=0;j<ds_dim;j++) {
				// construct vector[j] randomly as true or false			
				RandomElement re = getRandomElement(); 
				double random = re.uniform(0,1);				
				if (random > 0.5) {
					vector[j] = true;
				} else {
					vector[j] = false;
				}
			}
			individuum = new Individual_BS(i, ds_dim, problem.getObjectiveSpaceDimension(), vector, problem.getFitness(vector));
			ind.add(i, individuum);
		}
	}
	
	public void print() {
		for (Individual individual : ind) {
			System.out.println(((Individual_BS)individual).toString());
		}		
	}
	
}
