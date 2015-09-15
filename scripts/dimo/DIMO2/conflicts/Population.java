
package conflicts;

import java.util.LinkedList;
import cern.jet.random.engine.MersenneTwister;
import edu.cornell.lassp.houle.RngPack.RandomElement;

public class Population {
	int ds_dim = 2; // dimension of the decision space
	int mu = 10;    // number of individuals in the population

	LinkedList<Individual> ind = new LinkedList<Individual>(); // the individuals in the population

	Problem problem;
	RandomElement seed;


	private void initRandomGenerator() {
		seed = new MersenneTwister((int)System.currentTimeMillis()%1000);		 
	}
	
	public Population() {
		initRandomGenerator();
	}
	
	protected RandomElement getRandomElement() {
		return new MersenneTwister(seed.choose(Integer.MIN_VALUE, Integer.MAX_VALUE));
	}

	public int getDecisionSpaceDimension() {
		return ds_dim;
	}

	public LinkedList<Individual> getPopulation() {
		return ind;
	}

	public int getPopulationsize() {
		return mu;
	}

	public Problem getProblem() {
		return problem;
	}

	public void addIndividual(Individual ind) {
		this.ind.add(ind);
	}
	
}
