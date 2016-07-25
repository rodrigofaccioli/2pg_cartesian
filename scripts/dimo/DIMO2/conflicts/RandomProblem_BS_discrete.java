
package conflicts;

import edu.cornell.lassp.houle.RngPack.RandomElement;

public class RandomProblem_BS_discrete extends Problem_BS {

	RandomElement re;

	/**
	 * @param k>0 the dimension of the objective space
	 */
	public RandomProblem_BS_discrete(int k, RandomElement re) {
		// init RandomElement:
		this.re = re;
		if (k < 1) {
			this.os_dim=1;
		} else {
			this.os_dim=k;
		}
	}
	
	/** 
	 * Returns a randomly chosen objective vector of type double the
	 * entries of which are uniformly distributed in {0,1}.
	 */
	public double[] getFitness(boolean[] decisionVector) {
		double[] objectiveVector = new double[os_dim];		 
		for (int i=0;i<os_dim;i++) {			
			objectiveVector[i] = re.choose(0,1);        // uniform number out of {0,1}
		}
		return objectiveVector;
	}

}
