
package conflicts;

public class TestProblem_BS extends Problem_BS{
	static int a=0;
//	static int b=1;
	
	/**
	 * @param k>0 the dimension of the objective space
	 */
	public TestProblem_BS(int k) {		
		if (k < 1) {
			this.os_dim=1;
		} else {
			this.os_dim=k;
		}		
	}
	
	/** 
	 * Returns an objective vector of type double the
	 * entries as a test problem.
	 */
	public double[] getFitness(boolean[] decisionVector) {
		double[] objectiveVector = new double[os_dim];
		for (int i=0;i<os_dim;i++) {
			if (i%3 == 0 || i%3 == 1) {
				objectiveVector[i]=(i+1)*a;
				a++;
			} else {
				objectiveVector[i]=(1000-i*a);
				a--;
			}
		}
		return objectiveVector;
	}

}
