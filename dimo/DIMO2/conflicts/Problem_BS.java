
package conflicts;

public abstract class Problem_BS extends Problem {

	public abstract double[] getFitness(boolean[] decisionvector);
	
	public double[] getFitness(Object genotype) {
		return getFitness((boolean[])genotype);		
	}

}
