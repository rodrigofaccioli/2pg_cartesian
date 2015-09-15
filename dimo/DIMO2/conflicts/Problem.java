
package conflicts;

public abstract class Problem {
	protected int os_dim;

	public abstract double[] getFitness(Object genotype); 

	public int getObjectiveSpaceDimension() {
		return os_dim;
	}
}
