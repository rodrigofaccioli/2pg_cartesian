
package conflicts;

public class Individual {
	
	protected int ds_dim = 1;
	protected int os_dim = 1;
	protected double[] ov;    // the individual's objective vector
	protected int id;
	
	public int getID() {
		return id;
	}

	public double[] getObjectiveVector() {
		return ov;
	}

	public void setObjectiveVector(double[] ov) {
		this.ov = ov;
	}
	

}
