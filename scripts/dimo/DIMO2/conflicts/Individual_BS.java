
package conflicts;

public class Individual_BS extends Individual {

	boolean[] dv;   // the individual's decision vector
		
	public Individual_BS(int id, int ds_dim, int os_dim, boolean[] dv, double[] ov) {
		this.id = id;
		this.ds_dim = ds_dim;
		this.os_dim = os_dim;
		this.dv = dv;
		this.ov = ov;
	}
	
	public Individual_BS(int id, int ds_dim, int os_dim) {
		this.id = id;
		this.ds_dim = ds_dim;
		this.os_dim = os_dim;
		this.dv = new boolean[ds_dim];
		this.ov = new double[os_dim];
	}

	public boolean[] getDecisionVector() {
		return dv;
	}

	public void setDecisionVector(boolean[] dv) {
		this.dv = dv;
	}
	public String toString() {
		String str = id + ": ";
		for (int i=0;i<ds_dim;i++) {
			if (dv[i]) {
				str = str + 1;
			} else {
				str = str + 0;
			}
		}
		str = str + " - ";
		for (int i=0;i<os_dim;i++) {
			str = str + ov[i] + ", ";
		}
		return str;
	}
}
