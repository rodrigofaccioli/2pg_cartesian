
package conflicts;

import java.util.LinkedList;

public class FilePopulation extends Population {	

	public FilePopulation(FileProblem prob) {
		this.problem = prob;
		
		this.ds_dim = 32;
		this.ind = new LinkedList<Individual>();
		this.mu = prob.getNumberOfDifferentPoints();
		
		double[][] p = prob.getPoints();
		/* Inserting all Individuals from prob into the LinkedList ind*/
		for (int i=0; i<p.length; i++) {
			int id = new Double(p[i][0]).intValue();
			double[] obj_vector = new double[p[i].length-1];
			for (int j=0; j<p[i].length-1; j++) {
				obj_vector[j] = p[i][j+1];
			}
			/* Computing the individual's decision vector as binary representation of its id: */
			boolean[] dec_vector = new boolean[32];
			int b = id;
			for (int c=0; c<32; c++) {
				dec_vector[c] = ((b % 2)==0);
				b %= 2;
			}
			
			ind.add((Individual)(new Individual_BS(id, 32, prob.os_dim, dec_vector, obj_vector)));
		}

		this.seed = null;
		
	}
	
}
