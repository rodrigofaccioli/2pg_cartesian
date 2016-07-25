package conflicts.sets;

/**
 * Aggregation represents an objective aggregation with length(objectiveWeights)
 * objectives that are composed of k objectives according to the weights
 * 0<=objectiveWeight<=1 for each objective.
 */

public class Aggregation {
	
	private double[][] objectiveWeights;
	private int numberOfObjectives;
	private double maxError;
	private double averageError;
	
	public Aggregation(int k) {
		numberOfObjectives = k;
		objectiveWeights = new double[numberOfObjectives][numberOfObjectives];
		for (int i=0; i<numberOfObjectives; i++) {
			for (int j=0; j<numberOfObjectives; j++) {
				if (i==j) {
					objectiveWeights[i][j] = 1;
				} else {
					objectiveWeights[i][j] = 0;
				}
			}
		}
		maxError = 0.0;
		averageError = 0.0;
	}
	
	/* initializes all weights for numberOfAggregatedObjectives aggregated objectives
	 * as 0.0 where numberOfOriginalObjectives original objectives are assumed */
	public Aggregation(int numberOfOriginalObjectives, int numberOfAggregatedObjectives) {
		numberOfObjectives = numberOfOriginalObjectives;
		objectiveWeights = new double[numberOfAggregatedObjectives][numberOfOriginalObjectives];
		for (int i=0; i<numberOfAggregatedObjectives; i++) {
			for (int j=0; j<numberOfOriginalObjectives; j++) {
				objectiveWeights[i][j] = 0.0;
			}
		}
		maxError = Double.MAX_VALUE;
		averageError = Double.MAX_VALUE;
	}
	
	public Aggregation(double[][] weights, double maxError, double avgError) {
		numberOfObjectives = weights[0].length;
		objectiveWeights = new double[weights.length][weights[0].length];
		for (int i=0; i<weights.length; i++) {
			for (int j=0; j<weights[0].length; j++) {
				objectiveWeights[i][j] = weights[i][j];
			}
		}
		this.maxError = maxError;
		this.averageError = avgError;
	}
	
	/**
	 * returns aggregation as double array (deep copy):
	 *    rows: aggregated objectives
	 *    columns: original objectives
	 *    entries: weights, such that sum over all weights per row equals 1
	 */
	public double[][] getAggregation() {
		double[][] ret = new double[objectiveWeights.length][numberOfObjectives];
		for (int i=0; i<objectiveWeights.length; i++) {
			for (int j=0; j<numberOfObjectives; j++) {
				ret[i][j] = objectiveWeights[i][j];
			}
		}
		return ret;
	}
	
	/**
	 * aggregates current objectives i and j, given by the rows i and j,
	 * to the new objective alpha*i + (1-alpha)*j.
	 * Note that not necessarily i=f_i and j=f_j as i and j can be
	 * already aggregated objectives.
	 */
	public void aggregate(int i, int j, double alpha) {
		if (i>=j || i<0 || j<0 || i>=objectiveWeights.length || j>=objectiveWeights.length) {
			System.err.println("Wrong access of objectives i and j in Aggregation.aggregate(...)");
		}
		if (alpha < 0 || alpha >1) {
			System.err.println("aggregation weight alpha not between 0 and 1 in Aggregation.aggregate(...)");
		}
		/* error handling done */
		
		/* copy not changing objectives into newObjectiveWeights and extract objective
		 * weights for current objectives i and j for later aggregation */
		double[][] newObjectiveWeights = new double[objectiveWeights.length-1][numberOfObjectives];
		double[] weightsOfObjectiveI = new double[numberOfObjectives];
		double[] weightsOfObjectiveJ = new double[numberOfObjectives];
		int numberOfNewObjectives = 0;
		for (int row=0; row<objectiveWeights.length; row++) {
			if (row == i) {
				for (int obj=0; obj<numberOfObjectives; obj++) {
					weightsOfObjectiveI[obj] = objectiveWeights[row][obj];
				}	
				//row++;
			} else if (row == j) {
				for (int obj=0; obj<numberOfObjectives; obj++) {
					weightsOfObjectiveJ[obj] = objectiveWeights[row][obj];
				}
				//row++;
			} else {
				for (int obj=0; obj<numberOfObjectives; obj++) {
					newObjectiveWeights[numberOfNewObjectives][obj] = objectiveWeights[row][obj];
				}
				numberOfNewObjectives++;
			}
		}
		
		/* aggregate objectives i and j and attach new objective to newObjectiveWeights */
		for (int obj=0; obj<this.numberOfObjectives; obj++){
			newObjectiveWeights[numberOfNewObjectives][obj] = alpha*weightsOfObjectiveI[obj] + (1-alpha)*weightsOfObjectiveJ[obj];
		}
		
		/* change this.objectiveWeights to new aggregation in newObjectiveWeights */
		this.objectiveWeights = newObjectiveWeights;
	}
	
	public void setMaxError(double err) {
		this.maxError = err;
	}
	
	public double getMaxError() {
		return this.maxError;
	}
	
	public void setAverageError(double err) {
		this.averageError = err;
	}
	
	public double getAverageError() {
		return this.averageError;
	}
	
	/** returns a String representation of the Aggregation 
	 *  with the values of the weights for all objectives */
	public String toString() {
		String str = "";
		
		for (int i=0; i<this.objectiveWeights.length; i++) {
			for (int k=0; k<this.objectiveWeights[i].length; k++) {
				str += " " + this.objectiveWeights[i][k]; 
			}
			str += "\n";
		}
		str += "(max error: " + this.maxError + " )\n";
		str += "(avg error: " + this.averageError + " )\n";

		return str;
	}
	
}
