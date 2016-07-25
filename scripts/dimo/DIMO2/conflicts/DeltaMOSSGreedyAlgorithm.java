package conflicts;

import java.util.LinkedList;
import java.util.Vector;
import java.util.Collections; 
import conflicts.sets.Aggregation;
import conflicts.sets.ObjectiveSet;
import conflicts.sets.TreeSet;

public class DeltaMOSSGreedyAlgorithm {
	
	// matrix of objective values (columns) for individuals (rows) 
	double[][] values = new double[1][1];
	int os_dim = values[0].length;
	
	static double failure = 0.0000000001; 
	
	Relation[] relations;		// an array of the relations $\preceq_i$
	Relation dominanceRelation; // the relation $\preceq$
	
	public DeltaMOSSGreedyAlgorithm(double[][] d) {
		this.values = d;
		this.os_dim = values[0].length;
		
		/* preparing for initialization of relations: */
		Population pop = new Population();
		for (int i=0; i<values.length; i++) {
			boolean[] dv = new boolean[values.length];
			double[] ov = values[i];
			pop.addIndividual(new Individual_BS(i, values.length, this.os_dim, dv, ov));
		}
		this.initRelations(pop);
	}
	
	/* Constructor: initializes the weak dominance relations iff rel=1 */
	public DeltaMOSSGreedyAlgorithm(Population pop, int rel) {
		LinkedList<Individual> individuals = pop.getPopulation();
		int size = individuals.size();
		this.os_dim = (individuals.get(0)).os_dim;
		values = new double[size][os_dim];
		for(int i=0; i<size; i++) {
			Individual ind = individuals.get(i);
			values[i] = ind.getObjectiveVector();
		}
		if (rel==1) {
			this.initRelations(pop);
		}
	}

	private void initRelations(Population pop) {
		this.relations = new Relation[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			this.relations[i] = Controller.computeWeakRelation(i, pop);
		}
		this.dominanceRelation = Controller.computeWholeWeakRelation(pop);
	}
	
	public ObjectiveSet performGreedyAlgorithmGivenK(int k) {
		boolean[] elements = new boolean[this.os_dim];
		ObjectiveSet chosen = new ObjectiveSet(elements, Double.MAX_VALUE);
		int numberOfChosenObjectives = 0;
		
		while (numberOfChosenObjectives < k && chosen.getDelta() > 0) {
			chosen = computeNewObjectiveSet(chosen);
			numberOfChosenObjectives++;
		}
		return chosen;
	}
	
	/** PART OF GREEDY ALGORITHM FOR GIVEN K
	  * Returns the chosen objectives in the next iteration of the greedy algorithm.
	  * Given the set of already chosen objectives in 'chosen', the returned
	  * ObjectiveSet contains all formerly chosen objectives plus an additional
	  * objective that minimizes the delta-error for this new set of |chosen|+1
	  * objectives. The delta value of this new ObjectiveSet is computed, too. 
	  */	
	private ObjectiveSet computeNewObjectiveSet(ObjectiveSet chosen) {
		double delta_min = Double.MAX_VALUE;
		int bestObjective = 0;
		
		ObjectiveSet chosen_next = null;
		ObjectiveSet notchosen_next = null;
		boolean[] chosenObjectives = chosen.getElements();
		boolean[] notChosenObjectives = new boolean[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			if (!chosenObjectives[i]) {
				notChosenObjectives[i] = true;
			}
		}
		for (int i=0; i<this.os_dim; i++) {
			/* test for all non-chosen objectives, if delta-error gets minimal and
			 * store index of this objective in bestObjective: */
			if (!chosenObjectives[i]) {
				/* add i to chosen objectives and remove i from set of not chosen ones */
				chosenObjectives[i] = true;
				notChosenObjectives[i] = false;
				chosen_next = new ObjectiveSet(chosenObjectives);
				notchosen_next = new ObjectiveSet(notChosenObjectives);
				double delta = getDeltaMinFor(chosen_next, notchosen_next);
				if (delta < delta_min) {
					delta_min = delta;
					bestObjective = i;
				}
				/* remove i from chosen objectives and add i from notchosen again
				 * (get invariant working, that chosen and notchosen do not vary from i to i+1)*/
				chosenObjectives[i] = false;
				notChosenObjectives[i] = true;
			}
		}
		
		// prepare return value:
		chosenObjectives[bestObjective] = true;
		chosen_next = new ObjectiveSet(chosenObjectives, delta_min);
		
		return chosen_next;
	}
	
	/** PART OF GREEDY ALGORITHM FOR GIVEN K
	  * returns the minimal delta value such that the objective subset set1
	  * is delta-nonconflicting with set2.
	  */
	public double getDeltaMinFor(ObjectiveSet set1, ObjectiveSet set2) {
		double delta = 0;
		int numOfInds = values.length;
		for (int i=0; i<numOfInds; i++) {
			for (int j = 0; j < numOfInds; j++) {
				if (weaklyDominates(i, j, set1)) {
					/* for all objectives in set2 adjust the delta value according
					   to the current pair (i,j) of individuals: */
					boolean[] obj = set2.getElements();
					for (int k=0; k<obj.length; k++) {
						if (obj[k]) {
							delta = Math.max(delta, values[i][k]-values[j][k]);
						}
					}
				}
			}
		}
		return delta;
	}
	
	/** PART OF GREEDY ALGORITHM FOR GIVEN K
	  * returns true iff individual i weakly dominates individual j w.r.t. objective
	  * subset ef.
	  */
	private boolean weaklyDominates(int i, int j, ObjectiveSet ef) {
		boolean[] objectives = ef.getElements();
		for (int k=0; k<objectives.length; k++) {
			if (objectives[k]) {
				if (values[i][k] > values[j][k]) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * performs the greedy algorithm to reduce the number of objectives such that
	 * the error is at most delta allowing the aggregation of the initial objectives.
	 * The delta error is computed according to the parameter 'a': As the delta error,
	 * we use
	 * (i)  the maximum delta error (if a=1) and
     * (ii) the delta error averaged over all solution pairs (if a=2)
	 * 
	 */
	public conflicts.sets.Aggregation performGreedyAggregationAlgorithmGivenDelta(double delta, int a) {
		/* compute delta error for each solution pair (i,j) for the case that
		 * we wrongly assume that i is better than j (initialization) */
		double[][] deltaErrors = computeDeltaErrorsForAllSolutionPairs();
		
		Aggregation aggregation = new Aggregation(this.os_dim);
		double maxError = 0.0;
		double averageError = 0.0;
		int numOfErrorsConsidered = 0;
		int current_k = this.os_dim;
		
		/* try to reduce the objective set size as much as possible: */
		while (current_k > 0) {
			/* find best aggregation between objectives bestAggregation[0] and bestAggregation[1]: */
			double[] bestAggregation = getBestAggregation(this.values, deltaErrors, aggregation, a);
			
			/* stop if best possible aggregation yields a too large delta error (according to a)
			 * otherwise aggregate and continue */
			if (bestAggregation[a+2] > delta) {
				break;
			} else {
				/* compute best aggregation between objectives bestAggregation[0] and
				 * bestAggregation[1] with weight alpha=bestAggregation[2]. */
				aggregation.aggregate((int)bestAggregation[0], (int)bestAggregation[1], bestAggregation[2]);

				if (bestAggregation[3] > maxError) {
					maxError = bestAggregation[3];
				}
				averageError = bestAggregation[4];
				numOfErrorsConsidered++;

				aggregation.setMaxError(maxError);
				aggregation.setAverageError(averageError);
			}
			
			current_k = current_k - 1;
		}


		return aggregation;
	}
	
	/**
	 * performs the greedy algorithm to reduce the number of objectives to k
	 * allowing aggregation of the initial objectives where the aggregation is
	 * incrementally constructed, i.e., objectives are added to an empty
	 * aggregation with weights=0 until no improvement of the delta-error is
	 * possible.
	 * 
	 * We distinguish two different variants of the algorithm:
	 * Variant 1 adds each of the original objectives only once whereas in
	 * Variant 2 all original objectives can be used multiple times.
	 * 
	 * The delta error is computed according to the parameter 'a': As the delta
	 * error, we use
	 * (i)  the maximum delta error (if a=1) and
     * (ii) the delta error averaged over all solution pairs (if a=2)
	 * 
	 */
	public conflicts.sets.Aggregation performGreedyIncrementalAggregationAlgorithmGivenK(int k, int variant, int a){
		/* compute delta error for each solution pair (i,j) for the case that
		 * we wrongly assume that i is better than j (initialization) */
		double[][] deltaErrors = computeDeltaErrorsForAllSolutionPairs();
		
		/* aggregation contains the best aggregation of the last iteration, in the
		 * beginning, this is the empty aggregation */
		Aggregation aggregation = new Aggregation(this.os_dim, k);
		
		/* improvementPossible = false means there is no possible aggregation that
		 * decreases the delta-error (according to the choice of 'a'). */
		boolean improvementPossible = true;
		/* used(i) = true if objective f_i is already used in the aggregation (only
		 * needed in the case of variant 1) */
		boolean[] used = new boolean[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			used[i] = false;
		}
		boolean allUsed = false; // only true if all original objectives are already used
		/* bestAggregation contains best aggregation found so far in current iteration */
		Aggregation bestAggregation = new Aggregation(this.os_dim);
		bestAggregation.setAverageError(Double.MAX_VALUE);
		bestAggregation.setMaxError(Double.MAX_VALUE);
		
		
		while (!allUsed && improvementPossible) {
			improvementPossible = false;
			int usedObjective = -1; // stores number of objective that is used in best aggregation
			
			for (int i=0; i<this.os_dim; i++) {
				if (variant == 1 && used[i]) {
					continue; // goto next objective if already used in variant 1
				}

				for (int j=0;j<k;j++) {
					
//System.out.println("i= " + i + " j= " + j);
					
					/* compute optimal aggregation and delta error if original
					 * objective i is aggregated with already aggregated objective j. */
					Aggregation currentAggregation = computeOptimalAggregationAndDeltaError(this.values, deltaErrors, i, j, aggregation, a);

					
//System.out.println(currentAggregation.getMaxError());
					

					if ((a==1) && (currentAggregation.getMaxError() < bestAggregation.getMaxError())) {
						bestAggregation = currentAggregation;
						improvementPossible = true;
						usedObjective = i;
					} else if ((a==2) && (currentAggregation.getAverageError() < bestAggregation.getAverageError())) {
						bestAggregation = currentAggregation;
						improvementPossible = true;
						usedObjective = i;
					}
				}
			}
			
			if (improvementPossible) {
				aggregation = bestAggregation;
				used[usedObjective] = true;
				
//System.out.println("Improvement...");
				
			} else {
				/* in case one aggregated objective does not contain any objective,
				 * i.e., all weights are zero, we can have an improvement even
				 * if the error itself is not improving: */
				boolean someAggregatedObjectiveEmpty;
				double[][] weights = aggregation.getAggregation();
				for (int p=0; p<weights.length; p++) {
					someAggregatedObjectiveEmpty = true;
					for (int q=0; q<weights[p].length; q++) {
						if (weights[p][q] != 0.0) {
							someAggregatedObjectiveEmpty = false;
							break;
						}
					}
					if (someAggregatedObjectiveEmpty) {
						/* find first not used objective */
						int freeObjective = 0; 
						for (freeObjective=0; freeObjective<used.length; freeObjective++) {
							if (used[freeObjective] == false) {
								break;
							}
						}
						
						/* now, we can assign non-used objective to empty aggregated one: */
						weights[p][freeObjective] = 1.0;
						
						/* compute delta-errors for new aggregation */
						double[][] aggregatedObjValues = new double[this.values.length][weights.length];
						for (int i=0; i<this.values.length; i++) {
							for (int j=0; j<weights.length; j++) {
								aggregatedObjValues[i][j] = 0;
								for (int l=0; l<weights[0].length; l++) {
									aggregatedObjValues[i][j] += weights[j][l]*this.values[i][l];
								}
							}
						}					
						double[] errors = computeDeltaError(aggregatedObjValues, this.values);
						
						double maxError = errors[0];
						double avgError = errors[1];
						aggregation = new Aggregation(weights, maxError, avgError);
						improvementPossible = true;
						used[freeObjective] = true;
						
//System.out.println("--> freeObjective used: " + freeObjective);
//System.out.println("    maxError= " + maxError + " avgError= " + avgError);
					
						break;
					}
				}
			}
			
			/* update allUsed */
			if (variant == 1) {
				allUsed = true;
				for (int i=0; i<this.os_dim; i++) {
					if (!used[i]) {
						allUsed = false;
					}
				}
			}
			
			
//System.out.println("weights...");
//printDoubleMatrix(aggregation.getAggregation());
			
			
		}
		
		return aggregation;
	}
	
	
	/** compute optimal aggregation and delta error if the already aggregated objective j
	 * is aggregated with the original objective i if the objectives are already
	 * aggregated according to aggregation. */
	private Aggregation computeOptimalAggregationAndDeltaError(double[][] objectiveValues, double[][] deltaErrors, int i, int j, Aggregation aggregation, int a) {
		/* get weights of already aggregated objectives */
		double[][] weights = aggregation.getAggregation();
		
		/* compute aggregated objectives */
		double[][] aggregatedObjValues = new double[objectiveValues.length][weights.length+1];
		for (int x=0; x<objectiveValues.length; x++) {
			for (int y=0; y<weights.length; y++) {
				aggregatedObjValues[x][y] = 0;
				for (int z=0; z<objectiveValues[0].length; z++) {
					aggregatedObjValues[x][y] += weights[y][z] * objectiveValues[x][z];
				}
			}
		}
		
		/* in case of j being the empty aggregation, i.e., all weights are zero, add
		 * objective i with weight 1 to aggregation */
		boolean allWeightsZero = true;
		for (int x=0; x<weights[j].length; x++) {
			if (weights[j][x] != 0) {
				allWeightsZero = false;
			}
		}
		if (allWeightsZero) {
			weights[j][i] = 1;
			
			/* update aggregated objectives for changed objective j */
			for (int x=0; x<objectiveValues.length; x++) {
				aggregatedObjValues[x][j] = objectiveValues[x][i];
			}
			
			double[] maxAndAvgError = computeDeltaError(aggregatedObjValues, objectiveValues);
			Aggregation ret = new Aggregation(weights, maxAndAvgError[0], maxAndAvgError[1]);
			
			return ret;
		}
		
		/* otherwise, compute best aggregation and corresponding delta-errors
		 * ------------------------------------------------------------------ */
		
		/* add original objective i as last aggregated objective */
		for (int x=0; x<objectiveValues.length; x++) {
			aggregatedObjValues[x][aggregatedObjValues[0].length-1] = objectiveValues[x][i];
		}

		/* compute optimal weight and error if aggregated objective j is aggregated
		 * with original objective i */
		double[] optimalWeightAndError = computeOptimalWeightAndErrorForObjectivePairForIncrementAggregation(j, aggregatedObjValues[0].length-1, deltaErrors, aggregatedObjValues, objectiveValues, a);
		
		/* compute new aggregation */
		double optimalWeight = optimalWeightAndError[0];
		for (int y=0; y<weights[0].length; y++) {
			if (y==i) {
				weights[j][y] = optimalWeight * weights[j][y] + (1-optimalWeight)*1;
			} else {
				weights[j][y] = optimalWeight * weights[j][y];
			}
		}
				
		double deltaErrorMax = optimalWeightAndError[1];
		double deltaErrorAvg = optimalWeightAndError[2];
		
		Aggregation ret = new Aggregation(weights, deltaErrorMax, deltaErrorAvg);
		
		return ret;
	}
	
	/**
	 * PART OF THE GREEDY AGGREGATION HEURISTIC incrementally increasing the number
	 * of aggregated objectives. 
	 *
	 * returns the best weight combination ret[0] of objectives objective1 and
	 * objective2 such that the resulting delta error is minimal wrt the original
	 * objectives when the objectives objective1 and objective2 within
	 * aggregatedObjValues are aggregated.
	 * If a=1, the maximum delta error is considered (always returned in ret[1]);
	 * if a=2, the delta error averaged over all solution pairs is considered
	 * (and always returned in ret[2]). 
	 *  
	 */
	double[] computeOptimalWeightAndErrorForObjectivePairForIncrementAggregation(int objective1, int objective2, double[][] deltaErrors, double[][] aggregatedObjValues, double[][] originalObjValues,int a) {
		int numberOfSolutions = aggregatedObjValues.length;

		/* compute the two weight intervals and corresponding delta errors for
		 * each solution pair i,j (i<j) and store the interval boundary (index 0) as well
		 * as the two errors (indices 1 and 2) in intervalBound_error1_error2 */

		double[][] intervalBound_error1_error2
		   = new double[numberOfSolutions*(numberOfSolutions-1)/2][3];	// this variable contains the interval bound for each solution pair (0) and the corresponding delta errors left to the bound (1) and right to the bound (3)
		
		int numberOfCurrentSolutionPair = 0;
		for (int solution1=0; solution1<numberOfSolutions-1; solution1++) {
			for (int solution2=solution1+1; solution2<numberOfSolutions; solution2++) {
				intervalBound_error1_error2[numberOfCurrentSolutionPair]
				   = computeIntervalBoundAndErrorsForSolutionPairForIncrementalAggregation(solution1,solution2,objective1, objective2, deltaErrors, aggregatedObjValues, originalObjValues);
				numberOfCurrentSolutionPair++;
			}
		}
		
		/* 
		 * Extract best weight and delta errors (max/average)
		 * from intervalBound_error1_error2
		 */
		double[] ret = extractBestWeightAndDeltaValues(intervalBound_error1_error2, a);
	
		return ret;
	}

	/**
	 * PART OF THE GREEDY AGGREGATION HEURISTIC incrementally increasing the number
	 * of aggregated objectives. 
	 * 
	 * computes the two intervals [0,a] and [a,1] (ret[0]=a) together with the
	 * delta errors error1=ret[1] and error2=ret[2] such that the delta error for
	 * individuals solution1 and solution2 is at most error1 if objective1
	 * and objective2 of aggregatedObjValues are aggregated with a weight in
	 * [0,a] and error2 if the weight is chosen from [a,1]. The corresponding
	 * delta errors are taken from the array deltaErrors.
	 * 
	 * pre condition (in current version): objective2 == aggregatedObjValues[0].length !!!
	 */	
	double[] computeIntervalBoundAndErrorsForSolutionPairForIncrementalAggregation(int solution1, int solution2, int objective1, int objective2, double[][] deltaErrors, double[][] aggregatedObjValues, double[][] originalObjValues) {
		double[] ret = new double[3];
		
		/* extract aggregated objective vectors of examined solution pair */
		double[] objectiveVector1 = aggregatedObjValues[solution1];
		double[] objectiveVector2 = aggregatedObjValues[solution2];
		
		/* extract original objective vectors of examined solution pair */
		double[][] originalTwoObjVectors = new double[2][originalObjValues[0].length];
		for (int i=0; i<originalTwoObjVectors[0].length; i++) {
			originalTwoObjVectors[0][i] = originalObjValues[solution1][i];
			originalTwoObjVectors[1][i] = originalObjValues[solution2][i];
		}
	
		/* if the two objective vectors are indifferent or one strictly dominating
		 * the other one wrt objectives objective1 and objective 2, the weight of
		 * the aggregation does not influence the delta-error since the aggregation
		 * does not change the dominance relation:
		 */
		if ((objectiveVector1[objective1] == objectiveVector2[objective1] &&
				objectiveVector1[objective2] == objectiveVector2[objective2]) ||
			(objectiveVector1[objective1] < objectiveVector2[objective1] &&
					objectiveVector1[objective2] < objectiveVector2[objective2]) ||
			(objectiveVector1[objective1] > objectiveVector2[objective1] &&
					objectiveVector1[objective2] > objectiveVector2[objective2])) {
			
			ret[0] = 0.5; // choose weight in the middle of the interval
			
		} 
		/* if the two objective vectors are incomparable wrt objective1 and objective2
		 * or indifferent wrt one of the objectives but comparable wrt the other,
		 * the critical weight alpha can be computed */
		else {
			double numerator = objectiveVector1[objective2] - objectiveVector2[objective2];
			ret[0] = (numerator) /
			         (numerator + objectiveVector2[objective1]-objectiveVector1[objective1]);
		}
		
		/* check if pre-condition is fulfilled: */
		assert(objective2 == aggregatedObjValues[0].length);
		
		/* compute delta error for left interval */
		double[][] newAggregatedObjValues = new double[2][aggregatedObjValues[0].length-1];
		for (int j=0; j<newAggregatedObjValues[0].length; j++) {
			if (j==objective1) {
				newAggregatedObjValues[0][j] = objectiveVector1[objective2];
				newAggregatedObjValues[1][j] = objectiveVector2[objective2];
			} else {
				newAggregatedObjValues[0][j] = objectiveVector1[j];
				newAggregatedObjValues[1][j] = objectiveVector2[j];
			}
		}
		double[] deltaErrorMaxAndAvg = computeDeltaError(newAggregatedObjValues, originalTwoObjVectors);
		ret[1] = deltaErrorMaxAndAvg[0];
	
		/* compute delta error for right interval */
		newAggregatedObjValues = new double[2][aggregatedObjValues[0].length-1];
		for (int j=0; j<newAggregatedObjValues[0].length; j++) {
			if (j==objective1) {
				newAggregatedObjValues[0][j] = objectiveVector1[objective1];
				newAggregatedObjValues[1][j] = objectiveVector2[objective1];
			} else {
				newAggregatedObjValues[0][j] = objectiveVector1[j];
				newAggregatedObjValues[1][j] = objectiveVector2[j];
			}
		}
		deltaErrorMaxAndAvg = computeDeltaError(newAggregatedObjValues, originalTwoObjVectors);
		ret[2] = deltaErrorMaxAndAvg[0];
	
		return ret;
	}

	
	/**
	 * performs the greedy algorithm to reduce the number of objectives to k
	 * allowing aggregation of the initial objectives. The delta error is
	 * computed according to the parameter 'a': As the delta error, we use
	 * (i)  the maximum delta error (if a=1) and
     * (ii) the delta error averaged over all solution pairs (if a=2)
	 * 
	 */
	public conflicts.sets.Aggregation performGreedyAggregationAlgorithmGivenK(int k, int a) {
		/* compute delta error for each solution pair (i,j) for the case that
		 * we wrongly assume that i is better than j (initialization) */
		double[][] deltaErrors = computeDeltaErrorsForAllSolutionPairs();
		
		Aggregation aggregation = new Aggregation(this.os_dim);
		double maxError = 0.0;
		double averageError = 0.0;
		int current_k = this.os_dim;
		while (current_k > k) {
			/* compute best aggregation between objectives bestAggregation[0] and
			 * bestAggregation[1] with weight alpha=bestAggregation[2]. */
			double[] bestAggregation = getBestAggregation(this.values, deltaErrors, aggregation, a);
			aggregation.aggregate((int)bestAggregation[0], (int)bestAggregation[1], bestAggregation[2]);

			if (bestAggregation[3] > maxError) {
				maxError = bestAggregation[3];
			}
			averageError = bestAggregation[4];
		
			aggregation.setMaxError(maxError);
			aggregation.setAverageError(averageError);
			
			current_k = current_k - 1;
		}


		return aggregation;
	}
	
	/** PART OF GREEDY AGGREGATION ALGORITHM FOR GIVEN K
	 * returns the delta errors for each pair of solutions (i,j) if we
	 * (wrongly) assume that i is better than j.
	 */
	public double[][] computeDeltaErrorsForAllSolutionPairs() {
		double[][] deltaErrors = new double[this.values.length][this.values.length];
		for (int i=0; i<this.values.length; i++) {
			for (int j=0; j<this.values.length; j++) {
				double maxdelta = 0;
				for (int k=0; k<this.os_dim; k++) {
					if (this.values[i][k]-this.values[j][k] > maxdelta) {
						maxdelta = this.values[i][k]-this.values[j][k];
					}
				}
				deltaErrors[i][j] = maxdelta;
			}
		}
		return deltaErrors;
	}
	
	/** PART OF GREEDY AGGREGATION ALGORITHM FOR GIVEN K
	  * returns the two objectives i=ret[0] and j=ret[1] as well as weight alpha=ret[2]
	  * such that the aggregation of i and j via weight alpha minimizes the error
	  * delta=ret[3] if the original objective vectors are given by objVals and the current
	  * aggregation of the original objectives is given in aggregation (if a=1).
	  * If a=2, the two objectives i=ret[0] and j=ret[1] as well as weight alpha=ret[2]
	  * are returned such that the aggregation of i and j via weight alpha minimizes
	  * average error delta=ret[4] if the original objective vectors are given
	  * by objVals and the current aggregation of the original objectives is given
	  * in aggregation. In both cases, the other delta value (maximum delta value
	  * in ret[3] for a=2 and average error in ret[4] for the case a=1) is returned as
	  * well.  
	  */
	private double[] getBestAggregation(double[][] originalObjValues, double[][] deltaErrors, Aggregation aggregation, int a) {
		int numberOfSolutions = originalObjValues.length;
		double[][] objectiveWeights = aggregation.getAggregation();
		int numberOfAggregatedObjectives = objectiveWeights.length;
		double[][] aggregatedObjValues = new double[originalObjValues.length][this.os_dim];
		
		/* compute the transformed (aggregated) objective vectors: */
		for (int sol=0; sol<numberOfSolutions; sol++) {
			for (int aggObj=0; aggObj<numberOfAggregatedObjectives; aggObj++) {
				aggregatedObjValues[sol][aggObj] = 0;
				for (int origObj=0; origObj<this.os_dim; origObj++) {
					aggregatedObjValues[sol][aggObj]
					   += objectiveWeights[aggObj][origObj]*originalObjValues[sol][origObj];
				}
			}
		} /* aggregated objective vectors of all solutions now in aggregatedObjValues */
		
		/* compute delta error for every aggregated objective pair (i,j) according to
		 * the new weighting alpha between the two objectives (i<j): if a=1, with
		 * respect to maximum delta error; if a=2, with respect to error averaged over
		 * all solution pairs. */
		int bestI=0;
		int bestJ=0;
		/* store best maximum delta error and corresponding average error (if a=1) 
		 *    or corresponding maximum delta error and best average error (if a=2) */
		double[] bestDelta = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
		double bestWeight=0;
		for (int i=0; i<numberOfAggregatedObjectives-1; i++) {
			for (int j=i+1; j<numberOfAggregatedObjectives; j++) {
				double[] optimalWeightAndError = computeOptimalWeightAndErrorForObjectivePair(i, j, deltaErrors, aggregatedObjValues, a);
				if (optimalWeightAndError[a] < bestDelta[a-1]) {
					bestI=i;
					bestJ=j;
					bestWeight=optimalWeightAndError[0];
					bestDelta[0]=optimalWeightAndError[1];
					bestDelta[1]=optimalWeightAndError[2];
				}
			}
		}

		double[] ret = new double[5];
		ret[0] = bestI;
		ret[1] = bestJ;
		ret[2] = bestWeight;
		ret[3] = bestDelta[0];
		ret[4] = bestDelta[1];
		
		return ret;
	}
	
	/** PART OF GREEDY AGGREGATION ALGORITHM FOR GIVEN K
	  * returns the best weight combination ret[0] of objectives objective1 and
	  * objective2 such that the  resulting delta error is minimal wrt the
	  * original objectives when the objectives objective1 and objective2 are
	  * aggregated. If a=1, the maximum delta error is considered (always
	  * returned in ret[1]); if a=2, the delta error averaged over all solution
	  * pairs is considered (and always returned in ret[2]). 
	  */
	public double[] computeOptimalWeightAndErrorForObjectivePair(int objective1, int objective2, double[][] deltaErrors, double[][] aggregatedObjValues, int a) {
		int numberOfSolutions = aggregatedObjValues.length;

		/* compute the two weight intervals and corresponding delta errors for
		 * each solution pair i,j (i<j) and store the interval boundary (0) as well
		 * as the two errors (1 and 2) in intervalBound_error1_error2 */

		double[][] intervalBound_error1_error2
		   = new double[numberOfSolutions*(numberOfSolutions-1)/2][3];	// this variable contains the interval bound for each solution pair (0) and the corresponding delta errors left to the bound (1) and right to the bound (3)
		
		int numberOfCurrentSolutionPair = 0;
		for (int i=0; i<numberOfSolutions-1; i++) {
			for (int j=i+1; j<numberOfSolutions; j++) {
				intervalBound_error1_error2[numberOfCurrentSolutionPair]
				   = computeIntervalBoundAndErrorsForSolutionPair(i,j,objective1, objective2, deltaErrors, aggregatedObjValues);
				numberOfCurrentSolutionPair++;
			}
		}
		
		/* 
		 * Extract best weight and delta errors (max/average)
		 * from intervalBound_error1_error2
		 */
		double[] ret = extractBestWeightAndDeltaValues(intervalBound_error1_error2, a);
		
		
//for (int x=0; x<intervalBound_error1_error2.length;x++) {
//	for (int y=0; y<intervalBound_error1_error2[0].length;y++) {
//		System.out.print(" " + intervalBound_error1_error2[x][y]);
//	}System.out.println();
//}
//		
//System.out.println();
//System.out.println("ret[0]= " + ret[0]);
//System.out.println("ret[1]= " + ret[1]);
//System.out.println("ret[2]= " + ret[2]);
		
		return ret;
	}
	
	/** PART OF GREEDY AGGREGATION ALGORITHM FOR GIVEN K
	  *
	  * computes the two intervals [0,a] and [a,1] (ret[0]=a) together with
	  * the delta errors error1=ret[1] and error2=ret[2] such that
	  * the delta error for individuals i and j is at most error1 if objective1
	  * and objective2 are aggregated with a weight in [0,a] and error2 if
	  * the weight is chosen from [a,1].
	  * 
	  */ 
	public double[] computeIntervalBoundAndErrorsForSolutionPair(int i, int j, int objective1, int objective2, double[][] deltaErrors, double[][] aggregatedObjValues) {
		/* extract objective vectors of examined solution pair */
		double[] objectiveVector1 = aggregatedObjValues[i];
		double[] objectiveVector2 = aggregatedObjValues[j];
		
		
//System.out.println(" --------- -  - - - ");
//for (int x=0; x<objectiveVector1.length; x++) {
//	System.out.print(" " + objectiveVector1[x]);
//} System.out.println();
//for (int x=0; x<objectiveVector2.length; x++) {
//	System.out.print(" " + objectiveVector2[x]);
//} System.out.println();
//System.out.println();
//System.out.println("i= " + i);
//System.out.println("j= " + j);
//System.out.println("objective1= " + objective1);
//System.out.println("objective2= " + objective2);

		
		
		/* no error if objective vectors of i and j are indifferent wrt objective1 and
		 * objective2 */
		if (objectiveVector1[objective1] == objectiveVector2[objective1]
		   && objectiveVector1[objective2] == objectiveVector2[objective2]) {
			double[] ret = {0.5, 0.0, 0.0};
			
//System.out.println("jump back because of indifference");
			
			return ret;
		}
		/* no error if objective vectors of i and j are comparable wrt objective1 and
		 * objective2 */
		if ((objectiveVector1[objective1] < objectiveVector2[objective1]
		      && objectiveVector1[objective2] < objectiveVector2[objective2])
		   || (objectiveVector1[objective1] > objectiveVector2[objective1]
              && objectiveVector1[objective2] > objectiveVector2[objective2])) {
			double[] ret = {0.5, 0.0, 0.0};
			
//System.out.println("jump back because of comparability");
			
			return ret;
		}
		
		/* for the remaining cases, compute the objectives (except objective1 and
		 * objective2) where i has a lower objective value than j (IBetterThanJ),
		 * where j's objective values are lower than the ones of i (JbetterThanI).
		 */
		Vector<double[]> vecIBetterThanJ = new Vector<double[]>();
		Vector<double[]> vecJBetterThanI = new Vector<double[]>();
		for (int o=0; o<objectiveVector1.length; o++) {
			double[] objValues = new double[2];
			objValues[0] = objectiveVector1[o];
			objValues[1] = objectiveVector2[o];
			if (objectiveVector1[o] < objectiveVector2[o]) {
				vecIBetterThanJ.add(objValues);
			} else if (objectiveVector1[o] > objectiveVector2[o]) {
				vecJBetterThanI.add(objValues);
			}
		}
		double[][] IBetterThanJ = new double[vecIBetterThanJ.size()][2];
		IBetterThanJ = vecIBetterThanJ.toArray(IBetterThanJ);
		double[][] JBetterThanI = new double[vecJBetterThanI.size()][2];
		JBetterThanI = vecJBetterThanI.toArray(JBetterThanI);
		
		/* if objective vectors of i and j are indifferent wrt objective1 but
		 * comparable wrt objective2, there is only a possible error for
		 * weight = 1.0 */
		if (objectiveVector1[objective1] == objectiveVector2[objective1]
		   && objectiveVector1[objective2] != objectiveVector2[objective2]) {
			
			/* depending on whether i dominates j or j dominates i wrt
			 * objective2, the delta error is computed on the remaining
			 * objectives */
			double delta = 0.0;
			if (objectiveVector1[objective2] < objectiveVector2[objective2]) {
				/* Error only if no other objective with 
				 * objectiveVector1[objective2] < objectiveVector2[objective2]
				 * exists */
				if (IBetterThanJ.length <= 1) {
					delta = deltaErrors[j][i]; // by aggregating, we assume that j is better than i
				}
			} else {
				/* Error only if no other objective with 
				 * objectiveVector1[objective2] > objectiveVector2[objective2]
				 * exists */
				if (JBetterThanI.length <= 1) {
					delta = deltaErrors[i][j]; // by aggregating, we assume that i is better than j
				}
			}
			double[] ret = {1.0, 0.0, delta};
			return ret;
		}
		
		/* similar to the case above, an error can only occur for weight = 0.0
		 * if i is comparable to j wrt objective1 but indifferent wrt objective2 */
		if ((objectiveVector1[objective1] != objectiveVector2[objective1])
		   && (objectiveVector1[objective2] == objectiveVector2[objective2])) {
			
			/* depending on whether i dominates j or j dominates i wrt
			 * objective1, the delta error is computed on the remaining
			 * objectives */
			double delta = 0.0;
			if (objectiveVector1[objective1] < objectiveVector2[objective1]) {
				/* Error only if no other objective with 
				 * objectiveVector1[objective1] < objectiveVector2[objective1]
				 * exists */
				if (IBetterThanJ.length <= 1) {
					delta = deltaErrors[j][i]; // by aggregating, we assume that j is better than i
				}
			} else {
				/* Error only if no other objective with 
				 * objectiveVector1[objective1] > objectiveVector2[objective1]
				 * exists */
				if (JBetterThanI.length <= 1) {
					delta = deltaErrors[i][j]; // by aggregating, we assume that i is better than j
				}
			}
			double[] ret = {0.0, delta, 0.0};
			return ret;
		}
		
		/* Now, only the two cases i<j wrt objective1 but i>j wrt objective2
		 *                     and i>j wrt objective1 but i<j wrt objective2
		 * remain to be solved. */
		if (objectiveVector1[objective1] < objectiveVector2[objective1]
           && objectiveVector1[objective2] > objectiveVector2[objective2]) {
		
			/* weight for which the behavior changes*/
			double weight = (objectiveVector1[objective2] - objectiveVector2[objective2])
			   / (
			         objectiveVector2[objective1] - objectiveVector1[objective1]
			         + objectiveVector1[objective2] - objectiveVector2[objective2]
			     );
			/* case 0 <= \alpha <= weight, i.e., i>j wrt aggregated objective:
			 * error only if no other objective with i<j exists */
			double deltaLeft = 0.0;
			if (IBetterThanJ.length <= 1) {
				/* no other objective with i<j except objective1?
				 * then, we make an error wrt objective1 */
				deltaLeft = deltaErrors[j][i]; // by aggregating, we assume that j is better than i;
				
			}
			/* case weight < alpha <= 1, i.e., i<j wrt aggregated objective:
			 * error only if no other objective with i>j exists */
			double deltaRight = 0.0;
			if (JBetterThanI.length <= 1) {
				/* no other objective with i>j except objective2?
				 * then, we make an error wrt objective2 */
				deltaRight = deltaErrors[i][j]; // by aggregating, we assume that i is better than j
			}
			double[] ret = {weight, deltaLeft, deltaRight};
			return ret;
		}
		
		if (objectiveVector1[objective1] > objectiveVector2[objective1]
           && objectiveVector1[objective2] < objectiveVector2[objective2]) {
		
			/* weight for which the behavior changes*/
			double weight = (objectiveVector2[objective2] - objectiveVector1[objective2])
			   / (
			         objectiveVector1[objective1] - objectiveVector2[objective1]
			         + objectiveVector2[objective2] - objectiveVector1[objective2]
			     );
			/* case 0 <= \alpha <= weight, i.e., i<j wrt aggregated objective:
			 * error only if no other objective with i>j exists */
			double deltaLeft = 0.0;
			if (JBetterThanI.length <= 1) {
				/* no other objective with i>j except objective1?
				 * then, we make an error wrt objective1 */
				deltaLeft = deltaErrors[i][j]; // by aggregating, we assume that i is better than j
			}
			/* case weight < alpha <= 1, i.e., i>j wrt aggregated objective:
			 * error only if no other objective with i<j exists */
			double deltaRight = 0.0;
			if (IBetterThanJ.length <= 1) {
				/* no other objective with i<j except objective2?
				 * then, we make an error wrt objective 2: */
				deltaRight = deltaErrors[j][i]; // by aggregating, we assume that j is better than i
			}
			double[] ret = {weight, deltaLeft, deltaRight};
			return ret;
		}
		
		/* Hopefully, this part will never be reached... */
		System.out.println("Error occured in method");
		System.out.println("   computeIntervalBoundAndErrorsForSolutionPair(int i, int j, int objective1, int objective2, double[][] aggregatedObjValues)");
		System.out.println();
		return null;
	}
		
	
	/** PART OF GREEDY AGGREGATION ALGORITHM FOR GIVEN K
	  *
	  * extracts the best aggregation weight (ret[0]) and delta errors (ret[1] = maximum
	  * delta error, ret[2]=average delta error) for all solutions within
	  * intervalBound_error1_error2.
	  * 
	  * If a=1, the weight is determined by the maximum delta error, we make;
	  * if a=2, the weight is determined by the average delta error, we make.
	  * 
	  * intervalBound_error1_error2 contains the interval bounds
	  * b=intervalBound_error1_error2[i][0] and the corresponding delta errors of
	  * the intervals [0,b] (in intervalBound_error1_error2[i][1]) and
	  * [b,1] (in intervalBound_error1_error2[i][2]) for each solution pair (i.e.
	  * each row i of intervalBound_error1_error2).
	  */
	private double[] extractBestWeightAndDeltaValues(double[][] intervalBound_error1_error2, int a) {
		
		/* extract intervals as sorted array with second column for errors */
		Vector<Double> intervalBounds = new Vector<Double>();
		for (int i=0; i<intervalBound_error1_error2.length; i++) {
			/* insert the new interval [0,intervalBound] iff these bound was not inserted before */
			boolean included = false;
			for (Double bound : intervalBounds) {
				if (bound.doubleValue() == intervalBound_error1_error2[i][0]) {
					included = true;
				}
			}
			if (!included) {
				intervalBounds.add(new Double(intervalBound_error1_error2[i][0]));
			}
		}	
	
		/* insert interval [?,1] as well */
		intervalBounds.add(new Double(1.0));
		
		Collections.sort(intervalBounds);
		
		/* store intervals with maximum delta error (in intervalsWithErrors[j][1])
		 * and average delta error (in intervalsWithErrors[j][2]) for solution
		 * pair j */
		double[][] intervalsWithErrors = new double[intervalBounds.size()][3];
		int solutionPair=0;
		for (Double intervalBound : intervalBounds) {
			intervalsWithErrors[solutionPair][0] = intervalBound.doubleValue();
			intervalsWithErrors[solutionPair][1] = 0;
			intervalsWithErrors[solutionPair][2] = 0;
			solutionPair++;
		}
		
		/* compute maximal delta error over all solution pairs for all intervals */
		for (int i = 0; i<intervalBound_error1_error2.length; i++) {
			for (int j = 0; j<intervalsWithErrors.length; j++) {
				/* left from interval bound? --> take first error */
				if (intervalBound_error1_error2[i][0] >= intervalsWithErrors[j][0]) {
					/* test if current error is larger than worst error found so far */
					if (intervalBound_error1_error2[i][1] > intervalsWithErrors[j][1]) {
						intervalsWithErrors[j][1] = intervalBound_error1_error2[i][1]; 
					}
					intervalsWithErrors[j][2] += intervalBound_error1_error2[i][1];
				} else { /* right from interval bound? --> take second error */
					/* test if current error is larger than worst error found so far */
					if (intervalBound_error1_error2[i][2] > intervalsWithErrors[j][1]) {
						intervalsWithErrors[j][1] = intervalBound_error1_error2[i][2]; 
					}
					intervalsWithErrors[j][2] += intervalBound_error1_error2[i][2];
				}
			}
		}
		
		/* average error over number of solution pairs */
		for (int j = 0;  j<intervalsWithErrors.length; j++) {
			intervalsWithErrors[j][2] = intervalsWithErrors[j][2]/intervalBound_error1_error2.length;
		}

		/* 
		 * merge connected intervals that have the same (maximum for a=1;
		 * average if a=2) delta error:
		 */
		boolean connectedIntervalsWithSameError = false;
		for (int i=0; i<intervalsWithErrors.length-1; i++) {
			if (intervalsWithErrors[i][a]==intervalsWithErrors[i+1][a]) {
				connectedIntervalsWithSameError = true;
				break;
			}
		}
		if (connectedIntervalsWithSameError) {
			intervalsWithErrors = mergeConnectedIntervalsWithSameErrors(intervalsWithErrors, a);
		}
		
		/* Find interval with lowest error and compute weight as
		 * mid point of this interval. */
		double bestLeftIntervalBound = 0.0;
		double bestRightIntervalBound = 1.0;
		double leftIntervalBound = 0.0;
		double rightIntervalBound = 0.0;
		/* save current best delta error in minimumDelta[a-1] and the
		 * remaining delta error in minimumDelta[a mod 2] */
		double[] minimumDelta = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
		for (int i=0; i<intervalsWithErrors.length; i++) {
			rightIntervalBound = intervalsWithErrors[i][0];
			if (intervalsWithErrors[i][a] < minimumDelta[a-1]) {
				minimumDelta[a-1] = intervalsWithErrors[i][a];
				minimumDelta[a%2] = intervalsWithErrors[i][a%2 + 1];
				bestLeftIntervalBound = leftIntervalBound;
				bestRightIntervalBound = rightIntervalBound;
			}
			leftIntervalBound = rightIntervalBound;
		}
		
		/* return optimal weight (ret[0]), maximum delta error (ret[1]),
		 * and average delta error (ret[2]) */
		double[] ret = {(bestLeftIntervalBound+bestRightIntervalBound)/2, minimumDelta[0], minimumDelta[1]};
		/* TODO: one could use different strategies here to set the optimal weight:
		 *       e.g., - use mid point of best interval (as implemented)
		 *             - use always the left interval bound
		 *             - use always the right interval bound
		 *             - sample uniformly in the best interval 
		 */
		
		return ret;
	}
	
	/** PART OF GREEDY AGGREGATION ALGORITHM FOR GIVEN K
	 *
	 * merges connected intervals with the same (maximum for a=1;
	 * average if a=2) delta error into one and returns the
	 * intervals with the corresponding delta errors.
	 */	
	private double[][] mergeConnectedIntervalsWithSameErrors(double[][] intervalsWithErrors, int a) {
		/* go through intervals and add interval to LinkedList only if delta error changes */
		LinkedList<double[]> intervalList = new LinkedList<double[]>();
		for (int i=0; i<intervalsWithErrors.length-1; i++) {
			double[] lastInterval = intervalsWithErrors[i];
			double[] currentInterval = intervalsWithErrors[i+1];

			/* insert interval only if errors are not the same 
			 * otherwise explicit merging */
			if (lastInterval[a] != currentInterval[a]) {
				intervalList.add(lastInterval);
			} else {
				currentInterval[1] = Math.max(currentInterval[1], lastInterval[1]);
				currentInterval[2] = Math.max(currentInterval[2], lastInterval[2]);
			}
		}
		/* last interval is always inserted */
		intervalList.add(intervalsWithErrors[intervalsWithErrors.length-1]);
		
		double[][] returnIntervals = new double[intervalList.size()][3];
		intervalList.toArray(returnIntervals);
		
		return returnIntervals;
	}
	
	/** returns true iff individual i epsilon-dominates individual j w.r.t. objective
	  * subset ef, where the additive epsilon-dominance relation for minimization is used:
	  * 
	  * i epsilon-dominates j 
	  *    :\Leftrightarrow \forall k\in ef: f_k(i) - epsilon \leq f_k(j)
	  */
	private boolean epsilonDominates(int i, int j, ObjectiveSet ef, double epsilon) {
		boolean[] objectives = ef.getElements();		
		for (int k=0; k<objectives.length; k++) {
			if (objectives[k]) {				
				/* use static variable failure to handle arithmetic errors: */
				if ((values[i][k] - epsilon - failure) > values[j][k]) {
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Performs the greedy delta-MOSS algorithm as described in bz2006d.
	 *
	 */
	public ObjectiveSet performGreedyAlgorithmGivenDelta(double delta) {
		int a = values.length;
		boolean[][] rel = new boolean[a][a];
		for (int i=0;i<a; i++) {
			for (int j=0;j<a; j++) {
				rel[i][j] = true;
			}
		}
		Relation R = new Relation(0, rel);
		R = R.minus(0, this.dominanceRelation);
		boolean[] chosen = new boolean[this.os_dim];

		while (R.numberOfRelatedPairs != 0) {
			Relation smallest = new Relation(System.currentTimeMillis(), rel);
			int smallestSize = Integer.MAX_VALUE;
			int theSmallest = 0;
			for (int i=0; i<this.os_dim; i++) {
				if (!chosen[i]) {
					chosen[i] = true; // set only temporarily
					Relation zeroDeltaDominance = computeZeroDeltaDominance(chosen, delta);
					Relation current = this.relations[i].intersect(System.currentTimeMillis(), R);
					current = current.minus(System.currentTimeMillis(), zeroDeltaDominance);
					int currentSize = current.getNumberOfRelatedPairs();
					if (currentSize < smallestSize) {
						smallest = current;
						smallestSize = currentSize;
						theSmallest = i;
					}
					chosen[i] = false; // reset the temporarily chsen objective
				}
			}
			R = smallest;
			chosen[theSmallest] = true;
		}
		
		/* compute correct delta value of chosen objective set: */
		ObjectiveSet chosenObj = new ObjectiveSet(chosen);
		boolean[] notchosen = new boolean[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			notchosen[i] = !chosen[i];
		}
		ObjectiveSet notchosenObj = new ObjectiveSet(notchosen);
		ObjectiveSet ret = new ObjectiveSet(chosen, this.getDeltaMinFor(chosenObj, notchosenObj));
		
		return ret;
	}
	
	
	/**
	 *  computes the dominance relation $\preceq^{0,\delta}_{chosen, notchosen}$.
	 */
	private Relation computeZeroDeltaDominance(boolean[] chosen, double delta) {
		int ds_dim = this.values.length;
		Relation zerodelta = new Relation(ds_dim);
		ObjectiveSet chosenSet = new ObjectiveSet(chosen);
		boolean[] notchosen = new boolean[chosen.length];
		for (int i=0; i<chosen.length; i++) {
			notchosen[i] = !chosen[i];
		}
		ObjectiveSet notchosenSet = new ObjectiveSet(notchosen);
		/* compute the relation between the decision vectors p and q */
		for (int p=0; p<ds_dim; p++) {
			for (int q=0; q<ds_dim; q++) {
				if (this.weaklyDominates(p,q, chosenSet) &&
										this.epsilonDominates(p,q,notchosenSet, delta)) {
					zerodelta.setinrelation(p,q,true);
				}
			}
		}
		return zerodelta;
	}
	
	public ObjectiveSet performGreedyAlgorithm() {
		return performGreedyAlgorithmGivenDelta(0);
	}
	
	
	/** Alternative algorithm using the greedy algorithm for given K !!! */
	public ObjectiveSet performGreedyAlgorithmGivenDelta2(double delta) {
		int[] ints = {1};
		ObjectiveSet os = new ObjectiveSet(ints, 2, Double.MAX_VALUE);
		int k = 1;
		while (os.getDelta() > delta) {
			os = this.performGreedyAlgorithmGivenK(k);
			k++;
		}
		return os;
	}
	
	
	
	//*********************************************************************************
	
	/** Here, the implementation of an alternative greedy algorithm starts. It is based
	  * on a kind of hierarchical clustering of the objectives, cf. brockho2007a.
	  * 
	  * if method == 0, the choice which objective to choose in the single steps is determined
	  * by the delta-errors between the chosen objectives, e.g. if
	  *     A={f_1,f_2,f_3}, chosen obj. f_2 and
	  *     B={f_4,f_5}, chosen obj. f_4
	  * are considered, the delta error is computed only for the pair {f_2,f_4} and the new
	  * chosen objective can be either f_2 or f_4.
	  * 
	  * if method == 1, each combination of objectives in the considered sets are considered, 
	  * e.g. if
	  *     A={f_1,f_2,f_3}, chosen obj. f_2 and
	  *     B={f_4,f_5}, chosen obj. f_4
	  * are considered, the delta error is computed as the minimum delta error of
	  *     f_1 wrt. {f_2,f_3,f_4,f_5},
	  *     f_2 wrt. {f_1,f_3,f_4,f_5},
	  *     ...
	  *     f_5 wrt. {f_1,f_2,f_3,f_4}
	  * and each objective within the two sets can become the new chosen objective.  
	  * 
	  * The argument outputfilename indicates whether the output will be written to the file
	  * 'outputfilename' or to standard out (if 'outputfilename' is empty)
	  *   
	  */
	
	public void computeTree(int method, String outputfilename) {
		if (method == 0) {
			computeTreeWithPairwiseComparison(outputfilename);			
		} else if (method == 1 || method == 2) {
			computeTreeWithAllCombinations(method, outputfilename);
		}
	}

	public void computeTree2PG(int method, String outputfilename, String outputfilename2PG) {
                computeTreeWithPairwiseComparison2PG(outputfilename,outputfilename2PG);
                /*
		if (method == 0) {
			computeTreeWithPairwiseComparison(outputfilename);			
		} 
                else if (method == 1 || method == 2) {
			computeTreeWithAllCombinations(method, outputfilename);
		}*/
	}

	private void computeTreeWithPairwiseComparison2PG(String outputfilename, String outputfilename2PG) {
		Vector<String> toPrint = new Vector<String>(); // storing output
                Vector<String> toPrint2PG = new Vector<String>(); // storing output for 2PG
                
		/* preprocessing: compute minimal delta-errors between all
		 * objective pairs (i,j) and corresponding objective, with
		 * which the smaller error can be achieved. 
		 * 
		 * delta-error of objective pair (i,j) stands in data[i][j],
		 * the corresponding chosen objective in data[j][i].
		 */
		double[][] data = new double[this.os_dim][this.os_dim];
				
		/* temporary variables: */
		double deltaij;
		double deltaji;
		
		/* initialize objective sets in tree as TreeSets with one objective,
		 * delta=0 (because no objective is omitted) and the single objective
		 * as its chosen objective */
		int[] elements = new int[1];
		TreeSet[] objsets = new TreeSet[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			elements[0] = i;
			objsets[i] = new TreeSet(elements, this.os_dim, 0);
				// constructor of TreeSet chooses already the correct objective!			
		}		
			
		/* compute deltas and chosenobjectives: */
		for (int i=0; i<this.os_dim; i++) {
			for (int j=i+1; j<this.os_dim; j++) {
				deltaij = getDeltaMinFor(objsets[i], objsets[j]);
				deltaji = getDeltaMinFor(objsets[j], objsets[i]);
				
				if (deltaij < deltaji) {
					/* choose objective i because of better delta error */
					data[i][j] = deltaij;
					data[j][i] = i;
				} else {
					/* choose objective j because of better delta error */
					data[i][j] = deltaji;
					data[j][i] = j;
				}
			}
		} // end of preprocessing
		
		/* print current objective sets: */
		for (int i=0; i<objsets.length; i++) {
			toPrint.add(objsets[i].toString());
		}
				
		/* perform 'os_dim' cycles of merging two objective sets together: */
		TreeSet[] currentObjectiveSets = objsets;
		for (int k=this.os_dim; k>1; k--) {
			TreeSet[] newObjectiveSets = new TreeSet[k-1];
			
			int besti = 0;
			int bestj = 1;
			double bestDelta = Double.MAX_VALUE;
			/* find objective sets with smallest delta for merging: */
			for (int i=0; i<currentObjectiveSets.length; i++) {
				for (int j=i+1; j<currentObjectiveSets.length; j++) {
					if (data[currentObjectiveSets[i].getChosenObjective()][currentObjectiveSets[j].getChosenObjective()] < bestDelta) {
						besti = i;
						bestj = j;
						bestDelta = data[currentObjectiveSets[besti].getChosenObjective()][currentObjectiveSets[bestj].getChosenObjective()];
					}
				}
			}
			for (int i=0; i<k; i++) {
				if (i < besti) {
					newObjectiveSets[i] = currentObjectiveSets[i];
				} else if (i == besti) {
					/* compute merged objective set with new delta error as
					 * maximum of the errors of smaller sets and
					 * the error between obj_besti and obj_bestj
					 */
					
					int chosenObj1 = currentObjectiveSets[besti].getChosenObjective();
					int chosenObj2 = currentObjectiveSets[bestj].getChosenObjective();
					/* guarantee that chosenObj1 < chosenObj2 for accessing
					 * data[chosenObj1][chosneObj2]*/
					if (chosenObj1 > chosenObj2) {
						int temp = chosenObj1;
						chosenObj1 = chosenObj2;
						chosenObj2 = temp;						
					}
					newObjectiveSets[i] = new TreeSet(currentObjectiveSets[besti],
													currentObjectiveSets[bestj],
													Math.max(
															Math.max(data[chosenObj1][chosenObj2],
																	currentObjectiveSets[besti].getDelta()),
																	currentObjectiveSets[bestj].getDelta()),
																	new Double(data[chosenObj2][chosenObj1]).intValue());
				} else if (i < bestj) {
					newObjectiveSets[i] = currentObjectiveSets[i];
				} else if (i > bestj) {
					newObjectiveSets[i-1] = currentObjectiveSets[i];
				}
			}
			
			/* print new objective set: */
			toPrint.add("-------------------------------------");
                        toPrint2PG.add("-------------------------------------");
			for (int i=0; i<newObjectiveSets.length; i++) {
				toPrint.add(newObjectiveSets[i].toString());
                                toPrint2PG.add(newObjectiveSets[i].toString2PG());
                                
			}
			currentObjectiveSets = newObjectiveSets;
		}
		
		/* finally output everything */
		Output.print(toPrint, outputfilename);            
                Output.print2PG(toPrint2PG, outputfilename2PG);
	}
	
	private void computeTreeWithPairwiseComparison(String outputfilename) {
		Vector<String> toPrint = new Vector<String>(); // storing output
                
		/* preprocessing: compute minimal delta-errors between all
		 * objective pairs (i,j) and corresponding objective, with
		 * which the smaller error can be achieved. 
		 * 
		 * delta-error of objective pair (i,j) stands in data[i][j],
		 * the corresponding chosen objective in data[j][i].
		 */
		double[][] data = new double[this.os_dim][this.os_dim];
				
		/* temporary variables: */
		double deltaij;
		double deltaji;
		
		/* initialize objective sets in tree as TreeSets with one objective,
		 * delta=0 (because no objective is omitted) and the single objective
		 * as its chosen objective */
		int[] elements = new int[1];
		TreeSet[] objsets = new TreeSet[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			elements[0] = i;
			objsets[i] = new TreeSet(elements, this.os_dim, 0);
				// constructor of TreeSet chooses already the correct objective!			
		}		
			
		/* compute deltas and chosenobjectives: */
		for (int i=0; i<this.os_dim; i++) {
			for (int j=i+1; j<this.os_dim; j++) {
				deltaij = getDeltaMinFor(objsets[i], objsets[j]);
				deltaji = getDeltaMinFor(objsets[j], objsets[i]);
				
				if (deltaij < deltaji) {
					/* choose objective i because of better delta error */
					data[i][j] = deltaij;
					data[j][i] = i;
				} else {
					/* choose objective j because of better delta error */
					data[i][j] = deltaji;
					data[j][i] = j;
				}
			}
		} // end of preprocessing
		
		/* print current objective sets: */
		for (int i=0; i<objsets.length; i++) {
			toPrint.add(objsets[i].toString());
		}
				
		/* perform 'os_dim' cycles of merging two objective sets together: */
		TreeSet[] currentObjectiveSets = objsets;
		for (int k=this.os_dim; k>1; k--) {
			TreeSet[] newObjectiveSets = new TreeSet[k-1];
			
			int besti = 0;
			int bestj = 1;
			double bestDelta = Double.MAX_VALUE;
			/* find objective sets with smallest delta for merging: */
			for (int i=0; i<currentObjectiveSets.length; i++) {
				for (int j=i+1; j<currentObjectiveSets.length; j++) {
					if (data[currentObjectiveSets[i].getChosenObjective()][currentObjectiveSets[j].getChosenObjective()] < bestDelta) {
						besti = i;
						bestj = j;
						bestDelta = data[currentObjectiveSets[besti].getChosenObjective()][currentObjectiveSets[bestj].getChosenObjective()];
					}
				}
			}
			for (int i=0; i<k; i++) {
				if (i < besti) {
					newObjectiveSets[i] = currentObjectiveSets[i];
				} else if (i == besti) {
					/* compute merged objective set with new delta error as
					 * maximum of the errors of smaller sets and
					 * the error between obj_besti and obj_bestj
					 */
					
					int chosenObj1 = currentObjectiveSets[besti].getChosenObjective();
					int chosenObj2 = currentObjectiveSets[bestj].getChosenObjective();
					/* guarantee that chosenObj1 < chosenObj2 for accessing
					 * data[chosenObj1][chosneObj2]*/
					if (chosenObj1 > chosenObj2) {
						int temp = chosenObj1;
						chosenObj1 = chosenObj2;
						chosenObj2 = temp;						
					}
					newObjectiveSets[i] = new TreeSet(currentObjectiveSets[besti],
													currentObjectiveSets[bestj],
													Math.max(
															Math.max(data[chosenObj1][chosenObj2],
																	currentObjectiveSets[besti].getDelta()),
																	currentObjectiveSets[bestj].getDelta()),
																	new Double(data[chosenObj2][chosenObj1]).intValue());
				} else if (i < bestj) {
					newObjectiveSets[i] = currentObjectiveSets[i];
				} else if (i > bestj) {
					newObjectiveSets[i-1] = currentObjectiveSets[i];
				}
			}
			
			/* print new objective set: */
			toPrint.add("-------------------------------------");
			for (int i=0; i<newObjectiveSets.length; i++) {
				toPrint.add(newObjectiveSets[i].toString());
			}
			currentObjectiveSets = newObjectiveSets;
		}
		
		/* finally output everything */
		Output.print(toPrint, outputfilename);                
	}

	/* if method == 1, the delta error between the single objectives and
	 *                 all objectives within the corresponding subsets
	 * if method == 2, the delta error between the single objectives and
	 *                 the entire objective set
	 */
	private void computeTreeWithAllCombinations(int method, String outputfilename) {
		Vector<String> toPrint = new Vector<String>(); // storing output
		
		/* initialize objective sets in tree as TreeSets with one objective,
		 * delta=0 (because no objective is omitted) and the single objective
		 * as its chosen objective */
		int[] elements = new int[1];
		TreeSet[] objsets = new TreeSet[this.os_dim];
		for (int i=0; i<this.os_dim; i++) {
			elements[0] = i;
			objsets[i] = new TreeSet(elements, this.os_dim, 0);
			/* Note, that constructor of TreeSet already chooses the correct objective! */			
		}		
		
		/* print current objective sets: */
		for (int i=0; i<objsets.length; i++) {
			toPrint.add(objsets[i].toString());
		}
				
		/* perform 'os_dim' cycles of merging two objective sets together: */
		TreeSet[] currentObjectiveSets = objsets;
		for (int k=this.os_dim; k>1; k--) {
			TreeSet[] newObjectiveSets = new TreeSet[k-1];
			
			int besti = 0;
			int bestj = 1;
			double bestDelta = Double.MAX_VALUE;
			int bestChosen = -1; // this objective should be chosen to achieve the best delta
			                     // NOTE: this objective can be in {0,...,os_dim -1} !!!
			
			/* find objective sets with smallest delta for merging: */
			for (int i=0; i<currentObjectiveSets.length; i++) {
				for (int j=i+1; j<currentObjectiveSets.length; j++) {
					TreeSet currentBest = getCurrentBest(currentObjectiveSets[i], currentObjectiveSets[j], method);
					if (currentBest.getDelta() < bestDelta) {
						besti = i;
						bestj = j;
						bestDelta = currentBest.getDelta();
						bestChosen = currentBest.getChosenObjective();
					}
				}
			}
			for (int i=0; i<k; i++) {
				if (i < besti) {
					newObjectiveSets[i] = currentObjectiveSets[i];
				} else if (i == besti) {
					/* compute merged objective set with new delta error as
					 * maximum of the errors of smaller sets and
					 * the error between obj_besti and obj_bestj
					 */
					int chosenObj1 = currentObjectiveSets[besti].getChosenObjective();
					int chosenObj2 = currentObjectiveSets[bestj].getChosenObjective();
					/* guarantee that chosenObj1 < chosenObj2 for accessing
					 * data[chosenObj1][chosneObj2]*/
					if (chosenObj1 > chosenObj2) {
						int temp = chosenObj1;
						chosenObj1 = chosenObj2;
						chosenObj2 = temp;						
					}
					newObjectiveSets[i] = new TreeSet(currentObjectiveSets[besti],
													currentObjectiveSets[bestj],
													bestDelta,
													bestChosen);
				} else if (i < bestj) {
					newObjectiveSets[i] = currentObjectiveSets[i];
				} else if (i > bestj) {
					newObjectiveSets[i-1] = currentObjectiveSets[i];
				}
			}
			
			/* print new objective set: */
			toPrint.add("-------------------------------------");
			for (int i=0; i<newObjectiveSets.length; i++) {
				toPrint.add(newObjectiveSets[i].toString());
			}
			currentObjectiveSets = newObjectiveSets;
		}
		
		/* finally output everything */
		Output.print(toPrint, outputfilename);
	}
	
	/**
	 * PART OF TREE-BASED GREEDY ALGORITHM
	 * returns a TreeSet with error delta and chosen objective chosenObj, where
	 * delta is the minimum error over all possible single objectives wrt. the remaining
	 * objectives in ObjSet1 \cup objSet2. 
	 * 
	 * if method == 1, the delta error between the single objectives and
	 *                 all objectives within the corresponding subsets
	 * if method == 2, TODO THIS IS NOT CORRECTLY IMPLEMENTED YET !!!
	 *                 the delta error between the single objectives and
	 *                 the entire objective set 
	 */
	private TreeSet getCurrentBest(TreeSet objSet1, TreeSet objSet2, int method) {
		double currentBestDelta = Double.MAX_VALUE;
		int currentBestChosen = -1;
		TreeSet ret = new TreeSet();
		
		/* store all elements in objSet1 \cup objSet2 within variable 'element': */
		ObjectiveSet objSet = new ObjectiveSet(objSet1);
		objSet.addAll(objSet2);
		boolean[] element = objSet.getElements();
		
		/* instantiate an ObjectiveSet with all objectives: */
		boolean[] a = new boolean[this.os_dim];
		for (int s=0; s<a.length; s++) {
			a[s] = true;					
		}
		ObjectiveSet allObjectives = new ObjectiveSet(a);
		
		for (int i=0; i<element.length; i++) {
			if (element[i]) {
				/* temporary variables (NOT FAST AND CLEVERLY IMPLEMENTED): */
				boolean[] e = new boolean[this.os_dim];
				boolean[] rest = objSet.getComplement();
				for (int s=0; s<e.length; s++) {
					e[s] = false;					
				}
				e[i] = true;
				rest[i] = true;
				double delta = -1;
				/* compute delta error between objective i and rest according to method: */
				if (method == 1) {
					delta = getDeltaMinFor(new ObjectiveSet(e), objSet);
				} else if (method == 2) {
					delta = getDeltaMinFor(new ObjectiveSet(rest), allObjectives);
				}
				
				System.out.println("i=" + i + " and delta=" + delta);
				
				/* decide whether objective i is best objective */
				if (delta < currentBestDelta) {
					currentBestDelta = delta;
					currentBestChosen = i;
				}
			}
		}
	
		ret.setDelta(currentBestDelta);
		ret.setChosenObjective(currentBestChosen);
		
		return ret;
	}
	
	
	/**
	 * Computes the maximum and average delta-error between two solution sets given by
	 * their objective vectors.
	 * 
	 * @param objectiveVectorsSet1 first set of objective vectors
	 * @param objêctiveVectorsSet2 second set of objective vectors
	 * @return the maxmimum (in ret[0]) and average (in ret[1]) delta error
	 *         between the two solution sets, the objective vectors of which
	 *         are given in objectiveVectorsSet1 and objectiveVectorsSet2.
	 *         
	 * pre condition1: both objectiveVectorsSet1 and objectiveVectorsSet2 need to have
	 *                 the same number of solutions, i.e., the length of the arrays in
	 *                 the first dimension have to be the same. The number of objectives
	 *                 in both objective Vector sets do not need to be the same.
	 * pre condition2: the set 1 should always induce an order that is a refinement of the
	 *                 order induced by set 2, i.e. if a weakly dominates b wrt set 1
	 *                 a should weakly dominate b also wrt set 2 or a and b are incomparable
	 *                 wrt set 2 and if a and b are incomparable wrt set 1, they should also
	 *                 be incomparable wrt set 2.
	 */
	public double[] computeDeltaError(double[][] objectiveVectorsSet1, double[][] objectiveVectorsSet2) {
		
//System.out.println("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD");
//System.out.println("   computeDeltaError(...)");
//System.out.println("objectiveVector 1 wrt ObjectiveSet1: ");
//System.out.println("---------------------");
//for (int i=0; i<objectiveVectorsSet1.length; i++) {
//	for (int j=0; j<objectiveVectorsSet1[0].length; j++) {
//		System.out.print(" " + objectiveVectorsSet1[i][j]);
//	}
//	System.out.println();
//}
//System.out.println("objectiveVector 1 wrt ObjectiveSet2: ");
//System.out.println("---------------------");
//for (int i=0; i<objectiveVectorsSet2.length; i++) {
//	for (int j=0; j<objectiveVectorsSet2[0].length; j++) {
//		System.out.print(" " + objectiveVectorsSet2[i][j]);
//	}
//	System.out.println();
//}


		
		
		double maxError = 0;
		double avgError = 0;
		
		int relationshipInSet1 = 0; // relationship between the current solution pair in set 1:
   		                            // can be 0: indifferent
		                            //        1: a weakly dominates b but not indifferent
		                            //        2: a weakly dominated by b but not indifferent
		                            //        3: a and b are incomparable
		int relationshipInSet2 = 0; // relationship between the current solution pair in set 2:
                                    // can be 0: indifferent
                                    //        1: a weakly dominates b but not indifferent
                                    //        2: a weakly dominated by b but not indifferent
                                    //        3: a and b are incomparable
		double maxErrorBetterSet2 = 0; // maximum error max(f_i(a) - f_i(b), 0) if we wrongly
                                       // assume that a is better than b wrt set 1 (i\in set 2)
		double maxErrorWorseSet2 = 0;  // maximum error max(f_i(b) - f_i(a), 0) if we wrongly
                                       // assume that b is better than a wrt set 1 (i\in set 2)
		
		/* compute delta error for each solution pair*/
		for (int i = 0; i<objectiveVectorsSet1.length; i++) {
			for (int j = i+1; j<objectiveVectorsSet1.length; j++) {
				relationshipInSet1 = 0;
				relationshipInSet2 = 0;
				maxErrorBetterSet2 = 0;
				maxErrorWorseSet2 = 0;
				for (int k=0; k<objectiveVectorsSet1[0].length; k++) {
					/* compute relationship in set 1 */
					if (objectiveVectorsSet1[i][k] < objectiveVectorsSet1[j][k]) {
						if (relationshipInSet1 == 2 || relationshipInSet1 == 3) {
							relationshipInSet1 = 3;
						} else {
							relationshipInSet1 = 1;
						}
					} else if (objectiveVectorsSet1[i][k] > objectiveVectorsSet1[j][k]) {
						if (relationshipInSet1 == 1 || relationshipInSet1 == 3) {
							relationshipInSet1 = 3;
						} else {
							relationshipInSet1 = 2;
						}
					}
				}
				for (int k=0; k<objectiveVectorsSet2[0].length; k++) {
					/* compute relationship in set 2 */
					if (objectiveVectorsSet2[i][k] < objectiveVectorsSet2[j][k]) {
						if (relationshipInSet2 == 2 || relationshipInSet2 == 3) {
							relationshipInSet2 = 3;
						} else {
							relationshipInSet2 = 1;
						}
						maxErrorWorseSet2 = Math.max(maxErrorWorseSet2, objectiveVectorsSet2[j][k] - objectiveVectorsSet2[i][k]);
					} else if (objectiveVectorsSet2[i][k] > objectiveVectorsSet2[j][k]) {
						if (relationshipInSet2 == 1 || relationshipInSet2 == 3) {
							relationshipInSet2 = 3;
						} else {
							relationshipInSet2 = 2;
						}
						maxErrorBetterSet2 = Math.max(maxErrorBetterSet2, objectiveVectorsSet2[i][k] - objectiveVectorsSet2[j][k]);
					}
				}
							
				/* update delta errors maxError and avgError if relationship in set 1
				 * and set 2 are different */
				double error = 0; 
				if (relationshipInSet1 != relationshipInSet2) {
					if ((relationshipInSet1 == 0 && relationshipInSet2 == 1) ||
							(relationshipInSet1 == 2 && relationshipInSet2 == 3)) {
						error = maxErrorWorseSet2;
					}
					if ((relationshipInSet1 == 0 && relationshipInSet2 == 2) ||
							(relationshipInSet1 == 1 && relationshipInSet2 == 3)) {
						error = maxErrorBetterSet2;
						
					}
					if (relationshipInSet1 == 0 && relationshipInSet2 == 3) {
						error = Math.max(maxErrorWorseSet2, maxErrorBetterSet2);
					}
				}
				
				maxError = Math.max(maxError, error);
				avgError = avgError + error;
			}
		}
		
		double[] ret = new double[2];
		ret[0] = maxError;
		/* error averaged over all solution pairs = accumulated error / #(solution pairs) */
		ret[1] = avgError*2/(objectiveVectorsSet1.length*(objectiveVectorsSet1.length-1));
		
//System.out.println();
//System.out.println("maxError= " + ret[0]);
//System.out.println("avgError= " + ret[1]);
//System.out.println("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD");

		return ret;
	}
	
	private void printDoubleMatrix(double[][] d) {
		for (int i=0; i<d.length; i++) {
			for(int j=0; j<d[0].length; j++) {
				System.out.print(d[i][j] + " ");
			}System.out.println();
		}
	}
}
