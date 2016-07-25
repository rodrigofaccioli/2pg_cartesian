
package conflicts;

import conflicts.sets.IntSet;
import conflicts.sets.SetOfIntSets;

public class MOSSExactAlgorithm {

	Population pop;			// set of solutions
	int os_dim;				// number of objectives
	Relation[] relations;	// an array of the relations $\preceq_i$
	
	public MOSSExactAlgorithm(Population pop, int os_dim, Relation[] rels) {
		this.pop = pop;
		this.os_dim = os_dim;
		this.relations = rels;
	}

	
	/**
	 * Performs exact algorithm for MOSS problem as described in bz2007d,
	 * given the population pop and the dominance relations $\preceq_i$
	 * stored in relations.
	 */
	public IntSet performExactAlgorithm() {
		/* Possible sets of objectives are saved in possibleSets; 
		 * at the end, one set out of them is chosen.
		 */
		SetOfIntSets possibleSets = new SetOfIntSets();
		SetOfIntSets currentSet;

		/* A minimum set is only defined for more than 1 individual: */
		if (this.pop.mu < 2) {
			return null;
		}
		/* Compute the possible minimum objective sets for all pairs
		 * of individuals: */
		for (int i=0; i<this.pop.mu-1;i++) {
			for (int j=i+1;j<this.pop.mu;j++) {
				currentSet = computeSetOfSetsFor(i,j);
				if (currentSet != null) {
					possibleSets.union_ExactAlgo(currentSet);
				}
			}
		}
		return possibleSets.getSmallestIntSet();
	}
	
	
	/**
     * @param p
     * @param q
     * @return the set of all objective sets non-conflicting with the
     * entire objective set wrt a given pair (p,q) of individuals.
     * 
     * If p and q are incomparable, the method will return all objective sets
     * which explain the incomparability of p and q corresponding to the
     * current dominance relation. If p and q are comparable, the method will
     * return the set of objectives, where p and q are not equal. If p and q
     * are indifferent all single objective sets are handed back. 
     */
    private SetOfIntSets computeSetOfSetsFor(int p, int q) {
    	SetOfIntSets sois = new SetOfIntSets();
    	
    	SetOfIntSets pSet = new SetOfIntSets(); // pSet = {i|p >_i q and q \not>_i p}
    	SetOfIntSets qSet = new SetOfIntSets(); // qSet = {i|q >_i p and p \not>_i q}
    	
    	for (int i=0; i<this.os_dim; i++) {
    		if (relations[i].inrelation(p,q) && !(relations[i].inrelation(q,p))) {
    			int singleobjective[] = new int[1];
    			singleobjective[0] = i;
    			IntSet newset = new IntSet(singleobjective, this.os_dim);
    			pSet.add(newset);
    		} 
    		if (relations[i].inrelation(q,p) && !(relations[i].inrelation(p,q))) {
    			int singleobjective[] = new int[1];
    			singleobjective[0] = i;
    			IntSet newset = new IntSet(singleobjective, this.os_dim);
    			qSet.add(newset);
    		}        		
    	}
    	if (pSet.size() > 0 && qSet.size() > 0) {
    		pSet.union_ExactAlgo(qSet);
    		sois = pSet;        		
    	} else {
       		if (pSet.size() == 0 && qSet.size() == 0) {        			
    			/* S_{pq}:= {1,...,k} */
        		int singleobjective[] = new int[1];
        		for (int i=0; i<this.os_dim; i++) {
        			singleobjective[0] = i;
        			IntSet newset = new IntSet(singleobjective, this.os_dim);
        			sois.add(newset);
        		}
       		}
    		if (pSet.size() > 0){
    			sois = pSet;
    		}
    		if (qSet.size() > 0) {
    			sois = qSet;
    		}
    	}
    	        	
    	return sois;
    }
	
}
