
package conflicts;

import java.util.Arrays;
import java.util.Vector;
import conflicts.sets.IntSet;

public class MOSSGreedyAlgorithm {

	int os_dim;					// the number of objectives
	Relation[] relations;		// an array of the relations $\preceq_i$
	Relation dominanceRelation; // the relation $\preceq$
	
	public MOSSGreedyAlgorithm(int numberOfObjectives, Relation[] rels, Relation rel) {
		this.os_dim = numberOfObjectives;
		this.relations = rels;
		this.dominanceRelation = rel;		
	}
	
	/**
	 * Performs greedy algorithm for MOSS problem as described in bz2007d,
	 * given the population pop and the dominance relations $\preceq_i$
	 * stored in relations.
	 */
	public IntSet performGreedyAlgorithm() {
		Vector<Relation> temprel = new Vector<Relation>(); // temporary relation vector
		
		/* Copy the relations for manipulating them during the algorithm: */
		Vector<Relation> rel = new Vector<Relation>();
		for (int i=0; i<os_dim; i++) {			
			Relation r = this.relations[i].getExactCopy();
			r.setId(i); // relation's ID becomes the number of the corresponding objective
			rel.add(r);			
		}

		/* this is the index which indicates which objectives are in the non-redundant set */
		boolean[] index = new boolean[this.os_dim];
		for (int i=0;i<index.length;i++) {
			index[i] = false;     // initialized with false's
		}

		while (!(rel.isEmpty() || isIntersection(index, this.dominanceRelation, this.relations))) {
			/* Get a relation with minimum number of edges: */
			Relation smallest = getSmallest(rel);
			/* Use first as member in the intersection: */ 
			index[(new Long(smallest.getId())).intValue() ] = true;
			/* Delete first element in rel: */
			rel.remove(smallest);
			/* get all other relations in rel and update them: */
			for (Relation next : rel) {
				update(this.relations[(new Long(smallest.id)).intValue()], next);
				temprel.add(next);
			}
			Relation[] temp = (Relation[])temprel.toArray();
			Arrays.sort(temp);
			/* reset for next iteration: */
			temprel = new Vector<Relation>();
			rel = new Vector<Relation>();
			for (int i=0; i<temp.length;i++) {
				rel.add((Relation)temp[i]);
			}			
		}
		
		return new IntSet(index);
	}
	
	/**
	 * @return a relation from rel with the smallest number of related pairs, i.e.
	 *         with the smallest number of edges in the graph representation
	 *         
	 */
	public Relation getSmallest(Vector<Relation> rel) {
		Relation smallest = null;
		for (Relation next : rel) {
			if (smallest == null) {
				smallest = next;
			} else if (next.numberOfRelatedPairs < smallest.numberOfRelatedPairs) {
				smallest = next;
			}	
		}
		return smallest;
	}
		
	/**
	 * Part of the greedy algorithm: Updates the relations after the ID of deletedRel
	 * is registered in the index, i.e., set all relationships in toUpdate to false
	 * which are false in deletedRel.
	 * 
	 * @param deletedRel the relation which is deleted in the set of objectives
	 * @param toUpdate   the relation we want to update
	 */
	private void update(Relation deletedRel, Relation toUpdate) {
		for (int i=0;i<deletedRel.n;i++) {
			for (int j=0;j<deletedRel.n;j++) {
				if (!deletedRel.inrelation(i,j)) {
					toUpdate.setinrelation(i,j,false);
				}
			}
		}
	}

	/** 
	 * @param index 
	 * @param whole
	 * @param rel array of Relations where the intersection of all rel[i] should be 
	 *            the same than whole
	 * @return true iff the intersection of all rel[i], where index[i] is true, is the
	 *              relation whole
	 */
	private boolean isIntersection(boolean[] index, Relation whole, Relation[] rel) {
		if (rel == null) return (whole==null);
		
		/* initialize intersection as the relation {1,...n} x {1,...,n} */
		boolean[][] a= new boolean[whole.n][whole.n];
		for(int i=0;i<whole.n;i++) {
			for (int j=0;j<whole.n;j++) {
				a[i][j] = true;
			}
		}
		Relation intersection = new Relation(999999, a);
		for (int i=0; i<index.length; i++) {
			if (index[i]) {
				intersection = intersection.intersect(intersection.id, rel[i]);
			}
		}
				
		if (intersection.equal(whole)) {
			return true;
		} else {
			return false;
		}
	}



}
