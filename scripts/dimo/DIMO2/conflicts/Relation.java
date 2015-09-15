
package conflicts;

public class Relation implements Comparable<Relation>{
	
	long id;
	int n;
	boolean rel[][];
	int numberOfRelatedPairs; // the number of edges in the corresponding relation graph
	
	public Relation() {
		id = System.currentTimeMillis();
		n=1;
		rel = new boolean[n][n];
		numberOfRelatedPairs = 0;
	}

	public Relation(int n) {
		id = System.currentTimeMillis();
		this.n=n;
		rel = new boolean[n][n];
		numberOfRelatedPairs = 0;
	}
	
	public Relation(long id, int n) {
		this.id = id;
		this.n=n;
		rel = new boolean[n][n];
		numberOfRelatedPairs = 0;
	}
	
	
	/**
	 * Constructor 
	 * 
	 * pre: 2D array rel should be quadratic
	 */
	public Relation(long id, boolean[][] rel) {
		this.id = id;
		n = rel.length;
		numberOfRelatedPairs = 0;
		this.rel = new boolean[n][n];
		for (int i=0; i<n; i++) {
			for (int j=0;j<n;j++) {
				this.rel[i][j] = rel[i][j];
				if (rel[i][j]) {
					numberOfRelatedPairs++;
				}
			}			
		}
	}

	/**
	 * @param i element i of the set {1,...n}
	 * @param j element j of the set {1,...n}
	 * @return true if i is in relation with j / 
	 *                 iRj / 
	 *                 (i,j)\in R
	 *         false if i is not in relation with j
	 *           and if i or j is not in {1,...,n}
	 */
	public boolean inrelation(int i, int j) {
		if (i<n && i >= 0 && j<n && j>=0) {
			return rel[i][j]; 
		} else
			return false;
	}
	
		
	/**
	 * @param i index of first element 
	 * @param j index of second element
	 * @param b 
	 * 
	 * Sets i and j in relation iff b == true.
	 */
	public void setinrelation(int i, int j, boolean b) {
		if (i>=0 && i<= n && j>=0 && j<=n) {
			if (rel[i][j] == false && b) {
				this.numberOfRelatedPairs++;
			} else if (rel[i][j] == true && !b) {				
				this.numberOfRelatedPairs--;
			}		
			rel[i][j] = b;			
		} else {
			System.out.println("ERROR while setting relation: index out of bounds");
		}
			
	}
	
	
	/**
	 * @param id the ID of the resulting new relation
	 * @param r a relation on set {1,...,n} x {1,...,n} 
	 * @return the intersection of this relation and relation r
	 */
	public Relation intersect(long id, Relation r) {
		if (r.n != this.n) {
			return null;						
		}
		Relation s = new Relation(id, rel);
		for (int i=0;i<n;i++) {
			for (int j=0; j<n; j++) {
				s.setinrelation(i,j, (rel[i][j] & r.rel[i][j])); 
			}
		}		
	   return s;
	}
	
	/**
	 * Returns a new relation with ID id, containing all edges from 'this', except
	 * the edges in r.
	 */
	public Relation minus(long id, Relation r){
		Relation result = this.getExactCopy();
		for (int i=0; i<this.n; i++) {
			for (int j=0; j<this.n; j++) {
				if (r.inrelation(i,j)) {
					result.setinrelation(i,j, false);
				}
			}
		}
		return result;
	}
	
	
	/**
	 * @return the complementary relation
	 */
	public Relation getComplement() {
		Relation newRel = new Relation(this.n);
		for (int i=0; i<this.n; i++) {
			for (int j=0; j<this.n;j++) {
				newRel.setinrelation(i,j, !this.rel[i][j]);				
			}			
		}
		return newRel;
	}
	
	/* returns the number of pairs, which stand in relation*/
	public int getNumberOfRelatedPairs() {
		return numberOfRelatedPairs;
	}
	
	
	/**
	 * @param r a Relation
	 * @return true iff the relation r and this relation are the same
	 */
	public boolean equal(Relation r) {
		boolean equal = true;
		if (this.n != r.n) {
			return false;
		}
		for (int i=0; i<this.n; i++) {
			for(int j=0; j<this.n; j++) {
				if (this.rel[i][j] != r.inrelation(i,j)) {
					return false;
				}
			}
		}
		return equal;
	}

	/**
	 * Returns an exact (deep) copy of 'this'.
	 */
	public Relation getExactCopy() {
		Relation r = new Relation(this.id, this.n);
		for (int i=0;i<n;i++) {
			for(int j=0;j<n;j++) {
				r.setinrelation(i,j,this.rel[i][j]);
			}
		}
		return r;
	}
	
	public void print() {
		System.out.println("ID= " + id + " of " + this);
		for (int i=0;i<n;i++) {
			for (int j=0;j<n;j++) {
				if (this.inrelation(i,j)) {
					System.out.print("X");
				} else
					System.out.print("-");
			}
			System.out.println();
		}
	}
	
	
	/** 
	 * Compares this relation with 'o'. A relation is smaller than another
	 * if it has less true-entries than the other relation.
	 */
	public int compareTo(Relation o) {
		if (this.numberOfRelatedPairs < o.numberOfRelatedPairs) {
			return -1;
		} else if (this.numberOfRelatedPairs == o.numberOfRelatedPairs) {
			return 0;
		} else {
			return 1;
		}		
	}

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}

}
