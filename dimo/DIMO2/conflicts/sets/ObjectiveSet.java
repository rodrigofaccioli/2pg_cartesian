package conflicts.sets;

public class ObjectiveSet extends IntSet {
	double delta;
	
	public ObjectiveSet(int size) {
		super(size);
		this.delta = 0;
	}

	public ObjectiveSet(boolean[] elements) {
		super(elements);
		this.delta = 0;
	}

	public ObjectiveSet(IntSet is) {
		super(is);
		this.delta = 0;
	}

	public ObjectiveSet(int[] elements, int size) {
		super(elements, size);
		this.delta = 0;
	}
	
	public ObjectiveSet(boolean[] elements, double delta) {
		super(elements);
		this.delta = delta;
	}
	
	/**
	 * Creates a new Objective set as a subset of
	 * {0,..., size-1} with the elements in elements
	 * and error delta
	 *
	 */
	public ObjectiveSet(int[] elements, int size, double delta) {
		super(elements, size);
		this.delta = delta;
	}

	public double getDelta() {
		return this.delta;
	}	
	
	public boolean theSame(ObjectiveSet is) {
		boolean same = super.theSame(is);
		if (this.delta != is.delta) {
			same = false;
		}
		return same;
	}
	
	/**
	 * @return true iff the ObjectiveSet this is a superset
	 * of os and the delta value is greater or equal os's
	 * delta value.
	 */
	public boolean isSuperSetOf(ObjectiveSet os) {
		if (this.delta < os.delta) {
			return false;	
		} else if (!super.isSuperSetOf(os)) {
			return false;
		}
		return true;
	}
	
	public String toString() {
		String str = "{";		
		for (int i=0;i<this.set.length;i++) {
			if (this.set[i]) {
				str = str + " " + i;
			}
		}
		str += " } with error ";
		str += this.delta;
		return str;
	}
	
	public boolean[] getElements() {
		return this.set;
	}
	
	/**
	 * Returns a deep copy of this ObjectiveSet
	 */
	public ObjectiveSet deepCopy() {
		return new ObjectiveSet(this.set, this.delta);
	}
}
