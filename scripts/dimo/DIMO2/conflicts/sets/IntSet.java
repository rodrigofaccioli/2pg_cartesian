package conflicts.sets;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class IntSet implements Set<Integer> {
	boolean[] set;          // if set[i] == true, ith element is contained in set
	int numOfElements = 0;  // current number of elements in the set

	public IntSet(int size) {
		this.set = new boolean[size];
	}
	
	public IntSet(boolean[] elements) {
		this.set = new boolean[elements.length];
		for (int i=0; i<elements.length; i++) {
			if (elements[i]) {
				this.set[i] = true;
				this.numOfElements++;
			}
		}
	}

	public IntSet(IntSet is) {
		this.set = new boolean[is.set.length];
		for (int i=0; i<is.set.length; i++) {
			this.set[i] = is.set[i];
			if (is.set[i]) {
				this.numOfElements++;
			}
		}
	}

	/**
	 * @param elements the elements in the set
	 * @param size the elements in the set will be inside {0, size-1}
	 */
	public IntSet(int[] elements, int size) {
		this.set = new boolean[size];
		for (int i=0; i<elements.length; i++) {
			if (!this.set[elements[i]]) { // handle multiple entries in elements
				this.set[elements[i]] = true;
				this.numOfElements++;
			}
		}
	}

	public int size() {
		return this.numOfElements;
	}

	public boolean isEmpty() {
		return (this.numOfElements == 0);
	}

	public boolean contains(Object arg0) {
		// No error handling here!!!
		return this.set[((Integer)arg0).intValue()];
	}

	public boolean contains(int i) {
		return contains((new Integer(i)));
	}

	/**
	 * Returns an Iterator of type Integer, containing all integers
	 * in this set.
	 */
	public Iterator<Integer> iterator() {
		Vector<Integer> vec = new Vector<Integer>();
		int pos = 0;
		for (int i=0;i<this.set.length;i++) {
			if (this.set[i]) {
				vec.add(pos, (new Integer(i)));
				pos++;
			}
		}
		return vec.iterator();
	}

	/**
	 * Returns always an int[]-array of size this.size() with the elements of this.set
	 */
	public Object[] toArray() {
		int[] array = new int[this.size()];
		int pos = 0;
		for (int i=0;i<this.set.length;i++) {
			if (this.set[i]) {
				array[pos] = i;
				pos++;
			}
		}
		return null;
	}

	public <T> T[] toArray(T[] args) {
		// TODO Not supported yet
		return null;
	}

	public boolean add(Integer arg0) {
		int integer = arg0.intValue();
		boolean noEntry = !this.set[integer];
		if (noEntry) {
			this.numOfElements++;
			this.set[integer] = true;
		}
		return noEntry;
	}

	public boolean remove(Object arg0) {
		int integer = ((Integer)arg0).intValue();
		boolean entry = this.set[integer];
		if (entry) {
			this.numOfElements--;
			this.set[integer] = false;
		}
		return entry;
	}

	public boolean containsAll(Collection<?> arg0) {
		boolean all = true;
		Iterator<?> iter = arg0.iterator();
		while (iter.hasNext()) {
			if (!this.contains((Integer)iter.next())) {
				all = false;
				break;
			}
		}
		return all;
	}

	public boolean addAll(Collection<? extends Integer> arg0) {
		boolean b = true;
		Iterator<?> iter = arg0.iterator();
		while (iter.hasNext()) {
			if (!this.add((Integer)iter.next())) {
				b = false;
			}
		}
		return b;
	}

	/**
	 * pre: arg0 has to be an IntSet
	 * 
	 * post: this IntSet is the intersection between this IntSet and IntSet arg0
	 *
	 * @return false if problems during intersection occured and else true
	 */
	public boolean retainAll(Collection<?> arg0) {
		if (!arg0.getClass().isInstance(new IntSet(1))) {
			return false;
		}
		boolean[] intersection = new boolean[this.set.length];
		int newNumOfElements = 0;
		for (Object object : arg0) {
			int i = ((Integer)object).intValue();
			if (this.set[i]) {
				intersection[i] = true;
				newNumOfElements++;
			} else {
				intersection[i] = false;
			}
		}
		this.set = intersection;
		this.numOfElements = newNumOfElements;
		return true;
	}

	public boolean removeAll(Collection<?> arg0) {
		boolean b = true;
		for (Object object : arg0) {
			if (!this.remove((Integer)object)) {
				b = false;
			}
		}
		return b;
	}

	public void clear() {
		for (int i=0;i<this.set.length;i++) {
			this.set[i] = false;
		}
		this.numOfElements = 0;
	}
	
	public String toString() {
		String str = "";
		for (int i=0;i<this.set.length;i++) {
			if (this.set[i]) {
				str = str + " " + i;
			}
		}
		return str;
	}

	/**
	 * @param is
	 * @return true iff is and this IntSet contain the same elements
	 */
	public boolean theSame(IntSet is) {
		if (is.set.length != this.set.length) {
			return false;
		}
		boolean ret = true;
		int l = this.set.length;
		for (int i=0; i<l; i++) {
			if (this.set[i] != is.set[i]) {
				ret = false;
				break;
			}
		}
		return ret;
	}
	
	/**
	 * @param is another IntSet (has to be a subset of the same universal set {1,...,this.set.length)} 
	 * @return true iff this IntSet is a superset of is
	 */
	public boolean isSuperSetOf(IntSet is) {
		int l = is.set.length;
		/* to be a superset of is, 
		 * is has to have the same universal set as this.set! */
		if (l != this.set.length) {
			return false;
		}
		/* to be a superset of is, 
		 * this.set has to have at least as many elements than is! */
		if (this.numOfElements < is.numOfElements) {
			return false;
		}
		for (int i=0; i< l; i++) {
			if (is.set[i]) {
				if (!(this.set[i])) {
					return false;
				}
			}
		}
		return true;
	}
	
	public boolean[] getComplement() {
		boolean[] ret = new boolean[this.set.length];
		for (int i=0; i<this.set.length; i++) {
			ret[i] = !(this.set[i]);
		}		
		return ret;
	}

}
