package conflicts.sets;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class SetOfIntSets implements Set<IntSet> {

	Vector<IntSet> elements;
	
	public SetOfIntSets() {
		elements = new Vector<IntSet>();
	}
	
	public int size() {
		if (this.elements != null) {
			return this.elements.size();
		} else {
			return 0;
		}
	}

	public boolean isEmpty() {
		return this.elements.isEmpty();
	}

	public boolean contains(Object arg0) {
		boolean contained = false;
		for (IntSet is : this.elements) {
			if (is.theSame((IntSet)arg0)) {
				contained = true;
				break;
			}
		}
		return contained;
	}

	public Iterator<IntSet> iterator() {
		return this.elements.iterator();
	}

	public Object[] toArray() {
		return this.elements.toArray();
	}

	public <T> T[] toArray(T[] arg0) {
		return this.elements.toArray(arg0);
	}

	public boolean add(IntSet arg0) {
		boolean contained = false;
		for (IntSet is : this.elements) {
			if (is.theSame(arg0)) {
				contained = true;
				break;
			}
		}
		if (contained) {
			return false;
		} else {
			this.elements.add(arg0);
			return true;
		}		
	}

	public boolean remove(Object arg0) {
		IntSet is = (IntSet)arg0;
		boolean contained = false;
		Iterator<IntSet> iter = this.elements.iterator();
		IntSet next = new IntSet(is.set.length);
		while (iter.hasNext()) {
			next = (IntSet)(iter.next());
			if (next.theSame(is)) {
				contained = true;
				break;
			}
		}
		if (contained) {
			this.elements.remove(next);
			return true;
		} else {			
			return false;
		}
	}

	public boolean containsAll(Collection<?> arg0) {
		boolean ret = true;
		Iterator<?> iter = arg0.iterator();
		while (iter.hasNext()) {
			if (!(this.contains(iter.next()))) {
				ret = false;
			}
		}
		return ret;
	}

	public boolean addAll(Collection<? extends IntSet> arg0) {
		boolean ret = true;
		for (IntSet intSet : arg0) {
			if (!this.add(intSet)) {
				ret = false;
			}
		}
		return ret;
	}

	public boolean retainAll(Collection<?> arg0) {
		Vector<IntSet> newelements = new Vector<IntSet>();
		Iterator<?> iter = arg0.iterator();
		while (iter.hasNext()) {
			IntSet is = (IntSet)(iter.next());
			if (this.contains(is)) {
				newelements.add(is);
			}
		}
		this.elements = newelements;

		return true;
	}

	public boolean removeAll(Collection<?> arg0) {
		boolean ret = true;
		Iterator<?> iter = arg0.iterator();
		while (iter.hasNext()) {
			if (!this.remove(((IntSet)iter.next()))) {
				ret = false;
			}
		}
		return ret;
	}

	public void clear() {
		this.elements.clear();
	}
	
	public String toString() {
		String str ="";
		for (IntSet is : this.elements) {
			str = str + is.toString() + " | " ;
		}
		return str;
	}
	
	/**
	 * @return one IntSet in this SetOfIntSets with a minimal number of elements
	 */
	public IntSet getSmallestIntSet() {
		IntSet smallest = null;	
		for (IntSet is : this.elements) {
			if (smallest == null || is.size() < smallest.size()) {
				smallest = is;
			}
		}		
		if (smallest != null) {
			return smallest;
		} else {
			/* return an empty set with indefinite dimension, here for example dimension 1 */
			return new IntSet(1);
		}
	}
	
	
	/**
	 * This method is used during the exact algorithm for computing a minimal non-redundant
	 * set of objectives, not conflicting with the whole set of all objectives.
	 * 
	 * This SetOfIntSets is updated during this method with the information from currentSet.
	 * 
	 * post: The new list of possible sets in this.elements is the list we get if we union each
	 *       set of this.elements with each set in currentSet. After this computation all redundant
	 *       sets are deleted from this.elements.
	 * 
	 * @param currentSet
	 */
	public void union_ExactAlgo(SetOfIntSets currentSet) {		
		Vector<IntSet> newList = new Vector<IntSet>();
		/* nested while loop runs faster if smaller set is considered in inner loop: */
		SetOfIntSets smallset;
		SetOfIntSets largeset;
		if (this.size() > currentSet.size()) {
			smallset = currentSet;
			largeset = this;			
		} else {
			smallset = this;
			largeset = currentSet;			
		}
		if (smallset.isEmpty()) { 
			this.elements = largeset.elements; // union = largest set if smallest set empty
		} else {
			Iterator<IntSet> largeset_iter = largeset.elements.iterator();
			while (largeset_iter.hasNext()) {
				Iterator<IntSet> smallset_iter = smallset.elements.iterator();
				IntSet set1 = largeset_iter.next();
				while (smallset_iter.hasNext()) {
					IntSet set2 = smallset_iter.next();
					IntSet unionset = union(set1,set2);					
					addIfSubSet(unionset, newList);
				}
			}
			this.elements = newList;
		}
	}
	
	/**
	 * If IntSet is a subset of a set in list:
	 *    - include IntSet is into list
	 *    - delete all supersets of IntSet is from list
	 * If there is no subset of IntSet is in newList
	 *    - add IntSet is to list
	 *  
	 * @param is IntSet
	 * @param list a vector of IntSets 
	 */
	private void addIfSubSet(IntSet is, Vector<IntSet> list) {
		/* insert = true iff is has to be inserted into list: */
		boolean insert = true;
		/* a list of elements, we have to delete in list: */
		Vector<IntSet> delete = new Vector<IntSet>();
		
		Iterator<IntSet> iter = list.iterator();
		while (iter.hasNext()) {
			IntSet currentIs = (IntSet)(iter.next());
			if (currentIs.isSuperSetOf(is)) {
				delete.add(currentIs);
				insert = true;
			} else {
				if (is.isSuperSetOf(currentIs)) {
					insert = false;
				}
			}
		}
		if (insert) {
			/* delete elements in list */
			list.removeAll(delete);
			list.add(is);
		}
	}
	
	private IntSet union(IntSet s1, IntSet s2) {
		IntSet union;
		if (s1.size() > s2.size()) {
			union = new IntSet(s1); 
			union.addAll(s2);
		} else {
			union = new IntSet(s2); 
			union.addAll(s1);
		}			
		return union;
	}

}
