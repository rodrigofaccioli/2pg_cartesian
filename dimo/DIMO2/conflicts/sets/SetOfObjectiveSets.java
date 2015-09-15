package conflicts.sets;

import java.util.Collection;
import java.util.Iterator;
import java.util.Vector;

public class SetOfObjectiveSets {
	Vector<ObjectiveSet> elements;
	
	public SetOfObjectiveSets() {
		elements = new Vector<ObjectiveSet>();	
	}	
	
	public boolean contains(Object arg0) {
		ObjectiveSet os = (ObjectiveSet)arg0;
		int size = this.elements.size();
		for (int i=0; i<size; i++) {
			ObjectiveSet next = (ObjectiveSet)(this.elements.get(i));
			if (next.theSame(os)) {
				return true;				
			}
		}
		return false;
	}
	
	public boolean add(Object arg0) {
		if (this.contains(arg0)) {
			return false;
		} else {
			this.elements.add((ObjectiveSet)arg0);
			return true;
		}		
	}
	
	public boolean remove(Object arg0) {
		ObjectiveSet os = (ObjectiveSet)arg0;
		boolean contained = false;
		Iterator<ObjectiveSet> iter = this.elements.iterator();
		IntSet next = new IntSet(os.set.length);
		while (iter.hasNext()) {
			next = (IntSet)(iter.next());
			if (next.theSame(os)) {
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
	
	public boolean retainAll(Collection<ObjectiveSet> arg0) {
		Vector<ObjectiveSet> newelements = new Vector<ObjectiveSet>();
		for (ObjectiveSet os : arg0) {
			if (this.contains(os)) {
				newelements.add(os);
			}
		}
		this.elements = newelements;

		return true;
	}
	
	/**
	 * Returns true iff all objective sets in arg0 are removed,
	 * i.e., if a set A in arg0 is not in this, returns false.
	 */
	public boolean removeAll(Collection<ObjectiveSet> arg0) {
		boolean ret = true;
		for (ObjectiveSet objectiveSet : arg0) {
			if (!this.remove(objectiveSet)) {
				ret = false;
			}
		}
		return ret;
	}
		
	public String toString() {
		String str ="";
		for (ObjectiveSet objectiveSet : this.elements) {
			str = str + (objectiveSet).toString() + " | " ;
		}
		return str;
	}
	
	/**
	 * @return one IntSet in this SetOfIntSets with a minimal number of elements
	 * the delta value of which is less or equal to delta.
	 */
	public ObjectiveSet getSmallestObjectiveSetWithFailureLessOrEqualTo(double delta) {
		ObjectiveSet smallest = null;
		for (ObjectiveSet currentSet : this.elements) {
			if ((((smallest == null) && (currentSet != null)) && currentSet.getDelta() <= delta) ||
					(smallest != null && (currentSet.size() < smallest.size()) && (currentSet.getDelta() <= delta))) {
				smallest = currentSet;
			}
		}
		if (smallest != null) {
			return smallest;
		} else {
			// return an empty set with indefinite dimension, here for example dimension 1
			return new ObjectiveSet(1);
		}
	}
	
	/**
	 * @return an objective subset of size less or equal k (i.e. a set with <=k objectives)
	 * the delta value of which is the smallest among all sets of size k
	 */
	public ObjectiveSet getBestObjectiveSetOfSizeLessOrEqual(int k) {
		ObjectiveSet best = null;
		for (ObjectiveSet objectiveSet : this.elements) {
			if (best == null) {
				if (objectiveSet.size() <= k) {
					best = objectiveSet;
				}
			} else if ((objectiveSet.getDelta() < best.getDelta()) && (objectiveSet.size() <= k)) {
				best = objectiveSet;
			}
		}
		return best;
	}
	
	public Vector<ObjectiveSet> getElements() {
		return this.elements;
	}
	
	/**  
	 * post: this.elements contains the same ObjectiveSets as elements as a deep copy
	 */
	public void setElements(Vector<ObjectiveSet> elements) {
		this.elements = new Vector<ObjectiveSet>();
		for (ObjectiveSet objectiveSet : elements) {
			this.add(objectiveSet.deepCopy());
		}
	}
	
	
	
	//----------------EVERYTHING COPIED FROM SetOfIntSets.java: -------------
	
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

	public Object[] toArray() {
		return this.elements.toArray();
	}

	public Object[] toArray(Object[] arg0) {
		return this.elements.toArray(arg0);
	}

	public boolean containsAll(Collection<ObjectiveSet> arg0) {
		boolean ret = true;
		for (ObjectiveSet objectiveSet : arg0) {
			if (!(this.contains(objectiveSet))) {
				ret = false;
			}
		}
		return ret;
	}

	public boolean addAll(Collection<ObjectiveSet> arg0) {
		boolean ret = true;
		for (ObjectiveSet objectiveSet : arg0) {
			if (!this.add(objectiveSet)) {
				ret = false;
			}
		}
		return ret;
	}

	public void clear() {
		this.elements.clear();
	}

}
