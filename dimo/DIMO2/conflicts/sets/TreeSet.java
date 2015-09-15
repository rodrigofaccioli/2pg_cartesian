package conflicts.sets;

public class TreeSet extends ObjectiveSet {
	int chosenObjective = 0;
	
	public TreeSet() {
		super(0);
		this.delta = 0;
		this.chosenObjective = 0;
		this.numOfElements = 0;
		this.set = null;
	}
	
	public TreeSet(int[] elements, int size, double delta) {
		super(elements, size, delta);
		if (elements.length >= 1) {
			this.chosenObjective = elements[0];
		}
	}
	
	/**
	 * Creates a new TreeSet with all elements from both os1 and os2, delta value
	 * delta and chosen objective chosenObjective
	 */
	public TreeSet(ObjectiveSet os1, ObjectiveSet os2, double delta, int chosenObjective) {
		super(os1.getElements(), delta);
		this.chosenObjective = chosenObjective;
		this.addAll(os2);
	}
	
	public int getChosenObjective() {
		return this.chosenObjective;
	}
	
	public void setChosenObjective(int i) {
		this.chosenObjective = i;
	}
	
	public void setDelta(double d) {
		this.delta = d;
	}

	public String toString() {
		String str = "{";		
		for (int i=0;i<this.set.length;i++) {
			if (this.set[i]) {
				str = str + " " + i;
			}
		}
		str += "} with failure ";
		str += this.delta;
		str += " and chosen objective ";
		str += this.chosenObjective; 
		return str;
	}

	public String toString2PG() {
                /*It is used at 2PG output*/
		String str = "{";		
		for (int i=0;i<this.set.length;i++) {
			if (this.set[i]) {
				str = str + " " + i;
			}
		}
		str += "} ";
		str += this.delta;
		str += " ";
		str += this.chosenObjective; 
		return str;
	}

	
}
