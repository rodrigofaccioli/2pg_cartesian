
package test;

import java.util.Vector;

import conflicts.sets.IntSet;
import conflicts.sets.SetOfIntSets;
import junit.framework.TestCase;

public class TestSetOfIntSets extends TestCase {
	
	private int[] myints1 = {3, 6, 1, 9};
	private int[] myints2 = {2, 6, 1, 4};
	private int[] myints3 = {1, 6, 9, 3, 3, 7};
	private int[] myints5 = {1, 6, 9, 3, 3, 6};	// same set as myIntSet1!
	private int[] myints6 = {2, 6, 1, 4, 4, 2, 1, 1};
	private int[] myints7 = {1, 6, 9, 3, 7};	// same set as myIntSet1!
	private IntSet myIntSet1 = new IntSet(myints1, 10);
	private IntSet myIntSet2 = new IntSet(myints2, 10);
	private IntSet myIntSet3 = new IntSet(myints3, 10);
	private IntSet myIntSet4 = new IntSet(myints2, 10);	// same set as myIntSet2!
	private IntSet myIntSet5 = new IntSet(myints5, 10);
	private IntSet myIntSet6 = new IntSet(myints6, 10);	
	private IntSet myIntSet7 = new IntSet(myints7, 10); 
	
	public void testAddAllAndSize() {
		SetOfIntSets sois = new SetOfIntSets();
		Vector<IntSet> vec = new Vector<IntSet>();
		vec.add(myIntSet1);
		vec.add(myIntSet2);
		vec.add(myIntSet3);
		vec.add(myIntSet4);				
		sois.addAll(vec);
		assertEquals(3, sois.size());
	}
	
	public void testAddAndSize() {
		SetOfIntSets sois = new SetOfIntSets();
		assertEquals(0, sois.size());
				
		sois.add(myIntSet1);
		assertEquals(1, sois.size());
		
		sois.add(myIntSet2);
		assertEquals(2, sois.size());
				
		sois.add(myIntSet5);
		assertEquals(2, sois.size());
		
		sois.add(myIntSet4);				
		assertEquals(2, sois.size());
	}
	
	public void testGetSmallestIntSet() {
		SetOfIntSets sois = new SetOfIntSets();
		
		sois.add(myIntSet1);
				
		sois.add(myIntSet6);
				
		sois.add(myIntSet7);
		
		IntSet smallest = sois.getSmallestIntSet();
		assertEquals(4, smallest.size());
		
		assertEquals(true, myIntSet1.theSame(smallest)||myIntSet2.theSame(smallest));		
	}
	
	public void testUnion_ExactAlgo() {
		/* initialize first SetOfIntSets */
		SetOfIntSets first = new SetOfIntSets();
		int[] myints1 = {1, 2, 3};
		IntSet myIntSet1 = new IntSet(myints1, 10);
		first.add(myIntSet1);				
		int[] myints2 = {2, 3, 4, 4};
		IntSet myIntSet2 = new IntSet(myints2, 10);
		first.add(myIntSet2);				
		int[] myints3 = {1, 2};	// same set as myIntSet1!
		IntSet myIntSet3 = new IntSet(myints3, 10); 
		first.add(myIntSet3);
		/* initialize second SetOfIntSets */
		SetOfIntSets second = new SetOfIntSets();
		int[] myints1b = {1, 2, 3};
		IntSet myIntSet1b = new IntSet(myints1b, 10);
		second.add(myIntSet1b);				
		int[] myints2b = {5, 6};
		IntSet myIntSet2b = new IntSet(myints2b, 10);
		second.add(myIntSet2b);
		
		/* perform the union (union now in first) */
		first.union_ExactAlgo(second);
		
		/* the sets in the union "first |_| second" */
		int[] ints1 = {1, 2, 3};
		IntSet set1 = new IntSet(ints1, 10);
		int[] ints2 = {2, 3, 4, 5, 6};
		IntSet set2 = new IntSet(ints2, 10);
		int[] ints3 = {1, 2, 5, 6};
		IntSet set3 = new IntSet(ints3, 10);
		
		
		assertEquals(true, first.contains(set1));
		first.remove(set1);
		assertEquals(true, first.contains(set2));
		first.remove(set2);
		assertEquals(true, first.contains(set3));
		first.remove(set3);
		assertEquals(0, first.size());
	}
	
	public void testContains() {
		SetOfIntSets sois = new SetOfIntSets();
		sois.add(myIntSet1);
		sois.add(myIntSet2);
		sois.add(myIntSet3);
		sois.add(myIntSet4);
		sois.add(myIntSet5);
		sois.add(myIntSet6);
		sois.add(myIntSet7);
		int[] myints = {1, 6, 7, 6, 9, 3, 7};	// same set as myIntSet1!
		IntSet myIntSet = new IntSet(myints, 10);
		assertTrue(sois.contains(myIntSet));
		
		int[] newints = {1, 2, 3};	// same set as myIntSet1!
		IntSet newIntSet = new IntSet(newints, 10);
		Vector<IntSet> vec = new Vector<IntSet>();
		vec.add(myIntSet);
		vec.add(newIntSet);
		assertFalse(sois.containsAll(vec));
	}
	
	public void testIsEmpty() {
		SetOfIntSets sois = new SetOfIntSets();
		assertTrue(sois.isEmpty());
		
		sois.add(myIntSet2);
		assertFalse(sois.isEmpty());
		
		sois.remove(myIntSet2);
		assertTrue(sois.isEmpty());
	}
	
	public static void main(String[] args) {
		junit.swingui.TestRunner.run(TestSetOfIntSets.class);
	}

}
