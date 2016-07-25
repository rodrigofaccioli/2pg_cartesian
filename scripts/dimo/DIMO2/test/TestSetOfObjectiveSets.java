
package test;

import java.util.Vector;
import conflicts.sets.ObjectiveSet;
import conflicts.sets.SetOfObjectiveSets;
import junit.framework.TestCase;

public class TestSetOfObjectiveSets extends TestCase {
	private int[] A = {3, 6, 1, 9};
	private int[] A_prime = {3, 9, 1, 3, 6, 9};
	private int[] B = {2, 6, 1, 4};
	private int[] C = {1, 2, 3};
	private ObjectiveSet mySet1 = new ObjectiveSet(A, 10, 2.1);
	private ObjectiveSet mySet2 = new ObjectiveSet(A_prime, 10, 2.1);
	private ObjectiveSet mySet3 = new ObjectiveSet(B, 10, 2.1);
	private ObjectiveSet mySet4 = new ObjectiveSet(A, 10, 1.9);
	private ObjectiveSet mySet5 = new ObjectiveSet(C, 10, 2.4);
	
	public static void main(String[] args) {
		junit.swingui.TestRunner.run(TestSetOfObjectiveSets.class);
	}

	/*
	 * Test method for 'conflicts.sets.SetOfObjectiveSets.contains(Object)'
	 */
	public void testContains() {
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		soos.add(mySet1);
		soos.add(mySet3);
		assertTrue(soos.contains(mySet1));
		assertTrue(soos.contains(mySet2));
		assertTrue(soos.contains(mySet3));
		assertFalse(soos.contains(mySet4));
	}

	/*
	 * Test method for 'conflicts.sets.SetOfObjectiveSets.add(Object)'
	 */
	public void testAdd() {		
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		assertEquals(0, soos.size());
		assertTrue(soos.add(mySet1));
		assertEquals(1, soos.size());
		assertFalse(soos.add(mySet2));
		assertEquals(1, soos.size());
		assertTrue(soos.add(mySet3));
		assertEquals(2, soos.size());
		assertTrue(soos.add(mySet4));
		assertEquals(3, soos.size());		
	}

	/*
	 * Test method for 'conflicts.sets.SetOfObjectiveSets.remove(Object)'
	 */
	public void testRemove() {
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		assertEquals("size not 0", 0, soos.size());
		assertFalse("removing possible", soos.remove(mySet1));
		assertEquals("size not 0", 0, soos.size());
		assertTrue("adding not possible", soos.add(mySet1));
		assertEquals("size not 1", 1, soos.size());
		assertFalse("removing possible", soos.remove(mySet3));
		assertEquals("size not 1", 1, soos.size());
		assertTrue("adding not possible", soos.add(mySet3));
		assertEquals("size not 3", 2, soos.size());
		assertTrue("removing not possible", soos.remove(mySet2));
		assertEquals("size not 1", 1, soos.size());
	}

	/*
	 * Test method for 'conflicts.sets.SetOfObjectiveSets.retainAll(Collection)'
	 */
	public void testRetainAll() {
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		soos.add(mySet1);		
		soos.add(mySet3);
		Vector<ObjectiveSet> vec = new Vector<ObjectiveSet>();
		vec.add(mySet3);
		vec.add(mySet4);
		assertTrue("retaining caused trouble", soos.retainAll(vec));
		assertEquals("size is not 1", 1, soos.size());
	}

	/*
	 * Test method for 'conflicts.sets.SetOfObjectiveSets.removeAll(Collection)'
	 */
	public void testRemoveAll() {
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		soos.add(mySet1);		
		soos.add(mySet3);
		Vector<ObjectiveSet> vec = new Vector<ObjectiveSet>();
		vec.add(mySet2);
		vec.add(mySet4);
		soos.removeAll(vec);
		assertEquals("size is not 1", 1, soos.size());
	}

	/*
	 * Test method for 'conflicts.sets.SetOfObjectiveSets.getSmallestIntSet()'
	 */
	public void testGetSmallestObjectiveSetWithFailureLessOrEqualTo() {
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		soos.add(mySet2);
		soos.add(mySet3);
		soos.add(mySet4);
		soos.add(mySet5);
		ObjectiveSet smallest2point2 = soos.getSmallestObjectiveSetWithFailureLessOrEqualTo(2.2);
		assertTrue(mySet1.theSame(smallest2point2) ||
				mySet3.theSame(smallest2point2) ||
				mySet4.theSame(smallest2point2));
		assertTrue(mySet4.theSame(soos.getSmallestObjectiveSetWithFailureLessOrEqualTo(1.9)));
		assertEquals(0, soos.getSmallestObjectiveSetWithFailureLessOrEqualTo(0.84567).size());
	}
	
	public void testGetBestObjectiveSetOfSize() {
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		soos.add(mySet2);
		soos.add(mySet3);
		soos.add(mySet4);
		soos.add(mySet5);
		assertTrue("mySet4 not the same as soos.getBestObjectiveSetOfSizeLessOrEqual(5))", mySet4.theSame(soos.getBestObjectiveSetOfSizeLessOrEqual(5)));
		assertTrue("mySet4 not the same as soos.getBestObjectiveSetOfSizeLessOrEqual(4))", mySet4.theSame(soos.getBestObjectiveSetOfSizeLessOrEqual(4)));
		assertTrue("mySet5 not the same as soos.getBestObjectiveSetOfSizeLessOrEqual(3))", mySet5.theSame(soos.getBestObjectiveSetOfSizeLessOrEqual(3)));
		assertEquals("soos.getBestObjectiveSetOfSizeLessOrEqual(2) not equal to 'null'", soos.getBestObjectiveSetOfSizeLessOrEqual(2), null);
		soos.remove(mySet4);
		assertTrue("mySet4 is not the best objective set with k<=5", mySet4.theSame(soos.getBestObjectiveSetOfSizeLessOrEqual(5)));
		
	}
	
	public void testSetElements() {
		Vector<ObjectiveSet> vec = new Vector<ObjectiveSet>();
		vec.add(mySet2);
		vec.add(mySet3);
		vec.add(mySet4);
		SetOfObjectiveSets soos = new SetOfObjectiveSets();
		soos.add(mySet5);
		soos.setElements(vec);
		assertFalse(vec == soos.getElements());
	}

}
