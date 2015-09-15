
package test;

import conflicts.sets.ObjectiveSet;
import junit.framework.TestCase;

public class TestObjectiveSet extends TestCase {

	public static void main(String[] args) {
		junit.swingui.TestRunner.run(TestObjectiveSet.class);
	}

	/*
	 * Test method for 'conflicts.sets.ObjectiveSet.theSame(ObjectiveSet)'
	 */
	public void testTheSameObjectiveSet() {
		int[] A = {3, 6, 1, 9};
		int[] A_prime = {3, 9, 1, 3, 6, 9};
		int[] B = {2, 6, 1, 4};
		ObjectiveSet mySet1 = new ObjectiveSet(A, 10, 2.1);
		ObjectiveSet mySet2 = new ObjectiveSet(A_prime, 10, 2.1);
		ObjectiveSet mySet3 = new ObjectiveSet(B, 10, 2.1);
		ObjectiveSet mySet4 = new ObjectiveSet(A, 18, 2.1);
		ObjectiveSet mySet5 = new ObjectiveSet(A, 10, 1.9);
		assertTrue("1 vs. 2", mySet1.theSame(mySet2));
		assertFalse("1 vs. 3", mySet1.theSame(mySet3));
		assertFalse("1 vs. 4", mySet1.theSame(mySet4));
		assertFalse("1 vs. 5", mySet1.theSame(mySet5));				
	}

	/*
	 * Test method for 'conflicts.sets.ObjectiveSet.isSuperSetOf(ObjectiveSet)'
	 */
	public void testIsSuperSetOfObjectiveSet() {
		int[] A = {3, 6, 1, 9, 17, 23, 16, 4, 2};
		int[] B = {2, 6, 1, 4, 6, 9, 2, 1};
		int[] B2 = {1, 2, 4, 6, 9};
		int[] D = {3, 6, 5, 10};		
		ObjectiveSet one = new ObjectiveSet(A, 25, 2.1);
		ObjectiveSet two = new ObjectiveSet(B, 25, 2.1);
		ObjectiveSet three = new ObjectiveSet(B2, 25, 2.1);
		ObjectiveSet one_low = new ObjectiveSet(A, 25, 0.8);
		ObjectiveSet two_high = new ObjectiveSet(B, 25, 2.100000001);
		ObjectiveSet four = new ObjectiveSet(D, 25, 2.1);
		
		assertTrue(one.isSuperSetOf(one));
		
		assertTrue(one.isSuperSetOf(two));
		assertFalse(two.isSuperSetOf(one));

		assertTrue(two.isSuperSetOf(three));
		assertTrue(three.isSuperSetOf(two));
		
		assertTrue(one.isSuperSetOf(one_low));
		assertFalse(one_low.isSuperSetOf(one));
		
		assertFalse(one.isSuperSetOf(two_high));
		assertFalse(two_high.isSuperSetOf(one));
		
		assertFalse(one.isSuperSetOf(four));
		assertFalse(four.isSuperSetOf(one));
		
	}
	
	public void testDeepCopy() {
		int[] A = {3, 6, 1, 9, 17, 23, 16, 4, 2};
		ObjectiveSet os = new ObjectiveSet(A, 25, 2.1);
		ObjectiveSet deepCopy = os.deepCopy();
		assertTrue(os.theSame(deepCopy));
		assertTrue(deepCopy.theSame(os));
		os.add(new Integer(4));
		assertTrue(os.theSame(deepCopy));
		assertTrue(deepCopy.theSame(os));
		os.add(new Integer(5));
		assertFalse(os.theSame(deepCopy));
		assertFalse(deepCopy.theSame(os));		
	}

}
