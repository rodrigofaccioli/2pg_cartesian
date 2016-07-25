
package test;

import java.util.HashSet;
import conflicts.sets.IntSet;
import junit.framework.TestCase;

public class TestIntSet extends TestCase {

	public void testSize1() {
		
		for(int i=0; i<3; i++) {
			System.out.println(i);
		}
		
		
		boolean bool_array[] = {true, false, true, true, false};
		IntSet myIntSet = new IntSet(bool_array);
		assertEquals(3,myIntSet.size());
	}
	
	public void testSize2() {
		IntSet myIntSet = new IntSet(4);
		assertEquals(0, myIntSet.size());
	}
	
	public void testSize3() {
		int[] myints = {3, 6, 1, 9};
		IntSet myIntSet = new IntSet(myints, 10);
		assertEquals(4, myIntSet.size());
	}

	public void testRemove() {
		int[] myints = {3, 6, 1, 9, 17};
		IntSet myIntSet = new IntSet(myints, 20);
		myIntSet.remove(new Integer(9));
		assertEquals(4, myIntSet.size());
		myIntSet.remove(new Integer(9));
		assertEquals(4, myIntSet.size());
	}
	
	public void testAddAndContains() {
		IntSet myIntSet = new IntSet(20);
		myIntSet.add(new Integer(8));
		assertEquals(true, myIntSet.contains(8));
		myIntSet.add(new Integer(9));
		assertEquals(false, myIntSet.contains(10));
		assertEquals(2, myIntSet.size());
	}
	
	public void testTheSame() {
		int[] myints1 = {3, 6, 1, 9};
		IntSet myIntSet1 = new IntSet(myints1, 10);
		int[] myints2 = {2, 6, 1, 4};
		IntSet myIntSet2 = new IntSet(myints2, 10);
		int[] myints3 = {1, 6, 9, 3};
		IntSet myIntSet3 = new IntSet(myints3, 10);
		assertEquals(false, myIntSet1.theSame(myIntSet2));
		assertEquals(true, myIntSet1.theSame(myIntSet3));		
	}
	
	public void testIsSuperSetOf() {
		int[] myints1 = {1, 7, 2, 6, 1, 4};
		IntSet myIntSet1 = new IntSet(myints1, 10);
		int[] myints2 = {2, 6, 4, 1};
		IntSet myIntSet2 = new IntSet(myints2, 10);
		IntSet myIntSet3 = new IntSet(myints1, 20);
		assertEquals(true, myIntSet1.isSuperSetOf(myIntSet2));
		assertEquals(false, myIntSet2.isSuperSetOf(myIntSet1));
		assertEquals(false, myIntSet1.isSuperSetOf(myIntSet3));
	}
	
	public void testClear() {
		int[] myints1 = {1, 7, 2, 6, 1, 4};
		IntSet myIntSet1 = new IntSet(myints1, 10);
		myIntSet1.clear();
		assertEquals(0, myIntSet1.size());
		assertEquals(true, myIntSet1.theSame(new IntSet(10)));
	}
	
	public void testRetainAll() {
		int[] myints1 = {1, 7, 2, 6, 1, 4};
		IntSet myIntSet1 = new IntSet(myints1, 10);
		int[] myints2 = {2, 8, 6, 3, 5, 1};
		IntSet myIntSet2 = new IntSet(myints2, 10);
		int[] myints = {1, 2, 6};
		IntSet intersection = new IntSet(myints, 10);
		assertEquals(true, myIntSet1.retainAll(myIntSet2));
		assertEquals(true, intersection.theSame(myIntSet1));
		HashSet<Integer> myhashset = new HashSet<Integer>();
		myhashset.add(new Integer(1));
		myhashset.add(new Integer(2));
		myhashset.add(new Integer(6));
		assertEquals(false, myIntSet2.retainAll(myhashset));
		assertEquals(false, intersection.theSame(myIntSet2));
		
	}
	
	public static void main(String[] args) {
		junit.swingui.TestRunner.run(TestIntSet.class);
	}
}
