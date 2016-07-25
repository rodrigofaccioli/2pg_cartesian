
package test;

import conflicts.Relation;
import junit.framework.TestCase;

public class TestRelation extends TestCase {

	boolean[][] b1 = {{true, true}, {false, false}};
	Relation r1 = new Relation(System.currentTimeMillis(), b1);
	boolean[][] b2 = {{true, false}, {false, false}};
	Relation r2 = new Relation(System.currentTimeMillis(), b2);
	boolean[][] b3 = {{true, false}, {false, true}};
	Relation r3 = new Relation(System.currentTimeMillis(), b3);
	boolean[][] b4 = {{false, true}, {true, true}};
	Relation r4 = new Relation(System.currentTimeMillis(), b4);
	boolean[][] b5 = {{false, true}, {false, false}};
	Relation r5 = new Relation(System.currentTimeMillis(), b5);
	
	
	/*
	 * Test method for 'conflicts.Relation.inrelation(int, int)'
	 */
	public void testInrelation() {
		assertTrue(r1.inrelation(0,0));
		assertTrue(r1.inrelation(0,1));
		assertFalse(r1.inrelation(1,0));
		assertFalse(r1.inrelation(1,1));
	}

	/*
	 * Test method for 'conflicts.Relation.intersect(long, Relation)'
	 */
	public void testIntersect() {
		assertTrue(r5.equal(r1.intersect(0, r4)));
	}

	/*
	 * Test method for 'conflicts.Relation.minus(long, Relation)'
	 */
	public void testMinus() {
		assertTrue(r5.equal(r1.minus(0, r2)));		
	}

	/*
	 * Test method for 'conflicts.Relation.getComplement()'
	 */
	public void testGetComplement() {
		assertTrue(r2.equal(r4.getComplement()));
		assertFalse(r1.equal(r4.getComplement()));
		assertFalse(r1.equal(r3.getComplement()));
		assertFalse(r1.equal(r2.getComplement()));
		assertFalse(r2.equal(r3.getComplement()));
		assertFalse(r3.equal(r4.getComplement()));
	}

	/*
	 * Test method for 'conflicts.Relation.getNumberOfRelatedPairs()'
	 */
	public void testGetNumberOfRelatedPairs() {
		assertEquals(2, r1.getNumberOfRelatedPairs());
		assertEquals(1, r2.getNumberOfRelatedPairs());
		assertEquals(2, r3.getNumberOfRelatedPairs());
	}

	/*
	 * Test method for 'conflicts.Relation.getExactCopy()'
	 */
	public void testGetExactCopy() {
		Relation r = r1.getExactCopy();
		r.setinrelation(1,1,true);
		assertFalse(r.equal(r1));
	}

	/*
	 * Test method for 'conflicts.Relation.compareTo(Object)'
	 */
	public void testCompareTo() {
		assertEquals(1, r1.compareTo(r2));
		assertEquals(-1, r2.compareTo(r1));
		assertEquals(0, r1.compareTo(r3));
	}

}
