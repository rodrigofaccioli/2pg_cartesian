
package test;

import java.util.Vector;
import conflicts.MOSSGreedyAlgorithm;
import conflicts.Relation;
import junit.framework.TestCase;

public class TestMOSSGreedyAlgorithm extends TestCase {
	
	/*        1 1 1
	 *   A =  0 1 1
	 *        0 0 1  with |A|=6 */
	boolean[][] adjmatrixA = {{true, true, true}, {false, true, true}, {false, false, true}};
	Relation relA = new Relation(1, adjmatrixA);
	/*        1 1 1
	 *   B =  1 1 1
	 *        1 1 1  with |B|=9 */
	boolean[][] adjmatrixB = {{true, true, true}, {true, true, true}, {true, true, true}};
	Relation relB = new Relation(2, adjmatrixB);
	/*        1 1 1
	 *   C =  1 1 1
	 *        0 0 1  with |C|=7 */
	boolean[][] adjmatrixC = {{true, true, true}, {true, true, true}, {false, false, true}};
	Relation relC = new Relation(3, adjmatrixC);
	/*        1 0 0
	 *   D =  0 1 0
	 *        0 0 1  with |D|=3 */
	boolean[][] adjmatrixD = {{true, false, false}, {false, true, false}, {false, false, true}};
	Relation relD = new Relation(4, adjmatrixD);
	/*        1 1 1
	 *   E =  0 1 0
	 *        1 1 1  with |E|=7 */
	boolean[][] adjmatrixE = {{true, true, true}, {false, true, false}, {true, true, true}};
	Relation relE = new Relation(5, adjmatrixE);
	
	Vector<Relation> unsorted = new Vector<Relation>();

	
	protected void setUp() {
		unsorted.add(relA);
		unsorted.add(relB);
		unsorted.add(relC);
		unsorted.add(relD);
		unsorted.add(relE);		  
	  }

	public void testGetSmallest() {
		MOSSGreedyAlgorithm mga = new MOSSGreedyAlgorithm(0,null,null);
		Relation smallest = mga.getSmallest(unsorted);
		assertTrue(smallest.equal(relD));		
	}
	
	public static void main(String[] args) {
		junit.swingui.TestRunner.run(TestMOSSGreedyAlgorithm.class);
	}
	
}
