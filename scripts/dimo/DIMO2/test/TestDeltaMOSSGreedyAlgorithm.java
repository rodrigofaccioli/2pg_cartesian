
package test;

import conflicts.DeltaMOSSGreedyAlgorithm;
import conflicts.FileProblem;
import conflicts.sets.ObjectiveSet;
import junit.framework.TestCase;

public class TestDeltaMOSSGreedyAlgorithm extends TestCase {

	public static void main(String[] args) {
		junit.swingui.TestRunner.run(TestDeltaMOSSGreedyAlgorithm.class);
	}

	/*
	 * Test method for 'conflicts.DeltaMOSSGreedyAlgorithm.performGreedyAlgorithmGivenK(int)'
	 */
	public void testPerformGreedyAlgorithmGivenK() {
		/* preparing input for constructor of DeltaMOSSGreedyAlgorithm */
		FileProblem fp = new FileProblem("test/testExactAlgo1.txt");
		double[][] points = fp.getPoints();
		
		DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(points);		
		
		int[] set1 = {1};
		ObjectiveSet os1 = new ObjectiveSet(set1, 4, 2);
		int[] set2 = {0, 1};
		ObjectiveSet os2 = new ObjectiveSet(set2, 4, 0);
		assertTrue("givenK=1 failed,", os1.theSame(dmga.performGreedyAlgorithmGivenK(1)));
		assertTrue("givenK=2 failed,", os2.theSame(dmga.performGreedyAlgorithmGivenK(2)));
		assertTrue("givenK=3 failed,", os2.theSame(dmga.performGreedyAlgorithmGivenK(3)));
		assertTrue("givenK=4 failed,", os2.theSame(dmga.performGreedyAlgorithmGivenK(4)));
	}

	/*
	 * Test method for 'conflicts.DeltaMOSSGreedyAlgorithm.performGreedyAlgorithmGivenDelta(double)'
	 */
	public void testPerformGreedyAlgorithmGivenDelta() {
		// test for example with delta = 0.2
		boolean[] b = {true, true, false, true};
		ObjectiveSet expectedresult = new ObjectiveSet(b, 0.0); 
		
		double[][] ov = {{0.5, 0.1, 0.9, 0.6}, {0.7, 0.2, 0.5, 0.4}, {0.4, 0.7, 0.5, 0.2}, {0.3, 0.2, 0.5, 0.8}};
		DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(ov);
		
		ObjectiveSet result = dmga.performGreedyAlgorithmGivenDelta(0.2);
		assertTrue(expectedresult.theSame(result));
		
		// test same example with delta = 0.5
		boolean[] c = {false, true, false, true};
		ObjectiveSet expectedresult2 = new ObjectiveSet(c, 0.4); 
		ObjectiveSet result2 = dmga.performGreedyAlgorithmGivenDelta(0.5);
		assertTrue(expectedresult2.theSame(result2));
		
	}

	/*
	 * Test method for 'conflicts.DeltaMOSSGreedyAlgorithm.performGreedyAlgorithm()'
	 */
	public void testPerformGreedyAlgorithm() {
		boolean[] b = {true, true, false, true};
		ObjectiveSet expectedresult = new ObjectiveSet(b, 0.0); 
		
		double[][] ov = {{0.5, 0.1, 0.9, 0.6}, {0.7, 0.2, 0.5, 0.4}, {0.4, 0.7, 0.5, 0.2}, {0.3, 0.2, 0.5, 0.8}};
		DeltaMOSSGreedyAlgorithm dmga = new DeltaMOSSGreedyAlgorithm(ov);
		
		ObjectiveSet result = dmga.performGreedyAlgorithmGivenDelta(0.2);
		assertTrue(expectedresult.theSame(result));
	}

}
