
package test;

import junit.framework.Test;
import junit.framework.TestSuite;

public class AllTests {

	public static void main(String[] args) {
		junit.swingui.TestRunner.run(AllTests.class);
	}

	public static Test suite() {
		TestSuite suite = new TestSuite("Test for test");
		//$JUnit-BEGIN$
		suite.addTestSuite(TestIntSet.class);
		suite.addTestSuite(TestSetOfIntSets.class);
		suite.addTestSuite(TestMOSSGreedyAlgorithm.class);
		suite.addTestSuite(TestObjectiveSet.class);
		suite.addTestSuite(TestSetOfObjectiveSets.class);
		suite.addTestSuite(TestDeltaMOSSExactAlgorithm.class);
		suite.addTestSuite(TestDeltaMOSSGreedyAlgorithm.class);
		suite.addTestSuite(TestRelation.class);
		//$JUnit-END$
		return suite;
	}

}
