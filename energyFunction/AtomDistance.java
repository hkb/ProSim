package energyFunction;

import java.util.List;

import javax.vecmath.Point3d;

import dataStructure.ChainTree;

public class AtomDistance implements EnergyFunction {
	
	private ChainTree testing, target;	// the chain trees to compute the energy
	List<Point3d> targetPoints;
	
	/**
	 * 
	 * @param testing
	 * @param target
	 */
	public AtomDistance(ChainTree testing, ChainTree target) {
		this.testing = testing;
		this.target = target;
		
		// precompute target points
		this.targetPoints = this.target.getBackboneAtomPoints();
	}

	@Override
	public double compute() {
		double sum = 0;
		
		List<Point3d> testingPoints = this.testing.getBackboneAtomPoints();
		
		for (int i = 0, j = testingPoints.size(); i < j; i++) {	
			double diff = Math.abs(testingPoints.get(i).distance(this.targetPoints.get(i)));
			
			sum += diff * diff;
		}

		return Math.sqrt(sum / testingPoints.size());
	}

}
