package energyFunction;

import java.util.List;

import javax.vecmath.Point3d;

import dataStructure.ChainTree;

public class EndAtomDistance implements EnergyFunction {
	
	private ChainTree testing, target;	// the chain trees to compute the energy
	public Point3d targetPoint;
	
	/**
	 * 
	 * @param testing
	 * @param target
	 */
	public EndAtomDistance(ChainTree testing, ChainTree target) {
		this.testing = testing;
		this.target = target;
		
		// precompute target points
		this.targetPoint = this.target.getBackboneAtomPositions().get(this.target.length());
	}

	@Override
	public double compute() {
		return this.testing.getBackboneAtomPositions().get(this.testing.length()).distance(this.targetPoint);
	}

}
