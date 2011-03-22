package energyFunction;

import java.util.List;

import javax.vecmath.Point3d;

import dataStructure.ChainTree;

public class Compactness implements EnergyFunction {
	
	private ChainTree cTree;
	
	public Compactness(ChainTree cTree) {
		this.cTree = cTree;
	}
	

	@Override
	public double compute() {
		List<Point3d> points = this.cTree.getBackboneAtomPositions();
		
		return 0;
	}
	

}
