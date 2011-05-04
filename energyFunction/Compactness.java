package energyFunction;

import java.util.List;

import javax.vecmath.Point3d;

import math.Point3D;

import dataStructure.ChainTree;

public class Compactness implements EnergyFunction {
	
	private ChainTree cTree;
	
	public Compactness(ChainTree cTree) {
		this.cTree = cTree;
	}
	

	@Override
	public double compute() {
		List<Point3D> points = this.cTree.getBackboneAtomPositions();
		
		return 0;
	}
	

}
