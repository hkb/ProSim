package energyFunction;

import java.util.List;

import javax.vecmath.Point3d;

import math.Point3D;

import dataStructure.ChainTree;

public class LoopAtomDistance implements EnergyFunction {
	
	private ChainTree loop;	// the chain trees to compute the energy
	private int start, end;
	private Point3D[] targetPoints;
	
	/**
	 * 
	 * @param testing
	 * @param target
	 */
	public LoopAtomDistance(ChainTree loop, int start, int end) {
		this.loop = loop;
		this.start = start;
		this.end = end;
		
		// pre-compute target points
		List<Point3D> points = loop.getBackboneAtomPositions(start, end);
		this.targetPoints = new Point3D[points.size()];
		
		int i = 0;
		for(Point3D point : points) {
			this.targetPoints[i] = point;
			
			i++;
		}
	}

	@Override
	public double compute() {
		double sum = 0;
		int i = 0;
		
		for (Point3D point : this.loop.getBackboneAtomPositions(this.start, this.end)) {	
			double diff = point.distance(this.targetPoints[i]);
			
			sum += diff * diff;
			i++;
		}

		return Math.sqrt(sum / this.targetPoints.length);
	}

}
