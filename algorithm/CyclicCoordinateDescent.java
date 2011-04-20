package algorithm;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import tool.Tuple2;
import dataStructure.ChainTree;

/**
 * An implementation of the cyclic coordinate descent (CCD) algorithm
 * for choosing the rotation of bonds during loop closure.
 * 
 * @author hkb
 *
 */
public class CyclicCoordinateDescent {
	
	public static double TARGET_RMS_DISTANCE = 0.08; 
	public static int MAX_ITERATIONS = 5000;
	
	private ChainTree loop;						// the loop to close including the target anchor
	private Point3d[] target;					// the position of the target

	
	/**
	 * Creates a new CCD computer from the unfolded loop and the target
	 * residue.
	 * 
	 * @param loop The loop to close including the target residue.
	 * @param target The target residue.
	 */
	public CyclicCoordinateDescent(ChainTree loop, ChainTree target) {
		// fetch target position
		if (target.length() != 3)
			throw new IllegalArgumentException("Target must be a single residue!");
		
		this.target = new Point3d[3];
		
		int i = 0;
		for (Point3d position : target.getBackboneAtomPositions()) {
			this.target[i] = position;
			
			i++;
		}
		
		// store loop 
		this.loop = loop;
	}
	
	/**
	 * Returns the rotation around the given bond that brings the loop terminal
	 * residue closer to the target.
	 * 
	 * @return The rotation angle that minimises the loop terminals distance to the target.
	 */
	public double getRotationAngle(int bond) {
		/*
		 * Compute vectors.
		 */
		Vector3d O = new Vector3d(1,2,3);
		
		
		/*
		 * Compute alpha.
		 */
		int alpha = 0;
		
		/*
		 * Compute second derivative.
		 */
		//Math.atan2(arg0, arg1)
		
		return alpha;
	}
	
	/**
	 * Determines if the loop is closed.
	 */
	public boolean isLoopClosed() {
		return false;
	}
}
