package algorithm;

import java.util.List;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import math.Point3D;
import math.Tuple2;
import math.Vector3D;

import dataStructure.ChainTree;
import edu.math.Line;
import edu.math.Vector;

/**
 * An implementation of the cyclic coordinate descent (CCD) algorithm
 * for choosing the rotation of bonds during loop closure.
 * 
 * @author hkb
 *
 */
public class CyclicCoordinateDescent {
	
	public static double TARGET_RMS_DISTANCE = 0.08;
	public static int TARGET_LENGTH = 3;
	public static int MAX_ITERATIONS = 20;
	
	private ChainTree loop;						// the loop to close including the target anchor
	private Vector3D[] target;					// the position of the target

	
	/**
	 * Creates a new CCD computer from the unfolded loop and the target
	 * residue.
	 * 
	 * @param loop The loop to close including the target residue.
	 * @param target The target residue.
	 */
	public CyclicCoordinateDescent(ChainTree loop, ChainTree target) {
		// fetch target position
		if (target.length() != TARGET_LENGTH)
			throw new IllegalArgumentException("Target must be a single residue!");
		
		this.target = new Vector3D[3];
		
		int i = 0;
		for (Point3D position : target.getBackboneAtomPositions()) {
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
		// determine the rotation axis of the bond
		List<Point3D> bondAtoms = this.loop.getBackboneAtomPositions(bond, bond);
		Vector3D rotationAxis = bondAtoms.get(1).subtract(bondAtoms.get(0));
		
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length()-2, this.loop.length()-1);
		Vector3D[] moving = new Vector3D[TARGET_LENGTH];
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			moving[i] = movingTerminalAtoms.get(i);
		}
		
		// compute the vectors r, t, n
		Vector3D[] r = new Vector3D[TARGET_LENGTH];
		Vector3D[] f = new Vector3D[TARGET_LENGTH];
		Vector3D[] s = new Vector3D[TARGET_LENGTH];
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			Vector3D M = moving[i];
			Vector3D T = this.target[i];
			Vector3D O = M.norm().scale(rotationAxis.dot(M.norm()));
			
			r[i] = M.subtract(O);
			f[i] = T.subtract(O);
			s[i] = O.norm().cross(r[i].norm());
		}
		
		/*
		 * Compute alpha.
		 */
		double numerator = 0;
		double denominator = 0;
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			numerator += f[i].dot(s[i].norm()) * r[i].length();
			denominator += f[i].dot(r[i].norm()) * r[i].length();
		}
		
		/*
		 * Compute second derivative.
		 */
		//Math.atan2(arg0, arg1)
		
		return Math.atan(numerator / denominator);
	}
	
	/**
	 * Determines if the loop is closed.
	 */
	public boolean isLoopClosed() {
		return false;
	}
	
	/**
	 * Converts a point to a vector.
	 * 
	 * @param point
	 * @return
	 */
	private static Vector3D pointToVector(Point3d point) {
		return new Vector3D(point.x, point.y, point.z);
	}
}
