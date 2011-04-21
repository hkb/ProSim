package algorithm;

import java.util.List;

import javax.vecmath.Point3d;

import tool.Tuple2;
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
	public static int MAX_ITERATIONS = 5000;
	
	private ChainTree loop;						// the loop to close including the target anchor
	private Vector[] target;					// the position of the target

	
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
		
		this.target = new Vector[3];
		
		int i = 0;
		for (Point3d position : target.getBackboneAtomPositions()) {
			this.target[i] = pointToVector(position);
			
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
		List<Point3d> bondAtoms = this.loop.getBackboneAtomPositions(bond, bond);
		Line rotationAxis = new Line(pointToVector(bondAtoms.get(0)), pointToVector(bondAtoms.get(1)));
		
		// fetch the positions of the moving terminal residue
		List<Point3d> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length()-2, this.loop.length()-1);
		Vector[] moving = new Vector[3];
		
		for (int i = 0; i < 3; i++) {
			moving[i] = pointToVector(movingTerminalAtoms.get(i));
		}
		
		// compute the vectors r, t, n
		Vector[] r = new Vector[3];
		Vector[] t = new Vector[3];
		Vector[] n = new Vector[3];
		
		for (int i = 0; i < 3; i++) {
			Vector M = moving[i];
			Vector T = this.target[i];
			Vector O = rotationAxis.
		}
		
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
	
	/**
	 * Converts a point to a vector.
	 * 
	 * @param point
	 * @return
	 */
	private static Vector pointToVector(Point3d point) {
		return new Vector(point.x, point.y, point.z);
	}
}
