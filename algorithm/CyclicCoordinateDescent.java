package algorithm;

import java.awt.Color;
import java.util.List;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import tool.ChainTreeScene;

import math.Line3D;
import math.Point3D;
import math.Tuple2;
import math.Vector3D;

import dataStructure.ChainTree;
import edu.math.Line;
import edu.math.Vector;
import geom3d.Line3d;

/**
 * An implementation of the cyclic coordinate descent (CCD) algorithm
 * for choosing the rotation of bonds during loop closure.
 * 
 * @author hkb
 *
 */
public class CyclicCoordinateDescent {
	
	public static int TARGET_LENGTH;
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
		TARGET_LENGTH = target.length()+1;
		
		this.target = new Vector3D[TARGET_LENGTH];
		
		int i = 0;
		for (Point3D position : target.getBackboneAtomPositions()) {
			this.target[i] = new Vector3D(position);
			
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
	public double __getRotationAngle(int bond) {		
		// determine the rotation axis of the bond
		List<Point3D> bondAtoms = this.loop.getBackboneAtomPositions(bond, bond);
		Line3D rotationAxis = new Line3D(bondAtoms.get(0), bondAtoms.get(1));
		
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length()-2, this.loop.length()-1);
		Vector3D[] moving = new Vector3D[TARGET_LENGTH];
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			moving[i] = new Vector3D(movingTerminalAtoms.get(i));
		}
		
		// compute the vectors r, t, n
		Vector3D[] r = new Vector3D[TARGET_LENGTH];
		Vector3D[] f = new Vector3D[TARGET_LENGTH];
		Vector3D[] s = new Vector3D[TARGET_LENGTH];

		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			Vector3D M = moving[i];
			Vector3D F = this.target[i];
			Vector3D O = rotationAxis.projectOnto(new Point3D(M)).asVector();
			
			r[i] = O.vectorTo(M);		
			f[i] = O.vectorTo(F);
			s[i] = r[i].norm().cross(rotationAxis.x.vectorTo(rotationAxis.y).norm());
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
		
		double alpha = Math.atan(numerator / denominator);
		
		/*
		 * Compute second derivative.
		 */
		double secondDerivative = 0;
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			secondDerivative += Math.cos(alpha) * f[i].dot(r[i].norm()) * 2
								* r[i].length() + Math.sin(alpha) * f[i].dot(s[i].norm()) * 2
								* r[i].length();	
		}
		
		/*
		 * Compute resulting alpha.
		 */
		//Math.atan2(arg0, arg1)
		
		if (secondDerivative < 0) {
			if (alpha > 0) {
				alpha -= Math.PI;
			} else {
				alpha += Math.PI;
			}
		}
		
		return alpha;			
	}
	
	public double getRotationAngle(int bond) {
		// determine the rotation axis of the bond
		List<Point3D> bondAtoms = this.loop.getBackboneAtomPositions(bond, bond);
		Line3D rotationAxis = new Line3D(bondAtoms.get(0), bondAtoms.get(1));
		
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length()-TARGET_LENGTH+1, this.loop.length()-1);

		
		// compute the values b, c
		double b = 0;
		double c = 0;

		for (int i = 0; i < TARGET_LENGTH; i++) {
			Vector3D M = new Vector3D(movingTerminalAtoms.get(i));
			Vector3D F = this.target[i];
			Vector3D O = rotationAxis.projectOnto(new Point3D(M)).asVector();
			
			Vector3D r = O.vectorTo(M);
			Vector3D f = O.vectorTo(F);
			Vector3D s = r.norm().cross(rotationAxis.x.vectorTo(rotationAxis.y).norm());

			double r2 = 2 * r.length();
			b += r2 * f.dot(r.norm());
			c += r2 * f.dot(s.norm());
		}

		/*
		 * Compute alpha.
		 */
		double divisor = Math.sqrt(b*b+c*c);
		
		return Math.atan2(c / divisor, b / divisor);
	}
	
	/**
	 * Determines if the loop is closed.
	 */
	public double targetRMSDistance() {
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length()-TARGET_LENGTH+1, this.loop.length()-1);
		Point3D[] moving = new Point3D[TARGET_LENGTH];
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			moving[i] = movingTerminalAtoms.get(i);
		}
		
		//
		double rmsd = 0;
		
		for (int i = 0; i < TARGET_LENGTH; i++) {
			double tmp = moving[i].distance(new Point3D(this.target[i]));
			rmsd += tmp*tmp;
		}
		
		return Math.sqrt(rmsd / TARGET_LENGTH);
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
