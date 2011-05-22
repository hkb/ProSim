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
	
	private ChainTree loop;						// the loop to close including the target anchor
	private Vector3D[] target;					// the position of the target

	
	/**
	 * Creates a new CCD computer from the unfolded loop and the target
	 * residue.
	 * 
	 * @param cTree The backbone segment that contains the loop and ends in the target.
	 * @param target The target residue.
	 */
	public CyclicCoordinateDescent(ChainTree cTree, ChainTree target) {		
		this.target = new Vector3D[target.length()*3];
		
		int i = 0;
		for (Point3D position : target.getBackboneAtomPositions()) {
			this.target[i] = position.asVector();
			
			i++;
		}
		
		// store loop 
		this.loop = cTree;
	}

	
	public double getRotationAngle(int bond) {
		// determine the rotation axis of the bond
		int aminoAcid = this.loop.getAminoAcid(bond);
		List<Point3D> bondAtoms = this.loop.getBackboneAtomPositions(aminoAcid, aminoAcid);
		Line3D rotationAxis = new Line3D(bondAtoms.get(bond % 3), bondAtoms.get(bond % 3 + 1));
		
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length(), this.loop.length());
		
		// compute the values b, c
		double b = 0;
		double c = 0;
		
		for (int i = 0; i < this.target.length; i++) {
			Vector3D M = new Vector3D(movingTerminalAtoms.get(i));
			Vector3D F = this.target[i];
			Vector3D O = rotationAxis.projectOnto(new Point3D(M)).asVector();
			
			Vector3D r = O.vectorTo(M);
			Vector3D f = O.vectorTo(F);
			Vector3D o = rotationAxis.x.vectorTo(rotationAxis.y).norm();
			Vector3D s = r.norm().cross(o);

			double r2 = 2 * r.length();
			b += r2 * f.dot(r.norm());
			c += r2 * f.dot(s.norm());
		}
		
		/*
		 * Compute alpha.
		 */
		double divisor = Math.sqrt(b*b+c*c);
		
		// TODO why does the result come off in the opposit direction?
		return -Math.atan2(c / divisor, b / divisor);
	}
	
	/**
	 * Determines if the loop is closed.
	 */
	public double targetRMSDistance() {
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length(), this.loop.length());
		
		// calculate rmsd
		double rmsd = 0;
		
		for (int i = 0; i < this.target.length; i++) {
			double tmp = movingTerminalAtoms.get(i).distance(new Point3D(this.target[i]));
			rmsd += tmp*tmp;
		}
		
		return Math.sqrt(rmsd / this.target.length);
	}
	
	
	
	
	/* DO NOT USE!!! */

	/**
	 * Returns the rotation around the given bond that brings the loop terminal
	 * residue closer to the target.
	 * 
	 * @deprecated ONLY FOR COMPARISON WITH NEW IMPLEMENTATION!!!
	 * @return The rotation angle that minimises the loop terminals distance to the target.
	 */
	public double __getRotationAngle(int bond) {		
		// determine the rotation axis of the bond
		List<Point3D> bondAtoms = this.loop.getBackboneAtomPositions(bond, bond);
		Line3D rotationAxis = new Line3D(bondAtoms.get(0), bondAtoms.get(1));
		
		// fetch the positions of the moving terminal residue
		List<Point3D> movingTerminalAtoms = this.loop.getBackboneAtomPositions(this.loop.length()-2, this.loop.length()-1);
		Vector3D[] moving = new Vector3D[this.target.length];
		
		for (int i = 0; i < this.target.length; i++) {
			moving[i] = new Vector3D(movingTerminalAtoms.get(i));
		}
		
		// compute the vectors r, t, n
		Vector3D[] r = new Vector3D[this.target.length];
		Vector3D[] f = new Vector3D[this.target.length];
		Vector3D[] s = new Vector3D[this.target.length];

		
		for (int i = 0; i < this.target.length; i++) {
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
		
		for (int i = 0; i < this.target.length; i++) {
			numerator += f[i].dot(s[i].norm()) * r[i].length();
			denominator += f[i].dot(r[i].norm()) * r[i].length();
		}
		
		double alpha = Math.atan(numerator / denominator);
		
		/*
		 * Compute second derivative.
		 */
		double secondDerivative = 0;
		
		for (int i = 0; i < this.target.length; i++) {
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
}
