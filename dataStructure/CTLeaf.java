package dataStructure;

import math.Point3D;
import math.matrix.TransformationMatrix;


import boundingVolume.LinesegmentSweptSphere;


public class CTLeaf extends CTNode {
	
	public static double atomRadius = 1.7;	// the radius of the atom where this bond starts
	private double angle;					// the rotation angle of this bond
	
	
	
	/**
	 * Create a leaf node from its position and the position of the previous node.
	 *  
	 * @param relativePosition The position of this node relative to the previous node.
	 * @param i The index of the bond in the backbone.
	 */
	public CTLeaf(Point3D relativePosition, int i) {
		// the span of this sub chain (just this bond)
		super.low = super.high = i;
		
		// transformation matrix
		this.transformationMatrix = new TransformationMatrix(relativePosition.x, relativePosition.y, relativePosition.z);
		
		// bounding volume
		this.boundingVolume = new LinesegmentSweptSphere(relativePosition, atomRadius/2);
	}
	
	
	
	/**
	 * Rotates the bond by the given angle.
	 * 
	 * @param angle The angle to rotate with.
	 */
	public void rotate(double angle) {
		this.angle -= angle;
		
		transformationMatrix.rotate(this.angle);
	}
	
	@Override
	public String toString() {
		return "" + this.low;
	}
}
