package dataStructure;

import javax.vecmath.Point3d;

import boundingVolume.LinesegmentSweptSphere;

import math.matrix.TransformationMatrix;

public class CTLeaf extends CTNode {
	
	private double angle;	// the rotation angle of this bond
	
	/**
	 * Create a leaf node from its position and the position of the previous node.
	 *  
	 * @param relativePosition The position of this node relative to the previous node.
	 */
	public CTLeaf(Point3d relativePosition) {
		// transformation matrix
		this.transformationMatrix = new TransformationMatrix(relativePosition.x, relativePosition.y, relativePosition.z);
		
		// bounding volume
		this.boundingVolume = new LinesegmentSweptSphere(relativePosition, 1.0);
	}
	
	
	
	/**
	 * Returns the position of the leaf.
	 * 
	 * @return The position of the leaf.
	 */
	public Point3d getPosition() {
		TransformationMatrix m = this.transformationMatrix;
		return new Point3d(m.a14, m.a24, m.a34);
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
}
