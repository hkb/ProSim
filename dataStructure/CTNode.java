package dataStructure;

import math.matrix.TransformationMatrix;
import boundingVolume.BoundingVolume;
import boundingVolume.LinesegmentSweptSphere;

public class CTNode {
	
	public CTNode left, right, parent;					// connected nodes in the tree
	public BoundingVolume boundingVolume;				// the bounding volume of the node
	public TransformationMatrix transformationMatrix;	// the nodes transformation matrix
	public int height;									// the height of the nodes subtree
	public int low, high;								// the lowest and highest covered backbone bond

	public boolean isLocked = false;					// is this node locked?
	public boolean active = false;
	
	/**
	 * Creates a new node from its children.
	 * 
	 * @param left The left child node.
	 * @param right The right child node.
	 */
	public CTNode(CTNode left, CTNode right) {
		// create tree structure
		this.left = left;
		this.right = right;
		
		left.parent = this;
		right.parent = this;
		
		// compute height and backbone span
		this.height = Math.max(this.left.height, this.right.height) + 1;
		this.low = this.left.low;
		this.high = this.right.high;

		this.update();
	}
	
	/**
	 * Dummy constructor.
	 */
	protected CTNode() {
		// dummy constructor used by CTLeaf
	}
	
	
	
	/**
	 * Returns the left child of this node.
	 * 
	 * @return The left child.
	 */
	public CTNode getLeft() {
		return this.left;
	}
	
	/**
	 * Returns the right child of this node.
	 * 
	 * @return The right child.
	 */
	public CTNode getRight() {
		return this.right;
	}
	
	/**
	 * Returns the height of this node.
	 * 
	 * @return The height of the node.
	 */
	public int getHeight() {
		return this.height;		
	}
	
	/**
	 * Is this node a leaf?
	 */
	public boolean isLeaf() {
		return this.left == null || this.right == null;
	}
	
	/**
	 * Updates the information stored in the node.
	 */
	public void update() {
		// transformation matrix
		this.transformationMatrix = new TransformationMatrix(this.left.transformationMatrix, this.right.transformationMatrix);

		
		// bounding volume
		this.boundingVolume = this.left.boundingVolume.combine(this.right.boundingVolume.transform(this.left.transformationMatrix));
	}
	
	@Override
	public String toString() {
		return "";//+this.height;
	}

}
