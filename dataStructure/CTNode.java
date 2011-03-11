package dataStructure;

import java.util.Stack;

import energyFunction.EnergyFunction;

import math.matrix.TransformationMatrix;
import boundingVolume.BoundingVolume;

public class CTNode {
	
	protected CTNode left, right, parent;					// connected nodes in the tree
	private double energy;									// the energy of the current node
	protected BoundingVolume boundingVolume;				// the bounding volume of the node
	protected TransformationMatrix transformationMatrix;	// the nodes transformation matrix
	private int height;										// the height of the nodes subtree
	
	public boolean isLocked;								// is this node locked?
	
	
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
		
		this.update();
	}
	
	/**
	 * Dummy constructor.
	 */
	protected CTNode() {
		// dummy constructor used by CTLeaf
	}
	
	
	
	/**
	 * Updates the information stored in the node.
	 */
	public void update() {
		// compute height from the highest subtree
		this.height = (this.left.height > this.right.height ? this.left.height : this.right.height) + 1;
		
		// transformation matrix
		this.transformationMatrix = new TransformationMatrix(this.left.transformationMatrix, this.right.transformationMatrix);

		// bounding volume
		this.boundingVolume = this.left.boundingVolume.combine(this.left.boundingVolume, this.right.boundingVolume.transform(this.transformationMatrix));
	}
	
	/**
	 * Gets the left child.
	 * 
	 * @return The left child.
	 */
	public CTNode getLeft() {
		return this.left;
	}
	
	/**
	 * Gets the right child.
	 * 
	 * @return The right child.
	 */
	public CTNode getRight() {
		return this.right;
	}
	
	/**
	 * The height of the nodes subtree.
	 * 
	 * @return The height of the nodes subtree.
	 */
	public int getHeight() {
		return this.height;
	}
	
	/**
	 * Computes the energy of the node.
	 * 
	 * @param energyFunction
	 */
	public void computeEnergy(EnergyFunction energyFunction) {
		this.energy = energyFunction.compute(this);
	}
	
	@Override
	public String toString() {
		return this.height + "";
	}

}
