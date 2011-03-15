package dataStructure;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.vecmath.Point3d;

import tool.PDBParser;

import math.matrix.TransformationMatrix;

public class ChainTree {
	
	protected CTNode root;				// the root node of the tree
	protected Point3d position;			// the position of the left most node in the world
	protected CTLeaf[] backboneBonds;	// the leaf nodes of the tree (the bonds of the protein backbone) 
	
	protected Set<Integer> alphaHelix = new HashSet<Integer>();	// set of all bonds in alpha helices
	protected Set<Integer> betaSheet = new HashSet<Integer>();	// set of all bonds in beta sheets 
	
	private int lastRotatedBond;	// the last rotated bond
	
	

	/**
	 * Create a chain tree from its PDB id.
	 * 
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public ChainTree(String pdbId) {
		this(new PDBParser(pdbId).backbone);
		this.alphaHelix = new PDBParser(pdbId).alphaHelix;
		this.betaSheet = new PDBParser(pdbId).betaSheet;
	}
	
	/**
	 * Creates a chain tree from a list of 3D points.
	 * 
	 * @param points The points of the protein backbone atoms.
	 */
	public ChainTree(List<Point3d> points) {
		this.position = points.get(0);
				
		/*
		 * Create leaf nodes.
		 */
		this.backboneBonds = new CTLeaf[points.size()-1];
		
		Point3d current = this.position; // necessary to get relative positions of first element
		Point3d next;
		
		for (int i = 0, j = this.backboneBonds.length; i < j; i++) {
			next = points.get(i+1);
			
			// create new leaf from its relative position to the next leaf
			this.backboneBonds[i] = new CTLeaf(new Point3d(next.x-current.x, next.y-current.y, next.z-current.z), i);
			
			current = next;
		}
		
		/*
		 * Build tree structure.
		 * 
		 * Strategy: Maintain two levels of the tree; the current (initially the leafs)
		 * while building the next level by (greedy) paring nodes from left to right. If
		 * a single node is left at the rightmost node then just propagate it to the next
		 * level.
		 * This will end up with just one node in the current level, that is the root of 
		 * the tree.
		 */
		List<CTNode> currentLevel = new LinkedList<CTNode>();
		for (CTNode l : this.backboneBonds) { currentLevel.add(l); }
		
		List<CTNode> nextLevel; 
		
		do {
			nextLevel = new LinkedList<CTNode>();
			
			for (int i = 0, j = currentLevel.size(); i < j; i += 2) {
				if (i+1 < currentLevel.size()) { // two or more nodes left
					nextLevel.add(new CTNode(currentLevel.get(i), currentLevel.get(i+1)));
					
				} else if (i < currentLevel.size()) { // only one node left
					nextLevel.add(currentLevel.get(i));
				}
			}

			currentLevel = nextLevel;
		} while (currentLevel.size() > 1);
		
		// the root has been found
		this.root = currentLevel.get(0);
	}	
	
	
	
	/**
	 * Returns the root node of the tree.
	 * 
	 * @return The root node.
	 */
	public CTNode getRoot() {
		return this.root;
	}
	
	/**
	 * The length of the backbone.
	 */
	public int length() {
		return this.backboneBonds.length;
	}
	
	/**
	 * 
	 */
	public void unfold() {
		for (int i = 0, j = this.backboneBonds.length; i < j; i++) {
			this.changeRotationAngle(i, 180);
		}
	}
	
	/**
	 * Returns the points of the protein bonds.
	 * 
	 * @return The points of the protein bonds.
	 */
	public List<Point3d> getBackbonePoints () {
		TransformationMatrix transformationMatrix = new TransformationMatrix(this.position.x, this.position.y, this.position.z);
		
		List<Point3d> points = new ArrayList<Point3d>();
		
		for (int i = 0; i < this.backboneBonds.length; i++) {
			transformationMatrix.multR(this.backboneBonds[i].transformationMatrix);
			points.add(new Point3d(transformationMatrix.a14, transformationMatrix.a24, transformationMatrix.a34));
		}

		return points;
	}
	
	/**
	 * Is the bond part of a alpha helix.
	 * 
	 * @param i The index of the bond.
	 */
	public boolean isInAlphaHelix(int i) {
		return this.alphaHelix.contains(i);
	}
	
	/**
	 * Is the bond part of a beta sheet.
	 * 
	 * @param i The index of the bond.
	 */
	public boolean isInBetaSheet(int i) {
		return this.betaSheet.contains(i);
	}
	
	/**
	 * Tests the tree for a self-clash.
	 * 
	 * @return true if the tree clashes with it self.
	 */

	public boolean isClashing() {
		return isClashing(this.root, this.root);
	}
	
	/**
	 * Check if the sub-chains of two nodes does clash.
	 * 
	 * @param node1 A node to check for overlap.
	 * @param node2 A node to check for overlap.
	 * @return true if there is a clash else false
	 */
	private boolean isClashing(CTNode left, CTNode right) {
		
		// NOTE: This is a purely technical check to avoid double check of the same subtrees
		if (right.low < left.low)
			return false;
		
		// neighbouring atoms does not cause a clash
		if (left.isLeaf() && right.isLeaf() && left.low + 1 >= right.high)
			return false;
		
		// if no change has occurred between the trees then they have not changed position internal
		if (!this.hasChanged(left, right))
			return false;
		
		// check for overlap
		boolean overlap = left.boundingVolume.isOverlaping(right.boundingVolume.transform(this.getTransformationMatrix(left.low, right.low)));

		if (!overlap)
			return false;
		
		// if leaves then report clash
		if (left.isLeaf() && right.isLeaf())
			return true;
		
		// continue search
		if(left.isLeaf()) {
			return isClashing(left, right.left) || isClashing(left, right.right);
			
		} else if (right.isLeaf()) {
			return isClashing(left.left, right) || isClashing(left.right, right);
			
		} else {
			// only split the larger volume
			// TODO WHY DOES THIS WORK?!
			if (left.boundingVolume.volume() > right.boundingVolume.volume()) {
				return isClashing(left.left, right) || isClashing(left.right, right);
			} else {
				return isClashing(left, right.left) || isClashing(left, right.right);
			}		   
		}
	}
	
	/**
	 * Determines if the two sub chains has changed internally by the last update.
	 * 
	 * @param left The left most sub chain.
	 * @param right The rightmost sub chain.
	 * @require left.low <= right.low
	 * @return true if a bond affecting the subtrees has changed
	 */
	private boolean hasChanged(CTNode left, CTNode right) {
		if (right.low < left.low) {
			throw new IllegalArgumentException("Invalid argument order!");
		}
		
		int lastRotatedBond = this.lastRotatedBond;
		
		return (left.low <= lastRotatedBond && lastRotatedBond <= left.high) ||   // in left subtree
			   (right.low <= lastRotatedBond && lastRotatedBond <= right.high) || // in right subtree
			   (left.low <= lastRotatedBond && lastRotatedBond <= right.high);    // between left and right subtree
	}

	/**
	 * Calculates the matrix to transform the j-th coordinate system
	 * into the i-th coordinate system.
	 * 
	 * @param i The i-th bond.
	 * @param j The j-th bond.
	 * @require i < j
	 * @return The matrix to transform j into i.
	 */
	public TransformationMatrix getTransformationMatrix(int i, int j) {
		if (j < i) {
			throw new IllegalArgumentException("i ("+i+") must be less than or equal to j ("+j+")");
		}
		
		// unit transformation
		if (i == j) {
			return new TransformationMatrix();
		}
		
		// find nearest common ancestor
		CTNode ancestor = this.backboneBonds[i].parent;
		
		while (ancestor.high < j) {
			ancestor = ancestor.parent;
		}
		
		// if the ancestor exactly covers the nodes then just return it
		if (ancestor.low == i && ancestor.high == j) {
			return new TransformationMatrix(ancestor.transformationMatrix);
		}
		
		// go down from ancestor to the bonds while building a transformation matrix
		// TODO understand this section!!
		TransformationMatrix transformationMatrix = new TransformationMatrix();
		CTNode node;
		
		// traverse left subtree of ancestor
		node = ancestor.left;
		while (node != null) {
			if (node.low == i) {
				transformationMatrix.multL(node.transformationMatrix);
				break;
			} else {
				// witch subtree to choose
				if (i <= node.left.high) { // left
					transformationMatrix.multL(node.right.transformationMatrix);
					node = node.left;
				} else { // right
					node = node.right;
				}
			}
		}
		
		// traverse right subtree of ancestor
		node = ancestor.right;
		while (node != null) {
			if (node.high == j) {
				transformationMatrix.multR(node.transformationMatrix);
				break;
			} else {
				// witch subtree to choose
				if (j >= node.right.low) { // right
					transformationMatrix.multR(node.left.transformationMatrix);
					node = node.right;
				} else { // left
					node = node.left;
				}
			}
		}

		return transformationMatrix;
	}
	
	/**
	 * Changes the rotation angle of the i-th bond by the specified angle.
	 * 
	 * @param i The index of the bond.
	 * @param angle The angle to rotate the bond by.
	 */
	public void changeRotationAngle(int i, double angle) {
		CTLeaf bond = this.backboneBonds[i];
		
		// update the bonds transformation matrix
		bond.rotate(angle);
		
		// propagate changes up thought the tree
		CTNode node = (CTNode) bond;
		
		while(node.parent != null) {
			node = node.parent;
			node.update();
		}
		
		this.lastRotatedBond = i;
	}
}
