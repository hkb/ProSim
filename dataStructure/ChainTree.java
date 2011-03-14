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
	
	private CTNode root;		// the root node of the tree
	private Point3d position;	// the position of the left most node in the world
	private CTLeaf[] backboneBonds;	// the leaf nodes of the tree (the bonds of the protein backbone) 
	
	public Set<Integer> alphaHelix = new HashSet<Integer>();	// set of all bonds in alpha helices
	public Set<Integer> betaSheet = new HashSet<Integer>();		// set of all bonds in beta sheets 
	
	private int lastRotatedBond;
	
	

	/**
	 * Create a chain tree from its PDB id.
	 * Protein(pdbId, 2, true)
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public ChainTree(String pdbId) {
		this(new PDBParser(pdbId).backbone);
		this.alphaHelix = new PDBParser(pdbId).alphaHelix;
		this.betaSheet = new PDBParser(pdbId).betaSheet;
	}
	
	/**
	 * Creates a chain tree from a list of 3D points.
	 */
	public ChainTree(List<Point3d> points) {
		this(points, new Point3d(0, 0, 0));
	}
	
	/**
	 * Creates a chain tree from a list of 3D points with a given offset.
	 */
	public ChainTree(List<Point3d> points, Point3d offset) {
		this.position = offset;
				
		/*
		 * Create leaf nodes.
		 */
		this.backboneBonds = new CTLeaf[points.size()-1];
		
		Point3d current = points.get(0); // necessary to get relative positions of first element
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
	 * Is the bond an peptide bond?
	 * 
	 * @param i The index of the bond.
	 */
	public boolean isPeptide(int i) {
		return i % 3 == 1;
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

		// neighbouring atoms does not cause a clash
		if (left.isLeaf() && right.isLeaf() && left.low + 1 >= right.high) {
			return false;
		}
		
		// if no change has occurred between the trees then they have not changed position internal
		if (!(left.low <= this.lastRotatedBond && this.lastRotatedBond <= right.high)) {
			return false;
		}
		
		// overlap?
		boolean overlap = true;
		if (left.low != right.low) {			
			overlap = left.boundingVolume.isOverlaping(right.boundingVolume.transform(this.getTransformationMatrix(left.low, right.low)));
		}

		// if not then break
		if (!overlap) {
			return false;
		}
		
		// if leaves then report clash
		if (left.isLeaf() && right.isLeaf()) {
			return true;
		}
		
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
	 * Determines if the two nodes has changed internally by the last update,
	 * @param left
	 * @param right
	 * @return
	 */
	private boolean hasChanged(CTNode left, CTNode right) {
		int lastRotatedBond = this.lastRotatedBond;
		
		return left.low <= lastRotatedBond && lastRotatedBond <= left.high || 
			   right.low <= lastRotatedBond && lastRotatedBond <= right.high ||
			   left.low <= lastRotatedBond && lastRotatedBond <= right.high;
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
	private TransformationMatrix getTransformationMatrix(int i, int j) {
		if (j<=i) {
			throw new IllegalArgumentException(i+"<="+j);
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
		// update the bonds transformation matrix
		this.backboneBonds[i].rotate(angle);
		
		// propagate changes up thought the tree
		CTNode node = this.backboneBonds[i];
		
		while(node.parent != null) {
			node = node.parent;
			node.update();
		}
		
		this.lastRotatedBond = i;
	}
}
