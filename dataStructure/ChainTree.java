package dataStructure;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.vecmath.Point3d;

import tool.BinaryTreePainter;
import tool.PDBParser;

import math.matrix.TransformationMatrix;

public class ChainTree {
	
	private CTNode root;		// the root node of the tree
	private Point3d position;	// the position of the left most node in the world
	private CTLeaf[] backbone;	// the leaf nodes of the tree (the protein backbone) 
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
		this.backbone = new CTLeaf[points.size()];
		
		Point3d current = offset;
		Point3d next;
		
		for (int i = 0, j = this.backbone.length; i < j; i++) {
			next = points.get(i);
			
			// create new leaf from its relative position to the next leaf
			this.backbone[i] = new CTLeaf(new Point3d(next.x-current.x, next.y-current.y, next.z-current.z), i);
			
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
		for (CTNode l : this.backbone) { currentLevel.add(l); }
		
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
		return this.backbone.length;
	}
	
	/**
	 * Returns the points of the protein bonds.
	 * 
	 * @return The points of the protein bonds.
import java.util.Arrays;
	 */
	public List<Point3d> getBackbonePoints () {
		TransformationMatrix transformationMatrix = new TransformationMatrix(this.position.x, this.position.y, this.position.z);
		
		List<Point3d> points = new ArrayList<Point3d>();
		
		for (int i = 0; i < this.backbone.length; i++) {
			transformationMatrix.multR(this.backbone[i].transformationMatrix);
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
	
	/*
	 * Clash detection stuff.
	 */
	
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
	int indent = 0;
	private boolean isClashing(CTNode left, CTNode right) {
		indent++;
		print(left + " vs " + right);

		// neighbouring atoms may not collide
		if (left.low + 3 > right.high) {
			print("neighbours");
			indent--;
			return false;
		}
		
		// if no change has occurred between the trees then they have not changed position internal
		if (!(left.low <= this.lastRotatedBond && this.lastRotatedBond <= right.high)) {
			print("no change");
			indent--;
			return false;
		}
		
		// overlap?
		boolean overlap = true;
		if (left.low != right.low) {
			overlap = left.boundingVolume.isOverlaping(right.boundingVolume.transform(this.getTransformationMatrix(left.low, right.low)));
		}

		// if not then break
		if (!overlap) {
			print("no overlap");
			indent--;
			return false;
		}
		
		// if leaves then report clash
		if (left.isLeaf() && right.isLeaf()) {
			error("CLASH");
			indent--;
			return true;
		}
		
		// continue search
		boolean r;
		if(left.isLeaf()) {
			print("ll");
			r = isClashing(left, right.left) || isClashing(left, right.right);
			
		} else if (right.isLeaf()) {
			print("rl");
			r = isClashing(left.left, right) || isClashing(left.right, right);
			
		} else {
			print("llrl");
			r = isClashing(left, right.left) || isClashing(left, right.right) ||
				   isClashing(left.left, right) || isClashing(left.right, right);
		}
		
		indent--;
		return r;
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

		CTNode ancestor = this.backbone[i].parent;
		
		// find nearest common ancestor
		while (ancestor.high < j) {
			ancestor = ancestor.parent;
		}
		
		// if the ancestor exactly covers the nodes then just return it
		if (ancestor.low == i && ancestor.high == j) {
			return ancestor.transformationMatrix;
		}
		
		// go down from ancestor to the bonds while building a transformation matrix
		TransformationMatrix transformationMatrix = new TransformationMatrix();
		CTNode node;
			
		// traverse left subtree of ancestor
		node = ancestor.left;
		while (node != null) {
			if (node.low == i) {
				transformationMatrix.multL(node.transformationMatrix);
				node = null;
			} else {
				// witch subtree to choose
				if (i <= node.left.high) { // left
					transformationMatrix.multL(node.transformationMatrix);
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
				node = null;
			} else {
				// witch subtree to choose
				if (j >= node.right.low) { // right
					transformationMatrix.multR(node.transformationMatrix);
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
		this.backbone[i].rotate(angle);
		
		// propagate changes up thought the tree
		CTNode node = this.backbone[i];
		
		while(node.parent != null) {
			node = node.parent;
			node.update();
		}
		
		this.lastRotatedBond = i;
	}
	
	private  void print(String str) {
		String indent = "";
		for (int i = 0; i < this.indent; i++) indent += "  ";
		System.out.println(indent + str);
	}
	
	private  void error(String str) {
		System.err.println(str);
	}
}
