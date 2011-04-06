package dataStructure;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.Point3d;

import boundingVolume.BoundingVolume;

import edu.math.Vector;

import tool.PDBParser;

import math.matrix.TransformationMatrix;

public class ChainTree {
	
	public CTNode root;											// the root node of the tree
	public CTLeaf[] backboneBonds;								// the leaf nodes of the tree (the bonds of the protein backbone)
	
	private Point3d position;									// the position of the left most node in the world
	private double angle;										// the rotating angle of this backbone
	private TransformationMatrix worldTransformation;			// the transformation to transform from the proteins local coordinates to the world 
	
	protected Set<Integer> alphaHelix = new HashSet<Integer>();	// set of all bonds in alpha helices
	protected Set<Integer> betaSheet = new HashSet<Integer>();	// set of all bonds in beta sheets 
	protected Set<Integer> heteroAtoms = new HashSet<Integer>();// set of all hetero atoms
	
	private int lastRotatedBond;								// the last rotated bond
	
	

	/**
	 * Create a chain tree from its PDB id.
	 * 
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public ChainTree(String pdbId) {
		this(new PDBParser(pdbId).backbone);
		
		PDBParser parser = new PDBParser(pdbId); // we have to parse the file as java require that the first statement is a constructor
												 // call and i don't want to make a static method.
		this.alphaHelix = parser.alphaHelix;
		this.betaSheet = parser.betaSheet;
		this.heteroAtoms = parser.heteroAtoms;
	}
	
	/**
	 * Creates a chain tree from a list of 3D points.
	 * 
	 * @param points The points of the protein backbone atoms.
	 */
	public ChainTree(List<Point3d> points) {
		// store the absolute position of the chain
		Point3d first = points.get(0);
		this.position = new Point3d(first.x, first.y, first.z);
		
		this.worldTransformation = new TransformationMatrix(this.position.x, this.position.y, this.position.z);
				
		/*
		 * Create leaf nodes for each bond.
		 */
		this.backboneBonds = new CTLeaf[points.size()-1];
		
		Point3d start = this.position;	// the start of the current bond
		Point3d end;					// the end of the current bond
		
		for (int i = 0, j = this.backboneBonds.length; i < j; i++) {
			end = points.get(i+1);
			
			// check validity of the backbone chain
			if (start.distance(end) > CTLeaf.atomRadius) {
				throw new IllegalArgumentException("The provided backbone is not intact!");
			}
			
			// create new leaf from its relative position to the next leaf
			this.backboneBonds[i] = new CTLeaf(new Point3d(end.x-start.x, end.y-start.y, end.z-start.z), i);
			
			start = end;
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
		List<CTNode> currentLevel = new ArrayList<CTNode>();
		List<CTNode> nextLevel; 
		
		// initialise current level with all the leaves
		for (CTNode l : this.backboneBonds) { 
			currentLevel.add(l); 
		}

		// if the current level contains more than one node then group them in a new level
		do {
			nextLevel = new ArrayList<CTNode>();
			
			for (int i = 0, j = currentLevel.size(); i < j; i += 2) {
				if (i+1 < currentLevel.size()) {
					// two or more nodes left
					nextLevel.add(new CTNode(currentLevel.get(i), currentLevel.get(i+1)));
					
				} else if (i < currentLevel.size()) {
					// only one node left
					nextLevel.add(currentLevel.get(i));
				}
			}

			currentLevel = nextLevel;
		} while (currentLevel.size() > 1);
		
		// when only on node in the current level then store it as the root
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
	 * Returns the absolute position of the protein backbone atoms.
	 * 
	 * @return The points of the atoms.
	 */
	public List<Point3d> getBackboneAtomPositions() {
		TransformationMatrix transformationMatrix = new TransformationMatrix(this.worldTransformation);
		
		List<Point3d> points = new ArrayList<Point3d>();
		
		points.add(this.position);
		
		for (int i = 0; i < this.backboneBonds.length; i++) {
			transformationMatrix.multR(this.backboneBonds[i].transformationMatrix);
			points.add(new Point3d(transformationMatrix.a14, transformationMatrix.a24, transformationMatrix.a34));
		}

		return points;
	}
	
	/**
	 * Returns a chain tree representing a sub chain of the backbone.
	 * 
	 * @param i The starting atom (included).
	 * @param j The ending atom (included).
	 * @require Neither i nor j may split a amino acid.
	 */
	public ChainTree getSubchain(int i, int j) {
		if (i % 3 != 0 || j % 3 != 2) {
			throw new IllegalArgumentException("You can't split an amino acid!");
		}
		
		ChainTree cTree = new ChainTree(this.getBackboneAtomPositions().subList(i, j+1));
		
		// copy secondary structure information
		for (int k = i; k <= j; k++) {
			if (this.alphaHelix.contains(k)) {
				cTree.alphaHelix.add(k-i);
			}
		}

		for (int k = i; k <= j; k++) {
			if (this.betaSheet.contains(k)) {
				cTree.betaSheet.add(k-i);
			}
		}
		
		for (int k = i; k <= j; k++) {
			if (this.heteroAtoms.contains(k)) {
				cTree.heteroAtoms.add(k-i);
			}
		}
		
		return cTree;
	}
	
	/**
	 * Returns the dihedral angles defined by the backbone bonds.
	 * 
	 * @return The dihedral angles.
	 */
	public List<Double> getDihedralAngles() {
		List<Point3d> points = this.getBackboneAtomPositions();
		
		List<Double> dihedralAngles = new ArrayList<Double>();
		
		// first angle is zero
		dihedralAngles.add(0.0);
		
		// init computation
		Point3d q = points.get(0);
		Point3d r = points.get(1);
		Point3d s = points.get(2);
		
		Vector qr = new Vector(q.x-r.x, q.y-r.y, q.z-r.z);
		Vector rs = new Vector(r.x-s.x, r.y-s.y, r.z-s.z);
		Vector pq;
		
		// compute all angles
		for (int i = 1, j = points.size()-2; i < j; i++) {
			q = r;
			r = s;
			pq = qr;
			qr = rs;
			s = points.get(i+2);
			rs = new Vector(r.x-s.x, r.y-s.y, r.z-s.z);
			
			dihedralAngles.add((double) Vector.dihedralAngle(pq, qr, rs));
		}

		// last angle is zero
		dihedralAngles.add(0.0);
		
		return dihedralAngles;
	}
	
	/**
	 * Returns a transformation matrix from the given backbone coordinate system
	 * to the world coordinate system.
	 * 
	 * @param j The bond to transform from.
	 * @return The transformation from the j-th coordinate system into the world coordinate system
	 */
	public TransformationMatrix getWorldTransformation(int j) {
		TransformationMatrix transformationMatrix = this.getTransformationMatrix(0, j);
		
		transformationMatrix.multL(this.worldTransformation);
		
		return transformationMatrix;
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
		if (j < i) {
			throw new IllegalArgumentException("i ("+i+") must be less than or equal to j ("+j+")");
		}

		// unit transformation
		if (i == j) {
			return new TransformationMatrix();
		}
		
		// NOTE: from this point onwards we look for the transformation up to, but not including, the j-th coordinate system
		j--;  
		
		// find nearest common ancestor
		CTNode ancestor = this.root;
		
		while (ancestor.left != null) {
			if (j <= ancestor.left.high) {
				ancestor = ancestor.left;
			} else if (ancestor.right.low <= i) {
				ancestor = ancestor.right;
			} else {
				break;
			}
		}
		
		// if the ancestor exactly covers the nodes then just return it
		if (ancestor.low == i && ancestor.high == j) {
			return new TransformationMatrix(ancestor.transformationMatrix);
		}
		
		// go down from ancestor to the bonds while building the transformation matrix
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
	 * Is the bond to the right of the a hetero atom.
	 * 
	 * @param i
	 * @return
	 */
	public boolean isHeteroAtomBond(int i) {
		return this.heteroAtoms.contains(i);
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
		// neither does the neighbours neighbour since the bounding box of the i-th bond completely covers the i+1-th atom
		if (left.isLeaf() && right.isLeaf() && left.low + 2 >= right.high)
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
			// only split the larger volume to avoid future repeated checks
			// TODO prove that this works
			if (left.boundingVolume.volume() > right.boundingVolume.volume()) {
				return isClashing(left.left, right) || isClashing(left.right, right);
			} else {
				return isClashing(left, right.left) || isClashing(left, right.right);
			}		   
		}
	}
	
	/**
	 * Determines if the two sub chains has changed internally by the last rotation.
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
	 * Determines if this chain tree clashes with the other.
	 * 
	 * @param other The other chain tree.
	 * @return true if a clash occurs else false.
	 */
	public boolean areClashing(ChainTree other) {
		return areClashing(this.root, other.root, other);
	}
	
	/**
	 * Determines if the sub chain in this tree represented by this node clashes with the
	 * sub chain represented by the node in the other tree.
	 * 
	 * @param thisNode The node to test in this tree.
	 * @param otherNode The node to test in the other tree.
	 * @param other The other tree.
	 * @return true if the trees clash else false.
	 */
	private boolean areClashing(CTNode thisNode, CTNode otherNode, ChainTree other) {
		// has this node been moved in the world?
		// this works because only the chain segment to the right of the last rotation is moved in the world
		if (thisNode.high < this.lastRotatedBond)
			return false;
		
		// check for overlap
		BoundingVolume thisVolume = thisNode.boundingVolume.transform(this.getWorldTransformation(thisNode.low));
		BoundingVolume otherVolume = otherNode.boundingVolume.transform(other.getWorldTransformation(otherNode.low));

		// if no overlap then stop
		if(!thisVolume.isOverlaping(otherVolume))
			return false;
		
		// leaves are clashing
		if(thisNode.isLeaf() && otherNode.isLeaf())
			return true;
		
		// continue search
		if(thisNode.isLeaf()) {
			return this.areClashing(thisNode, otherNode.left, other) ||
			       this.areClashing(thisNode, otherNode.right, other);
			
		} else if(otherNode.isLeaf()) {
			return this.areClashing(thisNode.left, otherNode, other) ||
		       	   this.areClashing(thisNode.right, otherNode, other);
			
		} else {
			// only split the larger volume to avoid future repeated checks
			if (thisNode.boundingVolume.volume() > otherNode.boundingVolume.volume()) {
				return this.areClashing(thisNode.left, otherNode, other)  ||
					   this.areClashing(thisNode.right, otherNode, other);
			} else {
				return this.areClashing(thisNode, otherNode.left, other)  ||
					   this.areClashing(thisNode, otherNode.right, other);
			}
		}
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
		
		// remember the last rotated bond for checking algorithms
		this.lastRotatedBond = i;
	}
	
	/**
	 * Unfolds the protein into some, non clashing, confirmation.
	 */
	public void unfold() {
		List<Double> dihedralAngles = this.getDihedralAngles();
		
		for (int i = 0, j = this.backboneBonds.length; i < j; i++) {
			int d = 180;
			do {
				this.changeRotationAngle(i, d-dihedralAngles.get(i));
				d--;
			} while(this.isClashing());
		}
	}
	
	/**
	 * Move the entire protein.
	 * 
	 * @param move The vector that defines the movement.
	 */
	public void move(Point3d move) {
		this.position.add(move);
		
		this.worldTransformation = new TransformationMatrix(this.position.x, this.position.y, this.position.z);
		this.worldTransformation.rotate(this.angle);
	}
	
	/**
	 * Rotate the entire protein by the given angle.
	 * 
	 * @param The angle to rotate with.
	 */
	public void rotate(double angle) {
		this.angle -= angle;
		this.worldTransformation.rotate(this.angle);
	}
	
	@Override
	public String toString() {
		StringBuilder tmpString = new StringBuilder();
		
		int i = 0;
		for (CTLeaf bond : this.backboneBonds) {
			tmpString.append(bond);
			
			if (this.isInAlphaHelix(i)) {
				tmpString.append("(HELIX)");
			} else if (this.isInBetaSheet(i)) {
				tmpString.append("(SHEET)");
			} else if (this.isHeteroAtomBond(i)) {
				tmpString.append("(HETERO)");
			}
			
			tmpString.append("\n");
			i++;
		}
		
		return tmpString.toString();
	}
}
