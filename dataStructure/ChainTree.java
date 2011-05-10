package dataStructure;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import chemestry.AminoAcid;

import math.Point3D;
import math.Tuple2;
import math.Vector3D;
import math.matrix.TransformationMatrix;

import boundingVolume.BoundingVolume;

import edu.math.Vector;

import tool.PDBParser;


public class ChainTree {
	
	public CTNode root;												// the root node of the tree
	public CTLeaf[] backboneBonds;									// the leaf nodes of the tree (the bonds of the protein backbone)
	
	private Point3D position;										// the position of the left most node in the world
	private double angle;											// the rotating angle of this backbone
	public TransformationMatrix worldTransformation;				// the transformation to transform from the proteins local coordinates to the world 
	
	public List<AminoAcid.Type> primaryStructure = new ArrayList<AminoAcid.Type>();
	public List<Tuple2<Integer,Integer>> helixes = new ArrayList<Tuple2<Integer,Integer>>();		// set of all bonds in alpha helices
	public List<Tuple2<Integer,Integer>> sheets = new ArrayList<Tuple2<Integer,Integer>>();		// set of all bonds in beta sheets 
	
	private Set<Integer> rotatedBonds = new HashSet<Integer>();		// the last rotated bond
	private int lowestRotatedBond = Integer.MAX_VALUE;				// the index of the leftmost rotated bond
	
	

	/**
	 * Create a chain tree from its PDB id.
	 * 
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public ChainTree(String pdbId) {
		this(new PDBParser(pdbId).backboneAtomPositions);
		
		PDBParser parser = new PDBParser(pdbId); // we have to parse the file as java require that the first statement is a constructor
												 // call and i don't want to make a static method.
		this.helixes = parser.helixes;
		this.sheets = parser.sheets;
		this.primaryStructure = parser.primaryStructure;
	}
	
	/**
	 * Creates a chain tree from an array of chain trees.
	 * 
	 * @param cTrees The trees to create the chain tree from.
	 */
	public ChainTree(ChainTree[] cTrees) {
		this(getChainTreesCombinedBackboneAtomPositions(cTrees));
		
		// stuff to copy protein data
		int i = 0;
		for (ChainTree cTree : cTrees) {
			
			// copy secondary structure information
			for (Tuple2<Integer,Integer> helix : cTree.helixes) {
				this.helixes.add(new Tuple2<Integer,Integer>(helix.x+i, helix.y+i));
			}

			for (Tuple2<Integer,Integer> sheet : cTree.sheets) {
				this.sheets.add(new Tuple2<Integer,Integer>(sheet.x+i, sheet.y+i));
			}
			
			i += cTree.length();
		}
	}
	
	/**
	 * Creates a chain tree from a list of 3D points.
	 * 
	 * @param points The points of the protein backbone atoms.
	 */
	public ChainTree(List<Point3D> points) {
		// store the absolute position of the chain
		Point3D first = points.get(0);
		this.position = new Point3D(first.x, first.y, first.z);
		
		this.worldTransformation = new TransformationMatrix(this.position.x, this.position.y, this.position.z);
				
		/*
		 * Create leaf nodes for each bond.
		 */
		this.backboneBonds = new CTLeaf[points.size()-1];
		
		Point3D start = this.position;	// the start of the current bond
		Point3D end;					// the end of the current bond
		
		for (int i = 0, j = this.backboneBonds.length; i < j; i++) {
			end = points.get(i+1);
			
			// check validity of the backbone chain
			if (start.distance(end) > CTLeaf.atomRadius) {
				throw new IllegalArgumentException("The provided backbone is not intact!");
			}
			
			// create new leaf from its relative position to the next leaf
			this.backboneBonds[i] = new CTLeaf(new Point3D(end.x-start.x, end.y-start.y, end.z-start.z), i);
			
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
	 * The number of amino acids in the backbone.
	 * 
	 * @return Number of amino acids.
	 */
	public int length() {
		return (this.backboneBonds.length+1) / 3; // there are 3*n-1 bonds in a backbone of n amino acids
	}
	
	/**
	 * Returns the absolute position of the protein backbone atoms.
	 * 
	 * @return The points of the atoms.
	 */
	public List<Point3D> getBackboneAtomPositions() {
		return this.getBackboneAtomPositions(1, this.length());
	}

	/**
	 * Returns the absolute positions of the atoms in a subsegment of the protein backbone.
	 *  
	 * @param start The first amino acid of the segment.
	 * @param end The last amino acid of the segment.
	 * @require start <= end
	 * @return The points of the atoms in the segment.
	 */
	public List<Point3D> getBackboneAtomPositions(int start, int end) {
		int i = this.getPhi(start);
		int j = this.getPsi(end);
		
		TransformationMatrix transformationMatrix = new TransformationMatrix(this.worldTransformation);
		
		// if we aren't starting from the first bond then modify the transformation  
		if (i > 0) {
			transformationMatrix.multR(this.getTransformationMatrix(0, i));
		}
		
		
		List<Point3D> points = new ArrayList<Point3D>();
		
		// place all atoms related to bonds
		for (; i <= j; i++) {
			points.add(new Point3D(transformationMatrix.a14, transformationMatrix.a24, transformationMatrix.a34));
			
			transformationMatrix.multR(this.backboneBonds[i].transformationMatrix);
		}
		
		// place an extra atom after the last bond
		points.add(new Point3D(transformationMatrix.a14, transformationMatrix.a24, transformationMatrix.a34));
		
		return points;
	}

	/**
	 * Returns a list of all the rotatable (non locked) bonds ofprivate
	 * the backbone.
	 * 
	 * @return The rotatable bonds of the backbone.
	 */
	public List<Integer> rotatableBonds() {
		ArrayList<Integer> rotatableBonds = new ArrayList<Integer>();
		
		for (CTLeaf node : this.backboneBonds) {
			if (!node.isLocked) {
				rotatableBonds.add(node.low);
			}
		}
		
		return rotatableBonds;
	}	
	
	/**
	 * Returns a chain tree representing a sub chain of the backbone.
	 * 
	 * @param i The starting amino acid (included).
	 * @param j The ending amino acid (included).
	 * @require start <= end
	 */
	public ChainTree getSubchain(int start, int end) {
		ChainTree cTree = new ChainTree(this.getBackboneAtomPositions(start, end));
		
		// copy secondary structure information
		for (Tuple2<Integer, Integer> helix : this.helixes) {
			if (start <= helix.x && helix.x <= end || start <= helix.y && helix.y <= end) {
				cTree.helixes.add(new Tuple2<Integer,Integer>(Math.max(start, helix.x)-start+1, Math.min(end, helix.y)-start+1));
			}
		}

		for (Tuple2<Integer, Integer> sheet : this.sheets) {
			if (start <= sheet.x && sheet.x <= end || start <= sheet.y && sheet.y <= end) {
				cTree.sheets.add(new Tuple2<Integer,Integer>(Math.max(start, sheet.x)-start+1, Math.min(end, sheet.y)-start+1));
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
		List<Point3D> points = this.getBackboneAtomPositions();
		
		List<Double> dihedralAngles = new ArrayList<Double>();
		
		// first angle is zero
		dihedralAngles.add(0.0);
		
		// init computation
		Point3D q = points.get(0);
		Point3D r = points.get(1);
		Point3D s = points.get(2);
		
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
		
		// find nearest common ancestor from root
		CTNode ancestor = this.root;
		
		while (!ancestor.isLeaf()) {
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
	 * Is the amino acid part of a helix.
	 * 
	 * @param i The sequence number of the amino acid.
	 */
	public boolean isInHelix(int aminoAcid) {
		for (Tuple2<Integer,Integer> helix : this.helixes) {
			if (helix.x <= aminoAcid && aminoAcid <= helix.y)
				return true;
		}
		
		return false;
	}
	
	/**
	 * Is the amino acid part of a sheet.
	 * 
	 * @param i The sequence number of the amino acid.
	 */
	public boolean isInSheet(int aminoAcid) {
		for (Tuple2<Integer,Integer> sheet : this.sheets) {
			if (sheet.x <= aminoAcid && aminoAcid <= sheet.y)
				return true;
		}
		
		return false;
	}
	
	/**
	 * Tests the tree for a self-clash.
	 * 
	 * @return true if the tree clashes with it self.
	 */
	public boolean isClashing() {
		boolean isClashing = isClashing(this.root, this.root);
		
		this.rotatedBonds.clear();
		
		return isClashing;
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
		// neither does the neighbours neighbour since the bounding box of the i-th bond covers the i+1-th atom
		if (left.low + 2 >= right.high)
			return false;
		
		// if no change has occurred between the trees then they have not changed position internaly
		boolean hasChanged = false;
		
		for(int bond : this.rotatedBonds) {
			hasChanged = (left.low <= bond && bond <= left.high) ||   // in left subtree
		  				 (right.low <= bond && bond <= right.high) || // in right subtree
		  				 (left.low <= bond && bond <= right.high);    // between left and right subtree
			
			if(hasChanged)
				break;
		}
		
		if(!hasChanged)
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
			return this.isClashing(left, right.left) || 
				   this.isClashing(left, right.right);
			
		} else if (right.isLeaf()) {
			return this.isClashing(left.left, right) || 
				   this.isClashing(left.right, right);
			
		} else {
			// only split the larger volume to avoid future repeated checks
			// this still works as it is depth first into the largest volumes
			if (left.boundingVolume.volume() > right.boundingVolume.volume()) {
				return this.isClashing(left.left, right) || 
					   this.isClashing(left.right, right);
			} else {
				return this.isClashing(left, right.left) || 
					   this.isClashing(left, right.right);
			}		   
		}
	}
	
	/**
	 * Determines if this chain tree clashes with the other.
	 * 
	 * @param other The other chain tree.
	 * @return true if a clash occurs else false.
	 */
	public boolean areClashing(ChainTree other) {
		boolean areClashing = areClashing(this.root, other.root, other);
		
		this.lowestRotatedBond = Integer.MAX_VALUE;
		
		return areClashing;
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
		if (thisNode.high < this.lowestRotatedBond)
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
	 * @param angle The angle to rotate the bond by in radians.
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
		this.rotatedBonds.add(i);
		this.lowestRotatedBond = (i < this.lowestRotatedBond) ? i : this.lowestRotatedBond;
	}

	/**
	 * Unfolds the protein into some, non clashing, confirmation.
	 */
	public void unfold() {
		List<Double> dihedralAngles = new ArrayList<Double>();
		
		// filter phi and psi angles
		int i = 0;
		for (double angle : this.getDihedralAngles()) {
			if (angle != 0.0 && i % 3 != 2) {
				dihedralAngles.add(angle);
			}
			
			i++;
		}
		
		// unfold
		do {
			for (int bond : this.rotatableBonds()) {
				double angle = dihedralAngles.get((int) (Math.random() * dihedralAngles.size()));
				
				this.changeRotationAngle(bond, angle);
			}
		} while(this.isClashing());
	}
	
	/**
	 * Move the entire protein.
	 * 
	 * @param move The vector that defines the movement.
	 */
	public void move(Vector3D move) {
		this.position = new Point3D(new Vector3D(this.position).add(move));
		
		this.worldTransformation.multR(new TransformationMatrix(this.position.x, this.position.y, this.position.z));
	}
	
	/**
	 * Rotate the entire protein by the given angle.
	 * 
	 * @param The angle to rotate with in radians.
	 * @warning THIS METHOD ONLY ROTATES ABOUT THE AXIS FROM ORIGO TO THE FIRST ATOM!
	 */
	public void rotate(double angle) {
		this.angle += angle;
		this.worldTransformation.rotate(this.angle);
	}
	
	@Override
	public String toString() {
		StringBuilder tmpString = new StringBuilder();
		
		int i = 0;
		for (CTLeaf bond : this.backboneBonds) {
			tmpString.append(bond);
			
			if (this.isInHelix(i/3 + 1)) {
				tmpString.append("(HELIX)");
			} else if (this.isInSheet(i/3 + 1)) {
				tmpString.append("(SHEET)");
			}
			
			i++;
			
			if (i < this.backboneBonds.length)
				tmpString.append(", ");
		}
		
		return tmpString.toString();
	}
	
	
	/**
	 * The the phi bond from the amino acid.
	 * 
	 * @param aminoAcid The amino acid sequence number (1-indexed).
	 * @return The index of the phi bond (0-indexed).
	 */
	public int getPhi(int aminoAcid) {
		return aminoAcid * 3 - 3;
	}

	/**
	 * The the psi bond from the amino acid.
	 * 
	 * @param aminoAcid The amino acid sequence number (1-indexed).
	 * @return The index of the psi bond (0-indexed).
	 */
	public int getPsi(int aminoAcid) {
		return aminoAcid * 3 - 2;
	}

	/**
	 * The the omega bond from the amino acid.
	 * 
	 * @param aminoAcid The amino acid sequence number (1-indexed).
	 * @return The index of the omega bond (0-indexed).
	 */
	public int getOmega(int aminoAcid) {
		return aminoAcid * 3 - 1;
	}
	
	/**
	 * Get the sequence number of the alpha helix the bond is a part of. 
	 * 
	 * @param bond
	 * @return
	 */
	public int getAminoAcid(int bond) {
		return bond/3+1;
	}
	
	/*
	 * Static methods.
	 */
	
	/**
	 * Combines the absolute position of all backbone atoms in an array of ChainTrees.
	 * 
	 * The order of the trees determines the order of the backbone atoms.
	 * 
	 * @param cTrees The trees to calculate the positions from.
	 * @return A list of the combined backbone atom positions.
	 */
	private static List<Point3D> getChainTreesCombinedBackboneAtomPositions(ChainTree[] cTrees) {
		List<Point3D> points = new LinkedList<Point3D>();
		
		for (ChainTree cTree : cTrees) {
			points.addAll(cTree.getBackboneAtomPositions());
		}
		
		return points;
	}
}
