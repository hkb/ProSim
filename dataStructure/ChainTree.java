package dataStructure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import chemestry.AminoAcid;
import chemestry.AminoAcid.BondType;
import chemestry.AminoAcid.SecondaryStructure;
import chemestry.AminoAcid.Type;

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
	
	protected List<Tuple2<Type,SecondaryStructure>> proteinInformation;// information about the properties of the protein
	
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

		this.proteinInformation = parser.proteinInformation;
	}
	
	/**
	 * Creates a chain tree from an array of chain trees.
	 * 
	 * @param cTrees The trees to create the chain tree from.
	 */
	public ChainTree(ChainTree[] cTrees) {
		this(getChainTreesCombinedBackboneAtomPositions(cTrees));
		
		// stuff to copy protein data
		for (ChainTree cTree : cTrees) {
			// copy primary structure information
			this.proteinInformation.addAll(cTree.proteinInformation);

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
		
		// lock the end leafs as rotation about them is nonsense
		this.backboneBonds[0].isLocked = true;
		this.backboneBonds[this.backboneBonds.length-1].isLocked = true;
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
		return (this.backboneBonds.length / 3)+1;
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
		
		// copy protein information
		cTree.proteinInformation = this.proteinInformation.subList(start-1, end);
		
		return cTree;
	}
	
	/**
	 * Returns the dihedral angles defined by the backbone bonds.
	 * 
	 * @return The dihedral angles.
	 */
	public List<Double> getDihedralAngles() {
		return this.getDihedralAngles(1, this.length());
	}
	
	public List<Double> getDihedralAngles(int start, int end) {
		List<Point3D> points = this.getBackboneAtomPositions((start != 1) ? start-1 : start, 
															 (end != this.length()) ? end+1 : end);

		List<Double> dihedralAngles = new ArrayList<Double>();
		
		// first angle is zero
		if(start == 1)
			dihedralAngles.add(0.0);
		 
		// init computation
		int j = (end != this.length()) ? points.size()-2 : points.size();
		int i = (start == 1) ? 3 : 5;
		
		Vector3D p1 = points.get(i-3).asVector();
		Vector3D p2 = points.get(i-2).asVector();
		Vector3D p3 = points.get(i-1).asVector();
		Vector3D p4;

		Vector3D v1 = p1.vectorTo(p2);
		Vector3D v2 = p2.vectorTo(p3);
		
		// compute all angles			
		for (; i < j; i++) {
			p4 = points.get(i).asVector();
			Vector3D v3 = p3.vectorTo(p4);			
			
			double a = v2.length() * v1.dot(v2.cross(v3));
			double b = v1.cross(v2).dot(v2.cross(v3));
			
			dihedralAngles.add(Math.atan2(a,b));
			
			p3 = p4;
			v1 = v2;
			v2 = v3;
		}

		// last angle is zero
		if(end == this.length())
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
	public CTNode l1 = null;
	public CTNode l2 = null;
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
		if (left.isLeaf() && right.isLeaf()) {
			this.l1 = left;
			this.l2 = right;
			return true;
		}
		
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
	 * Determines if this chain tree clashes with the others.
	 * 
	 * @param other The other chain trees.
	 * @return true if a clash occurs else false.
	 */
	public boolean areClashing(ChainTree[] others) {
		boolean areClashing = false;
		
		for(ChainTree other : others) {
			areClashing = areClashing(this.root, other.root, other);
			
			if(areClashing)
				break;
		}
		
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
	 * Sets the dihedral angle around the bond to the specified angle.
	 * 
	 * @param i The bond to rotate about.
	 * @param angle The desired dihedral angle.
	 */
	public void setRotationAngle(int i, double angle) {
		int aminoAcid = this.getAminoAcid(i);
		BondType bondType = this.getBondType(i);
		
		if(aminoAcid == 1 && bondType == BondType.PHI || aminoAcid == this.length() && bondType != BondType.PHI) {
			throw new IllegalArgumentException("Dihedral angles as nonsense for this bond!");
		}
		
		List<Double> angles = this.getDihedralAngles(aminoAcid, (aminoAcid == this.length()) ? aminoAcid : aminoAcid+1);
		
		double currentAngle = Double.NaN;
		
		switch(bondType) {
			case PHI:
				currentAngle = angles.get(0);
				break;
			case PSI:
				currentAngle = angles.get(1);
				break;
			case OMEGA:
				currentAngle = angles.get(2);
				break;
		}
		
		this.changeRotationAngle(i, angle-currentAngle);
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
				
				this.setRotationAngle(bond, angle);
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
		for (Tuple2<Type,SecondaryStructure> aminoAcid : this.proteinInformation) {
			tmpString.append(aminoAcid.x);
			
			if(aminoAcid.y == SecondaryStructure.HELIX) {
				tmpString.append(" (HELIX)");
			} else if(aminoAcid.y == SecondaryStructure.SHEET) {
				tmpString.append(" (SHEET)");
			}
			
			i++;
			
			if (i < this.proteinInformation.size())
				tmpString.append(", ");
		}
		
		return tmpString.toString();
	}
	
	
	/*
	 * Protein specific information.
	 */
	
	/**
	 * Is the amino acid part of a helix.
	 * 
	 * @param i The sequence number of the amino acid.
	 */
	public boolean isInHelix(int aminoAcid) {
		if(this.proteinInformation == null)
			return false;
		
		return this.proteinInformation.get(aminoAcid-1).y == SecondaryStructure.HELIX;
	}
	
	/**
	 * Is the amino acid part of a sheet.
	 * 
	 * @param i The sequence number of the amino acid.
	 */
	public boolean isInSheet(int aminoAcid) {
		if(this.proteinInformation == null)
			return false;
		
		return this.proteinInformation.get(aminoAcid-1).y == SecondaryStructure.SHEET;
	}
	
	/**
	 * Returns segments of all known helixes in the protein.
	 * 
	 * @return A sorted, none overlapping list of pairs indexes of the first and last 
	 * 		   residues (both included) of the known helixes. 
	 */
	public List<Tuple2<Integer,Integer>> getHelixSegments() {
		return this.getSecondaryStructureSegments(SecondaryStructure.HELIX);
	}

	/**
	 * Returns segments of all known sheets in the protein.
	 * 
	 * @return A sorted, none overlapping list of pairs indexes of the first and last 
	 * 		   residues (both included) of the known sheets. 
	 */
	public List<Tuple2<Integer,Integer>> getSheetSegments() {
		return this.getSecondaryStructureSegments(SecondaryStructure.SHEET);
	}

	/**
	 * Returns segments of all known segments not in a known secondary structure in the protein.
	 * 
	 * @return A sorted, none overlapping list of pairs indexes of the first and last 
	 * 		   residues (both included) of the segments not in a known secondary structure. 
	 */
	public List<Tuple2<Integer,Integer>> getIntermediateSegments() {
		return this.getSecondaryStructureSegments(SecondaryStructure.NONE);
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
	
	/**
	 * Gives the type of the amino acid.
	 * 
	 * @param aminoAcid The index of the amino acid.
	 * @return the type of the amino acid.
	 */
	public Type getAminoAcidType(int aminoAcid) {
		return this.proteinInformation.get(aminoAcid-1).x;
	}
	
	/**
	 * Returns the type of the bond.
	 * 
	 * @param bond The bond to get the type of.
	 * @return The type of the bond.
	 */
	public BondType getBondType(int bond) {
		switch(bond % 3) {
			case 0: return BondType.PHI;
			case 1: return BondType.PSI;
			case 2: return BondType.OMEGA;
		}
		
		// this should never happen
		throw new IllegalArgumentException("Unknown bond value!");
	}
	
	/**
	 * Analyses the backbone for segments of secondary structures of a specific type. 
	 * 
	 * @param type The type of the secondary structure.
	 * @return A sorted, none overlapping list of pairs indexes of the first and last 
	 * 		   residues (both included) of the segment. 
	 */
	private List<Tuple2<Integer,Integer>> getSecondaryStructureSegments(SecondaryStructure type) {
		int start = (this.proteinInformation.get(0).y == type) ? 1 : -1;
		int i = 1;
		
		List<Tuple2<Integer,Integer>> structure = new ArrayList<Tuple2<Integer,Integer>>();
		
		for(Tuple2<Type,SecondaryStructure> aminoAcid : this.proteinInformation) {
			if(start == -1 && aminoAcid.y == type) {
				start = i;
				
			} else if(start != -1 && aminoAcid.y != type) {
				structure.add(new Tuple2<Integer,Integer>(start, i-1));
				start = -1;
			}
			
			i++;
		}
		
		if(start != -1)
			structure.add(new Tuple2<Integer,Integer>(start, i-1));
		
		return structure;
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
