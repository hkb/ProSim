package dataStructure;

import geom3d.PointSet3d;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3d;

import tool.PDBParser;

import math.matrix.TransformationMatrix;
import molecule.Protein;

public class ChainTree {
	
	private CTNode root;		// the root node of the tree
	private Point3d position;	// the position of the left most node in the world
	private CTLeaf[] backbone;	// the leaf nodes of the tree (the protein backbone) 
	
	

	/**
	 * Create a chain tree from its PDB id.
	 * Protein(pdbId, 2, true)
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public ChainTree(String pdbId) {
		this(new PDBParser(pdbId).backbone);
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
			this.backbone[i] = new CTLeaf(new Point3d(next.x-current.x, next.y-current.y, next.z-current.z));
			
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
		return false;
	}
	
	/**
	 * Is the bond part of a alpha helix.
	 * 
	 * @param i The index of the bond.
	 */
	public boolean isInAlphaHelix(int i) {
		return false;
	}
	
	/**
	 * Is the bond part of a beta sheet.
	 * 
	 * @param i The index of the bond.
	 */
	public boolean isInBetaSheet(int i) {
		return false;
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
	}
	
	/**
	 * Returns a set of the points in the protein.
	 * 
	 * THIS IS NECCESERY ONLY BECAUSE JAVA REQUIRES THE super CALL TO
	 * BE THE FIRST IN THE CONSTRUCTOR!
	 * 
	 * @param protein The protein to calculate the points for.
	 * @return The points of the protein.
	 */
	private static List<Point3d> extractProteinPoints(Protein protein) {
		List<Point3d> points = new LinkedList<Point3d>();
		
		for (geom3d.Point3d point : protein.getPointSet()) {
			points.add(new Point3d(point.x, point.y, point.z));
		}
		
		return points;
	}
}
