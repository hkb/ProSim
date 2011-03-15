package dataStructure;

import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3d;

import boundingVolume.LinesegmentSweptSphere;

import edu.geom3D.Capsule;
import edu.math.Vector;

import math.matrix.TransformationMatrix;

import tool.BinaryTreePainter;

public class AdjustableChainTree extends ChainTree {
	
	/**
	 * Create a chain tree from its PDB id.
	 * 
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public AdjustableChainTree(String pdbId) {
		super(pdbId);

		this.lockAndGroupPeptidePlanes();
		this.lockAndGroupAlphaHelices();
		this.lockAndGroupBetaSheets();
	}
	
	/**
	 * Creates a chain tree from a list of 3D points.
	 * 
	 * @param points The points of the protein backbone atoms.
	 */
	public AdjustableChainTree(List<Point3d> points) {
		super(points);

		this.lockAndGroupPeptidePlanes();
		this.lockAndGroupAlphaHelices();
		this.lockAndGroupBetaSheets();
	}
	
	@Override
	public void changeRotationAngle(int i, double angle) {
		if (super.backboneBonds[i].isLocked) {
			throw new IllegalArgumentException("Can't rotate locked bonds!");
		}
		
		super.changeRotationAngle(i, angle);
	}
	
	/**
	 * Returns a list of all the rotatable (non locked) bonds of
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
	 * Locks and groups peptide planes.
	 */
	private void lockAndGroupPeptidePlanes() {
		for (int i = 2, j = super.backboneBonds.length; i < j; i = i+3) {
			this.backboneBonds[i].isLocked = true;
		}
	}
	
	/**
	 * Locks and groups alpha helices.
	 */
	private void lockAndGroupAlphaHelices() {
		if (super.alphaHelix.isEmpty())
			return;
		
		List<Integer> helixBonds = new ArrayList<Integer>();
		
		// lock alpha helices
		for (int i : super.alphaHelix) {
			this.backboneBonds[i].isLocked = true;
			helixBonds.add(i);
		}
		
		// detect alpha helix subsections
		Collections.sort(helixBonds);
		
		int start = helixBonds.get(0);
		int last = start;
				
		for (int i : helixBonds) {
			if (i > last+1) {
				this.groupAndLockSecondaryStructure(start, last);
				start = i;
			}
			
			last = i;
		}
		
		if (last != start) {
			this.groupAndLockSecondaryStructure(start, last);
		}
	}
	
	/**
	 * Locks and groups beta sheets. 
	 */
	private void lockAndGroupBetaSheets() {
		if (super.betaSheet.isEmpty())
			return;
		
		List<Integer> betaSheets = new ArrayList<Integer>();
		
		// lock beta sheets
		for (int i : super.betaSheet) {
			this.backboneBonds[i].isLocked = true;
			betaSheets.add(i);
		}
		
		// detect beta sheet subsections
		Collections.sort(betaSheets);
		
		int start = betaSheets.get(0);
		int last = start;
				
		for (int i : betaSheets) {
			if (i > last+1) {
				this.groupAndLockSecondaryStructure(start, last);
				start = i;
			}
			
			last = i;
		}
		
		if (last != start) {
			this.groupAndLockSecondaryStructure(start, last);
		}
	}
	
	/**
	 * Groups and locks any given secondary structure covering the i-th to
	 * j-th backbone atom.
	 * 
	 * @param i The first backbone atom of the structure.
	 * @param j The last backbone atom of the structure. 
	 */
	private void groupAndLockSecondaryStructure(int i, int j) {
		CTNode node = this.group(i, j);

		this.lockSubtree(node);		
		this.computeTightBoundingVolume(node);
	}
	
	/**
	 * Locks an entire subtree.
	 */
	private void lockSubtree(CTNode node) {
		if (node != null) {
			node.isLocked = true;
			this.lockSubtree(node.left);
			this.lockSubtree(node.right);
		}
	}
	
	/**
	 * Computes a tight bounding box for the leafs in the subtree.
	 * 
	 * @param node The node to compute the tight bounding box for.
	 */
	private void computeTightBoundingVolume(CTNode node) {		
		List<Point3d> points = new LinkedList<Point3d>();

		// translate all points into origo
		TransformationMatrix offset = this.backboneBonds[node.low].transformationMatrix;
		TransformationMatrix transformationMatrix = new TransformationMatrix(-1*offset.a14, -1*offset.a24, -1*offset.a34);
		
		// compute all points in the sub chain
		for (int i = node.low; i <= node.high; i++) {
			transformationMatrix.multR(this.backboneBonds[i].transformationMatrix);
			points.add(new Point3d(transformationMatrix.a14, transformationMatrix.a24, transformationMatrix.a34));
		}
		
		// update
		node.boundingVolume = new LinesegmentSweptSphere(points);
		
		// update subtrees
		if (!node.isLeaf()) {
			this.computeTightBoundingVolume(node.left);
			this.computeTightBoundingVolume(node.right);
		}
	}
	
	/**
	 * Groups the leafs i to j (both included) in their own subtree.
	 * 
	 * @param i The left most leaf.
	 * @param j The right most leaf.
	 * @require i < j
	 * @return The root node of the new subtree.
	 */
	private CTNode group(int i, int j) {
		CTLeaf left = super.backboneBonds[i];
		CTLeaf right = super.backboneBonds[j];
		
		// first go up from left to find the nearest common ancestor 
		CTNode node = left;
		
		while (node.high < j) {
			// if node is the left leaf of its parent then just continue up
			if (node == node.parent.left) {
				node = node.parent;
				
			// else node is the right leaf and we must perform a rotation
			} else {
				if (node.parent.parent.right == node.parent) {
					this.rotateLeft(node.parent.parent);
				} else {
					this.rotateRight(node.parent.parent);
				}
			}
		}
		
		// go down right to rotate the tree to only contain the specified sub chain
		node = right;
		
		while (node.low > i) {
			// if node is the left leaf of its parent then just continue up
			if (node == node.parent.right) {
				node = node.parent;
				
			// else node is the right leaf and we must perform a rotation
			} else {
				if (node.parent.parent.left == node.parent) {
					this.rotateRight(node.parent.parent);
				} else {
					this.rotateLeft(node.parent.parent);
				}
			}
		}

		return node;
	}
	
	/**
	 * Performs a left rotation rooted at the given node in the tree.
	 * http://en.wikipedia.org/wiki/Tree_rotation
	 * 
	 * @param node The node to rotate.
	 */
	private void rotateLeft(CTNode node) {
		CTNode a = node.right;
		CTNode b = node;
		CTNode d = a.left;
	
		// update sub chain information
		a.low = b.low;
		b.high = d.high;
	
		// compute new rotation matrix
		a.transformationMatrix = b.transformationMatrix;
		b.transformationMatrix = new TransformationMatrix(b.left.transformationMatrix, d.transformationMatrix);
		
		// rotate tree nodes
		a.left = node;
		a.parent = node.parent;
		
		if (a.parent == null) {
			super.root = a;
		} else {
			if (b.parent.left == b) {
				b.parent.left = a;
			} else {
				b.parent.right = a;
			}
		}
		
		b.right = d;
		b.parent = a;
		d.parent = b;
		
		// recompute heights
		b.height = Math.max(b.left.height, b.right.height) + 1;
		a.height = Math.max(a.left.height, a.right.height) + 1;
		
		// update nodes above the rotated subtree
		node = a.parent;
		while (node != null) {
			node.height = Math.max(node.left.height, node.right.height) + 1;
			node.update();
			node = node.parent;
		}
	}

	/**
	 * Performs a right rotation rooted at the given node in the tree.
	 * http://en.wikipedia.org/wiki/Tree_rotation
	 * 
	 * @param node The node to rotate.
	 */
	private void rotateRight(CTNode node) {
		CTNode a = node;
		CTNode b = node.left;
		CTNode d = b.right;
	
		// update sub chain information
		b.high = a.high;
		a.low = d.low;
	
		// compute new rotation matrix
		b.transformationMatrix = a.transformationMatrix;
		a.transformationMatrix = new TransformationMatrix(d.transformationMatrix, a.right.transformationMatrix);
		
		// rotate tree nodes
		b.right = a;
		b.parent = a.parent;
		
		if (b.parent == null) {
			super.root = b;
		} else {
			if (node.parent.left == node) {
				node.parent.left = b;
			} else {
				node.parent.right = b;
			}
		}
		
		a.left = d;
		a.parent = b;
		d.parent = a;
		
		// recompute heights
		a.height = Math.max(b.left.height, b.right.height) + 1;
		b.height = Math.max(b.left.height, b.right.height) + 1;
		
		// update nodes above the rotated subtree
		node = b.parent;
		while (node != null) {
			node.height = Math.max(node.left.height, node.right.height) + 1;
			node.update();
			node = node.parent;
		}
	}
}
