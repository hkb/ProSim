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
	
	private List<CTNode> lockedSubtrees = new LinkedList<CTNode>();
	
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
		this.rebalance();
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
		this.rebalance();
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
		
		this.lockedSubtrees.add(node);
	}
	
	/**
	 * 
	 */
	private void toggleLockedSubtreeVisibility(boolean visible) {
		for (CTNode node : this.lockedSubtrees) {
			// toggle this trees hight
			node.height = (visible) ? Math.max(node.left.height, node.right.height) + 1: 0; 
			
			// update ancestors
			node = node.parent;
			while (node != null) {
				node.height = Math.max(node.left.height, node.right.height) + 1;
				node = node.parent;
			}
		}
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

		// compute all points in the sub chain in the coordinate system of node.low
		TransformationMatrix transformationMatrix = new TransformationMatrix();

		for (int i = node.low; i <= node.high; i++) {
			points.add(new Point3d(transformationMatrix.a14, transformationMatrix.a24, transformationMatrix.a34));
			transformationMatrix.multR(this.backboneBonds[i].transformationMatrix);
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
	 * Re-balances the entire tree.
	 */
	private void rebalance() {
		this.toggleLockedSubtreeVisibility(false);
		this.rebalanceSubtree(this.root);
		this.toggleLockedSubtreeVisibility(true);
	}
	
	/**
	 * Re-balances some subtree.
	 * 
	 * @param node The root node of the subtree to be re-balanced.
	 * @require node.height >= 3
	 */
	private void rebalanceSubtree(CTNode node) {
		// only re-balance trees with minimum height at 3
		if (node.height >= 3) {
		
			// subtrees must be re-balanced first
			this.rebalanceSubtree(node.left);
			this.rebalanceSubtree(node.right);
		
			// then re-balance root
			this.rebalanceNode(node);
		}
	}
	
	private void rebalanceNode(CTNode node) {
		CTNode a = node;
		CTNode b = a.left;
		CTNode c = a.right;
		
		// left tree (b) is more that 1 higher so we must re-balance
		if (b.height - c.height > 1) {
			CTNode d = b.left;
			CTNode e = b.right;
			
			if (d.height >= e.height) {
				// perform double rotation
				this.rotateRight(a);
				
				this.rebalanceNode(a);
				this.rebalanceNode(b);
			} else {
				// perform single rotation
				this.rotateLeft(b);
				this.rotateRight(a);
				
				this.rebalanceNode(a);
				this.rebalanceNode(e);
			}
		}
		
		// right tree (c) is more than 1 higher so we must re-balance
		if (c.height - b.height > 1) {
			CTNode d = c.left;
			CTNode e = c.right;
			
			if (e.height >= d.height) {
				// perform double rotation
				this.rotateLeft(a);
				
				this.rebalanceNode(a);
				this.rebalanceNode(c);
			} else {
				// perform single rotation
				this.rotateRight(c);
				this.rotateLeft(a);
				
				this.rebalanceNode(a);
				this.rebalanceNode(d);
			}
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
		a.left = b;
		a.parent = b.parent;
		
		if (a.parent == null) {
			super.root = a;
		} else {
			if (b.parent.left == b) {
				b.parent.left = a;
			} else {
				b.parent.right = a;
			}node.update();
		}
		
		b.right = d;
		b.parent = a;
		d.parent = b;
		
		// update nodes above the rotated subtree
		// starts from d as both a and b are ancestors to d
		node = d.parent;
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
		
		// update nodes above the rotated subtree
		// starts from d as both a and b are ancestors to d
		node = d.parent;
		while (node != null) {
			node.height = Math.max(node.left.height, node.right.height) + 1;
			node.update();
			node = node.parent;
		}
	}
}
