package dataStructure;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

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
		// lock alpha helices
		for (int i : super.alphaHelix) {
			this.backboneBonds[i].isLocked = true;
		}
	}
	
	/**
	 * Locks and groups beta sheets. 
	 */
	private void lockAndGroupBetaSheets() {
		// lock beta sheets
		for (int i : super.betaSheet) {
			this.backboneBonds[i].isLocked = true;
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
	 * Groups the leafs i to j (both included) in their own subtree.
	 * 
	 * @param i The left most leaf.
	 * @param j The right most leaf.
	 * @require i < j
	 * @return The root node of the new subtree.
	 */
	public CTNode group(int i, int j) {
		CTLeaf left = super.backboneBonds[i];
		CTLeaf right = super.backboneBonds[j];
		
		// first go up from left to find the nearest common ancestor 
		CTNode node = left.parent;
		
		while (node.high < j) {
			System.out.println(node);
			
			// if node is the left leaf of its parent then just continue up
			if (node == node.parent.left) {
				node = node.parent;
				
			// else node is the right leaf and we must perform a rotation
			} else {
				if (node.parent.right == node) {
					System.out.println("left");
					this.rotateLeft(node);
				} else {
					System.out.println("right");
					this.rotateRight(node);
				}
			}
		}
		
		return node;
	}
	
	/**
	 * Performs a left rotation on the given node in the tree.
	 * http://en.wikipedia.org/wiki/Tree_rotation
	 * 
	 * @param node The node to rotate.
	 */
	private void rotateLeft(CTNode node) {
		CTNode Q = node;
		CTNode R = Q.parent.parent;
		CTNode P = Q.parent;
		CTNode B = Q.left;
		
		// update root
		if (R == null) {
			super.root = Q;
			Q.parent = null;
		} else if (R.left == P) {
			R.setLeft(Q);
		} else {
			R.setRight(Q);
		}

		// update leafs
		P.setRight(B);
		Q.setLeft(P);
	}

	/**
	 * Performs a right rotation on the given node in the tree.
	 * http://en.wikipedia.org/wiki/Tree_rotation
	 * 
	 * @param node The node to rotate.
	 */
	private void rotateRight(CTNode node) {
		CTNode P = node;
		CTNode R = P.parent.parent;
		CTNode Q = P.parent;
		CTNode B = P.right;
		
		// update root
		if (R == null) {
			super.root = P;
			P.parent = null;
		} else if (R.left == Q) {
			R.setLeft(P);
		} else {
			R.setRight(P);
		}
		
		// update leafs
		P.setRight(Q);
		Q.setLeft(B);
	}
}
