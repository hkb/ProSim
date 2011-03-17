package test;

import geom3d.Capsule3d;
import geom3d.Point3d;
import geom3d.Shape3d;

import java.awt.Color;
import java.io.IOException;
import java.util.LinkedList;

import math.matrix.TransformationMatrix;

import boundingVolume.LinesegmentSweptSphere;

import tool.BinaryTreePainter;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.CTNode;
import dataStructure.ChainTree;
import edu.geom3D.Capsule;

public class AdjustableChainTreeDebugger {
	public static void main(String[] args) throws InterruptedException, IOException {

		String pdbId = "1PUX";	// protein id
		int height = 0;			// paint bounding volumes at height x
		boolean locked = true;	// painting bounding volumes of locked trees
		
		
		
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		
		new BinaryTreePainter(cTree);

		
		ChainTreeScene scene = new ChainTreeScene(cTree);
		paintBoundingVolume(cTree, cTree.getRoot(), scene, height, locked);
	}
	
	private static void paintBoundingVolume(ChainTree cTree, CTNode node, ChainTreeScene scene, int height, boolean locked) {
		if (!node.isLeaf()) {
			if (node.height == height || height == -1 || (locked && node.isLocked)) {
				TransformationMatrix w = new TransformationMatrix(cTree.position.x, cTree.position.y, cTree.position.z);
				TransformationMatrix t = cTree.getTransformationMatrix(0, node.low);
				t.multR(w);
				
				Capsule vol = ((LinesegmentSweptSphere) node.boundingVolume.transform(t)).volume;

				Color color = (node.isLocked) ? new Color(255,0,0,100) : new Color(0,0,255,100);
				scene.scene.addShape(new Capsule3d(new Point3d(vol.p1.x(), vol.p1.y(), vol.p1.z()), new Point3d(vol.p2.x(), vol.p2.y(), vol.p2.z()), vol.rad), color);
			}

			if (!(locked && node.isLocked)) {
				paintBoundingVolume(cTree, node.left, scene, height, locked);
				paintBoundingVolume(cTree, node.right, scene, height, locked);
			}
		}
	}
}