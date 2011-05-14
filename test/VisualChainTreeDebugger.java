package test;

import java.awt.Color;
import java.io.File;

import javax.media.j3d.BadTransformException;

import boundingVolume.LinesegmentSweptSphere;
import geom3d.Capsule3d;
import geom3d.Point3d;
import tool.BinaryTreePainter;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.CTNode;
import dataStructure.ChainTree;
import edu.geom3D.Capsule;

public class VisualChainTreeDebugger {
	public static void main(String[] args) {
		File[] listOfFiles = new File("/home/hkb/.workspace/ProSim/pdb_files").listFiles();
		String pdbId = listOfFiles[(int)(Math.random() * listOfFiles.length)].getName().substring(0,4);
		System.out.println(pdbId);
		//String pdbId = "1JB0";
		
		
		//
		ChainTree cTree = new AdjustableChainTree(pdbId);
		
		
		//showBinaryTreeStructure(cTree);
		showProtein(cTree);
		//showBoundingVolumes(cTree);
		
		// paint single volume
		//paintBoundingVolume(cTree, cTree.backboneBonds[1].parent, new ChainTreeScene(cTree));
	}
	
	public static void showBinaryTreeStructure(ChainTree cTree) {
		new BinaryTreePainter(cTree);	
	}
	
	public static void showProtein(ChainTree cTree) {
		new ChainTreeScene(cTree);
	}
	
	public static void showBoundingVolumes(ChainTree cTree) {
		ChainTreeScene scene = new ChainTreeScene(cTree);
		
		// paint leaf bounding volumes
		for (CTNode node : cTree.backboneBonds) {
			paintBoundingVolume(cTree, node, scene);
		}
		
		// secondary structures
		paintSecondaryStructureBoundingVolumes(cTree, cTree.root, scene);
		
		// paint root
		//paintBoundingVolume(cTree, cTree.root, scene);
	}
	
	
	private static void paintSecondaryStructureBoundingVolumes(ChainTree cTree, CTNode node, ChainTreeScene scene) {
		if (node.isLeaf())
			return;
		
		if (node.isLocked) {
			paintBoundingVolume(cTree, node, scene);
			return;
		}
		
		paintSecondaryStructureBoundingVolumes(cTree, node.left, scene);
		paintSecondaryStructureBoundingVolumes(cTree, node.right, scene);
	}
	
	public static void paintBoundingVolume(ChainTree cTree, CTNode node, ChainTreeScene scene) {
		try {
			Capsule vol = ((LinesegmentSweptSphere) node.boundingVolume.transform(cTree.getWorldTransformation(node.low))).volume;

			Color color = (node.isLocked) ? new Color(255,0,0, 150) : new Color(0,0,255, 150);
			scene.scene.addShape(new Capsule3d(new Point3d(vol.p1.x(), vol.p1.y(), vol.p1.z()), new Point3d(vol.p2.x(), vol.p2.y(), vol.p2.z()), vol.rad), color);
		} catch (Exception e) {
			System.err.println(e);
		}
	}
}
