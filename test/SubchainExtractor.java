package test;


import java.util.List;

import javax.vecmath.Point3d;

import math.Vector3D;


import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;

public class SubchainExtractor {
	public static void main(String[] args) throws InterruptedException {
		
		String pdbId = "1PUX";	// protein id
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		
		ChainTreeScene scene = new ChainTreeScene();
		
		AdjustableChainTree t1 = cTree.getSubchain(0,7);
		AdjustableChainTree t2 = cTree.getSubchain(9, 16);
		AdjustableChainTree t12 = cTree.getSubchain(0, 16);

		t1.move(new Vector3D(0.0,0.0,-3.0));
		
		scene.add(t1);
		scene.add(t2);
		
		t12.move(new Vector3D(3.0,3.0,3.0));
		scene.add(t12);
		
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		
		// toy rotation
		List<Integer> rotateableBonds = t1.rotatableBonds();
		int direction = 1;
		int bond = 0;
		
		while(true) {
			t1.changeRotationAngle(rotateableBonds.get(bond), direction);
			scene.repaint();
			Thread.sleep(500);

			if (t1.areClashing(t12)) {
				if (bond + 1 == rotateableBonds.size()) {
					direction *= -1;
				}
				
				bond += direction;
			}
		}
	}
}
