package test;


import javax.vecmath.Point3d;


import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;

public class SubchainExtractor {
	public static void main(String[] args) throws InterruptedException {
		
		String pdbId = "1PUX";	// protein id
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		
		ChainTreeScene scene = new ChainTreeScene();
		
		AdjustableChainTree t1 = cTree.getSubchain(0,8);
		AdjustableChainTree t2 = cTree.getSubchain(9, 17);
		AdjustableChainTree t12 = cTree.getSubchain(0, 17);

		t1.move(new Point3d(0,0,-3));
		
		scene.add(t1);
		scene.add(t2);
		
		t12.move(new Point3d(3,3,3));
		scene.add(t12);
		
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		
		while(true) {
			t1.changeRotationAngle(0, 1);
			scene.repaint();

			if (t1.areClashing(t12)) {
				AdjustableChainTreeDebugger.paintVolume(t1, t1.c1, scene);
				AdjustableChainTreeDebugger.paintVolume(t12, t1.c2, scene);
				break;
			}
		}
	}
}
