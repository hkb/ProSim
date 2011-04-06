package test;

import geom3d.Point3d;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;

public class ProteinViewer {
	public static void main(String[] args) {
		String pdbId = "2B7T";
		
		new ChainTreeScene(new AdjustableChainTree(pdbId));
	}
}
