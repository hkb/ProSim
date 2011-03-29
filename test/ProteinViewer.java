package test;

import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;

public class ProteinViewer {
	public static void main(String[] args) {
		String pdbId = "1JN1";
		
		new ChainTreeScene(new AdjustableChainTree(pdbId));
	}
}
