package test;

import java.io.IOException;

import tool.BinaryTreePainter;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;

public class AdjustableChainTreeProperties {
	public static void main(String[] args) throws InterruptedException, IOException {

		String pdbId = "1TMZ";
		
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		
		cTree.group(5, 10);

		new BinaryTreePainter(cTree);
	}
}
