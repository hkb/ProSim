package test;

import java.io.IOException;
import java.util.List;

import tool.ChainTreeScene;

import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;

public class BuildChainTree {

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws InterruptedException, IOException {

		AdjustableChainTree cTree = new AdjustableChainTree("1PUX");
		List<Integer> rotateableBonds = cTree.rotatableBonds(); 

		ChainTreeScene scene = new ChainTreeScene(cTree);
		
		
		while(true) {
			int i = rotateableBonds.get((int) (Math.random() * rotateableBonds.size()));
			double angle = (Math.random()-0.5)*15*(Math.PI/180);
			
			System.out.println(i + " - " + angle);
			
			cTree.changeRotationAngle(i, angle);
			
			if(cTree.isClashing()) {
				System.err.println("Self clash!");
				cTree.changeRotationAngle(i, -angle);				
			} else {
				scene.repaint();
			}
		}
	}
}
