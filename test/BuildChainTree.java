package test;

import java.io.IOException;

import tool.ChainTreeScene;

import dataStructure.ChainTree;

public class BuildChainTree {

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws InterruptedException, IOException {

		ChainTree cTree = new ChainTree("1PUX");

		ChainTreeScene scene = new ChainTreeScene(cTree);
				
		while(true) {
			int i = (int) (Math.random() * cTree.length());
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
