package test;

import java.util.LinkedList;

import javax.vecmath.Point3d;

import tool.BinaryTreePainter;
import tool.ChainTreeScene;

import dataStructure.ChainTree;

public class BuildChainTree {

	/**
	 * @param args
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws InterruptedException {

		ChainTree cTree = new ChainTree("1PUX");

		//ChainTreeScene scene = new ChainTreeScene(cTree);
				
		while(true) {
			Thread.sleep(500);

			int i = (int) (Math.random() * cTree.length());
			double angle = (Math.random()-0.5)*15*(Math.PI/180);
			
			System.out.println(i + " - " + angle);
			
			cTree.changeRotationAngle(i, angle);
			
			if(cTree.isClashing()) {
				System.err.println("Self clash!");
				System.exit(0);
				cTree.changeRotationAngle(i, -angle);
				
			} else {
				//scene.repaint();
			}
		}
	}
}
