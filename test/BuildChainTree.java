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

		System.out.println("Creating molecule");
		ChainTree cTree = new ChainTree("1PUX");
		System.out.println("Painting molecule");
		new BinaryTreePainter(cTree);
/*
		ChainTreeScene scene = new ChainTreeScene(cTree);
		
		while(true) {
			Thread.sleep(100000);
			cTree.changeRotationAngle(3, Math.random()*180);
			scene.repaint();
		}
		*/
	}
}
