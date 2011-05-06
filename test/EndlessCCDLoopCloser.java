package test;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import algorithm.CyclicCoordinateDescent;

import math.Tuple2;
import math.Vector3D;
import tool.BackboneSegmentAnalyser;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;

public class EndlessCCDLoopCloser {
	public static void main(String[] args) throws Exception {
		/*
		 * Configuration.
		 */
		String pdbId = "1SIS"; // 1PUX, 1RKI, 1T0G, 1F3U, 1XJH, 1JN1, 1X6J, 2B7T, 1SIS, 1E2B
		int segmentNo = 1;
		double targetRMSDistance = 0.04;
		int maxIterations = 10000;
		

		/*
		 * Setup.
		 */
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		cTree.rotate(3);
		
		Tuple2<Integer, Integer> segment = BackboneSegmentAnalyser.extractIntermediateSegments(cTree).get(segmentNo);
		int start = segment.x;
		int end = segment.y;

		// create subtrees
		AdjustableChainTree t0 = cTree.getSubchain(0, start);
		AdjustableChainTree t1 = cTree.getSubchain(0, end-1);
		AdjustableChainTree t2 = cTree.getSubchain(end+1, cTree.length()-1);
		
		ChainTree[] cTrees = {t0, t2}; 
		ChainTreeScene scene = new ChainTreeScene(cTrees);
		
		// find rotateable bonds in the segment
		List<Integer> rotateableBonds = new ArrayList<Integer>();
		
		for (int i : t1.rotatableBonds()) {
			if (start <= i && i <= end) {
				rotateableBonds.add(i);
			}
		}
		
		// set up CCD algorithm
		ChainTree target = cTree.getSubchain(end+1, end+2);		
		ChainTree testLoop = cTree.getSubchain(start-1, end+2);
		testLoop.unfold();
		scene.add(testLoop, 100);
		
		CyclicCoordinateDescent anglePredictor = new CyclicCoordinateDescent(testLoop, target);

		/*
		 * Close loop.
		 */
		int itt = 0;
		
		while (true) {			
			for (int bond : rotateableBonds) {
				double angle = anglePredictor.getRotationAngle(bond-start+1);
				
				testLoop.changeRotationAngle(bond-start+1, angle);
			
				scene.repaint(testLoop);
			}
			
			if (anglePredictor.targetRMSDistance() < targetRMSDistance) { 
				System.out.println("Loop closed in " + itt + " iterations!");
				itt = 0;
				Thread.sleep(500);
				testLoop.unfold();
				scene.repaint(testLoop);
				Thread.sleep(500);
			} else {
				itt++;
			
				if (itt >= maxIterations) {
					itt = 0;
					System.err.println("TO MANY ITERATIONS!");
					testLoop.unfold();
				}
			}
		}
	}
}
