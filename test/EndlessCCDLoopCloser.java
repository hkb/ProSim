package test;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import boundingVolume.Empty;

import algorithm.CyclicCoordinateDescent;

import math.Tuple2;
import math.Vector3D;
import tool.BackboneSegmentAnalyser;
import tool.ChainTreeScene;
import tool.NeighborDependentRamachandranDistribution;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;

public class EndlessCCDLoopCloser {
	public static void main(String[] args) throws Exception {
		/*
		 * Configuration.
		 */
		String pdbId = "1SIS"; // 1PUX, 1RKI, 1T0G, 1F3U, 1XJH, 1JN1, 1X6J, 2B7T, 1SIS, 1E2B
		int segmentNo = 1;
		boolean closeTheHelix = true;
		double targetRMSDistance = 0.04;
		int maxIterations = 10000;


		/*
		 * Setup.
		 */
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		cTree.rotate(3);
		
		// find segment endpoints
		int start;
		int end;
		if(closeTheHelix) {
			Tuple2<Integer, Integer> segment1 = BackboneSegmentAnalyser.extractIntermediateSegments(cTree).get(segmentNo);
			Tuple2<Integer, Integer> segment2 = BackboneSegmentAnalyser.extractIntermediateSegments(cTree).get(segmentNo+1);
			start = segment1.x;
			end = segment2.y;
		} else {
			Tuple2<Integer, Integer> segment = BackboneSegmentAnalyser.extractIntermediateSegments(cTree).get(segmentNo);
			start = segment.x;
			end = segment.y;
		}

		// create subtrees
		AdjustableChainTree t1 = cTree.getSubchain(0, end-1+3);
		AdjustableChainTree t2 = cTree.getSubchain(end+1+3, cTree.length()-1);
		
		ChainTree[] cTrees = {t1, t2}; 
		ChainTreeScene scene = new ChainTreeScene(cTrees);
		
		// phi, psi angles
		List<Double> dihedralAngles = new ArrayList<Double>();

		int i = 0;
		for (double angle : cTree.getDihedralAngles()) {
			if (angle != 0.0 && i % 3 != 2) {
				dihedralAngles.add(angle);
			}
			
			i++;
		}
		
		// avoid false posetive self clash
		t1.backboneBonds[t1.backboneBonds.length-1].boundingVolume = new Empty();

		// find rotateable bonds in the segment
		List<Integer> rotateableBonds = new ArrayList<Integer>();
		
		for (int bond : t1.rotatableBonds()) {
			if (start <= bond && bond <= end) {
				rotateableBonds.add(bond);
			}
		}
		
		// set up CCD algorithm
		CyclicCoordinateDescent anglePredictor = new CyclicCoordinateDescent(t1, t1.getSubchain(t1.length()-2, t1.length()-1));
		scene.add(cTree.getSubchain(end+1, end+5), 20);

		// ready
		Thread.sleep(5000);
		unfold(t1, dihedralAngles, start, end);		
		scene.repaint();
		Thread.sleep(2000);

		/*
		 * Close loop.
		 */
		int itt = 0;
		
		while (true) {			
			for (int bond : rotateableBonds) {				
				double angle = anglePredictor.getRotationAngle(bond);

				t1.changeRotationAngle(bond, angle);
			
				scene.repaint(t1);

				if (anglePredictor.targetRMSDistance() < targetRMSDistance) { 
					System.out.println("Loop closed in " + itt + " iterations!");
					
					if (t1.isClashing() || t1.areClashing(t2)) {
						Thread.sleep(1);
						System.err.println("Conformation is clashing!");
					}
					
					itt = 0;
					Thread.sleep(500);
					unfold(t1, dihedralAngles, start, end);
					scene.repaint(t1);
					Thread.sleep(500);
				}	
			}
			
			itt++;
			
			if (itt > maxIterations) {
				System.err.println("TO MANY ITERATIONS!");
				
				itt = 0;
				Thread.sleep(500);
				unfold(t1, dihedralAngles, start, end);
				scene.repaint(t1);
				Thread.sleep(500);
			}
		}
	}
	
	/**
	 * 
	 * @param cTree
	 */
	private static void unfold(ChainTree cTree, List<Double> angles, int start, int end) {
		do {
			for (int bond : cTree.rotatableBonds()) {
				if (start <= bond && bond <= end) {
					double angle = angles.get((int) (Math.random() * angles.size()));
					cTree.changeRotationAngle(bond, angle);			
				}
			}
		} while (cTree.isClashing());
	}
}
