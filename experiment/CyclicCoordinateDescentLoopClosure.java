package experiment;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import boundingVolume.Empty;

import algorithm.CyclicCoordinateDescent;

import math.Tuple2;
import math.Vector3D;
import tool.BackboneSegmentAnalyser;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;

public class CyclicCoordinateDescentLoopClosure {
	public static void main(String[] args) throws Exception {
		/*
		 * Configuration.
		 */
		String pdbId = "1SIS"; // 1PUX, 1RKI, 1T0G, 1F3U, 1XJH, 1JN1, 1X6J, 2B7T, 1SIS, 1E2B
		int segmentNo = 1;
		double targetRMSDistance = 0.08;
		int maxIterations = 5000;
		

		/*
		 * Setup.
		 */
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		cTree.rotate(3);
		
		Tuple2<Integer, Integer> segment = BackboneSegmentAnalyser.extractIntermediateSegments(cTree).get(segmentNo);
		int start = segment.x;
		int end = segment.y;

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

		unfold(t1, dihedralAngles, start, end);		
		scene.repaint();

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
		for (int bond : cTree.rotatableBonds()) {
			if (start <= bond && bond <= end) {
				double angle = angles.get((int) (Math.random() * angles.size()));
				
				do {
					cTree.changeRotationAngle(bond, angle);

					angle -= Math.PI / 100;
				} while(cTree.isClashing());
			}
		}
	}
}
