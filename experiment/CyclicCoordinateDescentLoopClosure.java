package experiment;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import boundingVolume.Empty;

import algorithm.CyclicCoordinateDescent;

import math.Tuple2;
import math.Tuple3;
import math.Vector3D;
import test.VisualChainTreeDebugger;
import tool.BackboneSegmentAnalyser;
import tool.ChainTreeScene;
import tool.RamachandranDistribution;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;
import energyFunction.EnergyFunction;
import energyFunction.LoopAtomDistance;

public class CyclicCoordinateDescentLoopClosure {
	
	private static enum Restriction {NONE, NEIGHBOUR_INDEPENDENT, NEIGHBOUR_DEPENDENT};
	private static RamachandranDistribution ramachandranDistribution = new RamachandranDistribution();
	
	
	private static ChainTreeScene scene = new ChainTreeScene();
	private static boolean GUI = true;
	
	public static void main(String[] args) throws Exception {
		/*
		 * Configuration.
		 */
		String pdbId = "1PUX"; // 1PUX, 1RKI, 1T0G, 1F3U, 1XJH, 1JN1, 1X6J, 2B7T, 1SIS, 1E2B
		int startSegment = 1;
		int endSegment = 2;
		Restriction restriction = Restriction.NONE;
		double targetRMSDistance = 0.08;
		int maxIterationsPerClose = 5000;
		int numberOfClosedLoops = 10;
		
		
		
		
		/*
		 * Run experiments.
		 */
		closeAllLoops(pdbId, restriction, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
		
		
	}
	
	private static void closeAllSheetLoops (String pdbId, Restriction restriction, double targetRMSDistance, int maxIterationsPerClose, int numberOfClosedLoops) {
		ChainTree cTree = new ChainTree(pdbId);
		List<Tuple2<Integer, Integer>> segments = BackboneSegmentAnalyser.getSheetSegments(cTree);
		
		System.out.println(segments);
		
		for(int i = 0; i < segments.size()-1; i++) {
			Tuple2<Integer, Integer> sheet1 = segments.get(i);
			Tuple2<Integer, Integer> sheet2 = segments.get(i+1);
			
			System.out.println(sheet1 +" - "+ sheet2);
			
			if(sheet2.x - sheet1.y >= 4) {
				System.out.println("Closing " + new Tuple2(sheet1.y, sheet2.x));
				closeLoop(pdbId, sheet1.y, sheet2.x, restriction, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
			}
		}
	}
	
	private static void closeAllLoops(String pdbId, Restriction restriction, double targetRMSDistance, int maxIterationsPerClose, int numberOfClosedLoops) {
		ChainTree cTree = new ChainTree(pdbId);
		List<Tuple2<Integer, Integer>> segments = BackboneSegmentAnalyser.getIntermediateSegments(cTree);
		
		scene.add(cTree,20);
		
		for(int i = 0; i < segments.size(); i++) {
			Tuple2<Integer, Integer> segment = segments.get(i);
			
			if(segment.y - segment.x >= 4) {
				System.out.println("Closing " + segment);
				closeLoop(pdbId, segment.x, segment.y, restriction, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
			}
		}
	}


	private static Tuple3<Collection<Double>, Integer, Integer> closeLoop(String pdbId, int start, int end, Restriction restriction, double targetRMSDistance, int maxIterationsPerClose, int numberOfClosedLoops) {
		// setup chain trees
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		
		AdjustableChainTree cTreeLoop = cTree.getSubchain(1, end);
		AdjustableChainTree cTreeRemainder = cTree.getSubchain(end+1, cTree.length());

		
		// compute energy
		EnergyFunction energyFunction = new LoopAtomDistance(cTreeLoop, start, end);
		cTreeLoop.backboneBonds[cTreeLoop.backboneBonds.length-1].boundingVolume = new Empty();
		
		// find rotateable bonds in the segment
		List<Integer> rotateableBonds = new ArrayList<Integer>();
		
		for (int bond : cTreeLoop.rotatableBonds()) {
			int aminoAcid = cTreeLoop.getAminoAcid(bond);
			
			if (start <= aminoAcid && aminoAcid <= end) {
				rotateableBonds.add(bond);
			}
		}
		
		// phi, psi angles
		List<Double> dihedralAngles = new ArrayList<Double>();

		int k = 0;
		for (double angle : cTree.getDihedralAngles()) {
			if (angle != 0.0 && k % 3 != 2) {
				dihedralAngles.add(angle);
			}
			
			k++;
		}
		
		// init CCD
		CyclicCoordinateDescent anglePredictor = new CyclicCoordinateDescent(cTreeLoop, cTreeLoop.getSubchain(cTreeLoop.length(), cTreeLoop.length()));
		
		// GUI
		ChainTree guiLoop = null;
		if(GUI) {
			guiLoop = cTree.getSubchain(start-1, end+1);
			scene.add(guiLoop);
			scene.scene.centerCamera();
			scene.scene.autoZoom();
		}
		
		/*
		 * Close loops.
		 */
		int itterations = 0;
		int unclosed = 0;
		int clashes = 0;
		Collection<Double> energies = new ArrayList<Double>();
		double minEnergy = Double.MAX_VALUE;
		
		unfold(cTreeLoop, cTreeRemainder, dihedralAngles, start, end);
		
		while(energies.size() < numberOfClosedLoops) {
			for (int i = 0, j = rotateableBonds.size(); i < j; i += 2) { // loop over pairs of phi, psi bonds
				// basic residue information
				int bondPhi = rotateableBonds.get(i);
				int bondPsi = rotateableBonds.get(i+1);
				int aminoAcid = cTreeLoop.getAminoAcid(bondPhi);
				
				if(restriction == Restriction.NONE) {
					// compute new phi, psi angle angles
					double deltaPhi = anglePredictor.getRotationAngle(bondPhi);
					cTreeLoop.changeRotationAngle(bondPhi, deltaPhi);
				
					double deltaPsi = anglePredictor.getRotationAngle(bondPsi);
					cTreeLoop.changeRotationAngle(bondPsi, deltaPsi);
					
				} else {
					// get old phi, psi angles
					List<Double> oldAngles = cTreeLoop.getDihedralAngles(aminoAcid, aminoAcid);
					double oldPhi = oldAngles.get(0);
					double oldPsi = oldAngles.get(1);
				
					// compute new phi, psi angle proposals
					double deltaPhi = anglePredictor.getRotationAngle(bondPhi);
					cTreeLoop.changeRotationAngle(bondPhi, deltaPhi);
				
					double deltaPsi = anglePredictor.getRotationAngle(bondPsi);
					cTreeLoop.changeRotationAngle(bondPsi, deltaPsi);
				
					// get new phi, psi angles
					List<Double> newAngles = cTreeLoop.getDihedralAngles(aminoAcid, aminoAcid);
					double newPhi = newAngles.get(0);
					double newPsi = newAngles.get(1);
				
					// compute the probability of the new conformation
					double oldP;
					double newP;
				
					switch(restriction) {
						case NEIGHBOUR_INDEPENDENT:
							oldP = ramachandranDistribution.probability(cTreeLoop.aminoAcidType(aminoAcid), oldPhi, oldPsi);
							newP = ramachandranDistribution.probability(cTreeLoop.aminoAcidType(aminoAcid), newPhi, newPsi);
							break;
						case NEIGHBOUR_DEPENDENT: 
							oldP = ramachandranDistribution.probability(cTreeLoop.aminoAcidType(aminoAcid), cTreeLoop.aminoAcidType(aminoAcid-1), cTreeLoop.aminoAcidType(aminoAcid+1), oldPhi, oldPsi);
							newP = ramachandranDistribution.probability(cTreeLoop.aminoAcidType(aminoAcid), cTreeLoop.aminoAcidType(aminoAcid-1), cTreeLoop.aminoAcidType(aminoAcid+1), newPhi, newPsi);
							break;
						default:
							oldP = Double.NaN;
							newP = Double.NaN;
					}
				
					// reject the new conformation?
					if(Math.random() > newP / oldP) {
						cTreeLoop.changeRotationAngle(bondPhi, -deltaPhi);
						cTreeLoop.changeRotationAngle(bondPsi, -deltaPsi);
					}
				}
			}
			
			itterations++;
			
			// is loop closed?
			if (anglePredictor.targetRMSDistance() < targetRMSDistance) { 
				
				if (cTreeLoop.isClashing() || cTreeLoop.areClashing(cTreeRemainder)) {
					//System.err.println("Clashing!");
					clashes++;
				} else {
					double energy = energyFunction.compute();
					energies.add(energy);
					
					//System.out.println((int) (100*energies.size()/ ((float)numberOfClosedLoops)) + "%");
					
					if(energy < minEnergy) {
						System.out.println("New conformation "+energies.size()+": " + energy);
						
						if(GUI) {
							scene.remove(guiLoop);
							guiLoop = cTree.getSubchain(start-1, end+1);
							scene.add(guiLoop);
							scene.scene.centerCamera();
							scene.scene.autoZoom();
						}
						
						
						minEnergy = energy;
					}
				}
				
				unfold(cTreeLoop, cTreeRemainder, dihedralAngles, start, end);
				itterations = 0;
				continue;
			}
			
			// to many iterations?
			if(itterations >= maxIterationsPerClose) {
				//System.err.println("To many iterations!");
				unclosed++;
				unfold(cTreeLoop, cTreeRemainder, dihedralAngles, start, end);
				itterations = 0;
			}
		}
		
		System.out.println("Done, best energy: " + minEnergy);
		return new Tuple3<Collection<Double>, Integer, Integer>(energies, clashes, unclosed);
	}
	
	private static void unfold(AdjustableChainTree cTree, AdjustableChainTree other, List<Double> angles, int start, int end) {
		int startBond = cTree.getPhi(start);
		int endBond = cTree.getPsi(end);
		
		do {
			for (int bond : cTree.rotatableBonds()) {
				if (startBond <= bond && bond <= endBond) {
					double angle = angles.get((int) (Math.random() * angles.size()));
					cTree.changeRotationAngle(bond, angle);			
				}
			}
		} while (cTree.isClashing() || cTree.areClashing(other));
	}
}
