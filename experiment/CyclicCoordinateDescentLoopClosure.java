package experiment;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

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
	
	// this is sick!
	private static List<Tuple2<Tuple3<String,Integer,Integer>,Collection<Collection<Tuple2<Integer,Double>>>>> conformationCache = new LinkedList<Tuple2<Tuple3<String,Integer,Integer>,Collection<Collection<Tuple2<Integer,Double>>>>>();
	
	public static void main(String[] args) throws Exception {
		/*
		 * Configuration.
		 */
		String pdbId = "1PUX"; // 1PUX, 1RKI, 1T0G, 1F3U, 1XJH, 1JN1, 1X6J, 2B7T, 1SIS, 1E2B
		Restriction restriction = Restriction.NONE;
		double targetRMSDistance = 0.08;
		int maxIterationsPerClose = 5000;
		int numberOfClosedLoops = 1000;
		
		
		
		
		/*
		 * Run experiments.
		 */
		System.out.println("No constraint:");
		closeAllSheetLoops(pdbId, Restriction.NONE, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
		scene.clear();
		System.out.println("Nighbour independend:");
		closeAllSheetLoops(pdbId, Restriction.NEIGHBOUR_INDEPENDENT, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
		scene.clear();
		System.out.println("Nighbour dependend:");
		closeAllSheetLoops(pdbId, Restriction.NEIGHBOUR_DEPENDENT, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
		
	}
	
	private static void closeAllSheetLoops (String pdbId, Restriction restriction, double targetRMSDistance, int maxIterationsPerClose, int numberOfClosedLoops) {
		ChainTree cTree = new ChainTree(pdbId);
		List<Tuple2<Integer, Integer>> segments = BackboneSegmentAnalyser.getSheetSegments(cTree);
		
		scene.add(cTree.getSubchain(1, segments.get(0).x));
		for(Tuple2<Integer, Integer> segment : segments) {
			scene.add(cTree.getSubchain(segment.x, segment.y));
		}
		scene.add(cTree.getSubchain(segments.get(segments.size()-1).y, cTree.length()));
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		
		for(int i = 0; i < segments.size()-1; i++) {
			Tuple2<Integer, Integer> sheet1 = segments.get(i);
			Tuple2<Integer, Integer> sheet2 = segments.get(i+1);
						
			if(sheet2.x - sheet1.y >= 4) {
				closeLoop(pdbId, sheet1.y+1, sheet2.x-1, restriction, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
			}
		}
		
		scene.add(cTree,20);
	}
	
	private static void closeAllLoops(String pdbId, Restriction restriction, double targetRMSDistance, int maxIterationsPerClose, int numberOfClosedLoops) {
		ChainTree cTree = new ChainTree(pdbId);
		List<Tuple2<Integer, Integer>> segments = BackboneSegmentAnalyser.getIntermediateSegments(cTree);
		
		scene.add(cTree,20);
		
		for(int i = 0; i < segments.size(); i++) {
			Tuple2<Integer, Integer> segment = segments.get(i);
			
			if(segment.y - segment.x >= 4) {	
				closeLoop(pdbId, segment.x, segment.y, restriction, targetRMSDistance, maxIterationsPerClose, numberOfClosedLoops);
			}
		}
	}


	private static Tuple3<Collection<Double>, Integer, Integer> closeLoop(String pdbId, int start, int end, Restriction restriction, double targetRMSDistance, int maxIterationsPerClose, int numberOfClosedLoops) {
		// setup chain trees
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		
		AdjustableChainTree cTreeLoop = cTree.getSubchain(1, end+1);
		AdjustableChainTree cTreeRemainder = cTree.getSubchain(end+2, cTree.length());
		
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
		
		// generate unfold conformations
		Tuple3<String,Integer,Integer> conformationKey = new Tuple3<String,Integer,Integer>(pdbId, start, end);
		Collection<Collection<Tuple2<Integer,Double>>> conformations = null;
		
		for(Tuple2<Tuple3<String,Integer,Integer>,Collection<Collection<Tuple2<Integer,Double>>>> cs : conformationCache) {
			if(cs.x.equals(conformationKey))
				conformations = cs.y;
		}
		if(conformations == null) {
			conformations = generateConformations(cTreeLoop.getSubchain(1, cTreeLoop.length()), cTreeRemainder, dihedralAngles, start, end, numberOfClosedLoops);
			conformationCache.add(new Tuple2<Tuple3<String,Integer,Integer>,Collection<Collection<Tuple2<Integer,Double>>>>(conformationKey, conformations));
		}
		
		// init CCD
		CyclicCoordinateDescent anglePredictor = new CyclicCoordinateDescent(cTreeLoop, cTreeLoop.getSubchain(cTreeLoop.length(), cTreeLoop.length()));
		
		// GUI
		ChainTree guiLoop = null;
		
		/*
		 * Close loops.
		 */
		int itterations = 0;
		int unclosed = 0;
		int clashes = 0;
		Collection<Double> energies = new ArrayList<Double>();
		double minEnergy = Double.MAX_VALUE;
				
		for(Collection<Tuple2<Integer,Double>> conformation : conformations) {
			itterations = 0;
			unforldIntoConformation(cTreeLoop, conformation);
			
			while(true) {
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
						clashes++;
					} else {
						double energy = energyFunction.compute();
						energies.add(energy);
						
						if(energy < minEnergy) {
	
							if(guiLoop != null)
								scene.remove(guiLoop);
								guiLoop = cTreeLoop.getSubchain(start-1, end+1);
								scene.add(guiLoop);						
							
							minEnergy = energy;
						}
					}
					
					break;
				}
				
				// to many iterations?
				if(itterations >= maxIterationsPerClose) {
					unclosed++;
					break;
				}
			}
		}
		
		System.out.println(conformationKey + " done, best energy: " + minEnergy);
		return new Tuple3<Collection<Double>, Integer, Integer>(energies, clashes, unclosed);
	}
	
	private static Collection<Tuple2<Integer,Double>> unfold(AdjustableChainTree cTree, AdjustableChainTree other, List<Double> angles, int start, int end) {
		int startBond = cTree.getPhi(start);
		int endBond = cTree.getPsi(end);
		
		Collection<Tuple2<Integer,Double>> conformation;
		
		do {
			conformation = new LinkedList<Tuple2<Integer,Double>>();
			
			for (int bond : cTree.rotatableBonds()) {
				if (startBond <= bond && bond <= endBond) {
					double angle = angles.get((int) (Math.random() * angles.size()));
					cTree.changeRotationAngle(bond, angle);	
					
					conformation.add(new Tuple2<Integer,Double>(bond, angle));
				}
			}
		} while (cTree.isClashing() || cTree.areClashing(other));
		
		return conformation;
	}
	
	private static Collection<Collection<Tuple2<Integer,Double>>> generateConformations(AdjustableChainTree cTree, AdjustableChainTree other, List<Double> angles, int start, int end, int n) {
		List<Collection<Tuple2<Integer,Double>>> conformations = new ArrayList<Collection<Tuple2<Integer,Double>>>();
		
		for(int i = 0; i < n; i++) {
			conformations.add(unfold(cTree, other, angles, start, end));
		}
		
		return conformations;
	}
	
	private static void unforldIntoConformation(ChainTree cTree, Collection<Tuple2<Integer,Double>> conformation) {
		for(Tuple2<Integer,Double> bondInfo : conformation) {
			cTree.changeRotationAngle(bondInfo.x, bondInfo.y);
		}
	}
}
