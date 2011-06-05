package experiment;

import java.awt.Color;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import algorithm.CyclicCoordinateDescent;
import boundingVolume.Empty;

import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;
import energyFunction.EnergyFunction;
import energyFunction.LoopAtomDistance;

import math.Tuple2;

import tool.ChainTreeScene;
import tool.RamachandranDistribution;

public class CyclicCoordinateDescentAbInitioFolder {

	private static RamachandranDistribution ramachandranDistribution = new RamachandranDistribution();
	private static ChainTreeScene scene = new ChainTreeScene();

	private static double TARGET_RMSD = 0.08;
	private static double TARGET_LOOP_RMSD = 1.0;
	private static int MAX_ITERATIONS_PER_CLOSE = 5000;
	
	private static String[] PDBs = {"1PUX"};  
		
	public static void main(String[] args) throws Exception {
		
		
		for(String pdb : PDBs) {
			AdjustableChainTree cTree = new AdjustableChainTree(pdb);
			Set<Integer> rotateableBonds = new HashSet<Integer>(cTree.rotatableBonds());
			List<Tuple2<Integer,Integer>> segments = cTree.getSheetSegments();
			
			// GUI draw base
			scene.add(cTree, Color.DARK_GRAY, 100, 0.1);
			scene.add(cTree.getSubchain(1, segments.get(0).x));
			scene.add(cTree.getSubchain(segments.get(segments.size()-1).y, cTree.length()));
			
			for(Tuple2<Integer,Integer> segment : segments) {
				scene.add(cTree.getSubchain(segment.x, segment.y));
			}
			
			scene.scene.centerCamera();
			scene.scene.autoZoom();
			
			
			Tuple2<Integer,Integer> previousSegment = null;
			
			for(Tuple2<Integer,Integer> segment : segments) {
				if(previousSegment == null) {
					previousSegment = segment;
					continue;
				}
				
				// init tree
				int start = previousSegment.y+1;
				int end = segment.x-1;
				
				AdjustableChainTree cTreeLoop = cTree.getSubchain(1, end+1);
				AdjustableChainTree cTreeRemainder = cTree.getSubchain(end+2, cTree.length());
				cTreeLoop.backboneBonds[cTreeLoop.backboneBonds.length-1].boundingVolume = new Empty();
				
				// compute energy
				EnergyFunction energyFunction = new LoopAtomDistance(cTreeLoop, start, end);
				double minEnergy = Double.MAX_VALUE;
				
				// init ccd
				CyclicCoordinateDescent anglePredictor = new CyclicCoordinateDescent(cTreeLoop, cTreeLoop.getSubchain(cTreeLoop.length(), cTreeLoop.length()));
				
				// GUI
				ChainTree guiLoop = null;

				int i = 0;
				
				// close loop
				while(minEnergy > TARGET_LOOP_RMSD) {
					// unfold
					for (int residue = start; residue <= end; residue++) {
						if(rotateableBonds.contains(cTreeLoop.getPhi(residue))) {
							if(!rotateableBonds.contains(cTreeLoop.getPsi(residue))) {
								throw new IllegalArgumentException("Both phi and psi angles must be rotateable!");
							}
							
							Tuple2<Double,Double> angles = ramachandranDistribution.purposeAngle(cTreeLoop.getAminoAcidType(cTreeLoop.getAminoAcid(residue)));
							
							cTreeLoop.setRotationAngle(cTreeLoop.getPhi(residue), angles.x);
							cTreeLoop.setRotationAngle(cTreeLoop.getPsi(residue), angles.x);
						}
					}
					
					// close the loop
					for(int j = 0; j < MAX_ITERATIONS_PER_CLOSE; j++) {
						for (int residue = start; residue <= end; residue++) {
							// basic residue information
							int bondPhi = cTreeLoop.getPhi(residue);
							int bondPsi = cTreeLoop.getPsi(residue);
							int aminoAcid = cTreeLoop.getAminoAcid(bondPhi);
							
							if(rotateableBonds.contains(bondPhi)) {
								if(!rotateableBonds.contains(bondPsi)) {
									throw new IllegalArgumentException("Both phi and psi angles must be rotateable!");
								}

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
								double oldP = ramachandranDistribution.probability(cTreeLoop.getAminoAcidType(aminoAcid), oldPhi, oldPsi);
								double newP = ramachandranDistribution.probability(cTreeLoop.getAminoAcidType(aminoAcid), newPhi, newPsi);
								
								// reject the new conformation?
								if(Math.random() > newP / oldP) {
									cTreeLoop.changeRotationAngle(bondPhi, -deltaPhi);
									cTreeLoop.changeRotationAngle(bondPsi, -deltaPsi);
								}
							}
						}
						
						
						// have we reached our goal
						if (anglePredictor.targetRMSDistance() < TARGET_RMSD) { 
							
							if (!cTreeLoop.isClashing() && !cTreeLoop.areClashing(cTreeRemainder)) {
								double energy = energyFunction.compute();
								
								if(energy < minEnergy) {
									if(guiLoop != null)
										scene.remove(guiLoop);
									
									guiLoop = cTreeLoop.getSubchain(start-1, end+1);
									scene.add(guiLoop);		
									
									minEnergy = energy;
									System.out.println(pdb + " " +start+ "-" +end+ ": " + i + " " + minEnergy);
								}
							}
							
							break;
						}
					}
					
					i++;
				}

				previousSegment = segment;
			}
		}
		
		System.out.println("DONE!");
	}
}
