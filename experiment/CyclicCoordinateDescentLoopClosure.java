package experiment;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Date;
import java.util.EnumSet;
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
import tool.ChainTreeScene;
import tool.Profiler;
import tool.RamachandranDistribution;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;
import energyFunction.EnergyFunction;
import energyFunction.LoopAtomDistance;

public class CyclicCoordinateDescentLoopClosure {
	
	private static enum FoldingRestriction {NONE, NEIGHBOUR_INDEPENDENT, NEIGHBOUR_DEPENDENT};
	private static enum UnfoldingRestriction {NONE, NATURAL};
	private static RamachandranDistribution ramachandranDistribution = new RamachandranDistribution();
	private static ChainTreeScene scene = null;//new ChainTreeScene();
	private static Profiler profiler = new Profiler();
	private static Writer output;

	/*
	 * Configuration.
	 */
	private static String[] pdbIds = {"1PUX","1T0G","1E2B","2J3L","1M4J","1K7C","2YS4","1QDD","1O26","2OV0","1TI8","1UW0","1Z0J","1SPK","2CSK","2G1E","1VZI",
									  "2DNE","1X45","1TRE","2DAO","1H6X","1AYJ","1YO4","2OV0","2E2F","1WGV","2H41","2A7R","2BYE"};
	private static String[] sheetDominantPDBs = {"1PUX","1TRE"};
	private static double TARGET_RMSD = 0.08;
	private static int MAX_ITERATIONS_PER_CLOSE = 5000;
	private static int NUMBER_OF_TRIAL_LOOPS = 100;
	
	private static int NUMBER_OF_CONCURRENT_THREADS = 1;
	
	public static void main(String[] args) throws Exception {
		
		/**
		 * Thread for folding a single protein.
		 */
		final class ProteinFolder implements Runnable {
			private String pdbId;
			
			public ProteinFolder(String pdbId) {
				this.pdbId = pdbId;
			}
			
			@Override
			public void run() {
				closeAllLoops(this.pdbId);
			}
		}
		
		
		/*
		 * Run experiments.
		 */
		Date today = Calendar.getInstance().getTime();
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd_hh.mm.ss");
		output = new FileWriter("/home/hkb/Documents/Datalogi/Bachelor/Bachelorprojekt/data/experiments-"+formatter.format(today)+".tsv");
		
		
		for(int i = 0; i < pdbIds.length; i += NUMBER_OF_CONCURRENT_THREADS) {
			Thread[] threads = new Thread[NUMBER_OF_CONCURRENT_THREADS];
			
			for(int k = 0; k < NUMBER_OF_CONCURRENT_THREADS && i+k < pdbIds.length; k++) {
				System.out.println("starting thread for: " + pdbIds[i+k]);
				threads[k] = new Thread(new ProteinFolder(pdbIds[i+k]));
				threads[k].start();
			}
			
			for(int k = 0; k < NUMBER_OF_CONCURRENT_THREADS; k++) {
				threads[k].join();
			}
		}
		
		output.close();
		
		System.out.println("TOTALLY DONE!");
	}
	
	private static void closeAllSheetLoops (String pdbId) {
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		List<Tuple2<Integer, Integer>> segments = cTree.getSheetSegments();
		/*		
		scene.add(cTree.getSubchain(1, segments.get(0).x));
		for(Tuple2<Integer, Integer> segment : segments) {
			scene.add(cTree.getSubchain(segment.x, segment.y));
		}
		scene.add(cTree.getSubchain(segments.get(segments.size()-1).y, cTree.length()));
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		*/
		
		for(int i = 0; i < segments.size()-1; i++) {
			Tuple2<Integer, Integer> sheet1 = segments.get(i);
			Tuple2<Integer, Integer> sheet2 = segments.get(i+1);
						
			if(sheet2.x - sheet1.y >= 3) { // min 4 residues long
				closeLoop(pdbId, cTree, sheet1.y+1, sheet2.x-1);
			}
		}
		
		//scene.remove(cTree);
	}
	
	private static void closeAllLoops(String pdbId) {
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		List<Tuple2<Integer, Integer>> segments = cTree.getIntermediateSegments();
		/*
		scene.add(cTree.getSubchain(segments.get(0).x, segments.get(0).y+1));
		
		for(int i = 0; i < segments.size()-1; i++) {
			scene.add(cTree.getSubchain(segments.get(i).y+1, segments.get(i+1).x-1));
			
			if(i > 0 && segments.get(i).y - segments.get(i).x < 4) {
				scene.add(cTree.getSubchain(segments.get(i).x-1, segments.get(i).y+1));
			}
		}
		scene.add(cTree.getSubchain(segments.get(segments.size()-1).x-1, segments.get(segments.size()-1).y));
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		*/
		
		for(int i = 1; i < segments.size()-1; i++) { // min 4 residues long
			Tuple2<Integer, Integer> segment = segments.get(i);
			
			if(segment.y - segment.x >= 3) {	
				closeLoop(pdbId, cTree, segment.x, segment.y);
			}
		}
		
		//scene.remove(cTree);

	}


	private static void closeLoop(String pdbId, AdjustableChainTree cTree, int start, int end) {
		// setup chain trees
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

		// init CCD
		CyclicCoordinateDescent anglePredictor = new CyclicCoordinateDescent(cTreeLoop, cTreeLoop.getSubchain(cTreeLoop.length(), cTreeLoop.length()));
		
		// GUI
		ChainTree guiLoop = null;

		/*
		 * Close loops.
		 */
		
		// generate unfold conformations
		AdjustableChainTree cTreeLoopCopy = cTreeLoop.getSubchain(1, cTreeLoop.length());
		
		for(UnfoldingRestriction unfolding : new UnfoldingRestriction[] {UnfoldingRestriction.NONE, UnfoldingRestriction.NATURAL}) {
			Collection<Collection<Tuple2<Integer,Double>>> conformations = generateConformations(cTreeLoopCopy, cTreeRemainder, dihedralAngles, start, end, unfolding);
			
			for(FoldingRestriction restriction : new FoldingRestriction[]{FoldingRestriction.NONE, FoldingRestriction.NEIGHBOUR_INDEPENDENT, FoldingRestriction.NEIGHBOUR_DEPENDENT}) {
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
							
							if(restriction == FoldingRestriction.NONE) {
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
						if (anglePredictor.targetRMSDistance() < TARGET_RMSD) { 
							
							if (cTreeLoop.isClashing() || cTreeLoop.areClashing(cTreeRemainder)) {
								clashes++;
							} else {
								double energy = energyFunction.compute();
								energies.add(energy);
								
								if(energy < minEnergy) {
									/*
									if(guiLoop != null)
										scene.remove(guiLoop);
										guiLoop = cTreeLoop.getSubchain(start-1, end+1);
										scene.add(guiLoop);						
									*/
									minEnergy = energy;
								}
							}
							
							break;
						}
						
						// to many iterations?
						if(itterations >= MAX_ITERATIONS_PER_CLOSE) {
							unclosed++;
							break;
						}
					}
				}
				
				log(pdbId + "\t" + start + "\t" + end + "\t" + unfolding + "\t" + restriction + "\t" + clashes + "\t" +unclosed+"\t"+energies);				
				System.out.println(pdbId + " " + start + "-" + end + " " + unfolding+ " " + restriction + ": " + clashes + " " +unclosed+" "+minEnergy);
			}
		}
	}
	
	private static Collection<Tuple2<Integer,Double>> unfoldNatural(AdjustableChainTree cTree, AdjustableChainTree other, List<Double> angles, int start, int end) {
		int startBond = cTree.getPhi(start);
		int endBond = cTree.getPsi(end);
		
		Collection<Tuple2<Integer,Double>> conformation;
		
		do {
			conformation = new LinkedList<Tuple2<Integer,Double>>();
			
			for (int bond : cTree.rotatableBonds()) {
				if (startBond <= bond && bond <= endBond) {
					double angle = angles.get((int) (Math.random() * angles.size()));
					
					List<Double> currentAngles = cTree.getDihedralAngles(cTree.getAminoAcid(bond), cTree.getAminoAcid(bond));
					double currentAngle = currentAngles.get((bond % 3 == 0) ? 0 : 1);
					cTree.changeRotationAngle(bond, currentAngle-angle);
					
					conformation.add(new Tuple2<Integer,Double>(bond, currentAngle-angle));
				}
			}
		} while (cTree.isClashing() || cTree.areClashing(other));

		return conformation;
	}
	
	private static Collection<Tuple2<Integer,Double>> unfoldRandom(AdjustableChainTree cTree, int start, int end) {
		int startBond = cTree.getPhi(start);
		int endBond = cTree.getPsi(end);
		
		Collection<Tuple2<Integer,Double>> conformation = new LinkedList<Tuple2<Integer,Double>>();
			
		for (int bond : cTree.rotatableBonds()) {
			if (startBond <= bond && bond <= endBond) {
				conformation.add(new Tuple2<Integer,Double>(bond, Math.random()*Math.PI*2));
			}
		}

		return conformation;
	}
	
	private static Collection<Collection<Tuple2<Integer,Double>>> generateConformations(AdjustableChainTree cTree, AdjustableChainTree other, List<Double> angles, int start, int end, UnfoldingRestriction restriction) {
		List<Collection<Tuple2<Integer,Double>>> conformations = new ArrayList<Collection<Tuple2<Integer,Double>>>();
		
		for(int i = 0; i < NUMBER_OF_TRIAL_LOOPS; i++) {
			if(restriction == UnfoldingRestriction.NATURAL) {
				conformations.add(unfoldNatural(cTree, other, angles, start, end));
			} else {
				conformations.add(unfoldRandom(cTree, start, end));
			}
		}
		
		return conformations;
	}
	
	private static void unforldIntoConformation(ChainTree cTree, Collection<Tuple2<Integer,Double>> conformation) {
		for(Tuple2<Integer,Double> bondInfo : conformation) {
			cTree.changeRotationAngle(bondInfo.x, bondInfo.y);
		}
	}
	
	private static synchronized void log(String str) {
		try {
			output.write(str+"\n");
			output.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
