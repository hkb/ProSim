package experiment;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collection;
import java.util.Date;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
	private static enum UnfoldingRestriction {NONE, NEIGHBOUR_INDEPENDENT, NEIGHBOUR_DEPENDENT, NEIGHBOUR_INDEPENDENT_NO_CLASH, NEIGHBOUR_DEPENDENT_NO_CLASH, FROM_PROTEIN};
	private static RamachandranDistribution ramachandranDistribution = new RamachandranDistribution();
	private static ChainTreeScene scene = null;//new ChainTreeScene();
	private static Profiler profiler = new Profiler();
	private static Writer output;

	/*
	 * Configuration.
	 */
	private static String[] sheetDominantPDBs = {"1PUX","1E2B","2YS4","1QDD","1O26","2OV0","1Z0J","2G1E","1X45","1TRE","2GCX","2COF","1BUJ",
												 "1JFM","1MDC","2JPH","1UC6","1E2B","1ZRI","1MKB","1QA7","2PND","1HX7","1YH5","1P9K",
												 "1O26","2JZR","1WEY","2YZ0","1XHJ","1QTW","1GO0"};
	
	private static String[] nonSheetDominantPDBs = {"1T0G","2J3L","1M4J","1K7C","1TI8","1UW0","1SPK","2CSK","1VZI", "2DNE","2DAO","1H6X","1AYJ",
													"1YO4","2E2F","1WGV","2H41","2A7R","1LTG","1Z5F","1TIZ","1HUX","2A7R","1FSG","1B1A","2RDQ",
													"2V1N","2Z9H","1H6Q","1IXH","2PRF","1KUF","2BYE"};
	
	private static double TARGET_RMSD = 0.08;
	private static int MAX_ITERATIONS_PER_CLOSE = 5000;
	private static int NUMBER_OF_TRIAL_LOOPS = 100;
	
	
		
	public static void main(String[] args) throws Exception {
		
		
		/**
		 * Thread for folding all loops in a single protein.
		 */
		final class LoopCloser implements Runnable {
			private String pdbId;
			
			public LoopCloser(String pdbId) {
				this.pdbId = pdbId;
			}
			
			@Override
			public void run() {
				closeAllLoops(this.pdbId);
			}
		}
		
		/**
		 * Thread for folding a single protein.
		 */
		final class SheetFolder implements Runnable {
			private String pdbId;
			
			public SheetFolder(String pdbId) {
				this.pdbId = pdbId;
			}
			
			@Override
			public void run() {
				closeAllSheetLoops(this.pdbId);
			}
		}
		
		List<Thread> threads = new LinkedList<Thread>();
		
		/*
		 * Close all loops.
		 */
		
		output = new FileWriter("/home/hkb/Documents/Datalogi/Bachelor/Bachelorprojekt/data/latest-all-loop-experiment.tsv");
		
		List<String> pdbIds = new LinkedList<String>(Arrays.asList(sheetDominantPDBs));
		pdbIds.addAll(Arrays.asList(nonSheetDominantPDBs));
		
		for(String pdbId : pdbIds) {
			Thread thread = new Thread(new LoopCloser(pdbId));
			thread.start();
			threads.add(thread);
		}
		
		for(Thread thread : threads) {
			thread.join();
		}
		
		output.close();
		
		System.out.println("Done with all loops.");
		
		
		
		/*
		 * Close all central beta sheets.
		 */
		/*
		output = new FileWriter("/home/hkb/Documents/Datalogi/Bachelor/Bachelorprojekt/data/latest-central-betasheet-experiment.tsv");
		
		
		for(String pdbId : sheetDominantPDBs) {
			Thread thread = new Thread(new SheetFolder(pdbId));
			thread.start();
			threads.add(thread);
		}
		
		for(Thread thread : threads) {
			thread.join();
		}
		
		output.close();
		
		System.out.println("Done with sheet dominant proteins.");
		*/
		
		
		/*
		 * Close all non-central beta sheets.
		 */
		/*
		output = new FileWriter("/home/hkb/Documents/Datalogi/Bachelor/Bachelorprojekt/data/latest-non-central-betasheet-experiment.tsv");
		
		for(String pdbId : nonSheetDominantPDBs) {
			Thread thread = new Thread(new SheetFolder(pdbId));
			thread.start();
			threads.add(thread);
		}
		
		for(Thread thread : threads) {
			thread.join();
		}
		
		output.close();
		
		System.out.println("Done with non sheet dominant proteins.");
		*/
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
						
			if(segmentLength(sheet1.y+1, sheet2.x-1) >= 4) { // min 4 residues long
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
		
		for(int i = 1; i < segments.size()-1; i++) { // segment must be minimum 4 residues long
			Tuple2<Integer, Integer> segment = segments.get(i);
			
			if(segmentLength(segment.x, segment.y) >= 4) {	
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
		
		// phi, psi angles but not from loop and not in secondary structure
		List<Tuple2<Double,Double>> phiPsiPairs = new ArrayList<Tuple2<Double,Double>>();
		for (int aminoAcid = 2; aminoAcid < cTree.length()-1; aminoAcid++) {
			List<Double> angles = cTree.getDihedralAngles(aminoAcid, aminoAcid);
			
			if ((aminoAcid < start || end < aminoAcid) && !cTree.isInHelix(aminoAcid) && !cTree.isInSheet(aminoAcid)) {
				phiPsiPairs.add(new Tuple2<Double,Double>(angles.get(0), angles.get(1)));
			}
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

		for(UnfoldingRestriction unfolding : new UnfoldingRestriction[] {UnfoldingRestriction.NONE, UnfoldingRestriction.NEIGHBOUR_INDEPENDENT, UnfoldingRestriction.NEIGHBOUR_INDEPENDENT_NO_CLASH, UnfoldingRestriction.NEIGHBOUR_DEPENDENT, UnfoldingRestriction.NEIGHBOUR_DEPENDENT_NO_CLASH, UnfoldingRestriction.FROM_PROTEIN}) {
			Collection<List<Tuple2<Integer,Tuple2<Double,Double>>>> conformations = generateConformations(cTreeLoopCopy, cTreeRemainder, phiPsiPairs, start, end, unfolding);
			
			for(FoldingRestriction restriction : new FoldingRestriction[]{FoldingRestriction.NONE, FoldingRestriction.NEIGHBOUR_INDEPENDENT, FoldingRestriction.NEIGHBOUR_DEPENDENT}) {
				int itterations = 0;
				int unclosed = 0;
				int clashes = 0;
				Collection<Double> energies = new ArrayList<Double>();
				double minEnergy = Double.MAX_VALUE;
				
				for(Collection<Tuple2<Integer,Tuple2<Double,Double>>> conformation : conformations) {
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
										oldP = ramachandranDistribution.probability(cTreeLoop.getAminoAcidType(aminoAcid), oldPhi, oldPsi);
										newP = ramachandranDistribution.probability(cTreeLoop.getAminoAcidType(aminoAcid), newPhi, newPsi);
										break;
									case NEIGHBOUR_DEPENDENT: 
										oldP = ramachandranDistribution.probability(cTreeLoop.getAminoAcidType(aminoAcid), cTreeLoop.getAminoAcidType(aminoAcid-1), cTreeLoop.getAminoAcidType(aminoAcid+1), oldPhi, oldPsi);
										newP = ramachandranDistribution.probability(cTreeLoop.getAminoAcidType(aminoAcid), cTreeLoop.getAminoAcidType(aminoAcid-1), cTreeLoop.getAminoAcidType(aminoAcid+1), newPhi, newPsi);
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
	
	
	
	
	
	private static List<Tuple2<Integer,Tuple2<Double,Double>>> unfold(AdjustableChainTree cTree, AdjustableChainTree other, List<Tuple2<Double,Double>> phiPsiPairs, int start, int end, UnfoldingRestriction restriction) {
		Set rotateableBonds = new HashSet(cTree.rotatableBonds());
		
		List<Tuple2<Integer,Tuple2<Double,Double>>> conformation;
		
		do {
			conformation = new LinkedList<Tuple2<Integer,Tuple2<Double,Double>>>();
			
			for (int i = start; i <= end; i++) {
				if(rotateableBonds.contains(cTree.getPhi(i))) {
					if(!rotateableBonds.contains(cTree.getPsi(i))) {
						throw new IllegalArgumentException("Both phi and psi angles must be rotateable!");
					}
				
					Tuple2<Double,Double> angles = null;
					
					if(restriction == UnfoldingRestriction.NONE) {
						angles = new Tuple2<Double,Double>(Math.random()*Math.PI*2, Math.random()*Math.PI*2);
					} else if (restriction == UnfoldingRestriction.FROM_PROTEIN) {
						angles = phiPsiPairs.get((int) (Math.random() * phiPsiPairs.size()));
					} else if (restriction == UnfoldingRestriction.NEIGHBOUR_INDEPENDENT || restriction == UnfoldingRestriction.NEIGHBOUR_INDEPENDENT_NO_CLASH) {
						angles = ramachandranDistribution.purposeAngle(cTree.getAminoAcidType(cTree.getAminoAcid(i)));
					} else {
						angles = ramachandranDistribution.purposeAngle(cTree.getAminoAcidType(cTree.getAminoAcid(i)), cTree.getAminoAcidType(cTree.getAminoAcid(i-1)), cTree.getAminoAcidType(cTree.getAminoAcid(i+1)));
					}
				
					cTree.setRotationAngle(cTree.getPhi(i), angles.x);
					cTree.setRotationAngle(cTree.getPsi(i), angles.y);
						
					conformation.add(new Tuple2<Integer,Tuple2<Double,Double>>(i, angles));
				}
			}
		} while ((restriction == UnfoldingRestriction.NEIGHBOUR_INDEPENDENT_NO_CLASH || restriction == UnfoldingRestriction.NEIGHBOUR_DEPENDENT_NO_CLASH) && (cTree.isClashing() || cTree.areClashing(other)));
		
		return conformation;
	}

	
	private static Collection<List<Tuple2<Integer,Tuple2<Double,Double>>>> generateConformations(AdjustableChainTree cTree, AdjustableChainTree other, List<Tuple2<Double,Double>> phiPsiPairs, int start, int end, UnfoldingRestriction restriction) {
		List<List<Tuple2<Integer,Tuple2<Double,Double>>>> conformations = new ArrayList<List<Tuple2<Integer,Tuple2<Double,Double>>>>();
		
		for(int i = 0; i < NUMBER_OF_TRIAL_LOOPS; i++) {
			conformations.add(unfold(cTree, other, phiPsiPairs, start, end, restriction));
		}
		
		return conformations;
	}
	
	private static void unforldIntoConformation(ChainTree cTree, Collection<Tuple2<Integer,Tuple2<Double,Double>>> conformation) {
		for(Tuple2<Integer,Tuple2<Double,Double>> bondInfo : conformation) {
			cTree.setRotationAngle(cTree.getPhi(bondInfo.x), bondInfo.y.x);
			cTree.setRotationAngle(cTree.getPsi(bondInfo.x), bondInfo.y.y);
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
	
	private static int segmentLength(int x, int y) {
		return y-x+1;
	}
}
