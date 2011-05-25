package test;

import java.util.List;

import math.Vector3D;

import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;
import dataStructure.NaiveProteinRepresentation;
import energyFunction.AtomDistance;
import energyFunction.EnergyFunction;
import geom3d.Point3d;

public class VisualKarlKoderBeating {

	public static void main(String[] args) throws InterruptedException{
		/*
		 * Configuration.
		 */
		String pdbId = "2QMT"; // 2EQP, 2QMT, 1CC8, 2EEM
		double errorTolerance = 0.01;
		

		
		/*
		 * Setup.
		 */
		final class ProteinFolder implements Runnable {
			private String prefix;
			private ChainTree cTree;
			private double errorTolerance;
			private ChainTreeScene scene;
			
			public ProteinFolder(String prefix, ChainTree cTree, double errorTolerance, ChainTreeScene scene) {
				this.prefix = prefix;
				this.cTree = cTree;
				this.errorTolerance = errorTolerance;
				this.scene = scene;
			}
			
			@Override
			public void run() {
				try {
					foldProtein(this.prefix, this.cTree, this.errorTolerance, this.scene);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		ChainTree cTree = new AdjustableChainTree(pdbId);
		ChainTree naive = new NaiveProteinRepresentation(pdbId);
		
		cTree.move(new Vector3D(40.0,-30.0,-15.0));
		
		ChainTreeScene scene = new ChainTreeScene();
		scene.add(cTree);
		scene.add(naive);
		scene.scene.addText("Karl Koder", new Point3d(-5.0,35.0,0.0), 5.0);
		scene.scene.addText("ChainTree", new Point3d(35.0,35.0,0.0), 5.0);
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		
		Thread.sleep(5000);
		
		Thread thread1 = new Thread(new ProteinFolder("KarlKoder", naive, errorTolerance, scene));
		Thread thread2 = new Thread(new ProteinFolder("ChainTree", cTree, errorTolerance, scene));
		
		thread1.start();
		thread2.start();
		
		thread1.join();
		thread2.join();
		
	}
	
	private static void foldProtein(String prefix, ChainTree cTree, double errorTolerance, ChainTreeScene scene) throws InterruptedException {
		ChainTree target = new ChainTree(cTree.getBackboneAtomPositions());
		
		EnergyFunction energyFunction = new AtomDistance(cTree, target);
		List<Integer> rotateableBonds = cTree.rotatableBonds();
		
		errorTolerance++; // simpler computation
		
		cTree.unfold();
		scene.repaint(cTree);
		scene.scene.autoZoom();
		Thread.sleep(5000);
		
		// ready, set
		double energy = energyFunction.compute();
		double energyUpperBound = energy;
		int iterations = 0;
		
		// GO!
		while(true) {
			iterations++;
			
			int i = rotateableBonds.get((int) (Math.random() * rotateableBonds.size()));
			double angle = (Math.random()-0.5)*15*(Math.PI/180);
			
			cTree.changeRotationAngle(i, angle);
			
			if(cTree.isClashing()) {
				// undo move if tree is clashing
				cTree.changeRotationAngle(i, -angle);
				
			} else {
				// if not clashing then test for energy efficiency
				double tmpEnergy = energyFunction.compute();
				
				if (tmpEnergy >= energyUpperBound) {
					// if energy is higher than the upper bound then discard it
					cTree.changeRotationAngle(i, -angle);	
				
				} else {
					// try the new conformation
					energy = tmpEnergy;
					
					if (tmpEnergy * errorTolerance < energyUpperBound) {
						// permanently accept the new conformation
						energyUpperBound = tmpEnergy * errorTolerance; 
						
						System.out.println(prefix + ": " + iterations + ": " + energy);
						
						scene.repaint(cTree);
					}
				}
			}
		}
	}
}
