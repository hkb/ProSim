package experiment;

import java.util.List;

import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import energyFunction.AtomDistance;	
import energyFunction.DihedralAngles;
import energyFunction.EnergyFunction;

public class MonteCarloSimulation {
	public static void main(String[] args){
		/*
		 * Setup
		 */
		String pdbId = "1PUX"; // 1PUX <3
		double errorTolerance = 0.01;
		double targetEnergy = 0.01;
		
		
		
		/*
		 * Experiment.
		 */
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);
		AdjustableChainTree target = new AdjustableChainTree(pdbId);
		
		EnergyFunction energyFunction = new AtomDistance(cTree, target);
		List<Integer> rotateableBonds = cTree.rotatableBonds();
		
		errorTolerance++; // simpler computation
		
		cTree.unfold();

		ChainTreeScene scene = new ChainTreeScene(cTree);
		ChainTreeScene targetScene = new ChainTreeScene(target);
		
		scene.scene.setWindowTitle(pdbId + " - Folding");
		targetScene.scene.setWindowTitle(pdbId + " - Target");
		
		double energy = energyFunction.compute();
		double energyUpperBound = energy;
		int iterations = 0;
		
		while(energy > targetEnergy) {
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
						
						System.out.println(iterations + ": " + energy);
						scene.repaint();
					}
				}
			}
		}
		
		System.out.println("Simulation done.");
	}
}
