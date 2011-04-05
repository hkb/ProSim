package experiment;

import java.util.ArrayList;
import java.util.List;

import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.ChainTree;
import energyFunction.AtomDistance;
import energyFunction.EnergyFunction;

public class PrimitiveLoopClosure {
	public static void main(String[] args) throws InterruptedException {
		/*
		 * Configuration.
		 */
		String pdbId = "1XJH";
		int start = 118;
		int end = 137;
		
		/*
		 * Setup.
		 */
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);	

		// create subtrees
		AdjustableChainTree t1 = cTree.getSubchain(0, end);
		AdjustableChainTree t2 = cTree.getSubchain(end+1, cTree.length());
		
		ChainTree[] cTrees = {t1, t2}; 
		ChainTreeScene scene = new ChainTreeScene(cTrees);
		
		// define energy function
		EnergyFunction energyFunction = new AtomDistance(t1, new AdjustableChainTree(t1));
		
		// find rotateable bonds in the segment
		List<Integer> rotateableBonds = new ArrayList<Integer>();
		
		for (int i : t1.rotatableBonds()) {
			if (start <= i && i <= end) {
				rotateableBonds.add(i);
			}
		}
		
		// unfold segment
		for (int i : rotateableBonds) {
			t1.changeRotationAngle(i, Math.random()*180);
		}
		
		scene.repaint(t1);
		
		Thread.sleep(2000);
		
		// make ready simulation
		double energy = energyFunction.compute();
		
		// simulate
		while(true) {
			int bond = rotateableBonds.get((int) (Math.random() * rotateableBonds.size()));
			double angle = (Math.random()-0.5)*15*(Math.PI/180);
			
			Thread.sleep(10);
			
			t1.changeRotationAngle(bond, angle);
			
			if (!t1.isClashing()) {
			
				double tmpEnergy = energyFunction.compute();
			
				if (tmpEnergy < energy) {
					energy = tmpEnergy;
					scene.repaint(t1);
					continue;
				}
			}
			
			t1.changeRotationAngle(bond, -angle);
		}
	}
}
