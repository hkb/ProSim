package energyFunction;

import dataStructure.CTNode;

public interface EnergyFunction {

	/**
	 * Computes the energy of the particular node in the tree.
	 * 
	 * @param node Node to compute the energy for.
	 * @return The energy of the node.
	 */
	public double compute(CTNode node);
}
