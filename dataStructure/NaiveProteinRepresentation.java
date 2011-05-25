package dataStructure;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import chemestry.AminoAcid.SecondaryStructure;
import chemestry.AminoAcid.Type;

import math.Point3D;
import math.Tuple2;
import math.Vector3D;
import math.matrix.TransformationMatrix;
import tool.ChainTreeScene;
import tool.PDBParser;

/**
 * A naive representation of a protein.
 * 
 * @author Karl Koder
 */
public class NaiveProteinRepresentation extends AdjustableChainTree {
	
	private List<Point3D> backboneAtoms = new ArrayList<Point3D>();										// the bonds of the protein backbone
	
	
	/**
	 * Create a chain tree from its PDB id.
	 * 
	 * @param pdbId The PDB id to create a chain tree for.
	 */
	public NaiveProteinRepresentation(String pdbId) {
		super(pdbId);
		
		this.backboneAtoms = super.getBackboneAtomPositions();
	}
	
	/**
	 * Creates a chain tree from a list of 3D points.
	 * 
	 * @param points The points of the protein backbone atoms.
	 */
	public NaiveProteinRepresentation(List<Point3D> points) {
		super(points);
		this.backboneAtoms = points;
	}
	
	
	
	/**
	 * Returns the absolute position of the protein backbone atoms.
	 * 
	 * @return The points of the atoms.
	 */
	public List<Point3D> getBackboneAtomPositions() {
		return this.backboneAtoms;
	}
	
	
	/**
	 * Changes the rotation angle of the i-th bond by the specified angle.
	 * 
	 * @param i The index of the bond.
	 * @param angle The angle to rotate the bond by in radians.
	 */
	public void changeRotationAngle(int i, double angle) {
		// compute the change vector
		Vector3D start = new Vector3D(this.backboneAtoms.get(i));
		Vector3D end = new Vector3D(this.backboneAtoms.get(i+1));
		
		// 
		TransformationMatrix transform = new TransformationMatrix(start.vectorTo(end));
		transform.rotate(angle);
		transform.a14 = 0; transform.a24 = 0; transform.a34 = 0; // no translation
		
		Vector3D base = end;
				

		// update the rest of the chain
		i+=2;
		for(int j = this.backboneAtoms.size(); i < j; i++) {
			Vector3D atom = new Vector3D(this.backboneAtoms.get(i));

			atom = transform.transform(atom.subtract(base)).add(base);
			
			this.backboneAtoms.set(i, new Point3D(atom));
		}
	}
	
	public boolean isClashing() {
		for(Point3D point : this.backboneAtoms) {
			for(Point3D other : this.backboneAtoms) {
				if(point == other)
					continue;
				
				if(point.distance(other) < CTLeaf.atomRadius/2)
					return true;
			}
		}
		
		return false;
	}
}
