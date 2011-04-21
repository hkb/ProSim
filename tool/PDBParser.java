package tool;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.StringTokenizer;

import javax.vecmath.Point3d;

public class PDBParser {
	
	public List<Point3d> backbone = new ArrayList<Point3d>();
	public Set<Integer> alphaHelix = new HashSet<Integer>();
	public Set<Integer> betaSheet = new HashSet<Integer>();
	public Set<Integer> heteroAtoms = new HashSet<Integer>();
	
	// small state variables
	private boolean endOfBackbone;
	private int atomCount;
	
	
	/**
	 * Initialise the parser with the PDB directory.
	 * 
	 * @param directory The directory for the PDB files.
	 */
	public PDBParser(String pdbId) {
	    try {
			Scanner scanner = new Scanner(new File("pdb_files/" + pdbId + ".pdb"));
			
			while(scanner.hasNextLine()) {
				String record = scanner.nextLine();
				String type = recordType(record);
				
				if (type.equals("ENDMDL"))
					break;
				
				if (type.equals("ATOM")) {
					parseAtom(record);
				} else if (type.equals("HETATM")) {
					//parseHeteroAtom(record);
				} else if (type.equals("HELIX")) {
					parseHelix(record);
				} else if (type.equals("SHEET")) {
					parseSheet(record);
				} else if (type.equals("TER")) {
					parseTer(record);
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Parses a ATOM record
	 * 
	 * @param record The PDB record.
	 */
	private void parseAtom (String record) {
		String name = columns(record, 13, 16);
		
		// is backbone
		if(!this.endOfBackbone && (name.equals("N") || name.equals("CA") || name.equals("C"))) {
			double x = Double.parseDouble(columns(record, 31, 38));
			double y = Double.parseDouble(columns(record, 39, 46));
			double z = Double.parseDouble(columns(record, 47, 54));

			this.backbone.add(new Point3d(x, y, z));
			this.atomCount++;
		}
	}
	
	/**
	 * Parses a backbone hetero atom.
	 * 
	 * @param record The PDB record.
	 */
	private void parseHeteroAtom(String record) {
		this.parseAtom(record);
		
		// PDB files are 1-indexed but we need 0-indexing
		this.heteroAtoms.add(this.atomCount-1);
	}

	/**
	 * Parses a HELIX record
	 * 
	 * @param record The PDB record.
	 */
	private void parseHelix (String record) {
		int i = Integer.parseInt(columns(record, 22, 25)) * 3;
		int j = Integer.parseInt(columns(record, 34, 37)) * 3;
		
		while (i <= j) {
			this.alphaHelix.add(i);
			i++;
		}
	}

	/**
	 * Parses a SHEET record
	 * 
	 * @param record The PDB record.
	 */
	private void parseSheet (String record) {
		int i = Integer.parseInt(columns(record, 23, 26)) * 3;
		int j = Integer.parseInt(columns(record, 34, 37)) * 3;
		
		while (i <= j) {
			this.betaSheet.add(i);
			i++;
		}
	}
	
	/**
	 * Parses a TER record
	 * 
	 * @param record The PDB record.
	 */
	private void parseTer (String record) {
		this.endOfBackbone = true;
	}
	
	/**
	 * Returns the record type of the PDB record.
	 * 
	 * @param string
	 * @return
	 */
	private static String recordType(String string) {
		StringTokenizer tokens = new StringTokenizer(string);
		return tokens.nextToken();
	}
	
	/**
	 * Returns the specified columns of the record.
	 * The records are 1-indexed and both are included.
	 */
	private static String columns(String record, int i, int j) {
		return record.substring(i, j).trim();
	}
}
