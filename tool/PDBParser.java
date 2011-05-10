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

import chemestry.AminoAcid;

import math.Point3D;
import math.Tuple2;

public class PDBParser {
	
	public List<AminoAcid.Type> primaryStructure = new ArrayList<AminoAcid.Type>();
	public List<Point3D> backboneAtomPositions = new ArrayList<Point3D>();
	public List<Tuple2<Integer, Integer>> helixes = new ArrayList<Tuple2<Integer, Integer>>();
	public List<Tuple2<Integer, Integer>> sheets = new ArrayList<Tuple2<Integer, Integer>>();
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
					parsePrimaryStructure(record);
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

			this.backboneAtomPositions.add(new Point3D(x, y, z));
			this.atomCount++;
		}
	}
	
	private void parsePrimaryStructure(String record) {
		String name = columns(record, 13, 16);
		
		if(!this.endOfBackbone && name.equals("N")) {
			this.primaryStructure.add(AminoAcid.Type.valueOf(columns(record, 17, 20)));
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
		int i = Integer.parseInt(columns(record, 22, 25));
		int j = Integer.parseInt(columns(record, 34, 37));

		this.helixes.add(new Tuple2<Integer, Integer>(i, j));
	}

	/**
	 * Parses a SHEET record
	 * 
	 * @param record The PDB record.
	 */
	private void parseSheet (String record) {
		int i = Integer.parseInt(columns(record, 23, 26));
		int j = Integer.parseInt(columns(record, 34, 37));

		this.sheets.add(new Tuple2<Integer, Integer>(i, j));
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
