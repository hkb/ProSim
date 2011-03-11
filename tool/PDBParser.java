package tool;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
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
				
				if (type.equals("ATOM")) {
					parseAtom(record);
				} else if (type.equals("HELIX")) {
					parseHelix(record);
				} else if (type.equals("SHEET")) {
					parseSheet(record);
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Parses a ATOM record
	 * 
	 * @param tokens The PDB record.
	 */
	private void parseAtom (String record) {
		String name = columns(record, 13, 16);
		
		// is backbone
		if(name.charAt(0) == 'N' || name.charAt(0) == 'C' || name.substring(0, 2).equals("CA")) {
			double x = Double.parseDouble(columns(record, 31, 38));
			double y = Double.parseDouble(columns(record, 39, 46));
			double z = Double.parseDouble(columns(record, 47, 54));

			this.backbone.add(new Point3d(x, y, z));
		}
	}

	/**
	 * Parses a HELIX record
	 * 
	 * @param tokens The PDB record.
	 */
	private void parseHelix (String record) {
		
	}

	/**
	 * Parses a SHEET record
	 * 
	 * @param tokens The PDB record.
	 */
	private void parseSheet (String record) {
		
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
		return record.substring(i, j);
	}
}
