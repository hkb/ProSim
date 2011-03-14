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
	
	private boolean endOfBackbone;
	
	
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
	 * @param tokens The PDB record.
	 */
	private void parseAtom (String record) {
		String name = columns(record, 13, 16);
		
		// is backbone
		if(!endOfBackbone && (name.equals("C") || name.equals("N") || name.equals("CA"))) {
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
		int i = Integer.parseInt(columns(record, 22, 25));
		int j = Integer.parseInt(columns(record, 34, 37));
		
		while (i <= j) {
			this.alphaHelix.add(i);
			i++;
		}
	}

	/**
	 * Parses a SHEET record
	 * 
	 * @param tokens The PDB record.
	 */
	private void parseSheet (String record) {
		int i = Integer.parseInt(columns(record, 23, 26));
		int j = Integer.parseInt(columns(record, 34, 37));
		
		while (i <= j) {
			this.betaSheet.add(i);
			i++;
		}
	}
	
	/**
	 * Parses a TER record
	 * 
	 * @param tokens The PDB record.
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
