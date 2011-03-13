package test;

import tool.PDBParser;

public class ReadPDB {
	public static void main(String[] args) {

		PDBParser parser = new PDBParser("1X5R");
		
		System.out.print(parser.backbone.size());
	}
}
