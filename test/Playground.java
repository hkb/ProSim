package test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import math.Point3D;
import math.Tuple2;
import dataStructure.ChainTree;
import tool.ChainTreeScene;

public class Playground {
	private static String[] sheetDominantPDBs = {"1PUX","1E2B","2YS4","1QDD","1O26","2OV0","1Z0J","2G1E","1X45","1TRE","2GCX","2COF","1BUJ",
			 "1JFM","1MDC","2JPH","1UC6","1E2B","1ZRI","1MKB","1QA7","2PND","1HX7","1YH5","1P9K",
			 "1O26","2JZR","1WEY","2YZ0","1XHJ","1QTW","1GO0","1X5M","2E5R","1WGP","2IDA","2F6X","2DIV",
			 "1WE7","1ISM","1U2F","1C01","2CZQ","1CXW"};
	
	private static String[] nonSheetDominantPDBs = {"1T0G","2J3L","1M4J","1K7C","1TI8","1UW0","1SPK","2CSK","1VZI", "2DNE","2DAO","1H6X","1AYJ",
				"1YO4","2E2F","1WGV","2H41","2A7R","1LTG","1Z5F","1TIZ","1HUX","2A7R","1FSG","1B1A","2RDQ",
				"2V1N","2Z9H","1H6Q","1IXH","2PRF","1KUF","2BYE","1Q2Y","2QBU","1SGO","1ME5","1CLV","1YS1",
				"1VYX","1H2W","8ABP","1R0R","2E9H","2FQH","1B26","1COK","1PB5","1XAX","2BH1","1XN5","1RMD",
				"2QZF","1V5R","1O8R","1RKJ"};
	
	public static void main(String[] args) {
		ChainTreeScene scene = new ChainTreeScene();

		ChainTree cTree = new ChainTree("1SIS").getSubchain(1, 4);
		
		for(int bond : cTree.rotatableBonds()) {
			//cTree.setRotationAngle(bond, (bond % 2 == 0) ? Math.PI: 0);
		}
		
		scene.add(cTree);
		
		scene.scene.centerCamera();
		scene.scene.autoZoom();
		
		List<Point3D> atoms = cTree.getBackboneAtomPositions();
		double streetched = atoms.get(0).distance(atoms.get(atoms.size()-1));
		System.out.println(streetched);
		
		
		// collet bond lengths
		List<String> pdbIds = new LinkedList<String>(Arrays.asList(sheetDominantPDBs));
		pdbIds.addAll(Arrays.asList(nonSheetDominantPDBs));
		
		List<Double> lenghts = new ArrayList<Double>();
		List<Double> angles = new ArrayList<Double>();
		
		for(String pdb : pdbIds) {
			try {
				cTree = new ChainTree(pdb);
			
				for(Tuple2<Integer,Integer> segment : cTree.getIntermediateSegments()) {
					if (segment.y - segment.x == 3) {
						atoms = cTree.getBackboneAtomPositions(segment.x, segment.y);
						lenghts.add(atoms.get(0).distance(atoms.get(atoms.size()-1)));
						angles.addAll(cTree.getDihedralAngles(segment.x, segment.y));
					}
				}
			} catch (Exception e) {
				System.err.println(e);
			}
		}
		
		// find stretched bonds
		Collections.sort(lenghts);
		Collections.reverse(lenghts);
		System.out.println(lenghts);
		
		double value = lenghts.get(0);
		int i = 0;
		while(value > 9) {
			i++;
			value = lenghts.get(i);
		}
		
		System.out.println(i / ((double) lenghts.size()));
	}
}
