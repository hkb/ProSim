package tool;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point2i;
import javax.vecmath.Tuple2i;

import math.Tuple2;

import dataStructure.ChainTree;

public class BackboneSegmentAnalyser {
	
	public static List<Tuple2<Integer, Integer>> getHelixSegments(ChainTree cTree) {
		return cTree.helixes;
	}
	
	public static List<Tuple2<Integer, Integer>> getSheetSegments(ChainTree cTree) {
		List<Tuple2<Integer, Integer>> sortedSheets = new ArrayList<Tuple2<Integer, Integer>>();
		
		for (Tuple2<Integer, Integer> sheet : cTree.sheets) {
			int i = 0;
			while (i < sortedSheets.size() && sortedSheets.get(i).x < sheet.x) {
				i++;
			}
			
			sortedSheets.add(i, sheet);
		}
		return sortedSheets;
	}	

	public static List<Tuple2<Integer, Integer>> getIntermediateSegments(ChainTree cTree) {
		List<Tuple2<Integer, Integer>> secondaryStructures = new ArrayList<Tuple2<Integer, Integer>>();
		
		secondaryStructures.addAll(cTree.helixes);
		secondaryStructures.addAll(cTree.sheets);

		List<Integer> starts = new ArrayList<Integer>();
		List<Integer> ends = new ArrayList<Integer>();
		
		for (Tuple2<Integer, Integer> structure : secondaryStructures) {
			starts.add(structure.x);
			ends.add(structure.y);
		}
		
		Collections.sort(starts);
		Collections.sort(ends);
		
		starts.add(cTree.length());
		List<Tuple2<Integer, Integer>> intermediateSegments = new ArrayList<Tuple2<Integer, Integer>>();
		
		for (int i = 0, j = ends.size()-1; i < j; i++) {
			intermediateSegments.add(new Tuple2<Integer, Integer>(ends.get(i)+1, starts.get(i+1)-1));
		}
		
		return intermediateSegments;
	}
}
