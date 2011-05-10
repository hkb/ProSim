package tool;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point2i;
import javax.vecmath.Tuple2i;

import math.Tuple2;

import dataStructure.ChainTree;

public class BackboneSegmentAnalyser {

	public static List<Tuple2<Integer, Integer>> extractIntermediateSegments(ChainTree cTree) {
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
		
		intermediateSegments.add(new Tuple2<Integer, Integer>(1, starts.get(0)));
		for (int i = 0, j = ends.size(); i < j; i++) {
			intermediateSegments.add(new Tuple2<Integer, Integer>(ends.get(i)+1, starts.get(i+1)-1));
		}
		
		return intermediateSegments;
	}
}
