package tool;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point2i;
import javax.vecmath.Tuple2i;

import dataStructure.ChainTree;

public class BackboneSegmentAnalyser {

	public static List<Tuple2<Integer, Integer>> extractIntermediateSegments(ChainTree cTree) {
		List<Tuple2<Integer, Integer>> intermediateSegments = new ArrayList<Tuple2<Integer, Integer>>();
		
		int start = -1;		
		for (int i = 0, j = cTree.length(); i < j; i++) {
			if (cTree.isInAlphaHelix(i) || cTree.isInBetaSheet(i)) {
				if (start == -1)
					continue;
				
				intermediateSegments.add(new Tuple2<Integer, Integer>(start, i-1));
				start = -1;
			} else if (start == -1) {
				start = i;
			}
		}
		
		return intermediateSegments;
	}
}
