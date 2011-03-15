package test;

import geom3d.Capsule3d;
import geom3d.Point3d;
import geom3d.Shape3d;

import java.awt.Color;
import java.io.IOException;
import java.util.LinkedList;

import math.matrix.TransformationMatrix;

import boundingVolume.LinesegmentSweptSphere;

import tool.BinaryTreePainter;
import tool.ChainTreeScene;
import dataStructure.AdjustableChainTree;
import dataStructure.CTNode;
import dataStructure.ChainTree;

public class AdjustableChainTreeProperties {
	public static void main(String[] args) throws InterruptedException, IOException {

		String pdbId = "1PUX";
		
		AdjustableChainTree cTree = new AdjustableChainTree(pdbId);


		new BinaryTreePainter(cTree);
		
		ChainTreeScene scene = new ChainTreeScene(cTree);
	}
}
