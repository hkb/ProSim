package test;

import java.awt.Color;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3d;

import tool.BinaryTreePainter;
import tool.ChainTreeScene;

import dataStructure.ChainTree;
import edu.geom3D.Sphere;
import edu.math.Vector;

public class BuildChainTree {

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws InterruptedException, IOException {

		ChainTree cTree = new ChainTree("1PUX");

		ChainTreeScene scene = new ChainTreeScene(cTree);
				
		while(true) {
			InputStreamReader reader = new InputStreamReader(System.in);

			int i = (int) (Math.random() * cTree.length());
			double angle = (Math.random()-0.5)*15*(Math.PI/180);
			
			System.out.println(i + " - " + angle);
			
			cTree.changeRotationAngle(i, angle);
			scene.repaint();
			
			if(cTree.isClashing()) {
				System.err.println("Self clash!");
				//cTree.changeRotationAngle(i, -angle);
				
				// paint error
				List<Point3d> points = cTree.getBackbonePoints();
				Point3d point = points.get(cTree.cl1);
				scene.scene.addShape(new Sphere(new Vector(point.x, point.y, point.z), (float)1.7/2), new Color(255,0,0,250));
				point = points.get(cTree.cl2);
				scene.scene.addShape(new Sphere(new Vector(point.x, point.y, point.z), (float)1.7/2), new Color(255,0,0,250));
				
				reader.read();
				
			} else {
				//scene.repaint();
			}
		}
	}
}
