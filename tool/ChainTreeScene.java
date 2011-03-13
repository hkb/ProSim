package tool;

import java.awt.Color;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import dataStructure.ChainTree;
import edu.geom3D.Cylinder;
import edu.geom3D.Sphere;
import edu.j3dScene.J3DScene;
import edu.math.Vector;

/**
 * A 3D scene for painting multiple ChainTrees
 * 
 * @author hkb
 */
public class ChainTreeScene {
	
	private J3DScene scene = J3DScene.createJ3DSceneInFrame();
	private Map<ChainTree,GUINode[]> cTrees = new HashMap<ChainTree,GUINode[]>();
	
	/*
	 * Colours of the different part of the protein backbone.
	 */
	public Color colorAtom = new Color(0, 0, 255, 100);
	public Color colorBond = Color.CYAN;
	
	public Color colorPeptideBond = new Color(169, 169, 169, 210);
	public Color colorAlphaHelix = new Color(255, 255, 0, 255);
	public Color colourBetaSheet = new Color(0, 255, 0, 255);
	
	
	
	/**
	 * Create a new empty scene.
	 */
	public ChainTreeScene() {
		this.scene.setBackgroundColor(Color.DARK_GRAY);
	}
	
	/**
	 * Create a new scene with one tree.
	 * 
	 * @param cTree The tree to be in the scene.
	 */
	public ChainTreeScene(ChainTree cTree) {
		this();
		
		this.add(cTree);
		this.scene.centerCamera();
		this.scene.autoZoom();
	}
	
	/**
	 * Create a new scene with multiple trees.
	 * 
	 * @param cTrees The trees to be on the scene.
	 */
	public ChainTreeScene(Collection<ChainTree> cTrees) {
		this();
		
		for(ChainTree cTree : cTrees) {
			this.add(cTree);
		}
		this.scene.centerCamera();
		this.scene.autoZoom();
	}
	
	
	
	/**
	 * Add a chain tree to the scene.
	 * 
	 * @param cTree The chain tree to be added.
	 */
	public void add(ChainTree cTree) {
		if(!this.cTrees.containsKey(cTree)) {
			List<Point3d> points = cTree.getBackbonePoints();
			
			GUINode[] guiNodes = new GUINode[points.size()];
			
			for (int i = 0, j = points.size()-1; i < j; i++) {
				Point3d current = points.get(i);
				Point3d next = points.get(i+1);
				
				GUINode node = new GUINode(current, next);
				guiNodes[i] = node;
				
				// add GUI elements to scene
				Color bondColor;
				if (cTree.isInAlphaHelix(i)) {
					bondColor = this.colorAlphaHelix;
				} else if (cTree.isInBetaSheet(i)) {
					bondColor = this.colourBetaSheet;
				} else if (cTree.isPeptide(i)) {
					bondColor = this.colorPeptideBond;
				} else {
					bondColor = this.colorBond;
				}
				
				this.scene.addShape(node.sphere, this.colorAtom);
				this.scene.addShape(node.cylinder, bondColor);
			}
			
			// store the chain tree GUI info
			this.cTrees.put(cTree, guiNodes);
			
			this.repaint();
		}
	}
	
	/**
	 * Repaints the entire scene.
	 */
	public void repaint() {
		for (ChainTree cTree : this.cTrees.keySet()) {
			this.repaint(cTree);
		}
		
		this.scene.repaint();
	}
	
	/**
	 * Repaints only the specified tree.
	 * 
	 * @param cTree The tree to repaint.
	 */
	private void repaint(ChainTree cTree) {
		List<Point3d> points = cTree.getBackbonePoints();
		GUINode[] guiNodes = this.cTrees.get(cTree);
		
		for (int i = 0, j = points.size()-1; i < j; i++) {
			Point3d current = points.get(i);
			Point3d next = points.get(i+1);
			
			guiNodes[i].update(current, next);
		}
	}
	
	/**
	 * Private class to store information about the graphical representation 
	 * of an backbone element. 
	 */
	private class GUINode {
		public Sphere sphere = new Sphere(new Vector(0,0,0), 0.85f);
		public Cylinder cylinder = new Cylinder(new Vector(0,0,0), new Vector(0,0,0), 0.4f);
		
		public GUINode(Point3d current, Point3d next) {
			this.update(current, next);
		}
		
		public void update(Point3d current, Point3d next) {
			this.sphere.center = this.pointToVector(next);
			this.cylinder.p1 = this.pointToVector(current);
			this.cylinder.p2 = this.pointToVector(next);		
		}
		
		private Vector pointToVector(Point3d point) {
			return new Vector(point.x, point.y, point.z);
		}
	}
}
