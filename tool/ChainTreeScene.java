package tool;

import j3dScene.J3DScene;

import java.awt.Color;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import math.Point3D;
import math.Vector3D;

import dataStructure.CTLeaf;
import dataStructure.ChainTree;
import edu.geom3D.Cylinder;
import edu.geom3D.Sphere;
import edu.math.Vector;
import geom3d.Cylinder3d;
import geom3d.Sphere3d;

/**
 * A 3D scene for painting multiple ChainTrees
 * 
 * @author hkb
 */
public class ChainTreeScene {
	
	public static int MAX_REPAINTS_PER_SECOND = -1;
	
	public J3DScene scene = J3DScene.createJ3DSceneInFrame();
	private Map<ChainTree,GUINode[]> cTrees;
	private long lastRepaintTime;
	
	/*
	 * Colours of the different part of the protein backbone.
	 */
	public Color[] colorBackbone = {new Color(150,230,255),
									new Color(150,230,255),
							     	new Color(150,200,255)};

	public Color colorAlphaHelix = new Color(255, 255, 0);
	public Color colorBetaSheet = new Color(0, 255, 0);
	public Color colorHeteroAtom = new Color(255, 255, 255, 50);
	
	
	
	/**
	 * Create a new empty scene.
	 */
	public ChainTreeScene() {
		this.scene.setBackgroundColor(Color.WHITE);
		
		this.lastRepaintTime = System.currentTimeMillis();
		
		this.cTrees = new HashMap<ChainTree,GUINode[]>();
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
	public ChainTreeScene(ChainTree[] cTrees) {
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
		this.add(cTree, 100);
	}
	
	/**
	 * Add a chain tree to the scene.
	 * 
	 * @param cTree The chain tree to be added.
	 * @param visibility The overall visibility of the chain tree.
	 */
	public void add(ChainTree cTree, int visibility) {
		
		if(!this.cTrees.containsKey(cTree)) {
			List<Point3D> points = cTree.getBackboneAtomPositions();
			
			GUINode[] guiNodes = new GUINode[points.size()];
			
			for (int i = 0, j = points.size(); i < j; i++) {
				Point3D current = points.get(i);
				Point3D next = (i+1 < j) ? points.get(i+1) : null;
				
				GUINode node = new GUINode(current, next);
				guiNodes[i] = node;
				
				Color color;
				
				// add GUI elements to scene
				if (cTree.isInAlphaHelix(i)) {
					color = this.colorAlphaHelix;
				} else if (cTree.isInBetaSheet(i)) {
					color = this.colorBetaSheet;
				} else if (cTree.isHeteroAtomBond(i)) {
					color = this.colorHeteroAtom;
				} else {
					color = this.colorBackbone[i%3];
				}
				
				// draw shapes
				//this.scene.addShape(node.sphere, modifyTransparency(new Color(color.getRed(), color.getGreen(), color.getBlue(), 100), visibility));
				this.scene.addShape(node.sphere, modifyTransparency(color, visibility));
				
				if (next != null)
					this.scene.addShape(node.cylinder, modifyTransparency(color, visibility));
			}
			
			// store the chain tree GUI info
			this.cTrees.put(cTree, guiNodes);
		}
	}
	
	/**
	 * 
	 * @param start
	 * @param end
	 * @param radius
	 * @param color
	 */
	public void addCylinder(Vector3D start, Vector3D end, float radius, Color color) {
		this.scene.addShape(new Cylinder3d(vector3DToPoint3d(start), vector3DToPoint3d(end), radius), color);
	}
	
	public void addSphere(Vector3D center, float radius, Color color) {
		this.scene.addShape(new Sphere3d(vector3DToPoint3d(center), radius), color);
	}
	
	/**
	 * Repaints the entire scene.
	 */
	public void repaint() {
		for (ChainTree cTree : this.cTrees.keySet()) {
			this.repaint(cTree);
		}
	}
	
	/**
	 * Repaints only the specified tree.
	 * scene.scene.addShape(new);
	 * @param cTree The tree to repaint.
	 */
	public void repaint(ChainTree cTree) {
		if (MAX_REPAINTS_PER_SECOND == -1 || System.currentTimeMillis() - this.lastRepaintTime > 1000 / MAX_REPAINTS_PER_SECOND) {
			
			List<Point3D> points = cTree.getBackboneAtomPositions();
			GUINode[] guiNodes = this.cTrees.get(cTree);
	
			for (int i = 0, j = points.size(); i < j; i++) {
				Point3D current = points.get(i);
				Point3D next =(i+1 < j) ? points.get(i+1) : null;
				
				guiNodes[i].update(current, next);
			}
	
			this.lastRepaintTime = System.currentTimeMillis();
			this.scene.repaint();
		}
	}
	
	/**
	 * Modifies the transparency of a color.
	 * @param color The colour to modify.
	 * @param visibility The transparency.
	 * @return New color with modified transparency.
	 */
	private static Color modifyTransparency(Color color, int visibility) {
		int alpha = (int) Math.round(color.getAlpha() / 100.0 * visibility);
		
		return new Color(color.getRed(), color.getGreen(), color.getBlue(), alpha);
	}
	
	/**
	 * @return 
	 * 
	 */
	private static geom3d.Point3d vector3DToPoint3d(Vector3D vector) {
		return new geom3d.Point3d(vector.x, vector.y, vector.z);
	}
	
	/**
	 * Private class to store information about the graphical representation 
	 * of an backbone element. 
	 */
/*
	private class GUINode {
		public Sphere sphere = new Sphere(new Vector(0,0,0), (float) CTLeaf.atomRadius/2);
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
*/
	private class GUINode {
		public Sphere3d sphere = new Sphere3d(new geom3d.Point3d(0,0,0), 0.3f);
		public Cylinder3d cylinder = new Cylinder3d(new geom3d.Point3d(0,0,0), new geom3d.Point3d(0,0,0), 0.2f);
		
		public GUINode(Point3D current, Point3D next) {
			this.update(current, next);
		}
		
		public void update(Point3D current, Point3D next) {
			this.sphere.center = new geom3d.Point3d(current.x, current.y, current.z);
			
			if (next != null) {
				this.cylinder.getSegment().setA(new geom3d.Point3d(current.x, current.y, current.z));
				this.cylinder.getSegment().setB(new geom3d.Point3d(next.x, next.y, next.z));
			}
		}
		
		private Vector pointToVector(Point3D point) {
			return new Vector(point.x, point.y, point.z);
		}
	}
}
