package test;

import java.awt.Color;

import geom3d.Cylinder3d;
import geom3d.Point3d;
import math.Vector3D;
import j3dScene.J3DScene;

public class VisualVectorDebugger {
	public static void main(String[] args) {
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		
		Vector3D a = new Vector3D(2.0,0.0,0.0);
		Vector3D b = new Vector3D(0.0,1.0,0.0);
		Vector3D c = new Vector3D(1.0,1.0,1.0);
		Vector3D d = a.cross(b);
		Vector3D e = a.projectOnto(c);
		
		scene.addShape(new Cylinder3d(new Point3d(), vectorToPoint(a), 0.1f), Color.RED);
		scene.addShape(new Cylinder3d(new Point3d(), vectorToPoint(b), 0.1f), Color.GREEN);
		scene.addShape(new Cylinder3d(new Point3d(), vectorToPoint(c), 0.1f), Color.YELLOW);
		scene.addShape(new Cylinder3d(new Point3d(), vectorToPoint(d), 0.1f), Color.BLUE);
		scene.addShape(new Cylinder3d(new Point3d(), vectorToPoint(e), 0.2f), new Color(255, 255, 0, 100));
		
		scene.addShape(new Cylinder3d(vectorToPoint(c), vectorToPoint(e), 0.01f));
	}
	
	private static Point3d vectorToPoint(Vector3D vector) {
		return new Point3d(vector.x, vector.y, vector.z);
	}
}
