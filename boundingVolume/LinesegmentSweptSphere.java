package boundingVolume;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import dataStructure.CTLeaf;
import edu.geom3D.Capsule;
import edu.math.Vector;

import javax.vecmath.Point3d;

import math.matrix.TransformationMatrix;

public class LinesegmentSweptSphere implements BoundingVolume {
	
	public Capsule volume;	// the volume
	
	/**
	 * Create a new line segment swept sphere bounding volume.
	 * 
	 * @param line The line of the volume.
	 * @param radius The radius around the line.
	 */
	public LinesegmentSweptSphere(Point3d line, double radius) {
		 this(new Capsule(new Vector(0,0,0), new Vector(line.x, line.y, line.z), radius));
	}
	
	/**
	 * Creates a line segment swept sphere containing all the points.
	 * 
	 * @param points The points to be encapsulated.
	 */
	public LinesegmentSweptSphere(Collection<Point3d> points) {
		List<Vector> vectors = new ArrayList<Vector>();
		
		for (Point3d point : points) {
			vectors.add(new Vector(point.x, point.y, point.z));
		}
		
		this.volume = Capsule.createBoundingCapsule(vectors);
		this.volume.rad += CTLeaf.atomRadius/2;
	}
	
	/**
	 * Creates a new volume from an existing one.
	 */
	private LinesegmentSweptSphere(Capsule volume) {
		this.volume = volume;
	}

	@Override
	public boolean isOverlaping(BoundingVolume other) {
		return this.volume.overlaps(((LinesegmentSweptSphere) other).volume);
	}
	
	@Override
	public float volume() {
		return this.volume.volume();
	}
	
	@Override
	public BoundingVolume combine(BoundingVolume left, BoundingVolume right) {
		LinesegmentSweptSphere l = (LinesegmentSweptSphere) left;
		LinesegmentSweptSphere r = (LinesegmentSweptSphere) right;

		return new LinesegmentSweptSphere(Capsule.createBoundingCapsule_CovarianceFit(l.volume, r.volume));
	}

	@Override
	public BoundingVolume transform(TransformationMatrix transformationMatrix) {
		return new LinesegmentSweptSphere(new Capsule(transformationMatrix.transform(this.volume.p1),
													  transformationMatrix.transform(this.volume.p2),
				                                      this.volume.rad));
	}
	
	@Override
	public String toString() {
		return this.volume.toString();
	}
}
