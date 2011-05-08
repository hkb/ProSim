package boundingVolume;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import dataStructure.CTLeaf;
import edu.geom3D.Capsule;
import edu.math.Vector;
import geom3d.Capsule3d;


import math.Point3D;
import math.matrix.TransformationMatrix;



public class LinesegmentSweptSphere implements BoundingVolume {

	public Capsule volume;	// the volume

	/**
	 * Create a new line segment swept sphere bounding volume.
	 * 
	 * @param line The line of the volume.
	 * @param radius The radius around the line.
	 */
	public LinesegmentSweptSphere(Point3D line, double radius) {
		this(new Capsule(new Vector(0,0,0), new Vector(line.x, line.y, line.z), radius));
	}

	/**
	 * Creates a line segment swept sphere containing all the points.
	 * 
	 * @param points The points to be encapsulated.
	 */
	public LinesegmentSweptSphere(Collection<Point3D> points) {
		List<Vector> vectors = new ArrayList<Vector>();

		for (Point3D point : points) {
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
		if (other instanceof LinesegmentSweptSphere) {
			return this.volume.overlaps(((LinesegmentSweptSphere) other).volume);
			
		} else if (other instanceof Empty) {
			return false;

		} else {
			throw new IllegalArgumentException("Unsupported bounding volume!");
		}
	}

	@Override
	public float volume() {
		return this.volume.volume();
	}
	
	@Override
	public BoundingVolume combine(BoundingVolume other) {
		if (other instanceof LinesegmentSweptSphere) {
			LinesegmentSweptSphere otherSphere = (LinesegmentSweptSphere) other; 

			// intermediate stuff
			Vector p1LOrig = this.volume.p1;
			Vector p2LOrig = this.volume.p2;
			
			geom3d.Point3d p1L = new geom3d.Point3d(p1LOrig.x(), p1LOrig.y(), p1LOrig.z());
			geom3d.Point3d p2L = new geom3d.Point3d(p2LOrig.x(), p2LOrig.y(), p2LOrig.z());
	
			Vector p1ROrig = otherSphere.volume.p1;
			Vector p2ROrig = otherSphere.volume.p2;
			
			geom3d.Point3d p1R = new geom3d.Point3d(p1ROrig.x(), p1ROrig.y(), p1ROrig.z());
			geom3d.Point3d p2R = new geom3d.Point3d(p2ROrig.x(), p2ROrig.y(), p2ROrig.z());
	
			Capsule3d newBound = Capsule3d.createBoundingCapsule_MaxDist(new Capsule3d(p1L, p2L, this.volume.rad), new Capsule3d(p1R, p2R, otherSphere.volume.rad));
	
			geom3d.Point3d newA = newBound.segment.getA();
			geom3d.Point3d newB = newBound.segment.getB();
			
			Vector p1New = new Vector(newA.x, newA.y, newA.z);
			Vector p2New = new Vector(newB.x, newB.y, newB.z);
			
			return new LinesegmentSweptSphere(new Capsule(p1New, p2New, newBound.rad));

		} else if (other instanceof Empty) {
			return this;

		} else {
			throw new IllegalArgumentException("Unsupported bounding volume!");
		}
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
