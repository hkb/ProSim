package math;

/**
 * A class for an infinite line.
 * 
 * @author hkb
 */
public class Line3D extends Tuple2<Vector3D, Vector3D> {

	/**
	 * Creates a new line from two points on it.
	 *  
	 * @param x The first point.Point3D
	 * @param y The second point.
	 */
	public Line3D(Point3D x, Point3D y) {
		super(new Vector3D(x), new Vector3D(y));
	}

	/**
	 * Project a point onto the line.
	 * 
	 * @return The point projected onto the line.
	 */
	public Point3D projectOnto(Point3D point) {
		Vector3D xy = this.x.vectorTo(this.y);
		double t = xy.dot(this.x.vectorTo(point.asVector())) / Math.pow(xy.length(), 2);
		
		return new Point3D(this.x.add(xy.scale(t)));
	}
}
