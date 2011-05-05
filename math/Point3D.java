package math;

/**
 * A simple point in 3-dimensional space.
 * 
 * @author hkb
 */
public class Point3D extends Tuple3<Double,Double,Double> {

	/**
	 * Creates a new point.
	 * 
	 * @param x The first coordinate of the point.
	 * @param y The second coordinate of the point.
	 * @param z The third coordinate of the point.
	 */
	public Point3D(Double x, Double y, Double z) {
		super(x, y, z);
	}

	/**
	 * Creates a point from a 3-tuple.
	 * 
	 * @param tuple The 3-tuple to create the point from.
	 */
	public Point3D(Tuple3<Double,Double,Double> tuple) {
		this(tuple.x, tuple.y, tuple.z);
	}
	
	/**
	 * The distance between this point and the other.
	 * 
	 * @param other The point to find the distance to.
	 * @return The distance between the points.
	 */
	public double distance(Point3D other) {		
		return Math.abs(other.asVector().subtract(this.asVector()).length());
	}
	
	/**
	 * This point as a vector.
	 * 
	 * @return This point vector.
	 */
	public Vector3D asVector() {
		return new Vector3D(this);
	}
}
