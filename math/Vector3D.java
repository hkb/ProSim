package math;

/**
 * A basic 3-dimensional vector.
 * 
 * @author hkb
 */
public class Vector3D extends Tuple3<Double,Double,Double> {
	
	/**
	 * 
	 */
	public Vector3D() {
		super(0.0, 0.0, 0.0);
	}

	/**
	 * Creates a new 3-dimensional vector.
	 * 
	 * @param x The first coordinate of the vector.
	 * @param y The second coordinate of the vector.
	 * @param z The third coordinate of the vector.
	 */
	public Vector3D(Double x, Double y, Double z) {
		super(x, y, z);
	}
	
	/**
	 * Creates a vector from any 3-tuple of doubles.
	 * 
	 * @param tuple The tuple to create the vector from.
	 */
	public Vector3D(Tuple3<Double,Double,Double> tuple) {
		this(tuple.x, tuple.y, tuple.z);
	}

	/**
	 * Adds this vector with the other.
	 * 
	 * @param other The vector to add with this one.
	 * @return The vector as the result f the addition.
	 */
	public Vector3D add(Vector3D other) {
		return new Vector3D(this.x + other.x, this.y + other.y, this.z + other.z);
	}

	/**
	 * Subtracts the other vector from this.
	 * 
	 * @param other The vector to subtract from this.
	 * @return The vector as the result of the subtraction.
	 */
	public Vector3D subtract(Vector3D other) {
		return new Vector3D(this.x - other.x, this.y - other.y, this.z - other.z);
	}	
	
	/**
	 * Scales the vector by a given factor.
	 * 
	 * @param scale the factor to scale with.
	 * @return The scaled vector.
	 */
	public Vector3D scale(double scale) {
		return new Vector3D(this.x * scale, this.y * scale, this.z * scale);
	}
	
	/**
	 * Returns the length of this vector.
	 * 
	 * @return The length of this vector.
	 */
	public double length() {
		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
	}
	
	/**
	 * Normalises the vector.
	 * 
	 * @return The normalised version of this vector.
	 */
	public Vector3D norm() {
		double length = this.length();
		
		return new Vector3D(this.x / length, this.y / length, this.z / length);
	}
	
	/**
	 * Computes the dot product between this and the other vector.
	 * 
	 * @param other The other vector to compute the dot product with.
	 * @return The dot product.
	 */
	public double dot(Vector3D other) {
		return this.x * other.x + this.y * other.y + this.z * other.z;
	}
	
	/**
	 * The cross product between this and the other vector.
	 * 
	 * @param The other vector to compute the cross product with.
	 * @return The cross product vector.
	 */
	public Vector3D cross(Vector3D other) {
		return new Vector3D(this.y * other.z - this.z * other.y,
						    this.z * other.x - this.x * other.z,
						    this.x * other.y - this.y * other.x);
	}
	
	/**
	 * Computes a vector from this point to the other.
	 * 
	 * @param other The vector to compute a vector to.
	 * @return The vector from this vector to the other.
	 */
	public Vector3D vectorTo(Vector3D other) {
		return other.subtract(this);
	}
	
	/**
	 * Returns the other vector projected onto this vector.
	 * 
	 * @param other The vector to project onto this vector.
	 * @return The projected vector.
	 */
	public Vector3D projectOnto(Vector3D other) {
		Vector3D norm = this.norm();
		
		return norm.scale(other.dot(norm));
	}
}
