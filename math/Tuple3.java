package math;

/**
 * Generic 2-tuple class.
 * 
 * @author hkb
 *
 * @param <T1> The type of the tuples first element.
 * @param <T2> The type of the tuples second element.
 * @param <T3> The type of the tuples third element.
 */
public class Tuple3<T1, T2, T3> {
	
	public T1 x;
	public T2 y;
	public T3 z;

	/**
	 * Creates a new 3-tuple.
	 * 
	 * @param x First element in the tuple.
	 * @param y Second element of the tuple.
	 * @param z Third element of the tuple.
	 */
	public Tuple3 (T1 x, T2 y, T3 z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	/**
	 * Determines if the elements of this and the other tuple are equal.
	 * 
	 * @param other The other tuple.
	 * @return 
	 */
	public boolean equals(Tuple3<T1, T2, T3> other) {
		return this.x == null && other.y == null ||
		       this.y == null && other.y == null ||
		       this.z == null && other.z == null ||
		       this.x.equals(other.x) && this.y.equals(other.y) && this.z.equals(other.z);
	}
	
	@Override
	public int hashCode() {
		return ((37 * this.x.hashCode()) ^ this.y.hashCode()) ^ this.z.hashCode(); // 37 - my favourite prime number
	}
	
	@Override
	public String toString() {
		return "(" + this.x + ", " + this.y + ", " + this.z + ")";
	}
}