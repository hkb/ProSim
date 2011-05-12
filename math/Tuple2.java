package math;

/**
 * Generic 2-tuple class.
 * 
 * @author hkb
 *
 * @param <T1> The type of the tuples first element.
 * @param <T2> The type of the tuples second element.
 */
public class Tuple2<T1, T2> {
	
	public T1 x;
	public T2 y;

	/**
	 * Creates a new 2-tuple.
	 * 
	 * @param x First element in the tuple.
	 * @param y Second element of the tuple.
	 */
	public Tuple2 (T1 x, T2 y) {
		this.x = x;
		this.y = y;
	}
	
	/**
	 * Determines if the elements of this and the other tuple are equal.
	 * 
	 * @param other The other tuple.
	 * @return 
	 */
	public boolean equals(Tuple2<T1, T2> other) {
		return this.x.equals(other.x) && this.y.equals(other.y);
	}
	
	@Override
	public int hashCode() {
		return (37 * this.x.hashCode()) ^ this.y.hashCode(); // 37 - my favourite prime number
	}
	
	@Override
	public String toString() {
		return "(" + this.x + ", " + this.y + ")";
	}
}
