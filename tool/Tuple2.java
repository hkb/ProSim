package tool;

/**
 * Generic 2-tuple class.
 * 
 * @author hkb
 *
 * @param <E1> The type of the tuples first element.
 * @param <E2> The type of the tuples second element.
 */
public class Tuple2<E1, E2> {
	
	public E1 e1;
	public E2 e2;

	/**
	 * Creates a new 2-tuple.
	 * 
	 * @param e1 First element in the tuple.
	 * @param e2 Second element of the tuple.
	 */
	public Tuple2 (E1 e1, E2 e2) {
		this.e1 = e1;
		this.e2 = e2;
	}
	
	/**
	 * Determines if the elements of this and the other tuple are equal.
	 * 
	 * @param other The other tuple.
	 * @return 
	 */
	public boolean equals(Tuple2<E1, E2> other) {
		return this.e1.equals(other.e1) && this.e2.equals(other.e2);
	}
	
	@Override
	public int hashCode() {
		return (37 * this.e1.hashCode()) ^ this.e2.hashCode(); // 37 - my favourite prime number
	}
	
	@Override
	public String toString() {
		return "(" + this.e1 + ", " + this.e2 + ")";
	}
}
