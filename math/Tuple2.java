package math;

/**
 * Generic 2-tuple class.
 * 
 * @author hkb
 *
 * @param <T1> The type of the tuples first element.
 * @param <T2> The type of the tuples second element.
 */
public class Tuple2<T1, T2> implements Comparable<Tuple2<T1, T2>> {
	
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
	
	@Override
	public boolean equals(Object other) {
		if(other == null)
			return false;
		if(!(other instanceof Tuple2))
			return false;
		
		// java got generics all wrong so we don't really know what's going to happen. 
		try {
			Tuple2 otherTuple = (Tuple2) other;
		
			return this.x.equals(otherTuple.x) && this.y.equals(otherTuple.y);
		} catch (ClassCastException e) {
			return false;
		}
	}
	
	@Override
	public int compareTo(Tuple2<T1, T2> other) {
		try {
			Comparable x = (Comparable) this.x;
			Comparable y = (Comparable) this.y;
			
			int xs = x.compareTo(other.x);
			
			return (xs != 0) ? xs : y.compareTo(other.y);
		} catch (ClassCastException e) {
			throw new IllegalArgumentException(other + " isn't comparable to " + this);
		}
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
