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
public class Tuple3<T1, T2, T3> implements Comparable<Tuple3<T1, T2, T3>> {
	
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
	
	@Override
	public boolean equals(Object other) {
		if(other == null)
			return false;
		if(!(other instanceof Tuple3))
			return false;
		
		// java got generics all wrong so we don't really know what's going to happen. 
		try {
			Tuple3 otherTuple = (Tuple3) other;
		
			return this.x.equals(otherTuple.x) && this.y.equals(otherTuple.y) && this.z.equals(otherTuple.z);
		} catch (ClassCastException e) {
			return false;
		}
	}
	
	@Override
	public int compareTo(Tuple3<T1, T2, T3> other) {
		try {
			Comparable x = (Comparable) this.x;
			Comparable y = (Comparable) this.y;
			Comparable z = (Comparable) this.z;
			
			int xs = x.compareTo(other.x);
			
			if(xs != 0)
				return xs;
			
			int ys = y.compareTo(other.y);
			
			return (ys != 0) ? ys : z.compareTo(other.z);
		} catch (ClassCastException e) {
			throw new IllegalArgumentException(other + " isn't comparable to " + this);
		}
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