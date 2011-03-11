package boundingVolume;

import math.matrix.TransformationMatrix;

public interface BoundingVolume {

	/**
	 * Does the bounding volume overlap with the other.
	 * 
	 * @param other The volume to test for overlap with.
	 * @return true if the bounding volumes overlap else false
	 */
	public boolean isOverlaping(BoundingVolume other);
	
	/**
	 * Computes a new bounding volume containing the combined volumes of the 
	 * given volumes.
	 * 
	 * @param left A BoundingVolume.
	 * @param right A BoundingVolume.
	 * @return A combined bounding volume.
	 */
	public BoundingVolume combine(BoundingVolume left, BoundingVolume right);
	
	/**
	 * Compute a transformed version of the bounding volume given the transformation matrix.
	 * 
	 * @param transformationMatrix The matrix to transform the volume with.
	 * @return A version of the bounding volume transformed accordingly to the transformation matrix.
	 */
	public BoundingVolume transform(TransformationMatrix transformationMatrix);
}
