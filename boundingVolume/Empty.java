package boundingVolume;

import math.matrix.TransformationMatrix;

/**
 * An empty bounding volume.
 * 
 * @author hkb
 */
public class Empty implements BoundingVolume {

	@Override
	public BoundingVolume combine(BoundingVolume other) {
		return other;
	}

	@Override
	public boolean isOverlaping(BoundingVolume other) {
		return false;
	}

	@Override
	public BoundingVolume transform(TransformationMatrix transformationMatrix) {
		return this;
	}

	@Override
	public float volume() {
		return 0;
	}

}
