package math.matrix;

import edu.math.Vector;

/**
 * A matrix for tarnsforming coordinates between coordinate systems.
 * 
 * @author hkb
 */
public class TransformationMatrix {
	public double a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34; // the entries of the matrix
	
	/**
	 * Create unit transformation matrix.
	 */
	public TransformationMatrix() {
		a11 = 1; a12 = 0; a13 = 0; a14 = 0;
		a21 = 0; a22 = 1; a23 = 0; a24 = 0;
		a31 = 0; a32 = 0; a33 = 1; a34 = 0;
	}
	
	/**
	 * Create transformation matrix transformation by some position.
	 * 
	 * @param x
	 * @param y
	 * @param z
	 */
	public TransformationMatrix(double x, double y, double z) {
		this(1, 0, 0, x, 
		     0, 1, 0, y, 
		     0, 0, 1, z);
	}
	
	/**
	 * Create transformation matrix from the specifiend fields.
	 * 
	 * @param a11
	 * @param a12
	 * @param a13
	 * @param a14
	 * @param a21
	 * @param a22
	 * @param a23
	 * @param a24
	 * @param a31
	 * @param a32
	 * @param a33
	 * @param a34
	 */
	public TransformationMatrix(double a11, double a12, double a13, double a14,
								double a21, double a22, double a23, double a24,
								double a31, double a32, double a33, double a34) {
		this.a11 = a11; this.a12 = a12; this.a13 = a13; this.a14 = a14;
		this.a21 = a21; this.a22 = a22; this.a23 = a23; this.a24 = a24;
		this.a31 = a31; this.a32 = a32; this.a33 = a33; this.a34 = a34;
	}
	
	/**
	 * Copy other transformation matrix.
	 * 
	 * @param m The matrix to copy.
	 */
	public TransformationMatrix(TransformationMatrix m) {
		a11 = m.a11; a12 = m.a12; a13 = m.a13; a14 = m.a14;
		a21 = m.a21; a22 = m.a22; a23 = m.a23; a24 = m.a24;
		a31 = m.a31; a32 = m.a32; a33 = m.a33; a34 = m.a34;
	}
	
	/**
	 * Create transformation matrix from CTNode children transformations. TODO improve
	 * 
	 * @param l
	 * @param r
	 */
	public TransformationMatrix(TransformationMatrix l, TransformationMatrix r) {
		a11 = l.a11*r.a11 + l.a12*r.a21 + l.a13*r.a31;
		a12 = l.a11*r.a12 + l.a12*r.a22 + l.a13*r.a32;
		a13 = l.a11*r.a13 + l.a12*r.a23 + l.a13*r.a33;
		a14 = l.a11*r.a14 + l.a12*r.a24 + l.a13*r.a34 + l.a14;
		a21 = l.a21*r.a11 + l.a22*r.a21 + l.a23*r.a31;
		a22 = l.a21*r.a12 + l.a22*r.a22 + l.a23*r.a32;
		a23 = l.a21*r.a13 + l.a22*r.a23 + l.a23*r.a33;
		a24 = l.a21*r.a14 + l.a22*r.a24 + l.a23*r.a34 + l.a24;
		a31 = l.a31*r.a11 + l.a32*r.a21 + l.a33*r.a31;
		a32 = l.a31*r.a12 + l.a32*r.a22 + l.a33*r.a32;
		a33 = l.a31*r.a13 + l.a32*r.a23 + l.a33*r.a33;
		a34 = l.a31*r.a14 + l.a32*r.a24 + l.a33*r.a34 + l.a34;

	}


	/* 
	 * This rotation matrix is right multiplied by the rotation matrix m. 
	 * Original content of the matrix destroyed
	 */
	public void multR(TransformationMatrix m) {
		double b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33, b34;
		b11 = a11*m.a11 + a12*m.a21 + a13*m.a31;
		b12 = a11*m.a12 + a12*m.a22 + a13*m.a32;
		b13 = a11*m.a13 + a12*m.a23 + a13*m.a33;
		b14 = a11*m.a14 + a12*m.a24 + a13*m.a34 + a14;
		b21 = a21*m.a11 + a22*m.a21 + a23*m.a31;
		b22 = a21*m.a12 + a22*m.a22 + a23*m.a32;
		b23 = a21*m.a13 + a22*m.a23 + a23*m.a33;
		b24 = a21*m.a14 + a22*m.a24 + a23*m.a34 + a24;

		b31 = a31*m.a11 + a32*m.a21 + a33*m.a31;
		b32 = a31*m.a12 + a32*m.a22 + a33*m.a32;
		b33 = a31*m.a13 + a32*m.a23 + a33*m.a33;
		b34 = a31*m.a14 + a32*m.a24 + a33*m.a34 + a34;

		a11 = b11; a12 = b12; a13 = b13; a14 = b14;
		a21 = b21; a22 = b22; a23 = b23; a24 = b24;
		a31 = b31; a32 = b32; a33 = b33; a34 = b34;
	}

	/* 
	 * This rotation matrix is left multiplied by the rotation matrix m. 
	 * Original content of the matrix id destroyed
	 */
	public void multL(TransformationMatrix m) {
		double b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33, b34;
		b11 = m.a11*a11 + m.a12*a21 + m.a13*a31;
		b12 = m.a11*a12 + m.a12*a22 + m.a13*a32;
		b13 = m.a11*a13 + m.a12*a23 + m.a13*a33;
		b14 = m.a11*a14 + m.a12*a24 + m.a13*a34 + m.a14;
		b21 = m.a21*a11 + m.a22*a21 + m.a23*a31;
		b22 = m.a21*a12 + m.a22*a22 + m.a23*a32;
		b23 = m.a21*a13 + m.a22*a23 + m.a23*a33;
		b24 = m.a21*a14 + m.a22*a24 + m.a23*a34 + m.a24;

		b31 = m.a31*a11 + m.a32*a21 + m.a33*a31;
		b32 = m.a31*a12 + m.a32*a22 + m.a33*a32;
		b33 = m.a31*a13 + m.a32*a23 + m.a33*a33;
		b34 = m.a31*a14 + m.a32*a24 + m.a33*a34 + m.a34;

		a11 = b11; a12 = b12; a13 = b13; a14 = b14;
		a21 = b21; a22 = b22; a23 = b23; a24 = b24;
		a31 = b31; a32 = b32; a33 = b33; a34 = b34;
	}

	/**
	 * Apply the rotation by the given angle to the matrix.
	 * 
	 * @param angle The angle to change the rotation with.
	 */
	public void rotate(double angle) {		
		// precompute values
		Vector t = new Vector(this.a14, this.a24, this.a34);
		if (t.length() > 0) t = t.norm();
		
		double x = t.x();
		double y = t.y();
		double z = t.z();
		
		double s = Math.sin(angle);
		double c = Math.cos(angle);
		double d = 1 - c;
		
		// precompute to avoid double computations
		double dxy = d*x*y;
		double dxz = d*x*z;
		double dyz = d*y*z;
		double xs = x*s;
		double ys = y*s;
		double zs = z*s;
		
		// update matrix
		a11 = d*x*x+c; a12 = dxy-zs;  a13 = dxz+ys;
		a21 = dxy+zs;  a22 = d*y*y+c; a23 = dyz-xs;
		a31 = dxz-ys;  a32 = dyz+xs;  a33 = d*z*z+c;
	}
	
	/**
	 * Transform a vector. TODO improve
	 * @param v
	 * @return
	 */
	public Vector transform(Vector v) {
		float x = v.x();
		float y = v.y();
		float z = v.z();
		
		return new Vector(a11*x + a12*y + a13*z + a14, 
						  a21*x + a22*y + a23*z + a24,
						  a31*x + a32*y + a33*z + a34);
	}
	
	@Override
	public String toString() {
		return String.format("[%s, %s, %s, %s]\n[%s, %s, %s, %s]\n[%s, %s, %s, %s]", 
							  a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34);
	}
}
