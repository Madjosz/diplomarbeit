package de.mfabricius.diplom.data;

import static de.mfabricius.diplom.semilagrange.Problem.EPS;

import java.util.Map;

import de.mfabricius.diplom.util.IO;


/**
 * This class represent a point in a plane spanned by a {@link Triangle} and calculates the barycentric coordinates of
 * this point respective to the {@link Triangle}. Also there are methods provided to calculate the value of the linear
 * interpolant at the given point with values given on the vertices of the {@link Triangle}.
 * 
 * @author Martin Fabricius
 */
public class Barycentric {

    /** The triangle for which the coordinates are calculated. */
    public final Triangle t;
    /** The vector whose coordiantes are calculated. */
    public final Vector v;
    private final double[] bary = new double[3];
    /** Indicates that the given vector is a vertex of the triangle. */
    public final boolean isVertex;

    /**
     * Creates a new triple of barycentric coordinates for the given vector relative to the given triangle.
     * 
     * @param t
     *            the triangle
     * @param v
     *            the vector
     */
    public Barycentric(Triangle t, Vector v) {
        this(t, v, false);
    }

    /**
     * Creates a new triple of barycentric coordinates for the given {@link Vector} relative to the given
     * {@link Triangle}.
     * 
     * @param t
     *            the triangle
     * @param v
     *            the vector
     * @param isVertex
     *            is the vector a vertex of the triangle?
     */
    public Barycentric(Triangle t, Vector v, boolean isVertex) {
        this.t = t;
        this.v = v;
        barycentric();
        this.isVertex = isVertex;
    }

    /** @return <code>true</code> if the all barycentric coordinates are numerically non-negative */
    public boolean isInner() {
        return bary[0] > -EPS && bary[1] > -EPS && bary[2] > -EPS;
    }

    /** @return the vertex of the {@link Triangle} for which the corresponding coordinate has the lowest value */
    public Vector getLowest() {
        final double min = Math.min(bary[0], Math.min(bary[1], bary[2]));
        if (min == bary[0]) return t.v1();
        if (min == bary[1]) return t.v2();
        return t.v3();
    }

    /** Calculates the barycentric coordiantes. */
    private void barycentric() {
        final double area = t.area();
        bary[0] = signedArea(new Triangle(v, t.v2(), t.v3())) / area;
        bary[1] = signedArea(new Triangle(t.v1(), v, t.v3())) / area;
        bary[2] = signedArea(new Triangle(t.v1(), t.v2(), v)) / area;
    }

    /**
     * Calculates the area of the given {@link Triangle} in the plane of the reference triangle with sign, indicating
     * same or converse machining direction.
     * 
     * @param t
     *            the triangle to calculate the area
     * @return the signed area
     */
    private double signedArea(Triangle t) {
        return Vector.angle(this.t.normal(), t.normal()) > 1.5 ? -t.area() : t.area();
    }

    /**
     * Computes the linearly interpolated value of the function {@code u} at the {@link Vector} from the values of the
     * vertices of the {@link Triangle}.
     * 
     * @param u
     *            the function of the nodal values
     * @return the interpolated value at {@code v}
     */
    public double interpolate(Map<Vector, Double> u) {
        final Vector[] v = t.vertices();
        double value = 0.0;
        for (int i = 0; i < v.length; i++)
            value += bary[i] * u.get(v[i]);
        return value;
    }

    /*
     * object methods
     */

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj instanceof Barycentric) {
            final Barycentric b = (Barycentric) obj;
            return t.equals(b.t) && v.equals(b.v);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return t.hashCode() + v.hashCode();
    }

    @Override
    public String toString() {
        return "(" + IO.format(bary[0]) + " : " + IO.format(bary[1]) + " : " + IO.format(bary[2]) + ")";
    }
}
