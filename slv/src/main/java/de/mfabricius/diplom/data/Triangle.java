package de.mfabricius.diplom.data;

import static de.mfabricius.diplom.data.Vector.xprod;
import static java.lang.Math.max;
import static java.lang.Math.sqrt;


/**
 * This class represents a triangle in the 3-dimensional euclidean space R³ defined by its 3 vertices. Also there are
 * basic calculation methods provided for the following properties
 * <ul>
 * <li>area</li>
 * <li>circumfence</li>
 * <li>diameter</li>
 * <li>incircle radius</li>
 * <li>circumscribed circle radius</li>
 * <li>length of the sides</li>
 * <li>ratios of certain properties</li>
 * </ul>
 * Also there are methods for calculating the integral and the L1/L2 norm of the linear interpolant of a function with
 * values given at the vertices.
 * 
 * @author Martin Fabricius
 */
public class Triangle {

    private final Vector v1;
    private final Vector v2;
    private final Vector v3;
    private final Vector normal;

    /**
     * Creates a new {@link Triangle} with the given vertices
     * 
     * @param v1
     *            vertex 1
     * @param v2
     *            vertex 2
     * @param v3
     *            vertex 3
     */
    public Triangle(Vector v1, Vector v2, Vector v3) {
        this.v1 = v1;
        this.v2 = v2;
        this.v3 = v3;
        normal = xprod(v2.sub(v1), v3.sub(v1));// .unitVector();
    }

    /**
     * Creates a new {@link Triangle} copying the given triangle. The vertices will be copied by reference.
     * 
     * @param t
     *            the triangle
     */
    public Triangle(Triangle t) {
        this(t.v1, t.v2, t.v3);
    }

    /** @return vertex 1 */
    public Vector v1() {
        return v1;
    }

    /** @return vertex 2 */
    public Vector v2() {
        return v2;
    }

    /** @return vertex 3 */
    public Vector v3() {
        return v3;
    }

    /** @return an array consisting of the three vertices v1, v2 and v3 */
    public Vector[] vertices() {
        return new Vector[] { v1, v2, v3 };
    }

    /**
     * Returns the vertex, which is not on the given edge.
     * 
     * @param e
     *            the edge
     * @return the vertex not on the edge
     */
    public Vector not(Edge e) {
        if (e.equals(new Edge(v1, v2))) return v3;
        if (e.equals(new Edge(v2, v3))) return v1;
        if (e.equals(new Edge(v3, v1))) return v2;
        return null;
    }

    /** @return normal to this triangle */
    public Vector normal() {
        return normal;
    }

    /** @return length of edge from vertex 1 to vertex 2 */
    public double l12() {
        return v1.distance(v2);
    }

    /** @return length of edge from vertex 2 to vertex 3 */
    public double l23() {
        return v2.distance(v3);
    }

    /** @return length of edge from vertex 3 to vertex 1 */
    public double l31() {
        return v3.distance(v1);
    }

    /** @return length of the circumfence of this triangle */
    public double circumfence() {
        return circumfence(l12(), l23(), l31());
    }

    /** @return area of this triangle */
    public double area() {
        return area(l12(), l23(), l31());
    }

    /** @return the diameter of the triangle (the longest edge) */
    public double diameter() {
        return diameter(l12(), l23(), l31());
    }

    /** @return radius of the circumscribed circle */
    public double umkreis() {
        return umkreis(l12(), l23(), l31());
    }

    /** @return radius of the incircle */
    public double inkreis() {
        return inkreis(l12(), l23(), l31());
    }

    /** @return the ratio of the diameter and the radius of the incircle */
    public double diamikRatio() {
        return diamikRatio(l12(), l23(), l31());
    }

    /** @return the ratio of the radii of the circumscribed circle and the incircle */
    public double ukikRatio() {
        return ukikRatio(l12(), l23(), l31());
    }

    /**
     * Computes the integral of the linear interpolant with nodal values at the vertices of the triangle.
     * 
     * @param f1
     *            the value at v1
     * @param f2
     *            the value at v2
     * @param f3
     *            the value at v3
     * @return the value of the integral
     */
    public double integral(double f1, double f2, double f3) {
        final double transf = xprod(v2.sub(v1), v3.sub(v1)).norm2();
        return transf * (f1 + f2 + f3) / 6.;
    }

    /**
     * Computes the L1-Norm (integral of abs(f)) of the linear interpolant with nodal values at the vertices of the
     * triangle.
     * 
     * @param f1
     *            the value at v1
     * @param f2
     *            the value at v2
     * @param f3
     *            the value at v3
     * @return the L1-Norm
     */
    public double normL1(double f1, double f2, double f3) {
        final double transf = xprod(v2.sub(v1), v3.sub(v1)).norm2();
        if (f1 >= 0. && f2 >= 0. && f3 >= 0.) return transf * (f1 + f2 + f3) / 6.;
        if (f1 < 0. && f2 < 0. && f3 < 0.) return -transf * (f1 + f2 + f3) / 6.;
        if (f1 >= 0. && f2 >= 0. && f3 < 0.)
            return transf * (f1 + f2 + f3 * (1 - 2 * f3 * f3 / (f3 - f1) / (f3 - f2))) / 6.;
        if (f1 >= 0. && f2 < 0. && f3 >= 0.)
            return transf * (f1 + f2 * (1 - 2 * f2 * f2 / (f2 - f1) / (f2 - f3)) + f3) / 6.;
        if (f1 < 0. && f2 >= 0. && f3 >= 0.)
            return transf * (f1 * (1 - 2 * f1 * f1 / (f1 - f2) / (f1 - f3)) + f2 + f3) / 6.;
        if (f1 >= 0. && f2 < 0. && f3 < 0.)
            return -transf * (f1 * (1 - 2 * f1 * f1 / (f1 - f2) / (f1 - f3)) + f2 + f3) / 6.;
        if (f1 < 0. && f2 >= 0. && f3 < 0.)
            return -transf * (f1 + f2 * (1 - 2 * f2 * f2 / (f2 - f1) / (f2 - f3)) + f3) / 6.;
        else // if (f1 < 0. && f2 < 0. && f3 >= 0.)
            return -transf * (f1 + f2 + f3 * (1 - 2 * f3 * f3 / (f3 - f1) / (f3 - f2))) / 6.;
    }

    /**
     * Computes the square of the L²-Norm (integral of f²) of the linear interpolant with nodal values at the vertices
     * of the triangle.
     * 
     * @param f1
     *            the value at v1
     * @param f2
     *            the value at v2
     * @param f3
     *            the value at v3
     * @return the L2-Norm
     */
    public double normL22(double f1, double f2, double f3) {
        final double transf = xprod(v2.sub(v1), v3.sub(v1)).norm2();
        return transf * (f1 * f1 + f1 * f2 + f1 * f3 + f2 * f2 + f2 * f3 + f3 * f3) / 12.;
    }

    /**
     * Calculates the area of a {@link Triangle} with the given side lengths with the Heron's formula. If the inputs are
     * invalid or no Triangle can be constructed with this inputs, 0.0 will be returned.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the area of the {@link Triangle} or 0.0
     */
    private static double area(double a, double b, double c) {
        final double s = circumfence(a, b, c) / 2.;
        final double rad = s * (s - a) * (s - b) * (s - c);
        return rad < 0.0 ? 0.0 : sqrt(rad);
    }

    /**
     * Calculates the circumfence of a {@link Triangle} with the given side lengths. There are no checks if a real
     * {@link Triangle} can be constructed with this values.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the circumfence of the {@link Triangle}
     */
    private static double circumfence(double a, double b, double c) {
        return a + b + c;
    }

    /**
     * Calculates the diameter of a {@link Triangle} with the given side lengths, which is the length of the longest
     * side. There are no checks if a real {@link Triangle} can be constructed with this values.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the diameter of the {@link Triangle}
     */
    private static double diameter(double a, double b, double c) {
        return max(a, max(b, c));
    }

    /**
     * Calculates the radius of the circumscribed circle of a {@link Triangle} with the given side lengths. There are no
     * checks if a real {@link Triangle} can be constructed with this values.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the radius of the circumscribed circle of the {@link Triangle}
     */
    private static double umkreis(double a, double b, double c) {
        return a * b * c / (4. * area(a, b, c));
    }

    /**
     * Calculates the radius of the inscribed circle of a {@link Triangle} with the given side lengths. There are no
     * checks if a real {@link Triangle} can be constructed with this values.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the radius of the inscribed circle of the {@link Triangle}
     */
    private static double inkreis(double a, double b, double c) {
        return 2. * area(a, b, c) / circumfence(a, b, c);
    }

    /**
     * Calculates the ratio of the diameter and the radius of the inscribed circle of a {@link Triangle} with the given
     * side lengths. There are no checks if a real {@link Triangle} can be constructed with this values.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the ratio diameter / incircle
     */
    private static double diamikRatio(double a, double b, double c) {
        return diameter(a, b, c) / inkreis(a, b, c);
    }

    /**
     * Calculates the ratio of the radii of the circumscribed with the inscribed circle of a {@link Triangle} with the
     * given side lengths. There are no checks if a real {@link Triangle} can be constructed with this values.
     * 
     * @param a
     *            length of side a
     * @param b
     *            length of side b
     * @param c
     *            length of side c
     * @return the ratio circumscribed / incircle
     */
    private static double ukikRatio(double a, double b, double c) {
        final double s = circumfence(a, b, c) / 2.;
        return a * b * c / (4. * (s - a) * (s - b) * (s - c));
    }

    /*
     * object methods
     */

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj instanceof Triangle) {
            final Triangle t = (Triangle) obj;
            if (v1.equals(t.v1)) {
                return (v2.equals(t.v2) && v3.equals(t.v3)) || (v2.equals(t.v3) && v3.equals(t.v2));
            } else if (v1.equals(t.v2)) {
                return (v2.equals(t.v1) && v3.equals(t.v3)) || (v2.equals(t.v3) && v3.equals(t.v1));
            } else if (v1.equals(t.v3))
                return (v2.equals(t.v2) && v3.equals(t.v1)) || (v2.equals(t.v1) && v3.equals(t.v2));
        }
        return false;
    }

    @Override
    public int hashCode() {
        return v1.hashCode() + v2.hashCode() + v3.hashCode();
    }

    @Override
    public String toString() {
        return "{" + v1.toString() + "," + v2.toString() + "," + v3.toString() + "}";
    }
}
