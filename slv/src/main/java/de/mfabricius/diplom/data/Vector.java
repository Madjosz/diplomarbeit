package de.mfabricius.diplom.data;

import static de.mfabricius.diplom.semilagrange.Problem.EPS;
import static java.lang.Math.acos;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.Arrays;

import Jama.Matrix;


/**
 * This class represents a vector in the 3-dimensional euclidean space R³ defined by its 3 real cartesian coordiantes.
 * Also there are basic vector space operations provided.
 * 
 * @author Martin Fabricius
 */
public class Vector {

    /** The zero {@link Vector}. */
    public static final Vector ZERO = new Vector(0., 0., 0.);

    private final double[] coords;

    /**
     * Creates a new {@link Vector} representing the three coordinates.
     * 
     * @param x
     *            the x coordinate
     * @param y
     *            the y coordinate
     * @param z
     *            the z coordinate
     */
    public Vector(double x, double y, double z) {
        coords = new double[] { x, y, z };
    }

    /**
     * Creates a new {@link Vector} copying the given Vector.
     * 
     * @param v
     *            the Vector
     */
    public Vector(Vector v) {
        this(Arrays.copyOf(v.coords, 3));
    }

    private Vector(double[] coords) {
        this.coords = coords;
    }

    /** @return the x coordinate */
    public double x() {
        return coords[0];
    }

    /** @return the y coordinate */
    public double y() {
        return coords[1];
    }

    /** @return the z coordinate */
    public double z() {
        return coords[2];
    }

    /**
     * @param v
     *            the Vector to add
     * @return new Vector with the coordinate-wise sum
     */
    public Vector add(Vector v) {
        return new Vector(coords[0] + v.x(), coords[1] + v.y(), coords[2] + v.z());
    }

    /**
     * @param v
     *            the Vector to subtract
     * @return new Vector with the coordinate-wise difference
     */
    public Vector sub(Vector v) {
        return new Vector(coords[0] - v.x(), coords[1] - v.y(), coords[2] - v.z());
    }

    /**
     * @param factor
     *            the factor to multiply
     * @return new Vector with the coordinate-wise product
     */
    public Vector mult(double factor) {
        return new Vector(factor * coords[0], factor * coords[1], factor * coords[2]);
    }

    /**
     * @param m
     *            the Matrix to multiply with from left
     * @return new Vector as result of Matrix multiplication
     */
    public Vector mult(Matrix m) {
        final double[] v = new double[3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                v[i] += m.get(i, j) * coords[j];
        return new Vector(v);
    }

    /** @return the dot product of the two position {@link Vector}s */
    public static double skprod(Vector v1, Vector v2) {
        return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
    }

    /** @return the cross product of the two position {@link Vector}s */
    public static Vector xprod(Vector v1, Vector v2) {
        return new Vector(v1.y() * v2.z() - v1.z() * v2.y(), //
                v1.z() * v2.x() - v1.x() * v2.z(), //
                v1.x() * v2.y() - v1.y() * v2.x());
    }

    /** @return the square of the euclidean norm of the position {@link Vector} */
    public double norm22() {
        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
    }

    /** @return the euclidean norm of the position {@link Vector} */
    public double norm2() {
        return sqrt(norm22());
    }

    /** @return unit {@link Vector} pointing in the same direction */
    public Vector unitVector() {
        return mult(1.0 / norm2());
    }

    /** @return the euclidean distance to the given position {@link Vector} */
    public double distance(Vector v) {
        return sqrt(pow(coords[0] - v.x(), 2.) + pow(coords[1] - v.y(), 2.) + pow(coords[2] - v.z(), 2.));
    }

    /** @return the angle between the given {@link Vector}s, {@code 0.0} if one is the zero vector */
    public static double angle(Vector p1, Vector p2) {
        final double denom = p1.norm2() * p2.norm2();
        return denom > EPS ? acos(skprod(p1, p2) / denom) : 0.0;
    }

    /* object methods     */

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj instanceof Vector) {
            final Vector v = (Vector) obj;
            if (Double.compare(coords[0], v.x()) == 0 && Double.compare(coords[1], v.y()) == 0
                    && Double.compare(coords[2], v.z()) == 0)
                return true;
        }
        return false;
    }

    @Override
    public int hashCode() {
        long bits = 7L;
        bits = 31L * bits + Double.doubleToLongBits(coords[0]);
        bits = 31L * bits + Double.doubleToLongBits(coords[1]);
        bits = 31L * bits + Double.doubleToLongBits(coords[2]);
        return (int) (bits ^ (bits >> 32));
    }

    @Override
    public String toString() {
        return Arrays.toString(coords);
    }
}
