package de.mfabricius.diplom.util;

import static de.mfabricius.diplom.data.Vector.skprod;

import Jama.Matrix;
import de.mfabricius.diplom.data.Vector;


/**
 * Factory class, that creates transportation fields beta.
 * 
 * @author Martin Fabricius
 */
public class BetaFactory {

    /**
     * Gives a vector field with the tangential part of the given direction {@link Vector} when projecting on a sphere.
     * 
     * @param dir
     *            the direction of the constant vector field
     * @return tangential vector field
     */
    public static final R4toR3 ConstantForSphere(Vector dir) {
        return (t, v) -> {
            final Vector n = v.unitVector();
            return dir.sub(n.mult(skprod(n, dir)));
        };
    }

    /**
     * Gives the velocity field of a rotating sphere around the given axis with the given angular speed.
     * 
     * @param axis
     *            the axis of rotation
     * @param h
     *            the angular speed
     * @return tangential vector field
     */
    public static final R4toR3 RotationForSphere(Vector axis, double h) {
        return PhiFactory.Rotation(axis.unitVector(), h).velocity;
    }

    /**
     * Gives the velocity field of a rotating sphere around the z-axis with the given angular speed.
     * 
     * @param h
     *            the angular speed
     * @return tangential vector field
     */
    public static final R4toR3 RotationZForSphere(double h) {
        return PhiFactory.RotationZ(h).velocity;
    }

    /**
     * Creates a tangential field for a torus which goes around the small circle.
     * 
     * @param h
     *            the tangential speed
     * @return tangential vector field
     */
    public static final R4toR3 UpwindTorus(double h) {
        return (t, v) -> {
            final Vector c = new Vector(v.x(), v.y(), 0).unitVector();
            final Vector r = v.sub(c).unitVector().mult(h);
            final Matrix m = new Matrix(new double[][] { { c.y() * c.y(), -c.x() * c.y(), -c.x() },
                    { -c.x() * c.y(), c.x() * c.x(), -c.y() }, { c.x(), c.y(), 0. } });
            return r.mult(m);
        };
    }
}