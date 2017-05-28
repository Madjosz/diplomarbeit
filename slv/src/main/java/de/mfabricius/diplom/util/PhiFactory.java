package de.mfabricius.diplom.util;

import static de.mfabricius.diplom.data.Vector.skprod;
import static de.mfabricius.diplom.data.Vector.xprod;
import static de.mfabricius.diplom.util.R4toR3.ZERO;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import de.mfabricius.diplom.data.Vector;


/**
 * Factory class, that creates movements for surfaces with velocity.
 * 
 * @author Martin Fabricius
 */
public class PhiFactory {

    /** @returns constant rotation around the z-axis with angular speed 1. */
    public static final PHI RotationZ = RotationZ(1);

    /**
     * Constant rotation around the z-axis with given angular speed.
     * 
     * @param speed
     *            angular speed
     * @return constant rotation around the z-axis with given angular speed
     */
    public static final PHI RotationZ(double speed) {
        return new PHI((t, x) -> {
            final double s = speed * t;
            return new Vector(cos(s) * x.x() - sin(s) * x.y(), sin(s) * x.x() + cos(s) * x.y(), x.z());
        }, (t, x) -> new Vector(-speed * x.y(), speed * x.x(), 0.));
    }

    /**
     * Constant rotation around a given axis with given angular speed.
     * 
     * @param axis
     *            the Rotation axis (has to be a unit vector)
     * @param speed
     *            the angular speed
     * @return constant rotation around a given axis with given angular speed.
     */
    public static final PHI Rotation(Vector axis, double speed) {
        return new PHI((t, x) -> {
            final double s = speed * t;
            return axis.mult(skprod(axis, x)).add(xprod(xprod(axis, x), axis).mult(cos(s)))
                    .add(xprod(axis, x).mult(sin(s)));
        }, (t, x) -> xprod(axis, x).mult(speed));
    }

    /**
     * Centrical expansion with the factor {@code a} from (0,0,0).
     * 
     * @param a
     *            the expansion factor
     * @return centrical expansion with the factor {@code a} from (0,0,0)
     */
    public static final PHI Expand(double a) {
        return new PHI((t, x) -> x.mult(a * t), (t, x) -> x.mult(1. / t));
    }

    /** Oscillation (pulsation) with (0,0,0) as center. */
    public static final PHI Oscillate = new PHI((t, x) -> x.mult(1 + Math.sin(t) / 2.),
            (t, x) -> x.mult(Math.cos(t) / (2. + Math.sin(t))));

    /**
     * Constant translation in the given {@code direction}.
     * 
     * @param direction
     *            the {@link Vector} of the direction
     * @return constant translation in the given {@code direction}
     */
    public static final PHI Translation(Vector direction) {
        return new PHI((t, x) -> x.add(direction.mult(t)), (t, x) -> direction);
    }

    /** @return no movement */
    public static final PHI Zero() {
        return new PHI(ZERO, ZERO);
    }

    /**
     * Class that describes a movement of a surface with velocity field as stated in section 2.4 in the diploma thesis.
     * 
     * @author Martin Fabricius
     */
    public static final class PHI {

        private PHI(R4toR3 move, R4toR3 velocity) {
            this.move = move;
            this.velocity = velocity;
        }

        /** The movement at time t for the given initial {@link Vector} x on Γ₀. */
        public final R4toR3 move;
        /** The velocity field at {@link Vector} x at time t. */
        public final R4toR3 velocity;
    }
}