package de.mfabricius.diplom;

import static de.mfabricius.diplom.util.R4toR.circleAroundY;
import static de.mfabricius.diplom.util.R4toR.exponentialDistanceFromY;
import static de.mfabricius.diplom.util.Triangulation.SPHERE;
import static de.mfabricius.diplom.util.Triangulation.TORUS;

import java.io.FileWriter;
import java.io.IOException;

import de.mfabricius.diplom.data.Vector;
import de.mfabricius.diplom.semilagrange.Problem;
import de.mfabricius.diplom.util.BetaFactory;
import de.mfabricius.diplom.util.PhiFactory;
import de.mfabricius.diplom.util.R4toR3;
import de.mfabricius.diplom.util.Triangulation;


/**
 * This is the main entry class for this project. It holds several preset Problem configurations with PHIs, betas
 * initial value generator and in some cases the exact solution. To write error norms to a file instead of the
 * {@code System.out}, overwrite the static {@code Triangulation.out}.
 * 
 * @author Martin Fabricius
 * @version 1.0.0 2017-05-22
 *
 */
public class Main {

    /**
     * A sphere rotating around the z-axis and the tangential portion of a constant vector field.
     * 
     * @param h
     *            the approximation parameter
     * @return the Problem
     */
    static Problem rotatingSphere_constantField(double h) {
        return new Problem(SPHERE, PhiFactory.RotationZ, BetaFactory.ConstantForSphere(new Vector(2., 1., -.5)),
                circleAroundY, h, null);
    }

    /**
     * A torus rotating around the z-axis and the tangential field which goes around the small circle.
     * 
     * @param h
     *            the approximation parameter
     * @return the Problem
     */
    static Problem rotatingTorus_upwind(double h) {
        return new Problem(TORUS, PhiFactory.RotationZ, BetaFactory.UpwindTorus(1), circleAroundY, h, null);
    }

    /**
     * A sphere rotating around the z-axis and the transportation is reverse the rotation.
     * 
     * @param step
     *            the number of bisection steps
     * @return the Problem
     */
    static Problem rotatingSphere_antiZrotatingField(int step) {
        return new Problem(SPHERE, PhiFactory.RotationZ, BetaFactory.RotationZForSphere(-1), exponentialDistanceFromY,
                step, exponentialDistanceFromY);
    }

    /**
     * A sphere rotating around the z-axis and the transportation is the same as the rotation velocity.
     * 
     * @param step
     *            the number of bisection steps
     * @return the Problem
     */
    static Problem rotatingSphere_doubleZrotatingField(int step) {
        final R4toR3 revert = PhiFactory.RotationZ(-2).move;
        return new Problem(SPHERE, PhiFactory.RotationZ, BetaFactory.RotationZForSphere(1), exponentialDistanceFromY,
                step, (t, v) -> exponentialDistanceFromY.apply(t, revert.apply(t, v)));

    }

    /**
     * A sphere rotating around the given axis and the transportation is the reverse as the rotation velocity.
     * 
     * @param step
     *            the number of bisection steps
     * @param axis
     *            the axis of rotation (must be a unit Vector)
     * @return the Problem
     */
    static Problem rotatingSphere_antirotatingField(int step, Vector axis) {
        return new Problem(SPHERE, PhiFactory.Rotation(axis, 1), BetaFactory.RotationForSphere(axis, -1),
                exponentialDistanceFromY, step, exponentialDistanceFromY);
    }

    /**
     * A sphere rotating around the given axis and the transportation is the same as the rotation velocity.
     * 
     * @param step
     *            the number of bisection steps
     * @param axis
     *            the axis of rotation (must be a unit Vector)
     * @return the Problem
     */
    static Problem rotatingSphere_doublerotatingField(int step, Vector axis) {
        final R4toR3 revert = PhiFactory.Rotation(axis, -2).move;
        return new Problem(SPHERE, PhiFactory.Rotation(axis, 1.), BetaFactory.RotationForSphere(axis, 1.),
                exponentialDistanceFromY, step, (t, v) -> exponentialDistanceFromY.apply(t, revert.apply(t, v)));
    }

    public static void main(String[] args) throws IOException {
        Triangulation.out = new FileWriter("L2errors2.txt");
        for (int i = 0; i < 10; i++) {
            final Problem p = rotatingSphere_doubleZrotatingField(i);
            p.solve(0.125, 5, "solution" + i);
            Triangulation.out.write("\r\n\r\n");
        }
    }
}