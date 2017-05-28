package de.mfabricius.diplom.semilagrange;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import de.mfabricius.diplom.data.Vector;
import de.mfabricius.diplom.util.BiMap;
import de.mfabricius.diplom.util.IO;
import de.mfabricius.diplom.util.PhiFactory.PHI;
import de.mfabricius.diplom.util.R4toR;
import de.mfabricius.diplom.util.R4toR3;
import de.mfabricius.diplom.util.Triangulation;
import de.mfabricius.diplom.util.Triangulation.TriangulationUtils;


/**
 * This class is a model for the partial differential equation to solve. The equation is given in LaTeX
 * <pre>\dot u + \nabla_\Gamma \cdot \beta = 0</pre>
 * 
 * @author Martin Fabricius
 */
public class Problem {

    /** Search algorithms for the approximative foot point of the characteristic. */
    public static enum SEARCH {
        /** Compare the point with all vertices and take the closest. */
        NearestVertex,
        /** Search from the old vertex along a path of the closest neighbors. */
        NearestVertexPath,
        /**
         * Like {@link SEARCH#NearestVertexPath} and than project on the adjacent Triangles. Take the closest
         * projection.
         */
        NearestNeighbours
    }

    /** The movement of the surface given by a {@link PHI}. */
    public final PHI PHI;
    /** The transportation field tangent to the surface. */
    public final R4toR3 beta;
    /** The function to compute the initial values on Gamma_0. */
    public final R4toR initialValues;
    /** The exact solution to compare the result to (May be {@code null}). */
    public final R4toR exactSolution;
    /** The inital {@link Triangulation} Gamma_0. */
    public final Triangulation initialT3;
    /** A coarse {@link Triangulation} for visualising the transportation field beta. */
    public final Triangulation coarseT3;
    /** The search algorithm to use for the Semi-Lagrange-Scheme. */
    public final SEARCH search;
    /** The machine accuracy. */
    public static final double EPS = 1e-12;

    /**
     * Creates a new {@link Problem} for a given movement {@link PHI}, transportation fiel {@code beta} on the given
     * {@code surface} with the accuracy {@code h}.
     * 
     * @param surface
     *            the {@link TriangulationUtils} for creating the {@code Triangulation}
     * @param PHI
     *            the movement of the surface
     * @param beta
     *            the transportation field tangent to the surface
     * @param initialValues
     *            the function for setting up the initial values on Gamma_0
     * @param h
     *            the accuracy parameter
     * @param exactSolution
     *            the exact solution to compute the L2-Norm of the numerical error, may be {@code null}
     */
    public Problem(TriangulationUtils surface, PHI PHI, R4toR3 beta, R4toR initialValues, double h,
            R4toR exactSolution) {
        this.PHI = PHI;
        this.beta = beta;
        this.initialValues = initialValues;
        this.exactSolution = exactSolution;
        initialT3 = surface.triangulate.apply(h);
        final Triangulation coarse = surface.t3steps.apply(1);
        coarseT3 = coarse.getVertices().size() > initialT3.getVertices().size() ? initialT3 : coarse;
        search = SEARCH.NearestNeighbours;
    }

    /**
     * Creates a new {@link Problem} for a given movement {@link PHI}, transportation fiel {@code beta} on the given
     * {@code surface} with the accuracy {@code h}.
     * 
     * @param surface
     *            the {@link TriangulationUtils} for creating the {@code Triangulation}
     * @param PHI
     *            the movement of the surface
     * @param beta
     *            the transportation field tangent to the surface
     * @param initialValues
     *            the function for setting up the initial values on Gamma_0
     * @param step
     *            the number of bisection steps
     * @param exactSolution
     *            the exact solution to compute the L2-Norm of the numerical error, may be {@code null}
     */
    public Problem(TriangulationUtils surface, PHI PHI, R4toR3 beta, R4toR initialValues, int step,
            R4toR exactSolution) {
        this.PHI = PHI;
        this.beta = beta;
        this.initialValues = initialValues;
        this.exactSolution = exactSolution;
        coarseT3 = surface.t3steps.apply(step < 1 ? step : 1);
        initialT3 = surface.t3steps.apply(step);
        search = SEARCH.NearestNeighbours;
    }

    /**
     * Solve the {@link Problem} numerically according to the algorithm given in section 7.1 of the diploma thesis.
     * 
     * @param timestep
     *            the discretisation of the time interval
     * @param maxT
     *            the end of the time interval
     * @param fileName
     *            filename prefix for printing the visualisation, {@code null} if no output needed
     */
    public void solve(double timestep, double maxT, String fileName) {
        double currentTime = timestep;

        // create Gamma_0
        final Set<Vector> vertices = initialT3.getVertices();
        final Map<Vector, Vector> move = new ConcurrentHashMap<>(vertices.size());
        vertices.parallelStream().forEach(v -> move.put(v, PHI.move.apply(0., v)));
        Triangulation oldT3 = Triangulation.move(initialT3, move);

        // set up initial values
        Map<Vector, Double> u0 = new HashMap<>(vertices.size());
        for (final Vector v : oldT3.getVertices())
            u0.put(v, initialValues.apply(0.0, v));

        // print initial values
        oldT3.printProperties(0.);
        if (fileName != null) oldT3.print(u0, beta, coarseT3.move(move), fileName + "_0.000");

        while (currentTime <= maxT) {
            final double t = currentTime;

            // create Gamma_t
            final BiMap<Vector, Vector> step = new BiMap<>();
            vertices.forEach(p0 -> {
                final Vector p2 = PHI.move.apply(t, p0);
                step.put(move.put(p0, p2), p2);
            });
            final Triangulation newT3 = Triangulation.move(initialT3, move);

            final double h = newT3.printProperties(t);

            // solve equation with Semi Lagrange Method
            final Map<Vector, Double> u = SLV.semilagrange(this, timestep, t, oldT3, newT3, step, u0);

            // print the numerical error in the L²-Norm
            if (exactSolution != null) newT3.printNumericalErrors(h, t, u, exactSolution);

            // print solution
            if (fileName != null) newT3.print(u, fileName + "_" + IO.ff.format(t));

            // override old values
            oldT3 = newT3;
            u0 = u;
            currentTime += timestep;
        }
    }
}