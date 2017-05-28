package de.mfabricius.diplom.semilagrange;

import static de.mfabricius.diplom.data.Vector.skprod;
import static de.mfabricius.diplom.semilagrange.Problem.EPS;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import Jama.Matrix;
import de.mfabricius.diplom.data.Barycentric;
import de.mfabricius.diplom.data.Triangle;
import de.mfabricius.diplom.data.Vector;
import de.mfabricius.diplom.semilagrange.Problem.SEARCH;
import de.mfabricius.diplom.util.BiMap;
import de.mfabricius.diplom.util.IO;
import de.mfabricius.diplom.util.Triangulation;


/**
 * This class does the actual Semi-Lagrange-Scheme. It sets up the characteristical system, solves it, projects the foot
 * on the old {@link Triangulation} and interpolates the function value.
 * 
 * @author Martin Fabricius
 */
public class SLV {

    /**
     * This method solves a step of the Semi-Lagrange-Scheme for the whole {@link Triangualtion}
     * 
     * @param problem
     *            the {@link Problem} for the movement and the transportation field
     * @param timestep
     *            the timestep
     * @param time
     *            the current time
     * @param oldT3
     *            the Triangulation at t-1
     * @param newT3
     *            the Triangulation at t
     * @param step
     *            {@link BiMap} of the vertex movement before->after
     * @param u
     *            the values of the function u at the vertices in t-1
     * @return the values of the function u at the vertices in t
     */
    public static Map<Vector, Double> semilagrange(Problem problem, double timestep, double time, Triangulation oldT3,
            Triangulation newT3, BiMap<Vector, Vector> step, Map<Vector, Double> u) {

        final AtomicInteger points = new AtomicInteger();
        final AtomicInteger vertex = new AtomicInteger();
        final ConcurrentHashMap<Vector, Double> u2 = new ConcurrentHashMap<>(u.size());
        final long benchmark = -System.nanoTime();
        newT3.getVertices().parallelStream().forEach(v -> {
            points.incrementAndGet();
            // calculate vector field at v for characteristic system
            final Vector dir = v.sub(step.getInverse(v)).mult(1. / timestep).add(problem.beta.apply(time, v));
            // final Vector dir = problem.PHI.velocity.apply(time, v).add(problem.beta.apply(time, v));

            // solve characteristic system for v
            final Vector foot = charSystem(timestep, v, dir);

            // project foot onto oldT3
            final Barycentric bary = projectOnT3(foot, oldT3, step.getInverse(v), problem.search);

            if (bary.isVertex) vertex.incrementAndGet();

            // interpolate function value
            final double value = bary.interpolate(u);

            u2.put(v, value);
        });
        System.out.println(points + " Points\t" + vertex + " on Vertex\t"
                + IO.format((benchmark + System.nanoTime()) / 1e9) + " s");
        return u2;
    }

    /**
     * This method solves the charactristic system with the Euler method.
     * 
     * @param timestep
     *            the timestep
     * @param v
     *            the vertex where to set up
     * @param direction
     *            the direction of the characteristic
     * @return {@link Vector} foot of the characteristic one {@code timestep} before
     */
    private static Vector charSystem(double timestep, Vector v, Vector direction) {
        return v.sub(direction.mult(timestep));
    }

    /**
     * This method projects the given {@link Vector} othogonally onto the {@link Triangulation}.
     * 
     * @param p
     *            the {@link Vector} to project
     * @param t3
     *            the {@link Triangulation} to project onto
     * @param old
     *            the {@link Vector} of the foot of the characteristic to start the search from
     * @param search
     *            the search algorithm
     * @return the projection on the {@link Triangulation} as {@link Barycentric} tupel
     * @see SEARCH
     */
    private static Barycentric projectOnT3(Vector p, Triangulation t3, Vector old, SEARCH search) {

        switch (search) {
        case NearestNeighbours:
        case NearestVertex:
        case NearestVertexPath:
            // calculate nearest vertex
            final Vector nearest;
            if (search == SEARCH.NearestVertex) nearest = getNearestVertex(p, t3);
            else nearest = getNearestVertexPath(p, t3, old);

            if (search == SEARCH.NearestNeighbours) {
                // check adjacent triangles
                Barycentric projection = projectOnNearTriangle(p, t3, nearest);
                if (projection != null) return projection;
                projection = projectOnNearEdge(p, t3, nearest);
                if (projection != null) return projection;
            }
            // return nearest vertex
            return new Barycentric(t3.vNeighbors.get(nearest).peek(), nearest, true);
        default:
            return null;
        }

    }

    // private static Barycentric barycentricProjection(Vector p, Triangle t, Vector n) {
    // final Matrix solve = solveProjection(t, p, n);
    // final Vector proj = p.add(n.mult(solve.get(2, 0)));
    // return new Barycentric(t, proj);
    // }

    /**
     * Calculate the orthogonal projection of the given {@link Vector} onto each adjacent {@link Triangle} at the
     * {@link Vector} {@code nearest} in the given {@link Triangulation}. If any inner Points are met, the closest will
     * be taken.
     * 
     * @param v
     *            the {@link Vector} to project
     * @param t3
     *            the {@link Triangulation} to get the the {@link Triangle}s
     * @param nearest
     *            the nearest vertex of the {@link Triangulation}
     * @return the nearest inner projection, {@code null} if all failed.
     */
    private static Barycentric projectOnNearTriangle(Vector v, Triangulation t3, Vector nearest) {
        Barycentric projection = null;
        double mindist = Double.MAX_VALUE;

        for (final Triangle t : t3.vNeighbors.get(nearest)) {
            final Matrix solve = solveProjection(t, v);
            final double r = solve.get(0, 0);
            final double s = solve.get(1, 0);
            if (r >= -EPS && s >= -EPS && r + s <= 1 + EPS) {
                final Vector proj = v.add(t.normal().mult(solve.get(2, 0)));
                final double distance = v.distance(proj);
                if (distance < mindist) {
                    mindist = distance;
                    projection = new Barycentric(t, proj);
                }
            }
        }
        return projection;
    }

    /**
     * Calculate the orthogonal projection of the given {@link Vector} onto each adjacent {@link Edge} at the
     * {@link Vector} {@code nearest} in the given {@link Triangulation}. If any inner Points are met, the closest will
     * be taken.
     * 
     * @param v
     *            the {@link Vector} to project
     * @param t3
     *            the {@link Triangulation} to get the the {@link Edge}s
     * @param nearest
     *            the nearest vertex of the {@link Triangulation}
     * @return the nearest inner projection, {@code null} if all failed.
     */
    private static Barycentric projectOnNearEdge(Vector v, Triangulation t3, Vector nearest) {
        Barycentric projection = null;
        double mindist = Double.MAX_VALUE;
        final Set<Vector> done = new HashSet<>(8);
        done.add(nearest);

        for (final Triangle t : t3.vNeighbors.get(nearest)) {
            for (final Vector v2 : t.vertices()) {
                if (done.add(v2)) {
                    final Vector v1 = nearest.sub(v2);
                    final double v2n2 = v1.norm22();
                    final Vector proj = v2.add(v1.mult(skprod(v.sub(v2), v1) / v2n2));
                    if (proj.distance(nearest) + proj.distance(v2) < Math.sqrt(v2n2) * (1 + EPS)) {
                        final double distance = v.distance(proj);
                        if (distance < mindist) {
                            mindist = distance;
                            projection = new Barycentric(t, proj);
                        }
                    }
                }
            }
        }
        return projection;
    }

    /**
     * Look for the nearest vertex of the {@link Triangulation} for the given {@link Vector}. Calculates the distances
     * to all vertices in a parallel stream. May take very long.
     * 
     * @param v
     *            the {@link Vector} to look for the nearest vertex
     * @param t3
     *            the {@link Triangulation}
     * @return the nearest vertex to v
     */
    private static Vector getNearestVertex(Vector v, Triangulation t3) {
        final Nearest nearest = new Nearest();
        t3.getVertices().parallelStream().forEach(o -> {
            final double dist = v.distance(o);
            if (dist < nearest.min) {
                synchronized (nearest) {
                    if (dist < nearest.min) {
                        nearest.min = dist;
                        nearest.v = o;
                    }
                }
            }
        });
        return nearest.v;
    }

    /**
     * Wrapper class for calculating the nearest vertex in parallel forEach.
     * 
     * @author Martin Fabricius
     */
    private static class Nearest {

        double min = Double.MAX_VALUE;
        Vector v = new Vector(Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE);
    }

    /**
     * Look for the nearest vertex on the {@link Triangulation} for the given {@link Vector} following a path of the
     * nearest neighboring vertex, starting at the given {@code start}. {@code start} has to be a vertex of the
     * {@link Triangulation}.
     * 
     * @param v
     *            the {@link Vector} to search the nearest for
     * @param t3
     *            the {@link Triangulation} to search on
     * @param start
     *            the vertex to start the searching path
     * @return the nearest vertex
     * @throws NullPointerException
     *             if {@code start} is not a vertex of the {@link Triangulation}
     */
    private static Vector getNearestVertexPath(Vector v, Triangulation t3, Vector start) {
        double min = v.distance(start);
        Vector newNearest = start;
        do {
            start = newNearest;
            for (final Vector o : t3.getNeighborVertices(start)) {
                final double dist = v.distance(o);
                if (dist < min) {
                    min = dist;
                    newNearest = o;
                }
            }
        } while (start != newNearest);
        return start;
    }

    /**
     * Solves the linear equation system for the intersection of the given {@link Triangle} and the normal going through
     * the given {@link Vector}.
     * 
     * @param t
     *            the {@link Triangle}
     * @param v
     *            the {@link Vector}
     * @return the parameter of the solution
     */
    private static Matrix solveProjection(Triangle t, Vector v) {
        final Vector v1 = t.v1();
        final Vector v2 = t.v2();
        final Vector v3 = t.v3();
        final Vector n = t.normal();
        final Matrix A = new Matrix(3, 3);
        A.set(0, 0, v2.x() - v1.x());
        A.set(1, 0, v2.y() - v1.y());
        A.set(2, 0, v2.z() - v1.z());
        A.set(0, 1, v3.x() - v1.x());
        A.set(1, 1, v3.y() - v1.y());
        A.set(2, 1, v3.z() - v1.z());
        A.set(0, 2, -n.x());
        A.set(1, 2, -n.y());
        A.set(2, 2, -n.z());
        // final QRDecomposition qr = A.qr();
        final Matrix b = new Matrix(3, 1);
        b.set(0, 0, v.x() - v1.x());
        b.set(1, 0, v.y() - v1.y());
        b.set(2, 0, v.z() - v1.z());
        return A.solve(b);
    }
}