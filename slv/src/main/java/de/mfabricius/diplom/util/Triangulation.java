package de.mfabricius.diplom.util;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.DoubleFunction;
import java.util.function.IntFunction;
import java.util.function.Supplier;
import java.util.function.UnaryOperator;

import de.mfabricius.diplom.data.Edge;
import de.mfabricius.diplom.data.Triangle;
import de.mfabricius.diplom.data.Vector;


/**
 * This class represents a {@code Triangulation} of 2-dimensional surfaces in the euclidean space R³. The Triangulation
 * consists of a {@link List} of {@link Triangle}s and for each vertex of the {@code Triangualtion} the adjacent
 * {@link Triangle}s are saved in a {@link Map}.
 * 
 * Static methods for creating a {@code Triangulation} for a plane, a sphere, a torus and a cylinder are provided. These
 * are accessible through the {@link TriangulationUtils} and the constants {@link #PLANE}, {@link #SPHERE},
 * {@link #TORUS} and {@link #CYLINDER}.
 * 
 * Also there are methods for moving, printing and calaculating several properties of a {@code Triangulation}.
 * 
 * @author Martin Fabricius
 */
public class Triangulation {

    /**
     * The writer for writing the numerical errors to (@link printNumericalErrors()}). Overwrite this variable to output
     * to another {@link Writer}.
     */
    public static Writer out = new PrintWriter(System.out);

    /** The {@link List} of {@link Triangle}s of this {@link Triangulation}. */
    public final List<Triangle> triangles;
    /** The Map of adjacent {@link Triangle}s for each vertex of this {@link Triangulation}. */
    public final Map<Vector, Queue<Triangle>> vNeighbors;

    /**
     * Construct a {@link Triangulation} from the given {@link List} of {@link Triangle}s. This constructor is
     * {@code private} to avoid setting up irregular Triangulations.
     * 
     * @param triangles
     *            the List of Triangles
     */
    private Triangulation(List<Triangle> triangles) {
        this.triangles = triangles;
        vNeighbors = new ConcurrentHashMap<>();
        triangles.parallelStream().forEach(t -> {
            for (final Vector v : t.vertices())
                vNeighbors.computeIfAbsent(v, key -> new ConcurrentLinkedQueue<>()).add(t);
        });
    }

    /** @return the {@link Set} of {@link Vector}s which are the vertices of the {@link Triangulation}. */
    public Set<Vector> getVertices() {
        return Collections.unmodifiableSet(vNeighbors.keySet());
    }

    /**
     * Returns a new {@link Set} of {@link Vector}s, who are vetices of adjacent {@link Triangle}s to the given vertex.
     * 
     * @param v
     *            the vertex to get the neighbors for
     * @return Set of neighboring vetices
     */
    public Set<Vector> getNeighborVertices(Vector v) {
        final Set<Vector> adjacent = new HashSet<>();
        for (final Triangle t : vNeighbors.get(v))
            adjacent.addAll(Arrays.asList(t.vertices()));
        adjacent.remove(v);
        return adjacent;
    }

    /* SPHERE */

    /** Creates a initial {@link Triangulation} of the unit sphere sphere as an octahedron. */
    private static Triangulation initialSphere() {
        final Vector[] points = { new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(-1, 0, 0), new Vector(0, -1, 0),
                new Vector(0, 0, 1), new Vector(0, 0, -1) };
        final ArrayList<Triangle> triangles = new ArrayList<>(8);
        triangles.add(new Triangle(points[0], points[1], points[4]));
        triangles.add(new Triangle(points[1], points[2], points[4]));
        triangles.add(new Triangle(points[2], points[3], points[4]));
        triangles.add(new Triangle(points[3], points[0], points[4]));
        triangles.add(new Triangle(points[3], points[2], points[5]));
        triangles.add(new Triangle(points[2], points[1], points[5]));
        triangles.add(new Triangle(points[1], points[0], points[5]));
        triangles.add(new Triangle(points[0], points[3], points[5]));
        return new Triangulation(triangles);
    }

    /** The projector to project {@link Vector}s on a unit sphere. */
    private static final UnaryOperator<Vector> PROJECTOR_SPHERE = Vector::unitVector;

    /** The {@link TriangulationUtils} for setting up a {@link Triangulation} of the unit sphere. */
    public static TriangulationUtils SPHERE = new TriangulationUtils("sphere", Triangulation::initialSphere,
            PROJECTOR_SPHERE);

    /* PLANE */

    /** Creates a initial {@link Triangulation} of the unit square. */
    private static Triangulation initialPlane() {
        final Vector[] points = { new Vector(0, 0, 0), new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(1, 1, 0) };
        final ArrayList<Triangle> triangles = new ArrayList<>(2);
        triangles.add(new Triangle(points[0], points[1], points[2]));
        triangles.add(new Triangle(points[1], points[3], points[2]));
        return new Triangulation(triangles);
    }

    /** The {@link TriangulationUtils} for setting up a {@link Triangulation} of the unit square. */
    public static TriangulationUtils PLANE = new TriangulationUtils("plane", Triangulation::initialPlane,
            UnaryOperator.<Vector>identity());

    /* CYLINDER */

    /** Creates a initial {@link Triangulation} of a cylinder with height 1.2 and radius 1. */
    private static Triangulation initialCylinder() {
        final double sqrt2 = 1. / Math.sqrt(2);
        final Vector[] points = { new Vector(1, 0, 0.6), new Vector(0, 1, 0.6), new Vector(-1, 0, 0.6),
                new Vector(0, -1, 0.6), new Vector(sqrt2, sqrt2, -0.6), new Vector(-sqrt2, sqrt2, -0.6),
                new Vector(-sqrt2, -sqrt2, -0.6), new Vector(sqrt2, -sqrt2, -0.6) };
        final ArrayList<Triangle> triangles = new ArrayList<>(8);
        triangles.add(new Triangle(points[0], points[4], points[1]));
        triangles.add(new Triangle(points[1], points[5], points[2]));
        triangles.add(new Triangle(points[2], points[6], points[3]));
        triangles.add(new Triangle(points[3], points[7], points[0]));
        triangles.add(new Triangle(points[4], points[5], points[1]));
        triangles.add(new Triangle(points[5], points[6], points[2]));
        triangles.add(new Triangle(points[6], points[7], points[3]));
        triangles.add(new Triangle(points[7], points[4], points[0]));
        return new Triangulation(triangles);
    }

    /** The projector to project {@link Vector}s on the cylinder. */
    private static final UnaryOperator<Vector> PROJECTOR_CYLINDER = v -> {
        final double f = 1. / Math.sqrt(v.x() * v.x() + v.y() * v.y());
        return new Vector(f * v.x(), f * v.y(), v.z());
    };

    /**
     * The {@link TriangulationUtils} for setting up a {@link Triangulation} of a cylinder with radius 1 and height 1.2.
     */
    public static TriangulationUtils CYLINDER = new TriangulationUtils("cylinder", Triangulation::initialCylinder,
            PROJECTOR_CYLINDER);

    /* TORUS */

    /** Creates a initial {@link Triangulation} of a torus with great radius 1 and little radius 0.5. */
    private static Triangulation initialTorus() {
        final double sqrt2 = 1. / Math.sqrt(2);
        final double sqrt15 = 3. * Math.sqrt(3) / 4.;
        final Vector[] points = { new Vector(.5, 0, 0), new Vector(0, .5, 0), new Vector(-.5, 0, 0),
                new Vector(0, -.5, 0), new Vector(1, 0, .5), new Vector(sqrt2, sqrt2, .5), new Vector(0, 1, .5),
                new Vector(-sqrt2, sqrt2, .5), new Vector(-1, 0, .5), new Vector(-sqrt2, -sqrt2, .5),
                new Vector(0, -1, .5), new Vector(sqrt2, -sqrt2, .5), new Vector(1, 0, -.5),
                new Vector(sqrt2, sqrt2, -.5), new Vector(0, 1, -.5), new Vector(-sqrt2, sqrt2, -.5),
                new Vector(-1, 0, -.5), new Vector(-sqrt2, -sqrt2, -.5), new Vector(0, -1, -.5),
                new Vector(sqrt2, -sqrt2, -.5), new Vector(1.5, 0, 0), new Vector(sqrt15, .75, 0),
                new Vector(.75, sqrt15, 0), new Vector(0, 1.5, 0), new Vector(-.75, sqrt15, 0),
                new Vector(-sqrt15, .75, 0), new Vector(-1.5, 0, 0), new Vector(-sqrt15, -.75, 0),
                new Vector(-.75, -sqrt15, 0), new Vector(0, -1.5, 0), new Vector(.75, -sqrt15, 0),
                new Vector(sqrt15, -.75, 0) };
        final ArrayList<Triangle> triangles = new ArrayList<>(64);
        for (int i = 0; i < 4; i++) {
            triangles.add(new Triangle(points[0 + i], points[4 + 2 * i], points[5 + 2 * i]));
            triangles.add(new Triangle(points[0 + i], points[5 + 2 * i], points[(1 + i) % 4]));
            triangles.add(new Triangle(points[(1 + i) % 4], points[5 + 2 * i], points[i != 3 ? 6 + 2 * i : 4]));
            triangles.add(new Triangle(points[0 + i], points[13 + 2 * i], points[12 + 2 * i]));
            triangles.add(new Triangle(points[0 + i], points[(1 + i) % 4], points[13 + 2 * i]));
            triangles.add(new Triangle(points[(1 + i) % 4], points[i != 3 ? 14 + 2 * i : 12], points[13 + 2 * i]));
            triangles.add(new Triangle(points[4 + 2 * i], points[20 + 3 * i], points[21 + 3 * i]));
            triangles.add(new Triangle(points[4 + 2 * i], points[21 + 3 * i], points[5 + 2 * i]));
            triangles.add(new Triangle(points[5 + 2 * i], points[21 + 3 * i], points[22 + 3 * i]));
            triangles.add(new Triangle(points[5 + 2 * i], points[22 + 3 * i], points[i != 3 ? 6 + 2 * i : 4]));
            triangles.add(
                    new Triangle(points[i != 3 ? 6 + 2 * i : 4], points[22 + 3 * i], points[i != 3 ? 23 + 3 * i : 20]));
            triangles.add(new Triangle(points[12 + 2 * i], points[21 + 3 * i], points[20 + 3 * i]));
            triangles.add(new Triangle(points[12 + 2 * i], points[13 + 2 * i], points[21 + 3 * i]));
            triangles.add(new Triangle(points[13 + 2 * i], points[22 + 3 * i], points[21 + 3 * i]));
            triangles.add(new Triangle(points[13 + 2 * i], points[i != 3 ? 14 + 2 * i : 12], points[22 + 3 * i]));
            triangles.add(new Triangle(points[i != 3 ? 14 + 2 * i : 12], points[i != 3 ? 23 + 3 * i : 20],
                    points[22 + 3 * i]));
        }
        return new Triangulation(triangles);
    }

    /** The projector to project {@link Vector}s on the torus. */
    private static final UnaryOperator<Vector> PROJECTOR_TORUS = v -> {
        final Vector c = new Vector(v.x(), v.y(), 0).unitVector();
        final Vector r = v.sub(c).unitVector().mult(0.5);
        return c.add(r);
    };

    /**
     * The {@link TriangulationUtils} for setting up a {@link Triangulation} of a torus with great radius 1 and little
     * radius 0.5.
     */
    public static TriangulationUtils TORUS = new TriangulationUtils("torus", Triangulation::initialTorus,
            PROJECTOR_TORUS);

    /* MANIPULATION METHODS */

    /**
     * Bisects every {@link Edge} of every {@link Triangle} in this {@link Triangulation}, projects the new
     * {@link Vector} on the surface and builds up 4 new {@link Triangle}s out of every {@link Triangle}.
     * 
     * @param triangulation
     *            the {@link Triangulation} to bisect
     * @param projector
     *            the projector to project the new {@link Vector}s onto the surface
     * @return the new {@link Triangulation}
     */
    private static Triangulation bisect(Triangulation triangulation, UnaryOperator<Vector> projector) {
        final List<Triangle> triangles = triangulation.triangles;
        final List<Triangle> newTriangles = Collections.synchronizedList(new ArrayList<>(4 * triangles.size()));
        final Map<Edge, Vector> midpoints = new ConcurrentHashMap<>();
        triangles.parallelStream().forEach(t -> {
            final Vector[] v = t.vertices();
            final Vector[] m = new Vector[3];
            for (int i = 0; i < 3; i++)
                m[i] = midpoints.computeIfAbsent(new Edge(v[i], v[(i + 1) % 3]), e -> projector.apply(e.mid()));
            for (int i = 0; i < 3; i++)
                newTriangles.add(new Triangle(v[i], m[i], m[(i + 2) % 3]));
            newTriangles.add(new Triangle(m[0], m[1], m[2]));
        });
        return new Triangulation(newTriangles);
    }

    /**
     * Moves the vertices of this {@link Triangulation} to the new location given by a {@link Map}.
     * 
     * @param move
     *            the {@link Map} of the old and new location
     * @return the moved {@link Triangulation}
     * @throws NullPointerException
     *             if a mapping is missing for any vertex
     */
    public Triangulation move(Map<Vector, Vector> move) {
        return move(this, move);
    }

    /**
     * Moves the vertices of this {@link Triangulation} to the new location given by a {@link UnaryOperator}, so it
     * calculates the new postions from the old ones.
     * 
     * @param move
     *            the {@link UnaryOperator} to calculate the new positions
     * @return the moved {@link Triangulation}
     */
    public Triangulation move(UnaryOperator<Vector> move) {
        return move(this, move);
    }

    /**
     * Moves the vertices of the given {@link Triangulation} to the new location given by a {@link Map}.
     * 
     * @param triangulation
     *            the {@link Triangulation} to move
     * @param move
     *            the {@link Map} of the old and new location
     * @return the moved {@link Triangulation}
     * @throws NullPointerException
     *             if a mapping is missing for any vertex
     */
    public static Triangulation move(Triangulation triangulation, Map<Vector, Vector> move) {
        final List<Triangle> triangles = triangulation.triangles;
        final List<Triangle> newTriangles = Collections.synchronizedList(new ArrayList<>(triangles.size()));
        triangulation.triangles.parallelStream().forEach(t -> {
            final Triangle newT = new Triangle(move.get(t.v1()), move.get(t.v2()), move.get(t.v3()));
            if (newT.area() > 0.0) newTriangles.add(newT);
        });
        return new Triangulation(newTriangles);
    }

    /**
     * Moves the vertices of the given {@link Triangulation} to the new location given by a {@link UnaryOperator}, so it
     * calculates the new postions from the old ones.
     * 
     * @param triangulation
     *            the {@link Triangulation} to move
     * @param move
     *            the {@link UnaryOperator} to calculate the new positions
     * @return the moved {@link Triangulation}
     */
    @SuppressWarnings("serial")
    public static Triangulation move(Triangulation triangulation, UnaryOperator<Vector> move) {
        return move(triangulation, new ConcurrentHashMap<Vector, Vector>() {

            @Override
            public Vector get(Object key) {
                return computeIfAbsent((Vector) key, move);
            };
        });
    }

    /**
     * Bisects the given initial {@link Triangulation} until the maximal diameter of the {@link Triangle}s is less than
     * or equal to {@code h}.
     * 
     * @param triangles
     *            the initial {@link Triangulation}
     * @param projector
     *            the projector to project the new {@link Vector}s onto the surface
     * @param h
     *            the maximal diameter of the returned {@link Triangulation}
     * @return {@link Triangulation} with max diameter {@code h}.
     */
    private static Triangulation triangulate(Triangulation triangles, UnaryOperator<Vector> projector, double h) {
        while (getMaxDiameter(triangles) > h)
            triangles = bisect(triangles, projector);

        return triangles;
    }

    /**
     * Bisects the given initial {@link Triangulation} the given {@code step} number of times.
     * 
     * @param triangles
     *            the initial {@link Triangulation}
     * @param projector
     *            the projector to project the new {@link Vector}s onto the surface
     * @param steps
     *            the number of bisections
     * @return {@link Triangulation} with max diameter {@code h}.
     */
    private static Triangulation triangulate(Triangulation triangles, UnaryOperator<Vector> projector, int steps) {
        while (steps-- > 0)
            triangles = bisect(triangles, projector);

        return triangles;
    }

    /* COMPUTING METHODS */

    /**
     * Calculates the maximum radius of the circumscribed circles of the {@link Triangle}s for the given
     * {@link Triangulation}.
     * 
     * @param triangulation
     *            the {@link Triangulation} to calculate
     * @return the maximum radius of the circumscribed circles
     */
    public static double getMaxUkr(Triangulation triangulation) {
        return triangulation.triangles.parallelStream().mapToDouble(Triangle::umkreis).max().getAsDouble();
    }

    /**
     * Calculates the maximum ratio of the radii of the circumscribed circle and the incircle of the {@link Triangle}s
     * for the given {@link Triangulation}.
     * 
     * @param triangulation
     *            the {@link Triangulation} to calculate
     * @return the maximum ratio of the radii of the circumscribed circle and the incircle
     */
    public static double getMaxUkikRatio(Triangulation triangulation) {
        return triangulation.triangles.parallelStream().mapToDouble(Triangle::ukikRatio).max().getAsDouble();
    }

    /**
     * Calculates the maximum diameter of the {@link Triangle}s for the given {@link Triangulation}.
     * 
     * @param triangulation
     *            the {@link Triangulation} to calculate
     * @return the maximum diameter
     */
    public static double getMaxDiameter(Triangulation triangulation) {
        return triangulation.triangles.parallelStream().mapToDouble(Triangle::diameter).max().getAsDouble();
    }

    /**
     * Calculates the maximum ratio of the diameter and the radius of the incircle of the {@link Triangle}s for the
     * given {@link Triangulation}.
     * 
     * @param triangulation
     *            the {@link Triangulation} to calculate
     * @return the maximum ratio of the diameter and the radius of the incircle
     */
    public static double getMaxDiamikRatio(Triangulation triangulation) {
        return triangulation.triangles.parallelStream().mapToDouble(Triangle::diamikRatio).max().getAsDouble();
    }

    /**
     * Calculates the L²-Norm of the linear interpolant with nodal values at the vertices of the {@link Triangle}s.
     * 
     * @param triangulation
     *            the {@link Triangulation} to calculate
     * @param u
     *            the mapping the vertices to the values
     * @return the L²-Norm
     */
    public static double getL2Norm(Triangulation triangulation, Map<Vector, Double> u) {
        return Math.sqrt(triangulation.triangles.parallelStream()
                .mapToDouble(t -> t.normL22(u.get(t.v1()), u.get(t.v2()), u.get(t.v3()))).sum());
    }

    /* OUTPUT METHODS */

    /**
     * Prints the wireframe of this {@link Triangulation} to a file with the name {@code [fileName].py}.
     * 
     * @param fileName
     *            the name of the file
     */
    public void print(String fileName) {
        IO.print(this, fileName);
    }

    /**
     * Prints the 3d model of this {@link Triangulation} with the colorcoded function values to a file with the name
     * {@code [fileName].py}.
     * 
     * @param u
     *            the values at the vertices of the {@link Triangulation}
     * @param fileName
     *            the name of the file
     */
    public void print(Map<Vector, Double> u, String fileName) {
        IO.print(this, u, fileName);
    }

    /**
     * Prints the 3d model of this {@link Triangulation} with the colorcoded function values and the given field as
     * arrows to a file with the name {@code [fileName].py}.
     * 
     * @param u
     *            the values at the vertices of the {@link Triangulation}
     * @param field
     *            the transportation field
     * @param fileName
     *            the name of the file
     */
    public void print(Map<Vector, Double> u, R4toR3 field, String fileName) {
        IO.print(this, u, field, fileName);
    }

    /**
     * Prints the 3d model of this {@link Triangulation} with the colorcoded function values and the given field as
     * arrows to a file with the name {@code [fileName].py}.
     * 
     * @param u
     *            the values at the vertices of the {@link Triangulation}
     * @param field
     *            the transportation field
     * @param fieldT3
     *            the coarse {@link Triangulation} for the field
     * @param fileName
     *            the name of the file
     */
    public void print(Map<Vector, Double> u, R4toR3 field, Triangulation fieldT3, String fileName) {
        IO.print(this, u, field, fieldT3, fileName);
    }

    /**
     * Print the properties (max diameter, max diameter/incircle) of this {@link Triangulation} prefixed with the given
     * time to the System.out.
     * 
     * @param time
     *            the current timestep to prefix the message
     * @return the max diameter
     */
    public double printProperties(double time) {
        final double maxDiameter = getMaxDiameter(this);
        System.out.println(IO.df.format(time) + "\tØ " + maxDiameter + "\tσ " + getMaxDiamikRatio(this));
        return maxDiameter;
    }

    /**
     * Print the properties (max diameter, max diameter/incircle, number of {@link Triangle}s, number of vertices) of
     * this {@link Triangulation} to the System.out.
     * 
     * @param step
     *            the current bisection step
     */
    public void printProperties(int step) {
        System.out.println(step + "\tØ " + getMaxDiameter(this) + "\tσ " + getMaxDiamikRatio(this) + "\tΔ "
                + triangles.size() + "\tV " + vNeighbors.size());
    }

    /**
     * If an exact solution is given, the error of the numerical solution with the linear interpolant of the numerical
     * solution in the L²-Norm will be printed to {@link #out}.
     * 
     * @param h
     *            the maximal diameter of this {@link Triangulation}
     * @param time
     *            the current timestep
     * @param u
     *            the numerical solution
     * @param exactSolution
     *            the exact solution
     */
    @SuppressWarnings("serial")
    public void printNumericalErrors(double h, double time, Map<Vector, Double> u, R4toR exactSolution) {
        if (exactSolution == null) return;
        final Map<Vector, Double> err = new ConcurrentHashMap<Vector, Double>() {

            @Override
            public Double get(Object key) {
                return computeIfAbsent((Vector) key, v -> u.get(v) - exactSolution.apply(time, v));
            }
        };
        try {
            // out.write(IO.ff.format(time) + "\tØ " + h + "\tL² " + getL2Norm(this, err) + "\r\n");
            out.write(IO.ff.format(time) + "\t" + getL2Norm(this, err) + "\r\n");
            out.flush();
        } catch (final IOException e) {}
    }

    /**
     * This class provides utilities to set up a macro {@link Triangulation} and refine it via bisection of the
     * {@link Edge}s.
     * 
     * @author Martin Fabricius
     */
    public static class TriangulationUtils {

        /**
         * Create new {@link TriangulationUtils}.
         * 
         * @param name
         *            the name of the surface (mainly for output)
         * @param initial
         *            {@link Supplier} for the initial {@link Triangulation}
         * @param projector
         *            for projecting the new vertices onto the surface
         */
        private TriangulationUtils(String name, Supplier<Triangulation> initial, UnaryOperator<Vector> projector) {
            this.name = name;
            this.initial = initial;
            this.projector = projector;
            bisect = t -> bisect(t, projector);
            triangulate = h -> triangulate(initial.get(), projector, h);
            t3steps = step -> triangulate(initial.get(), projector, step);
        }

        /** The name of the {@link Triangulation}. */
        public final String name;
        /** The {@link Supplier} for the initial {@link Triangulation}. */
        public final Supplier<Triangulation> initial;
        /** The projector for projecting the new vertices onto the surface. */
        public final UnaryOperator<Vector> projector;
        /**
         * The bisection operation for this type of surface. Do not apply a {@link Triangulation} of another type. The
         * resulting {@link Triangulation} is useless then.
         */
        public final UnaryOperator<Triangulation> bisect;
        /** Returns a {@link Triangulation} with given maximum diameter. */
        public final DoubleFunction<Triangulation> triangulate;
        /** Returns a {@link Triangulation} after the given number of bisection steps. */
        public final IntFunction<Triangulation> t3steps;
    }
}