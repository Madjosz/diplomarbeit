package de.mfabricius.diplom.data;

/**
 * This class represents an edge of a triangle in the 3-dimensional euclidean space R³ defined by its 2 vertices.
 * 
 * @author Martin Fabricius
 */
public class Edge {

    private final Vector v1;
    private final Vector v2;
    private final Vector mid;

    /**
     * Creates a new {@link Edge} with the given vertices as endpoints.
     * 
     * @param v1
     *            vertex 1
     * @param v2
     *            vertex 2
     */
    public Edge(Vector v1, Vector v2) {
        this.v1 = v1;
        this.v2 = v2;
        mid = v1.add(v2).mult(.5);
    }

    /** @return vertex 1 */
    public Vector v1() {
        return v1;
    }

    /** @return vertex 2 */
    public Vector v2() {
        return v2;
    }

    /** @returns the mid point of this edge */
    public Vector mid() {
        return mid;
    }

    /*
     * object methods
     */

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj instanceof Edge) {
            final Edge e = (Edge) obj;
            return (v1.equals(e.v1) && v2.equals(e.v2)) || (v1.equals(e.v2) && v2.equals(e.v1));
        }
        return false;
    }

    @Override
    public int hashCode() {
        return v1.hashCode() + v2.hashCode();
    }

    @Override
    public String toString() {
        return "(" + v1.toString() + "," + v2.toString() + ")";
    }
}
