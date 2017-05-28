package de.mfabricius.diplom.util;

import java.util.HashMap;
import java.util.Map;
import java.util.function.BiFunction;

import de.mfabricius.diplom.data.Vector;


/**
 * This is a functional interface for functions R x R³ -> R.
 * 
 * @author Martin Fabricius
 */
public interface R4toR extends BiFunction<Double, Vector, Double> {

    /** All Vectors in a cylinder of radius 0.5 have value 1, else 0 (not continuous). */
    public static R4toR circleAroundY = (t, v) -> Math.pow(v.x(), 2.) + Math.pow(v.z(), 2.) < 0.25 ? 1. : 0.;

    /** Test function on the cylinder of radius 0.8 around the y-axis with max=1. */
    public static R4toR exponentialDistanceFromY = (t, v) -> {
        final double dist = v.x() * v.x() + v.z() * v.z();
        return dist < 0.64 ? Math.E * Math.exp(0.64 / (dist - 0.64)) : 0.;
    };

    /** A Map that returns 0.0 for each input. Does not support Map interface properly. */
    @SuppressWarnings("serial")
    public static final Map<Vector, Double> ZERO_MAP = new HashMap<Vector, Double>() {

        @Override
        public Double get(Object key) {
            return 0.0;
        }
    };
}