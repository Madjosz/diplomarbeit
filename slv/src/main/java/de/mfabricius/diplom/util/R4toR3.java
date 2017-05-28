package de.mfabricius.diplom.util;

import java.util.function.BiFunction;

import de.mfabricius.diplom.data.Vector;


/**
 * This is a functional interface for functions R x R³ -> R³.
 * 
 * @author Martin Fabricius
 */
public interface R4toR3 extends BiFunction<Double, Vector, Vector> {

    public static final R4toR3 ZERO = (t, x) -> Vector.ZERO;
}