package de.mfabricius.diplom.util;

import static de.mfabricius.diplom.util.IO.df;


/**
 * This is a utility class for color coding function values and output the color as VPython code.
 * 
 * @author Martin Fabricius
 */
public class Color {

    private final double min;
    private final double max;

    /** Creates a colorcoder with the range -1 (orange), 0 (green), 1 (blue). */
    public Color() {
        this(-1., 1.);
    }

    /** Creates a colorcoder with the range min (orange), avg(min,max) (green), max (blue). */
    public Color(double min, double max) {
        this.min = min;
        this.max = max;
    }

    /**
     * Creates a colorcoder with the range (orange) over (green) to (blue), where 0 is green if {@code force0} is
     * {@code true}.
     */
    public Color(double min, double max, boolean force0) {
        this(force0 ? -Math.max(Math.abs(min), Math.abs(max)) : min,
                force0 ? Math.max(Math.abs(min), Math.abs(max)) : max);
    }

    /**
     * Create the VPython color String for the given {@code value}.
     * 
     * @param value
     *            the value
     * @return VPython color String
     */
    public String color(double value) {
        double factor = (value - min) / (max - min);
        if (factor < 0) factor = 0;
        if (factor > 1) factor = 1;
        factor = 2. * factor - 1.;
        return factor < 0 ? (df.format(-factor) + "," + df.format(0.5 * factor + 1.) + ",0")
                : ("0," + df.format(1 - factor) + "," + df.format(factor));
    }

}
