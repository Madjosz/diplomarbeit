package de.mfabricius.diplom.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import de.mfabricius.diplom.data.Edge;
import de.mfabricius.diplom.data.Triangle;
import de.mfabricius.diplom.data.Vector;


/**
 * This class provides utility methods for printing {@link Triangulation}s to VPython code for visualisation.
 * 
 * @author Martin Fabricius
 */
public class IO {

    /** The {@link DecimalFormat} for {@code double}s. */
    public static final DecimalFormat df = new DecimalFormat("0.####", DecimalFormatSymbols.getInstance(Locale.US));

    /** The {@link DecimalFormat} for file names. */
    public static final DecimalFormat ff = new DecimalFormat("0.000", DecimalFormatSymbols.getInstance(Locale.US));

    /**
     * Convertes a {@link Vector} to a {@code String} for VPython.
     * 
     * @param v
     *            the Vector
     * @return VPython string
     */
    public static String print(Vector v) {
        return format(v.x()) + "," + format(v.y()) + "," + format(v.z());
    }

    /**
     * Formats a double with the format "0.####".
     * 
     * @param d
     *            the double
     * @return String representation
     */
    public static String format(double d) {
        return df.format(d);
    }

    /**
     * Prints VPython code visualising the wireframe of the given {@link Triangulation} to
     * "{@code output\[fileName].py}".
     * 
     * @param triangulation
     *            the Triangulation to print
     * @param fileName
     *            the name of the file without file extension
     */
    public static void print(Triangulation triangulation, String fileName) {
        ensureDirectory();
        try (FileWriter writer = new FileWriter("output\\" + fileName + ".py")) {
            writeHead(writer, fileName);
            final HashSet<Edge> edges = new HashSet<>();
            for (final Triangle t : triangulation.triangles) {
                for (final Edge e : Arrays.asList(new Edge(t.v1(), t.v2()), new Edge(t.v2(), t.v3()),
                        new Edge(t.v3(), t.v1())))
                    if (edges.add(e))
                        writer.write("curve(pos=[(" + print(e.v1()) + "),(" + print(e.v2()) + ")], color=(0,0,0))\r\n");
            }
        } catch (final IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Prints VPython code visualising the surface of the given {@link Triangulation} with color coded function values
     * to "{@code output\[fileName].py}".
     * 
     * @param t3
     *            the Triangulation
     * @param u
     *            the values of the function at the vertices
     * @param fileName
     *            the name of the output file
     * @see Color
     */
    public static void print(Triangulation t3, Map<Vector, Double> u, String fileName) {
        final Color c = new Color();
        ensureDirectory();
        try (FileWriter writer = new FileWriter("output\\" + fileName + ".py")) {
            writeHead(writer, fileName);
            writer.write("f = faces()\r\n");
            for (final Triangle t : t3.triangles) {
                writer.write("f.append(pos=(" + print(t.v1()) + "), normal=(0,0,0), color=(" + c.color(u.get(t.v1()))
                        + "))\r\n");
                writer.write("f.append(pos=(" + print(t.v2()) + "), normal=(0,0,0), color=(" + c.color(u.get(t.v2()))
                        + "))\r\n");
                writer.write("f.append(pos=(" + print(t.v3()) + "), normal=(0,0,0), color=(" + c.color(u.get(t.v3()))
                        + "))\r\n");
            }
            // writer.write("f.make_normals()\r\n");
            writer.write("scene.ambient=1\r\n");
            writer.write("scene.lights=[]\r\n");
            writer.write("f.make_twosided()\r\n");
            // writer.write("f.smooth()\r\n");
        } catch (final IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Prints VPython code visualising the surface of the given {@link Triangulation} with color coded function values
     * and the transportation field to "{@code output\[fileName].py}".
     * 
     * @param t3
     *            the Triangulation
     * @param u
     *            the values of the function at the vertices
     * @param field
     *            the transportation field
     * @param fileName
     *            the name of the output file
     * @see Color
     */
    public static void print(Triangulation t3, Map<Vector, Double> u, R4toR3 field, String fileName) {
        print(t3, u, field, t3, fileName);
    }

    /**
     * Prints VPython code visualising the surface of the given {@link Triangulation} with color coded function values
     * and the transportation field to "{@code output\[fileName].py}".
     * 
     * @param t3
     *            the Triangulation
     * @param u
     *            the values of the function at the vertices
     * @param field
     *            the transportation field
     * @param fieldT3
     *            the vertices where to draw the field
     * @param fileName
     *            the name of the output file
     * @see Color
     */
    public static void print(Triangulation t3, Map<Vector, Double> u, R4toR3 field, Triangulation fieldT3,
            String fileName) {
        print(t3, u, fileName);
        try (FileWriter writer = new FileWriter("output\\" + fileName + ".py", true)) {
            int counter = 0;
            final Set<Vector> vertices = fieldT3.getVertices();
            final int limit = vertices.size() / 200 + 1;
            for (final Vector p : vertices) {
                if (counter++ % limit == 0) writer.write("arrow(pos=(" + print(p) + "),axis=("
                        + print(field.apply(0., p)) + "),shaftwidth=0.05,color=color.gray(0.5))\r\n");
            }
        } catch (final IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Writes the head of the VPython script with window and viewport settings.
     * 
     * @param writer
     *            the {@link Writer} to write to
     * @param title
     *            the caption of the window
     * @throws IOException
     *             if an I/O error occured
     */
    private static void writeHead(Writer writer, String title) throws IOException {
        writer.write("from visual import *\r\n");
        writer.write("scene.title='" + title + "'\r\n");
        writer.write("scene.up=(0,0,1)\r\n");
        writer.write("scene.autoscale=False\r\n");
        writer.write("scene.width=600\r\n");
        writer.write("scene.height=600\r\n");
        writer.write("scene.background=(1,1,1)\r\n");
        if (title.toLowerCase().contains("torus")) {
            writer.write("scene.forward=(-1,-.5,-.7)\r\n");
            writer.write("scene.scale=(.55,.55,.55)\r\n");
        } else {
            writer.write("scene.forward=(-1,-.1,-.1)\r\n");
            writer.write("scene.scale=(.75,.75,.75)\r\n");
        }
        // writer.write("scene.stereo='passive'\r\n");
    }

    /** Ensures the existence of the output directory. */
    private static void ensureDirectory() {
        new File("output").mkdir();
    }
}