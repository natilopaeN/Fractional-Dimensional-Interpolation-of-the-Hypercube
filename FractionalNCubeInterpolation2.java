import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class FractionalNCubeInterpolation2 extends JPanel {
    private double dimension = 3.0;

    // rotation angles for 6 planes in 4D
    private double rotationX = 0.0, rotationY = 0.0, rotationZ = 0.0;
    private double rotationW = 0.0, rotationV = 0.0, rotationU = 0.0;

    // sliders for dimension / projection zoom
    private final JSlider dimSlider;
    private final JSlider projectionDistanceSlider;

    // distance used in stereographic projection (controls zoom)
    private double projectionDistance = 4.0;

    public FractionalNCubeInterpolation2() {
        setPreferredSize(new Dimension(800, 800));
        setBackground(Color.BLACK);
        setLayout(new BorderLayout());

        // Dimension slider: ranges from 1.00 to 4.00 (scaled by 100)
        dimSlider = new JSlider(100, 400, 300);
        dimSlider.setMajorTickSpacing(50);
        dimSlider.setPaintTicks(true);
        dimSlider.setPaintLabels(true);
        dimSlider.addChangeListener(e -> {
            dimension = dimSlider.getValue() / 100.0;
            repaint();
        });

        // Projection distance slider: affects stereographic scaling
        projectionDistanceSlider = new JSlider(100, 1000, 400);
        projectionDistanceSlider.setMajorTickSpacing(100);
        projectionDistanceSlider.setPaintTicks(true);
        projectionDistanceSlider.setPaintLabels(true);
        projectionDistanceSlider.addChangeListener(e -> {
            projectionDistance = projectionDistanceSlider.getValue() / 100.0;
            repaint();
        });

        JPanel controls = new JPanel(new GridLayout(4, 1));
        controls.add(new JLabel("Dimension (1.00 to 4.00)"));
        controls.add(dimSlider);
        controls.add(new JLabel("Projection Distance (Zoom)"));
        controls.add(projectionDistanceSlider);

        add(controls, BorderLayout.SOUTH);

        // animation timer: updates rotation angles continuously
        Timer timer = new Timer(40, e -> {
            rotationX += 0.02;
            rotationY += 0.015;
            rotationZ += 0.01;
            rotationW += 0.012;
            rotationV += 0.008;
            rotationU += 0.01;
            repaint();
        });
        timer.start();
    }

    // Lanczos approximation for Gamma function
    private double gamma(double x) {
        double[] p = {
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };
        int g = 7;

        if (x < 0.5) {
            return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
        } else {
            x -= 1;
            double a = 0.99999999999980993;
            for (int i = 0; i < p.length; i++) {
                a += p[i] / (x + i + 1);
            }
            double t = x + g + 0.5;
            return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
        }
    }

    // Numerical integration using trapezoidal rule over [a,b] for (d - t)^(d - 1)
    private double integrate(double a, double b, double d, int steps) {
        if (b <= a) return 0;

        double h = (b - a) / steps;
        double sum = 0.5 * (Math.pow(d - a, d - 1) + Math.pow(d - b, d - 1));
        for (int i = 1; i < steps; i++) {
            double t = a + i * h;
            sum += Math.pow(d - t, d - 1);
        }
        return sum * h;
    }

    // fractional weight for coordinate index i at fractional dimension d
    private double fractionalWeight(double d, int i) {
        double lower = i - 1;
        double upper = Math.min(i, d);
        if (upper <= lower) return 0; // no contribution if dimension less than coordinate segment start

        double numerator = integrate(lower, upper, d, 1000); // 1000 steps integration for smoothness
        double denominator = gamma(d);
        return numerator / denominator;
    }

    // applies fractional integral weighting to a coordinate
    private double weightedCoord(double coord, double d, int i) {
        double weight = fractionalWeight(d, i);
        return coord * weight;
    }

    // 4D cube vertices w/ fractional weights applied
    private double[][] generateFractionalCubeVertices(double edgeLength, double d) {
        int maxDimension = 4;
        int vertexCount = 1 << maxDimension; // 2^4 = 16 vertices
        double half = edgeLength / 2.0;
        double[][] vertices = new double[vertexCount][maxDimension];

        for (int v = 0; v < vertexCount; v++) {
            for (int i = 0; i < maxDimension; i++) {
                // assign coordinate based on bit value of vertex index
                double originalCoord = ((v & (1 << i)) != 0) ? half : -half;
                // apply dimension-based fractional integral weighting
                vertices[v][i] = weightedCoord(originalCoord, d, i + 1);
            }
        }
        return vertices;
    }

    // generate edge list for a 4D hypercube: edges connect vertices that differ in 1 bit
    private List<int[]> generateEdges(int dimension) {
        int vertexCount = 1 << dimension;
        List<int[]> edges = new ArrayList<>();
        for (int v = 0; v < vertexCount; v++) {
            for (int i = 0; i < dimension; i++) {
                int neighbor = v ^ (1 << i); // flip the ith bit
                if (v < neighbor) {
                    edges.add(new int[]{v, neighbor});
                }
            }
        }
        return edges;
    }

    // perform 6D rotation (4D space with 6 rotation planes)
    private double[] rotate6D(double[] v, double rx, double ry, double rz, double rw, double rv, double ru) {
        double[] r = v.clone();
        double cos, sin, x, y;

        // 2D plane rotations in 4D space
        // XY rotation
        cos = Math.cos(rx); sin = Math.sin(rx);
        x = r[0]; y = r[1];
        r[0] = x * cos - y * sin; r[1] = x * sin + y * cos;

        // YZ rotation
        cos = Math.cos(ry); sin = Math.sin(ry);
        x = r[1]; y = r[2];
        r[1] = x * cos - y * sin; r[2] = x * sin + y * cos;

        // ZX rotation
        cos = Math.cos(rz); sin = Math.sin(rz);
        x = r[2]; y = r[0];
        r[2] = x * cos - y * sin; r[0] = x * sin + y * cos;

        // XW rotation
        cos = Math.cos(rw); sin = Math.sin(rw);
        x = r[0]; y = r[3];
        r[0] = x * cos - y * sin; r[3] = x * sin + y * cos;

        // YW rotation
        cos = Math.cos(rv); sin = Math.sin(rv);
        x = r[1]; y = r[3];
        r[1] = x * cos - y * sin; r[3] = x * sin + y * cos;

        // ZW rotation
        cos = Math.cos(ru); sin = Math.sin(ru);
        x = r[2]; y = r[3];
        r[2] = x * cos - y * sin; r[3] = x * sin + y * cos;

        return r;
    }

    // 4D point -> 3D point stereographic projection
    private double[] stereographicProjection(double[] v) {
        double w = projectionDistance - v[3]; // denominator
        if (Math.abs(w) < 0.01) w = 0.01; // prevent near-infinite scaling
        return new double[]{v[0] / w, v[1] / w, v[2] / w};
    }

    // draw edges of projected fractional cube
    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        int w = getWidth() / 2, h = getHeight() / 2;

        // generate weighted vertices and edges for fractional 4-cube
        double[][] vertices = generateFractionalCubeVertices(2.0, dimension);
        List<int[]> edges = generateEdges(4); // edges of 4D hypercube

        g2d.setColor(Color.WHITE);

        // draw each edge after rotating and projecting its endpoints
        for (int[] edge : edges) {
            double[] p1 = rotate6D(vertices[edge[0]], rotationX, rotationY, rotationZ, rotationW, rotationV, rotationU);
            double[] p2 = rotate6D(vertices[edge[1]], rotationX, rotationY, rotationZ, rotationW, rotationV, rotationU);

            double[] proj1 = stereographicProjection(p1);
            double[] proj2 = stereographicProjection(p2);

            if (Double.isNaN(proj1[0]) || Double.isNaN(proj2[0])) continue;

            // convert 3D to 2D screen coordinates
            int x1 = (int) (w + proj1[0] * 300);
            int y1 = (int) (h - proj1[1] * 300);
            int x2 = (int) (w + proj2[0] * 300);
            int y2 = (int) (h - proj2[1] * 300);

            g2d.drawLine(x1, y1, x2, y2);
        }

        // Draw info text
        g2d.setColor(Color.GREEN);
        g2d.drawString(String.format("Dimension: %.2f", dimension), 10, 20);
        g2d.drawString(String.format("Projection Distance: %.2f", projectionDistance), 10, 40);
    }

    public static void main(String[] args) {
        JFrame frame = new JFrame("Fractional Dimensional 4D Cube Interpolation");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(new FractionalNCubeInterpolation2());
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}
