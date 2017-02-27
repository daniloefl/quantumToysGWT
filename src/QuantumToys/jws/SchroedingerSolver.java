/**
 *
 * Java Web Start application, which solves the Schroedinger equation
 * numerically in 1D.
 * <p>
 * It produces an interface, on which the user
 * can write the potential V(x) and it solves the equation:
 * $$\frac{d^2 \psi(x)}{dx^2} + 2 \Big (E - V(x) \Big ) \psi(x) = 0$$
 * <p>
 * This is the Schroedinger equation for a particle with mass 1 in Hartree
 * atomic units.
 * It assumes \(\lim_{x \to \pm \infty} \psi(x) = 0\) and uses the Shooting Method
 * to search for the energy eigenvalue \(E\), given the number of roots of \(\psi(x)\).
 * <p>
 * When logGrid is true, it
 * solves the equation as a function of \(z = \ln(x)\), allowing for more points to
 * sampled close to the origin. This is particularly helpful for the Coulomb
 * central potential, which is extremely strong near the origin and weak elsewhere.
 * <p>
 * As \(\textrm{dz} = \frac{1}{x} \textrm{dx}\):
 * <p>
 *    $$\begin{eqnarray}
 *    \frac{d^2 \psi}{dx^2}&=&\frac{d \Big ( \frac{1}{x} \frac{d \psi}{dz} \Big ) }{dx} \nonumber \\
 *                         &=&-\frac{1}{x^2} \frac{d \psi}{dz} + \frac{1}{x^2} \frac{d^2 \psi}{dz^2} \nonumber \\
 *                         &=&\frac{1}{x^2} \Big (\frac{d^2 \psi}{dz^2} - \frac{d \psi}{dz} \Big) \nonumber
 *    \end{eqnarray}$$
 * <p>
 * Therefore this change of variables generates a first derivative and one would
 * not be able to use the Numerov method. However, one can rewrite the equation
 * as a function of \(\xi(z) = x^{-1/2} \psi(z)\):
 * <p>
 *    $$\begin{eqnarray}
 *     \frac{d \psi(z)}{dz}&=&\frac{d \Big (x^{1/2} \xi(z) \Big )}{dz} \nonumber \\
 *                         &=&\frac{1}{2} x^{1/2} \xi(z) + x^{1/2} \frac{d \xi}{dz} \nonumber \\
 *     \frac{d^2 \psi(z)}{dz^2}&=&\frac{1}{4} x^{1/2} \xi(z) + \frac{1}{2} x^{1/2} \frac{d \xi}{dz} + \nonumber \\
 *                     &&+ \frac{1}{2} x^{1/2} \frac{d \xi}{dz} + x^{1/2} \frac{d^2 \xi}{dz^2} \nonumber \\
 *                     &=& \frac{1}{4} x^{1/2} \xi(z) + x^{1/2} \frac{d \xi}{dz} + x^{1/2} \frac{d^2 \xi}{dz^2} \nonumber
 *    \end{eqnarray}$$
 * <p>
 *   $$\frac{d^2 \psi}{dz^2} - \frac{d \psi}{dz} = - \frac{1}{4} x^{1/2} \xi(z) + x^{1/2} \frac{d^2 \xi}{dz^2}$$
 * <p>
 * That is:
 *   $$\frac{d^2 \psi}{dx^2} = x^{-3/2} \Big ( \frac{d^2 \xi}{dz^2} - \frac{1}{4} \xi \Big )$$ 
 * And:
 *   $$\psi(x) = x^{1/2} \xi(z)$$
 * <p>
 * Therefore one can rewrite the Schroedinger equation as:
 *   $$\frac{d^2 \psi}{dx^2} + 2 \Big (E - V(x)\Big ) \psi(x) = 0$$
 *   $$\frac{d^2 \xi}{dz^2} + 2 x^2 \Big [E - V(x) - \frac{1}{4} \frac{1}{x^2} \Big ] \xi(z) = 0$$
 * <p>
 * If <code>logGrid == true</code>, the equation above is solved
 * and the result from the <code>getDatasetPsi()</code> method is:
 *   $$\psi_{\textrm{normalized}}(x) = \frac{x^{1/2} \xi(z)}{\sqrt{N}}$$
 * Where:
 *   $$N = \int_0^{\infty} \Big (\psi(x) x \Big )^2 \textrm{dx} = \int_0^{\infty} \Big (\psi(x) x \Big )^2 x \textrm{dz}.$$
 * <p>
 * Therefore \(|\psi(x)|^2\) in this case is to be interpreted as the probability of
 * finding the particle in the region \([x, x+dx]\).
 * <p>
 * If <code>logGrid == false</code>, the original Schroedinger equation is solved
 * in a linear Grid and the result of <code>getDatasetPsi()</code> is normalised such that:
 *   $$\psi_{\textrm{normalized}}(x) = \frac{\psi(x)}{\sqrt{N}}$$
 * where:
 *   $$N = \int_{-\infty}^{\infty} \psi^2(x) \textrm{dx}.$$
 * <p>
 *
 * The actual calculation is done by the SchroedingerCalculator class.
 *
 * @author Danilo Enoque Ferreira de Lima <daniloefl@gmail.com>
 * @version %I%, %G%
 * @see QuantumToys.SchroedingerCalculator
 *
 */

package QuantumToys.jws;

import QuantumToys.shared.SchroedingerCalculator;

import javax.jnlp.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.net.*;

import java.util.*;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.*;
import org.jfree.data.xy.*;
import org.jfree.data.category.*;
import org.jfree.ui.*;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.annotations.*;
import org.jfree.chart.axis.NumberAxis;

/**
 * SchroedingerSolver implements the graphical interface for this Java Web application.
 * It generates a window and draws all the controls. All calculations are done by SchroedingerCalculator.
 * And sent to the JFreeChart as a dataset.
 *
 * @author Danilo Enoque Ferreira de Lima <daniloefl@gmail.com>
 * @version %I%, %G%
 * @see QuantumToys.shared.SchroedingerCalculator
 */
public class SchroedingerSolver implements ActionListener, ItemListener {

    static BasicService basicService = null;

    /**
     * Window.
     */
    private JFrame m_frame;

    /**
     * Label indicating number of nodes.
     */
    private JLabel m_labelN;

    /**
     * Label indicating potential text box.
     */
    private JLabel m_labelPot;

    /**
     * Number of nodes text box.
     */
    private JTextField m_n;

    /**
     * Potential text box.
     */
    private JTextField m_pot;

    /**
     * Minimum x value label.
     */
    private JLabel m_labelXmin;

    /**
     * Maximum x value label.
     */
    private JLabel m_labelXmax;

    /**
     * Number of points in the grid label.
     */
    private JLabel m_labelNpoints;

    /**
     * Minimum and maximum x value text box.
     */
    private JTextField m_xmin;
    private JTextField m_xmax;

    /**
     * Number of points text box.
     */
    private JTextField m_Npoints;


    /**
     * Update button to recalculate and redraw.
     */
    private JButton m_go;


    /**
     * Check box that automatically calculates x range.
     */
    private JCheckBox m_autoRange;

    /**
     * Check box that switches to the radial equation solver.
     */
    private JCheckBox m_logGrid;

    /**
     * Panel showing the actual plots.
     */
    private ChartPanel m_plot;

    /**
     * JFreeChart object holding plot panel.
     */
    private JFreeChart m_lineChart;

    /**
     * Line renderers for wave function and potential plots, used to select style.
     */
    private XYLineAndShapeRenderer m_rendererPsi;
    private XYLineAndShapeRenderer m_rendererEnergy;

    /**
     * SchroedingerCalculator object, which does the actual calculation.
     */
    private SchroedingerCalculator m_calc;

    /**
     * Text annotation added to the plot to indicate the numerical energy eigenvalue solution.
     */
    private JLabel m_eAnnotation;

    /**
     * Second Y axis, on the left, showing energy values.
     */
    private NumberAxis m_eAxis;

    /**
     * Labels for number of iterations and parameter used to scale energy steps.
     */
    private JLabel m_labelIter;
    private JLabel m_labelAlpha;


    /**
     * Text boxes for number of iterations and parameter used to scale energy steps.
     */
    private JTextField m_iter;
    private JTextField m_alpha;

    /**
     * Label shown when method does not converge to solution within the number of iterations requested.
     */
    private JLabel m_labelConverge;

    /**
     * Constructor.
     */
    public SchroedingerSolver() {
    }

    /**
     * Method called when a check box is selected.
     * @param e (required) Indicates which item was selected.
     */
    public void itemStateChanged(ItemEvent e) {
        Object source = e.getItemSelectable();
        // If the "Auto range" check box was clicked, calculate the
        // x range automatically using the autoRange() method
        // and update the x range and number of points text boxes.
        // Also disable those text boxes as we control those parameters automatically now.
        // If it has been de-selected, make the text boxes editable again.
        if (source == m_autoRange) {
            if (m_autoRange.isSelected()) {
                m_calc.autoRange();
                m_xmin.setText(Double.toString(m_calc.xmin));
                m_xmax.setText(Double.toString(m_calc.xmax));
                m_Npoints.setText(Integer.toString(m_calc.N));
                m_Npoints.setEditable(false);
                m_xmin.setEditable(false);
                m_xmax.setEditable(false);
            } else {
                m_Npoints.setEditable(true);
                m_xmin.setEditable(true);
                m_xmax.setEditable(true);
            }
            // And update the plot with the new ranges
            updatePlot();
        } else if (source == m_logGrid) {
            // If the "Radial equation" check box was selected,
            // inform the SchroedingerCalculator class so that
            // next time it looks for a solution, it solves the
            // appropriate equation.
            // Also update the x range and energy step text boxes, as
            // this option automatically updates those parameters
            // to avoid having an x range that includes negative values when
            // we are plotting in log scale.
            // It also chooses parameters more suitable for the Coulomb potential.
            if (m_logGrid.isSelected()) {
                m_calc.setLogGrid(true);
                m_xmin.setText(Double.toString(m_calc.xmin));
                m_xmax.setText(Double.toString(m_calc.xmax));
                m_alpha.setText(Double.toString(m_calc.m_alpha));
                m_pot.setText("-1/x+0.25/(2*x^2)"); // set potential to the Hydrogen potential
            } else {
                m_calc.setLogGrid(false);
                m_xmin.setText(Double.toString(m_calc.xmin));
                m_xmax.setText(Double.toString(m_calc.xmax));
                m_alpha.setText(Double.toString(m_calc.m_alpha));
                m_pot.setText(m_calc.def_pot); // set potential to the harmonic oscillator
            }
        }
    }

    /**
     * Get the JFreeChart dataset with plots to be shown.
     * @return An XYDataset with the $\psi$ and $|\psi^2|$ plots.
     */
    public XYDataset getDatasetPsi() {
        // Create plots
        // Remember that in the radial equation, the X axis will show
        // exp(grid).
        double [] x = new double[m_calc.N];
        double [] psi = new double[m_calc.N];
        m_calc.getDatasetPsi(x, psi);
        XYSeries p = new XYSeries("psi(x)");
	for (int i = 0; i < m_calc.N; ++i) {
	    p.add(x[i], psi[i]);
	}
        XYSeries p2 = new XYSeries("|psi(x)|^2");
	for (int i = 0; i < m_calc.N; ++i) {
	    p2.add(x[i], psi[i]*psi[i]);
	}

	XYSeriesCollection dataset = new XYSeriesCollection();
	dataset.addSeries(p);
	dataset.addSeries(p2);
        return dataset;
    }

    /**
     * Get a flat plot with the energy value and the potential energy curve as a JFreeChart plot.
     * @return A JFreeChart dataset with two plots: the energy and the potential energy.
     */
    public XYDataset getDatasetEnergy() {
        // Create plots
        // Remember that in the radial equation, the X axis is actually
        // exp(grid)
        double [] x = new double[m_calc.N];
        double [] potential = new double[m_calc.N];
        double E = m_calc.getDatasetEnergy(x, potential);
        XYSeries pot = new XYSeries("Potential(x)");
	for (int i = 0; i < m_calc.N; ++i) {
	    pot.add(x[i], potential[i]);
	}
        XYSeries energy = new XYSeries("Energy");
	for (int i = 0; i < m_calc.N; ++i) {
	    energy.add(x[i], E);
	}
	XYSeriesCollection dataset = new XYSeriesCollection();
	dataset.addSeries(pot);
	dataset.addSeries(energy);
        return dataset;
   }

    /**
     * Update the plot, by getting the JFreeChart datasets from the SchroedingerCalculator.
     * This will update the existing plot to the new one.
     * It also adjusts the energy axis range and writes the new energy value as an annotation.
     */
    public void updatePlot() {
	this.m_lineChart.getXYPlot().setDataset(0, getDatasetPsi());
	this.m_lineChart.getXYPlot().setDataset(1, getDatasetEnergy());
	double rminE = m_calc.E-5*Math.abs(m_calc.E);
	double rmaxE = m_calc.E+5*Math.abs(m_calc.E);
	if (rminE < m_calc.Vmin) {
	    rminE = m_calc.Vmin;
	}
	if (rmaxE > m_calc.Vmax) {
	    rmaxE = m_calc.Vmax;
	}
	m_lineChart.getXYPlot().getRangeAxis(1).setRange(rminE, rmaxE);
        this.m_eAnnotation.setText(String.format("Energy = %4f a.u. = %4f eV", m_calc.E, m_calc.E*m_calc.eV));
	this.m_eAnnotation.setVisible(true);
	this.m_lineChart.setNotify(true);
	this.m_plot.repaint();
	if (!m_calc.converged) {
	    m_labelConverge.setVisible(true);
	} else {
	    m_labelConverge.setVisible(false);
	}
    }

    /**
     * Method called when the "Update" button is clicked.
     * It will read all text boxes and check boxes, update the SchroedingerCalculator
     * parameters, and call recalculatePotential to update the internal potential array,
     * and then try to calculate the wave function and energy.
     * It calls updatePlot to redraw the results.
     * @param ae (optional) Indicates which action happened.
     */
    public void actionPerformed(ActionEvent ae) {
	try {
	    this.m_calc.setN(Integer.parseInt(m_n.getText()));
	    this.m_calc.setX(Double.parseDouble(m_xmin.getText()), Double.parseDouble(m_xmax.getText()), Integer.parseInt(m_Npoints.getText()));
	    this.m_calc.setIterStep(Integer.parseInt(m_iter.getText()), Double.parseDouble(m_alpha.getText()));
	} catch (NumberFormatException e) {
	    this.m_calc.setN(m_calc.def_n);
	    this.m_calc.setX(m_calc.def_xmin, m_calc.def_xmax, m_calc.def_N);
	    this.m_calc.setIterStep(m_calc.def_iter, m_calc.def_alpha);
	}
        if (m_autoRange.isSelected()) {
	    this.m_calc.recalculatePotential(m_pot.getText());
            this.m_calc.autoRange();
            m_xmin.setText(Double.toString(m_calc.xmin));
            m_xmax.setText(Double.toString(m_calc.xmax));
            m_Npoints.setText(Integer.toString(m_calc.N));
            m_xmin.setEditable(false);
            m_xmax.setEditable(false);
        } else {
	    this.m_calc.recalculatePotential(m_pot.getText());
	    this.m_calc.recalculate();
        }
        updatePlot();
    }

    /**
     * Method called when Java Web Start application is loaded.
     * It makes the window and creates all graphical objects.
     * It also creates a new instance of the SchroedingerCalculator class,
     * used for the calculations.
     */
    public void load() {

        // Create a new window.
        this.m_frame = new JFrame("Schroedinger equation solver");

        // Loads the Java Web Start basic service
        try {
          basicService = (BasicService)  ServiceManager.lookup("javax.jnlp.BasicService");
        } catch (UnavailableServiceException e) {
          System.err.println("Lookup failed: " + e);
        }

        // instantiate the calculator, to make the math
	this.m_calc = new SchroedingerCalculator(new JWSParser());
	this.m_calc.recalculate(); // recalculate the default function to show something when starting

        // Now load all the graphical objects in their positions
        this.m_logGrid = new JCheckBox("Radial equation");
        this.m_logGrid.setBounds(0, 0, 200, 25);
        this.m_frame.getContentPane().add(this.m_logGrid);

        this.m_labelConverge = new JLabel();
        this.m_labelConverge.setText("Warning: Did not converge!");
	this.m_labelConverge.setBounds(0, 600, 200, 25);
        this.m_frame.getContentPane().add(this.m_labelConverge);
	m_labelConverge.setVisible(false); // set this invisible so that it is only visible if calculations fail

        this.m_labelXmin = new JLabel();
        this.m_labelXmin.setText("min(x):");
	this.m_labelXmin.setBounds(200, 0, 50, 25);
        this.m_frame.getContentPane().add(this.m_labelXmin);

        this.m_xmin = new JTextField(4);
	this.m_xmin.setBounds(250, 0, 50, 25);
	this.m_xmin.setText(Double.toString(m_calc.def_xmin));
	this.m_frame.getContentPane().add(this.m_xmin);

        this.m_labelXmax = new JLabel();
        this.m_labelXmax.setText("max(x):");
	this.m_labelXmax.setBounds(300, 0, 50, 25);
        this.m_frame.getContentPane().add(this.m_labelXmax);

        this.m_xmax = new JTextField(4);
	this.m_xmax.setBounds(350, 0, 50, 25);
	this.m_xmax.setText(Double.toString(m_calc.def_xmax));
	this.m_frame.getContentPane().add(this.m_xmax);

        this.m_labelNpoints = new JLabel();
        this.m_labelNpoints.setText("# points:");
	this.m_labelNpoints.setBounds(400, 0, 100, 25);
        this.m_frame.getContentPane().add(this.m_labelNpoints);

        this.m_Npoints = new JTextField(4);
	this.m_Npoints.setBounds(500, 0, 50, 25);
	this.m_Npoints.setText(Integer.toString(m_calc.def_N));
	this.m_frame.getContentPane().add(this.m_Npoints);

        this.m_go = new JButton("Update");
	this.m_go.setBounds(550, 0, 100, 25);
        this.m_frame.getContentPane().add(this.m_go);

        this.m_autoRange = new JCheckBox("Auto range");
        this.m_autoRange.setBounds(650, 0, 150, 25);
        this.m_frame.getContentPane().add(this.m_autoRange);

        this.m_labelN = new JLabel();
        this.m_labelN.setText("# nodes:");
	this.m_labelN.setBounds(0, 25, 100, 25);
        this.m_frame.getContentPane().add(this.m_labelN);

        this.m_n = new JTextField(2);
	this.m_n.setBounds(100, 25, 50, 25);
	this.m_n.setText(Integer.toString(m_calc.def_n));
	this.m_frame.getContentPane().add(this.m_n);

        this.m_labelPot = new JLabel();
        this.m_labelPot.setText("Potential:");
	this.m_labelPot.setBounds(150, 25, 100, 25);
        this.m_frame.getContentPane().add(this.m_labelPot);

        this.m_pot = new JTextField(100);
	this.m_pot.setBounds(250, 25, 200, 25);
	this.m_pot.setText(m_calc.def_pot);
	this.m_frame.getContentPane().add(this.m_pot);

        this.m_labelIter = new JLabel();
        this.m_labelIter.setText("# iterations:");
	this.m_labelIter.setBounds(450, 25, 100, 25);
        this.m_frame.getContentPane().add(this.m_labelIter);

        this.m_iter = new JTextField(5);
	this.m_iter.setBounds(550, 25, 50, 25);
	this.m_iter.setText(Integer.toString(m_calc.def_iter));
	this.m_frame.getContentPane().add(this.m_iter);

        this.m_labelAlpha = new JLabel();
        this.m_labelAlpha.setText("Energy step:");
	this.m_labelAlpha.setBounds(600, 25, 100, 25);
        this.m_frame.getContentPane().add(this.m_labelAlpha);

        this.m_alpha = new JTextField(5);
	this.m_alpha.setBounds(700, 25, 50, 25);
	this.m_alpha.setText(Double.toString(m_calc.def_alpha));
	this.m_frame.getContentPane().add(this.m_alpha);

        // main controls are loaded

        // make the chart object, to show the result
        this.m_lineChart = ChartFactory.createXYLineChart("",
                                                          "x", "psi(x) or |psi(x)|^2",
                                                          getDatasetPsi(),
                                                          PlotOrientation.VERTICAL,
                                                          true,true,false);
        // make the actual panel graphical object
        this.m_plot = new ChartPanel(this.m_lineChart);
        this.m_plot.setPreferredSize(new Dimension(800, 600));
	this.m_plot.setLayout(null);
	this.m_plot.setBounds(0, 50, 800, 600);
        this.m_plot.setMouseWheelEnabled(true);
        this.m_plot.setHorizontalAxisTrace(true);
        this.m_plot.setVerticalAxisTrace(true);

        // this sets the line styles for the plots
	this.m_rendererPsi = new XYLineAndShapeRenderer();
	this.m_rendererPsi.setSeriesPaint(0, Color.BLUE);
	this.m_rendererPsi.setSeriesPaint(1, Color.BLACK);
	this.m_rendererPsi.setSeriesStroke(0, new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 1.0f));
	this.m_rendererPsi.setSeriesStroke(1, new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 1.0f));
	this.m_rendererPsi.setSeriesShapesVisible(0, false);
	this.m_rendererPsi.setSeriesShapesVisible(1, false);
        this.m_lineChart.getXYPlot().setRenderer(0, this.m_rendererPsi);

	this.m_rendererEnergy = new XYLineAndShapeRenderer();
	this.m_rendererEnergy.setSeriesPaint(0, Color.RED);
	this.m_rendererEnergy.setSeriesPaint(1, new Color(0, 100, 0));
	this.m_rendererEnergy.setSeriesStroke(0, new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 1.0f));
	this.m_rendererEnergy.setSeriesStroke(1, new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 1.0f));
	this.m_rendererEnergy.setSeriesShapesVisible(0, false);
	this.m_rendererEnergy.setSeriesShapesVisible(1, false);
        this.m_lineChart.getXYPlot().setRenderer(1, this.m_rendererEnergy);

        // set the axes font style
	this.m_lineChart.getXYPlot().getDomainAxis().setLabelFont(new Font("SansSerif", Font.BOLD, 12));
	this.m_lineChart.getXYPlot().getDomainAxis().setTickLabelFont(new Font("SansSerif", Font.BOLD, 10));

	this.m_lineChart.getXYPlot().getRangeAxis().setLabelFont(new Font("SansSerif", Font.BOLD, 12));
	this.m_lineChart.getXYPlot().getRangeAxis().setTickLabelFont(new Font("SansSerif", Font.BOLD, 10));

        // create a new axis on the left side of the plot
        // All energy values are to be read in the left axis
        // This allows one to change the axis ranges independently
	m_eAxis = new NumberAxis("Energy [Hartree a.u.]");
	m_eAxis.setLabelFont(new Font("SansSerif", Font.BOLD, 12));
	m_eAxis.setTickLabelFont(new Font("SansSerif", Font.BOLD, 10));
	this.m_lineChart.getXYPlot().setRangeAxis(1, m_eAxis);

        // inform JFreeChart that the energy dataset is to be mapped to this axis
	this.m_lineChart.getXYPlot().setDataset(1, getDatasetEnergy());
	this.m_lineChart.getXYPlot().mapDatasetToRangeAxis(1, 1);

        // create the energy annotation in the top of the plot
        this.m_eAnnotation = new JLabel();
        this.m_eAnnotation.setText(String.format("Energy = %4f a.u. = %4f eV", m_calc.E, m_calc.E*m_calc.eV));
	this.m_eAnnotation.setBounds(500, 600, 300, 25);
        this.m_frame.getContentPane().add(this.m_eAnnotation);

        // now indicate that if the "Update" button is clicked
        // or the check boxes are selected, call the actionEvent and itemEvent methods
        // of the current class
        this.m_go.addActionListener(this);
        this.m_autoRange.addItemListener(this);
        this.m_logGrid.addItemListener(this);

        // add the plot panel into the window
	this.m_frame.getContentPane().add(this.m_plot);

        // set the window size, background and default configuration
        this.m_frame.setResizable(false);
        this.m_frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	this.m_frame.setLayout(null);
	this.m_frame.setSize(800, 650);
	this.m_frame.getContentPane().setPreferredSize(new Dimension(800, 650));
	this.m_frame.setBackground(Color.WHITE);
	this.m_frame.setLocationRelativeTo(null);

        // update plot for the first time to make sure all parameters are in sync
	updatePlot();

        // show window
	this.m_frame.pack();
	this.m_frame.setVisible(true);

        // program exits when window is closed
    }

    /**
     * Main function, called when the application starts.
     * Just creates an instance of the SchroedingerSolver class and
     * calls load() to create all the graphical objects.
     * @param args (optional) Command line arguments, if any.
     */
    public static void main(String args[]) {
        SchroedingerSolver m_main = new SchroedingerSolver();
	m_main.load();
    }

}

