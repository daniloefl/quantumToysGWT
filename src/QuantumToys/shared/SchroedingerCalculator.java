package QuantumToys.shared;

import java.util.*;

/**
 * SchroedingerCalculator solves the Schroedinger equation, given by a generic potential.
 * The only restriction is that the wave function goes to zero at both end points.
 *
 * Numerically solves the Schroedinger equation in 1D using
 * the Numerov method to approximate derivatives.
 * It uses exp4j to read the potential V(x) from a String.
 * <p>
 * It solves the following equation by default.
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
 *
 *
 * @author Danilo Enoque Ferreira de Lima <daniloefl@gmail.com>
 * @version %I%, %G%
 */
public class SchroedingerCalculator {

    /**
     * One electron-volt in Hartree atomic units.
     */
    public double eV = 27.2113966413442;

    /**
     * Object that implements mathematical expression parser.
     */
    public Parser m_p;

    /**
     * Number of nodes in the solution.
     */
    public int m_n;

    /**
     * Energy eigenvalue.
     */
    public double E;

    /**
     * Minimum and maximum x value.
     */
    public double xmin;
    public double xmax;

    /**
     * Whether to solve the radial equation or the original one.
     */
    public boolean logGrid = false;

    /**
     * Static values with the default x ranges, number of points, potential shape,
     * number of nodes, number of iterations and energy step size.
     */
    public double def_xmin = -8;
    public double def_xmax = 8;
    public int def_N = 200;
    public String def_pot = "x*x*0.5";
    public int def_n = 0;
    public int def_iter = 100;
    public double def_alpha = 1;

    /**
     * Whether the numerical method has converged to a solution or not.
     */
    public boolean converged;

    /**
     * Values of x and the number of points in the Grid.
     */
    double [] x;
    public int N;

    /**
     * Potential expression in exp4j class and String formats.
     */
    public String exp_V_str;

    /**
     * Potential values in the Grid.
     * If the radial equation option is activated, this is V(z(x)). Otherwise, it is V(x).
     */
    double [] V;

    /**
     * Minimum and maximum values of the potential (to help drawer).
     */
    public double Vmin;
    public double Vmax;

    /**
     * Minimum and maximum energy values to scan.
     */
    double Emin = 0;
    double Emax = 0;

    /**
     * Value of (Emax - Emin) required for convergence.
     */
    double eps = 1e-12;

    /**
     * Number of iterations to perform and value of the energy step.
     */
    public int m_iter;
    public double m_alpha;

    /**
     * Arrays keeping the value of $\psi(x)$ (or $\psi(z(x))$ if using the log-grid).
     * Also keeping the wave functions when assuming only initial condition in xmin or only in xmax for debugging.
     */
    double [] psi;
    double [] psi_in;
    double [] psi_out;

    /**
     * Constructor.
     * Sets x range, number of points, energy step and potential to default configurations.
     */
    public SchroedingerCalculator(Parser p) {
        m_p = p;
	m_n = def_n;
	setX(def_xmin, def_xmax, def_N);
	m_iter = def_iter;
	m_alpha = def_alpha;
	recalculatePotential(def_pot);
    }

    /**
     * Setter for the log-grid.
     * Also sets the x range when selected, to avoid negative x values in a log-grid.
     * @param lg (required) Whether to solve the radial equation or the 1D linear equation.
     */
    public void setLogGrid(boolean lg) {
        logGrid = lg;
        if (logGrid) {
            xmin = Math.exp(-9);
            xmax = Math.exp(3);
            m_alpha = 0.1;
        } else {
            xmin = def_xmin;
            xmax = def_xmax;
            m_alpha = def_alpha;
        }
    }

    /**
     * Automatically calculates the x range.
     * It calculates the wave function and energy for a set of x range configurations,
     * until the energy eigenvalue does not change much.
     * The resulting configuration is set in xmin, xmax and N.
     */
    public void autoRange() {
        double modx = 2.0; // start with some random setup
	if (logGrid) modx = 25.0;
        double dmodx = 0.05; // increase it in steps of this
        double E_last = 0;
        E = 0;
        double dEdmodx = 0;
        double peaks = 1;
        // while the energy still changes significantly, increase the x range by some step
        // and recalculate the energy eigenvalue
        do {
            E_last = E;
            modx += dmodx;
            if (logGrid) {
                setX(1e-9, modx, (int) (2*modx*20*peaks));
            } else {
                setX(-modx, modx, (int) (2*modx*20*peaks));
            }
	    recalculatePotential(exp_V_str);
	    recalculate();
            peaks = 1;
            for (int k = 1; k < N-3; ++k) {
                if (Math.abs(psi[k]) < Math.abs(psi[k+1]) && Math.abs(psi[k+1]) > Math.abs(psi[k+2]))
                  peaks += 1;
            }
            //System.out.format("AutoRange: Range [%f, %f], with %d points -- Energy = %f, delta E = %f, delta E/E = %f, peaks = %f\n", xmin, xmax, N, E, E-E_last, Math.abs(E-E_last)/Math.abs(E), peaks);
            dEdmodx = (E - E_last)/dmodx;
            if (E >= 0.9*Emax) { // require that we are not so close to the maximum value of energy allowed
                continue;
            }
            // check change in energy eigenvalue
        } while (Math.abs(E-E_last)/Math.abs(E) > 1e-6 && Math.abs(E-E_last)/Math.abs(E) <= 1.0);

        // if suddenly the energy changed by a lot, stop and
        // move back, toward a safe point
        if (Math.abs(E-E_last)/Math.abs(E) > 1.0) {
            modx -= 2*dmodx;
        }
        modx = (double) Math.ceil(modx*10)/((int) 10);
        if (logGrid) {
            setX(1e-9, modx, (int) (2*modx*20*peaks));
        } else {
            setX(-modx, modx, (int) (2*modx*20*peaks));
        }

        // recalculate it all
	recalculatePotential(exp_V_str);
	recalculate();
    }

    /**
     * Set number of iterations and energy step size.
     * @param iter (required) Number of iterations.
     * @param alpha (required) Energy step.
     */
    public void setIterStep(int iter, double alpha) {
        m_iter = iter;
	m_alpha = alpha;
    }

    /**
     * Set x range and number of grid points.
     * If in the radial equation setup, store the grid as going from log(x0)
     * to log(xf), so that the X axis in the plot will show exp(x[k]).
     * @param x0 (required) minimum value of x
     * @param xf (required) maximum value of x
     * @param Np (required) number of grid points
     */
    public void setX(double x0, double xf, int Np) {
        xmin = x0;
	xmax = xf;
	N = Np;
	x = null;
	V = null;
	psi = null;
	psi_in = null;
	psi_out = null;

        x = new double[N];
        V = new double[N];
	psi = new double[N];
	psi_in = new double[N];
	psi_out = new double[N];
        if (logGrid) {
	    double dx = (Math.log(xmax) - Math.log(xmin))/N;
	    for (int i = 0; i < N; ++i) {
	        x[i] = Math.log(xmin) + i*dx;
	    }
        } else {
	    double dx = (xmax - xmin)/N;
	    for (int i = 0; i < N; ++i) {
	        x[i] = xmin + i*dx;
	    }
        }
	for (int i = 0; i < N; ++i) {
	    V[i] = 0;
	}
	for (int i = 0; i < N; ++i) {
	    psi[i] = 0;
	    psi_in[i] = 0;
	    psi_out[i] = 0;
	}
    }

    /**
     * Sets number of nodes desired.
     * @param n (required) number of nodes
     */
    public void setN(int n) {
        m_n = n;
    }


    /**
     * Recalculate the potential from a string.
     * This uses exp4j to get the analytical form of the potential in the string s
     * and calculate it for each grid point.
     * @param s (required) Potential analytical form.
     */
    public void recalculatePotential(String s) {
        exp_V_str = s;
	Vmin = 1e10;
	Vmax = -1e10;
        // for the radial equation, exp(x[k]) is to be interpreted as the
        // distance to the origin. The plot shows V(exp(grid)) vs exp(grid)
        // therefore the x values stored are not what is read as "x" in the graphical version.
        for (int i = 0; i < N; ++i) {
	    if (logGrid) {
	        V[i] = m_p.evaluate(s, "x", Math.exp(x[i]));
		if (V[i] == Double.NaN) {
		    V[i] = 0;
		}
	    } else {
	        V[i] = m_p.evaluate(s, "x", x[i]);
		if (V[i] == Double.NaN) {
		    V[i] = 0;
		}
	    }
            // also save the maximum and minimum values of V to help drawing it
	    if (V[i] < Vmin) Vmin = V[i];
	    if (V[i] > Vmax) Vmax = V[i];
	}
    }

    /**
     * Calculates the auxiliary vectors a and f, given some energy.
     * To solve the equation
     * $$\frac{g(x)}{dx} + a(x) g(x) = 0$$
     * using the Numerov method, this calculates $a(x)$
     * and sets the value of $f(x) = 1 + a(x) \frac{\textrm{dx}^2}{12}$,
     * where $\textrm{dx}$ is the grid step size in a linear grid.
     * If using the radial equation, one must adapt the equation so that
     * the equation form is as shown above and the grid steps are constant.
     * See the comment in the beginning of the document on how the radial
     * equation is transformed into an equation in a linear parameter $z$.
     * @param energy (required) Value of energy currently probed.
     * @param a (required) Array of doubles to be set with values of the function $a(x)$.
     * @param f (required) Array of doubles to be set with values of the auxiliary function $f(x)$.
     * @return Index on the grid, where the $a(x)$ changes sign (ie: energy = potential).
     */
    int getAuxVectors(double energy, double [] a, double [] f) {
        double m = 1.0; // mass
        int icl = -1;
	for (int i = 0; i < N; ++i) {
	    double dx = 0;
            if (i == N-1) {
                dx = (x[i] - x[i-1]);
            } else {
                dx = (x[i+1] - x[i]);
            }
	    if (logGrid) {
	        double r = Math.exp(x[i]); // this is the distance to the origin
	        a[i] = 2*m*r*r*(energy - V[i]); // see text in the beginning of the document
	    } else {
	        a[i] = 2*m*(energy - V[i]); // standard Schroedinger equation
	    }
	    f[i] = 1 + a[i]*dx*dx/12.0;
	    if (icl < 0 && i >= 1 && a[i]*a[i-1] < 0) {
	        icl = i;
	    }
	}
	return icl;
    }

    /**
     * Returns the wave function after implementing the Richardson extrapolation
     * to improve precision.
     * @param y (required) Wave function found with step size dx.
     * @param yk (required) Wave function found with step size dx/k.
     * @param n (required) Precision order (ie: y is O(dx^n)).
     * @param k (required) Factor of reduction of step size.
     * @return Value of y with precision O(dx^(n+1)).
     */
    double Richardson(double y, double yk, double n, double k) {
	return (Math.pow(k, n)*y - yk)/(Math.pow(k, n) - 1);
    }

    /**
     * Calculates the solution of the second order equation assuming that the
     * solution goes to zero in the first point of the grid.
     * @param energy (required) Energy being probed.
     * @param y (required) Solution of the equation to be returned by reference.
     * @param f (required) Auxiliary function $f(x)$ from getAuxVectors().
     */
    void getSolutionOutward(double energy, double [] y, double [] f) {
        // initial condition: y -> 0 as x -> +/- infinity
	// y''(x) + a(x) y(x) = 0
	double dx = Math.abs(x[1] - x[0]);
	if (logGrid) { // should be more or less the same as above, but this is a better approximation for the Coulomb potential
	    y[0] = (Math.exp(x[0]*0.5));
	    y[1] = (Math.exp(x[1]*0.5));
	} else {
	    y[0] = eps*1e-12;
	    y[1] = y[0] + 0 - dx*dx*0.5*2*(energy - V[0])*y[0];
	}
        // recursively solve diff. equation using Numerov method
	for (int i = 1; i < N - 1; ++i) {
	    if (f[i+1] != 0) {
	        y[i+1] = ((12 - f[i]*10)*y[i] - f[i-1]*y[i-1])/f[i+1];
	    } else { // avoid division by zero
	        y[i+1] = y[i];
	    }
        }

        // solve it again with twice the step size
        double [] yk = new double[N/2];
	if (logGrid) {
	    yk[0] = (Math.exp(x[0]*0.5));
	    yk[1] = (Math.exp(x[2]*0.5));
	} else {
	    yk[0] = eps*1e-12;
	    yk[1] = yk[0] + 0 - dx*dx*0.5*2*(energy - V[0])*yk[0];
	}
	for (int i = 1; i < N/2 - 1; ++i) {
	    yk[i+1] = ((12 - f[2*i]*10)*yk[i] - f[2*i-2]*yk[i-1])/f[2*i+2];
        }
        // use Richardson extrapolation to improve precision
	for (int i = 0; i < N; ++i) {
	    //y[i] = Richardson(y[i], yk[i/2], 5, 2);
	}
    }

    /**
     * Solves the differential equation assuming the solution goes to zero in the
     * last point of the grid.
     * @param energy (required) Value of energy being probed.
     * @param y (required) Solution to be returned.
     * @param f (required) Auxiliary vector $f(x)$ returned from getAuxVectors().
     */
    void getSolutionInward(double energy, double [] y, double [] f) {
        // initial condition: psi -> 0 as x -> +/- infinity
	// y''(x) + a(x) y(x) = 0
	double dx = Math.abs(x[N-1] - x[N-2]);
	if (logGrid) { // should be the same as above, but better suited if using the Coulomb potential
	    y[N-1] = Math.exp(-Math.sqrt(-2*energy)*Math.exp(x[N-1]));
	    y[N-2] = Math.exp(-Math.sqrt(-2*energy)*Math.exp(x[N-2]));
	} else {
	    y[N-1] = eps*1e-12;
	    y[N-2] = y[N-1] - 0 - dx*dx*0.5*2*(energy - V[N-1])*y[N-1];
	}
        // Use Numerov method to get solution iteratively.
	for (int i = N-2; i > 0; --i) {
	    if (f[i-1] != 0) {
	        y[i-1] = ((12 - f[i]*10)*y[i] - f[i+1]*y[i+1])/f[i-1];
            } else {
	        y[i-1] = y[i];
	    }
        }

        // solve again with twice the step size
        double [] yk = new double[N/2];
	if (logGrid) {
	    y[N/2-1] = Math.exp(-Math.sqrt(-2*energy)*Math.exp(x[N-1]));
	    y[N/2-2] = Math.exp(-Math.sqrt(-2*energy)*Math.exp(x[N-3]));
	} else {
	    yk[N/2-1] = eps*1e-12;
	    yk[N/2-2] = yk[N/2-1] - 0 - dx*dx*0.5*2*(energy - V[N/2-1])*yk[N/2-1];
	}
	for (int i = N/2-2; i > 0; --i) {
	    yk[i-1] = ((12 - f[2*i]*10)*yk[i] - f[2*i+2]*yk[i+1])/f[2*i-2];
        }

        // Using Richardson extrapolation to improve precision
	for (int i = 0; i < N; ++i) {
	    //y[i] = Richardson(y[i], yk[i/2], 5, 2);
	}
    }

    /**
     * Use solution from getSolutionOutward() until the point on which the energy = potential
     * and from there one, use the solution from getSolutionInward().
     * Normalises the two solutions so that they match in the transition point.
     * @param y_out (required) Solution assuming with initial condition zero at the first points.
     * @param y_in (required) Solution assuming with initial condition zero at the last points.
     * @param icl (required) Point at which to match the two solutions.
     * @param y (required) Matched solution.
     */
    void matchInOut(double [] y_out, double [] y_in, int icl, double [] y) {
        double rat = 1;
	if (icl >= 0 && y_in[icl] != 0 && y_out[icl] != 0) {
            // get ratio of the solutions at icl
	    rat = y_out[icl]/y_in[icl];
	    for (int i = 0; i < icl; ++i) {
	        y[i] = y_out[i]; // use solution with initial conditions set at 0 in the beginning
	    }
	    for (int i = icl; i < N; ++i) {
	        y[i] = rat*y_in[i]; // use solution with initial conditions set at infinity in the end, after some renormalisation
	    }
	} else {
	    for (int i = 0; i < N; ++i) {
	        y[i] = y_out[i];
	    }
	}
    }

    /**
     * Returns a value, which the Numerov predicts to be zero if the solution satisfies the equation.
     * The Numerov method predicts <code>(12 - 10*f[i])*y[i] - f[i-1]*y[i-1] - f[i+1]*y[i+1] = 0</code>
     * at each grid point.
     * @param icl (required) Point at which to get this constraint.
     * @param f (required) Auxiliary vector.
     * @param y (required) Solution vector.
     * @return Value of the constraint.
     */
    double getConstraint(int icl, double [] f, double [] y) {
        return (12 - 10*f[icl])*y[icl] - f[icl-1]*y[icl-1] - f[icl+1]*y[icl+1];
    }

    /**
     * Solve the differential equation with a test value of the energy.
     * @param energy (required) Energy probed.
     * @param y (required) Solution to be returned.
     * @param f (required) Auxiliary vector to be returned.
     * @param y_out (required) Solution to be returned imposing initial conditions at the origin.
     * @param y_in (required) Solution to be returned imposing initial conditions at infinity.
     * @return Index of the grid vector on which the energy = potential.
     */
    int solve(double energy, double [] y, double [] f, double [] y_out, double [] y_in) {
	int icl = -1;
	double [] a = new double[N];
        // get auxiliary vectors
	icl = getAuxVectors(energy, a, f);
	double [] y_outward = y_out; //new double[N];
        // solve the equation with the initial conditions at the origin
	getSolutionOutward(energy, y_outward, f);
	double [] y_inward = y_in; //new double[N];
        // solve the equation with the initial conditions at infinity
	getSolutionInward(energy, y_inward, f);
        // match both solutions at the (first) point on which the energy == potential
        // Could be any other point, but this generates maximum discontinuity in the first derivative
        // if the two solutions are not compatible (ie: we are at the wrong energy).
        // So we can check the first derivative at the merging point to check if we
        // are in the wrong energy
	matchInOut(y_outward, y_inward, icl, y);
	return icl;
    }

    /**
     * Get value of the wave function properly normalised.
     * For the radial equation, it converts the $\xi(z)$ solution to $\psi(x)$ as well.
     * @param y (required) Differential equation solution
     * @param psi (required) Normalised wave function to be returned.
     */
    void toPsi(double [] y, double [] psi) {
        double norm = 0;
        // first calculate normalisation
        for (int i = 0; i < N; ++i) {
	    double dx = 0;
            if (i == N-1) {
                dx = Math.abs(x[i] - x[i-1]);
            } else {
                dx = Math.abs(x[i+1] - x[i]);
            }
            // in the radial equation, we actually calculate $\xi(z)$
            // and $\xi(z) = x^{-1/2} \psi(x)$. Furthermore,
            // $x = \exp(z)$ and the grid stored in the vector x[] is actually
            // to be interpreted as the variable $z$.
            // Normalise so that int |psi*r|^2 dr = 1
	    if (logGrid) {
	        double r = Math.exp(x[i]);
	        norm += Math.pow(r*y[i]/Math.sqrt(r), 2)*(r*dx);
            } else {
                // if not in the radial equation
                // Normalize so that int |psi|^2 dx = 1
	        norm += y[i]*y[i]*dx;
	    }
	}
        // we want to present $\psi$, so get the square root of the normalisation
	norm = Math.sqrt(norm);
        for (int i = 0; i < N; ++i) {
	    if (logGrid) {
	        double r = Math.exp(x[i]);
	        psi[i] = y[i]/Math.sqrt(r)/norm;
	    } else {
	        psi[i] = y[i]/norm;
	    }
	}
    }

    /**
     * Does the actual calculation to get energy eigenvalue and wave function.
     */
    public void recalculate() {
        converged = false;
        Emin = 1e10;
	Emax = -1e10;
        for (int i = 0; i < N; ++i) {
	    if (V[i] < Emin) {
	        Emin = V[i];
	    }
	    if (V[i] > Emax) {
	        Emax = V[i];
	    }
	}
	if (logGrid) {
            // this prevents the energy from going over the asymptote
	    Emax = 1e-5; // in the radial equation, the potential should have an asymptote at zero
	}
	E = (Emin+Emax)*0.5; // guess some initial energy value
        // iterate until it converges
	for (int k = 0; k < m_iter; ++k) {
            // creates some vectors
	    double dE_probe = E*0.001;
	    double [] y1 = new double[N];
	    double [] y1_in = new double[N];
	    double [] y1_out = new double[N];
	    double [] f1 = new double[N];
	    double [] y2_in = new double[N];
	    double [] y2_out = new double[N];
	    double [] y2 = new double[N];
	    double [] f2 = new double[N];
            // try to get a solution for the current energy E
	    int icl1 = solve(E, y1, f1, y1_out, y1_in);
            // and try now to get solution for the energy E shifted by a small amount dE_probe
	    int icl2 = solve(E+dE_probe, y2, f2, y2_out, y2_in);
	    double F1 = 0;
	    double F2 = 0;
	    double dE = 0;
	    if (icl1 > 0) {
                // now get the value of the constraints given by the Numerov method
                // for each solution if energies E and E+dE_probe
                // In the correct energy value, these constraints should be zero
                // as both solutions will be the same, as they respect the initial condition
                // at zero and infinity simultaneously
                F1 = getConstraint(icl1, f1, y1);
                F2 = getConstraint(icl1, f1, y2);
		if (F1 != F2) {
                    // the constraint F is seen as an energy dependent function, F(E)
                    // we can expand it in Taylor series:
                    // F(E+dE) = F(E) + dE F'(E) + O(dE^2)
                    // If we want F(E+dE) = 0:
                    // dE = -F(E)/F'(E)
                    // with F'(E) = (F2-F1)/dE_probe, this becomes:
                    // dE = -F1*dE_probe/(F2-F1)
                    // Here we add a factor alpha so that the user can probe smaller steps of energy
                    // This is helpful if the disregarded O(dE^2) terms are significant
	            dE = -F1*m_alpha*dE_probe/(F2 - F1);
		} else {
                    // avoid division by zero
                    // if changing the energy changes nothing, the probed step is the best guess
		    dE = dE_probe;
		}
            } else {
                // no point at which energy == potential
                // energy is either too big or too small
                // Should not happen since we limit the energy step by the potential minimum and maxima values
                // However, perhaps the energy step was too big and some of the checks failed
	        dE = -0.1;
	    }
            // normalise solution
            // Also translates $\xi(z)$ solution to $\psi(x)$ solution if using the radial equation
	    toPsi(y1, psi);
	    toPsi(y1_out, psi_out);
	    toPsi(y1_in, psi_in);
            // get number of zeroes for the solution
	    int nodes = 0;
            for (int i = 3; i < N-3; ++i) {
	        if (psi[i]*psi[i+1] < 0) {
		    nodes++;
		}
	    }
            // for debugging:
            int i_max = -1;
            double psi_max = 0;
            for (int i = 0; i < N-1; ++i) {
              if (Math.abs(psi[i]) > psi_max) {
                i_max = i;
                psi_max = Math.abs(psi[i]);
              }
            }
	    //System.out.format("Iteration %d: E = %f, nodes = %d, icl = %d, Emin = %f, Emax = %f, dE = %f\n", k, E, nodes, icl1, Emin, Emax, dE);
            // now check number of zeroes
            // If we have the wrong number of zeroes as the desired one, we
            // are converging to the wrong eigenstate
            // Stop this, by limiting the minimum and maxima energy values
	    if (nodes != m_n) {
		if (nodes > m_n) { // we have too many zeroes: energy is too big
                    // do not allow the energy to be bigger than the current energy
		    Emax = E;
                    // move energy to (Emin + Emax)*0.5
		    dE = (E + Emin)*0.5 - E;
		} else if (nodes < m_n) { // we have too few zeroes: energy is too small
                    // do not allow the energy to be smaller than this
		    Emin = E;
                    // move energy to (Emin + Emax)*0.5
		    dE = (E + Emax)*0.5 - E;
		}
		E += dE; // shift it
	    } else {
                // we have the correct number of zeroes, so we are close to the correct eigenstate
                // but if there is no crossing point, the energy is too big or too small
                // should not happen, but some checks might fail
	        if (icl1 < 0) {
		    dE = -0.1;
		}
		double E_old = E;
                // we are now about to shift the energy
                // but before that, check if we are about to cross the allowed boundaries
                // if the energy is below the minimum or above the maximum,
                // we could get close to a different eigenstate and get a different number of zeroes
		if (logGrid && E + dE < Emin) { // will reduce energy too much! stop this!
		    dE = (Emax + Emin)*0.5 - E; // get some average energy
		} else if (logGrid && E + dE > Emax) { // will increase energy too much! stop this!
		    dE = (Emax + Emin)*0.5 - E; // get some average energy
		}
		E += dE; // step the energy

                // if we went into a positive energy direction, establish the
                // last energy as a minimum, otherwise establish the last energy as a maximum
	        if (icl1 < 0) { // should not happen
		    Emin = E_old;
                } else if (dE > 0) {
		    Emin = E_old;
		} else if (dE < 0) {
		    Emax = E_old;
		}
	    }
            // now go to the next iteration to check again if the solution with initial
            // conditions at -infinity and +infinity lead to a continuous solution
            // But first check if the allowed energy interval is still significant
            // if it is smaller than eps, we have converged: stop the iteration.
	    if (Math.abs(Emin - Emax) < eps && nodes == m_n) {
	        converged = true;
	        break;
	    }
	}
    }

    /**
     * Get the dataset with plots to be shown.
     * @param xv (required) X axis values.
     * @param yv (required) Y axis values.
     */
    public void getDatasetPsi(double [] xv, double [] yv) {
        // Create plots
        // Remember that in the radial equation, the X axis will show
        // exp(grid).
	for (int i = 0; i < N; ++i) {
	    if (logGrid) {
	        xv[i] = Math.exp(x[i]);
                yv[i] = psi[i];
	    } else {
	        xv[i] = x[i];
                yv[i] = psi[i];
	    }
	}
    }

    /**
     * Get energy value and the potential energy curve.
     * @param xv (required) X axis values.
     * @param yv (required) Y axis values.
     * @return Energy value.
     */
    public double getDatasetEnergy(double [] xv, double [] yv) {
        // Create plots
        // Remember that in the radial equation, the X axis is actually
        // exp(grid)
	for (int i = 0; i < N; ++i) {
	    if (logGrid) {
	        xv[i] = Math.exp(x[i]);
                yv[i] = V[i];
	    } else {
	        xv[i] = x[i];
                yv[i] = V[i];
	    }
	}
        return E;
   }
}

