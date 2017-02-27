package QuantumToys.client;

import QuantumToys.shared.SchroedingerCalculator;

import com.googlecode.gchart.client.GChart;

import java.util.*;

/**
 * Uses GChart to draw the results in the web page.
 * @author Danilo Ferreira de Lima <daniloefl@gmail.com>
 */
public class PsiChart extends GChart {

  // Constructor
  PsiChart() {
    setChartTitle("");
    setChartSize(800, 600);
    setLegendYShift(-150);

    addCurve();
    getCurve().setLegendLabel("<i>&Psi;(x)</i>");
    getCurve().setYAxis(Y_AXIS);
    getCurve().getSymbol().setSymbolType(SymbolType.BOX_CENTER); 
    getCurve().getSymbol().setWidth(5);
    getCurve().getSymbol().setHeight(5);
    getCurve().getSymbol().setBorderWidth(0);
    getCurve().getSymbol().setBackgroundColor("navy");
    getCurve().getSymbol().setFillThickness(2);
    getCurve().getSymbol().setFillSpacing(5);

    addCurve();
    getCurve().setLegendLabel("<i>|&Psi;(x)|<sup>2</sup></i>");
    getCurve().setYAxis(Y_AXIS);
    getCurve().getSymbol().setSymbolType(SymbolType.BOX_CENTER); 
    getCurve().getSymbol().setWidth(5);
    getCurve().getSymbol().setHeight(5);
    getCurve().getSymbol().setBorderWidth(0);
    getCurve().getSymbol().setBackgroundColor("red");
    getCurve().getSymbol().setFillThickness(2);
    getCurve().getSymbol().setFillSpacing(5);

    addCurve();
    getCurve().setLegendLabel("<i>Potential(x)</i>");
    getCurve().setYAxis(Y2_AXIS);
    getCurve().getSymbol().setSymbolType(SymbolType.BOX_CENTER); 
    getCurve().getSymbol().setWidth(5);
    getCurve().getSymbol().setHeight(5);
    getCurve().getSymbol().setBorderWidth(0);
    getCurve().getSymbol().setBackgroundColor("green");
    getCurve().getSymbol().setFillThickness(2);
    getCurve().getSymbol().setFillSpacing(5);

    addCurve();
    getCurve().setLegendLabel("<i>Energy</i>");
    getCurve().setYAxis(Y2_AXIS);
    getCurve().getSymbol().setSymbolType(SymbolType.BOX_CENTER); 
    getCurve().getSymbol().setWidth(5);
    getCurve().getSymbol().setHeight(5);
    getCurve().getSymbol().setBorderWidth(0);
    getCurve().getSymbol().setBackgroundColor("black");
    getCurve().getSymbol().setFillThickness(2);
    getCurve().getSymbol().setFillSpacing(5);

    getXAxis().setAxisLabel("x");
    getXAxis().setAxisLabelThickness(20);
    getXAxis().setTickCount(13);
    getXAxis().setTicksPerLabel(2);
    getXAxis().setHasGridlines(true);

    getYAxis().setAxisLabel("&Psi;(x)");
    getYAxis().setAxisLabelThickness(20);
    getYAxis().setTickCount(13);
    getYAxis().setTicksPerLabel(2);
    getYAxis().setHasGridlines(true);

    getY2Axis().setAxisLabel("Energy<p></p>[Hartree a.u.]");
    getY2Axis().setAxisLabelThickness(20);
    getY2Axis().setTickCount(13);
    getY2Axis().setTicksPerLabel(2);
    getY2Axis().setHasGridlines(true);
  }

  /**
   * Get wave function data and draw it.
   * @param calc (required) Instance of the equation solver.
   */
  public void getPsiData(SchroedingerCalculator calc) {
    double [] x = new double[calc.N];
    double [] psi = new double[calc.N];
    double [] psi2 = new double[calc.N];
    calc.getDatasetPsi(x, psi);
    for (int k = 0; k < calc.N; ++k) {
      psi2[k] = psi[k]*psi[k];
    }
    getCurve(0).clearPoints();
    getCurve(1).clearPoints();
    for (int k = 0; k < calc.N; ++k) {
      getCurve(0).addPoint(x[k], psi[k]);
      getCurve(1).addPoint(x[k], psi2[k]);
    }
  }

  /**
   * Get potential and total energy data and draw it.
   * @param calc (required) Instance of the equation solver.
   */
  public void getEnergyData(SchroedingerCalculator calc) {
    double [] x = new double[calc.N];
    double [] potential = new double[calc.N];
    double E = calc.getDatasetEnergy(x, potential);
    getCurve(2).clearPoints();
    getCurve(3).clearPoints();
    double Emax = E+10*Math.abs(E);
    double Emin = E-2*Math.abs(E);
    for (int k = 0; k < calc.N; ++k) {
      if (potential[k] < Emin) Emin = potential[k];
      //if (potential[k] > Emax) Emax = potential[k];
      if (potential[k] <= Emax)
        getCurve(2).addPoint(x[k], potential[k]);
      getCurve(3).addPoint(x[k], E);
    }

    getY2Axis().setAxisMax(Emax);
    getY2Axis().setAxisMin(Emin);
  }

}

