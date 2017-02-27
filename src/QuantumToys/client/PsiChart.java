package QuantumToys.client;

import QuantumToys.client.SchroedingerCalculator;

import com.googlecode.gchart.client.GChart;

import java.util.*;

public class PsiChart extends GChart {

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

  public void getEnergyData(SchroedingerCalculator calc) {
    double [] x = new double[calc.N];
    double [] potential = new double[calc.N];
    double E = calc.getDatasetEnergy(x, potential);
    getCurve(2).clearPoints();
    getCurve(3).clearPoints();
    for (int k = 0; k < calc.N; ++k) {
      getCurve(2).addPoint(x[k], potential[k]);
      getCurve(3).addPoint(x[k], E);
    }
  }

}

