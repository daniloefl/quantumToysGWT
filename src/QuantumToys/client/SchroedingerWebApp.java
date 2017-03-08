package QuantumToys.client;

import com.googlecode.gchart.client.GChart;

import QuantumToys.shared.SchroedingerCalculator;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.event.dom.client.KeyCodes;
import com.google.gwt.event.dom.client.KeyUpEvent;
import com.google.gwt.event.dom.client.KeyUpHandler;
import com.google.gwt.user.client.rpc.AsyncCallback;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.CheckBox;
import com.google.gwt.user.client.ui.VerticalPanel;

/**
 * This is the entry point for the web application.
 * It solves the Schroedinger equation and shows results on the web.
 * This is only the interface, while the actual equation solving code
 * is in QuantumToys.shared.SchroedingerCalculator.
 * @author Danilo Ferreira de Lima <daniloefl@gmail.com>
 * @see QuantumToys.shared.SchroedingerCalculator
 */
public class SchroedingerWebApp implements EntryPoint {

  /**
   * Labels and controls to be shown in the web page.
   */

  CheckBox logGrid;
  Label labelConverge;
  Label labelEnergy;

  Label labelXmin;
  Label labelXmax;

  TextBox xmin;
  TextBox xmax;

  Label labelNpoints;
  TextBox Npoints;

  Button go;

  CheckBox autoRange;

  Label labelN;
  TextBox n;

  Label labelPot;
  TextBox pot;

  Label labelIter;
  TextBox iter;

  Label labelAlpha;
  TextBox alpha;


  /**
   * Instance of the calculator class to solve equation
   */
  SchroedingerCalculator m_calc;

  /**
   * Instance of the class used to draw results.
   */
  PsiChart m_psiChart;

  /**
   * Updates the plot and the contents of all controls if needed.
   */
  public void updatePlot() {
    n.setText(Integer.toString(m_calc.m_n));
    pot.setText(m_calc.exp_V_str);
    iter.setText(Integer.toString(m_calc.m_iter));
    alpha.setText(Double.toString(m_calc.m_alpha));
    Npoints.setText(Integer.toString(m_calc.N));
    xmax.setText(Double.toString(m_calc.xmax));
    xmin.setText(Double.toString(m_calc.xmin));
    m_psiChart.getPsiData(m_calc);
    m_psiChart.getEnergyData(m_calc);
    m_psiChart.update();
    String text = "Energy = ";
    text = text + Double.toString(m_calc.E);
    text = text + " a.u. = ";
    text = text + Double.toString(m_calc.E*m_calc.eV);
    text = text + " eV";
    labelEnergy.setText(text);
  }

  /**
   * This is the entry point method.
   */
  public void onModuleLoad() {
    logGrid = new CheckBox("Radial equation");
    RootPanel.get("logGrid").add(logGrid);

    labelConverge = new Label();
    labelConverge.setText("Warning: Did not converge!");
    labelConverge.setVisible(false);
    RootPanel.get("labelConverge").add(labelConverge);

    labelEnergy = new Label();
    labelEnergy.setText("Energy = ");
    RootPanel.get("labelEnergy").add(labelEnergy);

    labelXmin = new Label();
    labelXmin.setText("min(x)");
    RootPanel.get("labelXmin").add(labelXmin);

    xmin = new TextBox();
    xmin.setWidth("1em");
    RootPanel.get("xmin").add(xmin);

    labelXmax = new Label();
    labelXmax.setText("max(x)");
    RootPanel.get("labelXmax").add(labelXmax);

    xmax = new TextBox();
    xmax.setWidth("1em");
    RootPanel.get("xmax").add(xmax);

    labelNpoints = new Label();
    labelNpoints.setText("# points");
    RootPanel.get("labelNpoints").add(labelNpoints);

    Npoints = new TextBox();
    Npoints.setWidth("4em");
    RootPanel.get("Npoints").add(Npoints);

    go = new Button("  Update  ");
    RootPanel.get("go").add(go);

    autoRange = new CheckBox("Auto range");
    RootPanel.get("autoRange").add(autoRange);

    labelN = new Label();
    labelN.setText("# nodes");
    RootPanel.get("labelN").add(labelN);

    n = new TextBox();
    n.setWidth("1em");
    RootPanel.get("n").add(n);

    labelPot = new Label();
    labelPot.setText("Potential");
    RootPanel.get("labelPot").add(labelPot);

    pot = new TextBox();
    pot.setWidth("10em");
    RootPanel.get("pot").add(pot);

    // Focus the cursor on the field when the app loads
    pot.setFocus(true);
    pot.selectAll();


    labelIter = new Label();
    labelIter.setText("# iterations");
    RootPanel.get("labelIter").add(labelIter);

    iter = new TextBox();
    iter.setWidth("4em");
    RootPanel.get("iter").add(iter);

    labelAlpha = new Label();
    labelAlpha.setText("Energy step");
    RootPanel.get("labelAlpha").add(labelAlpha);

    alpha = new TextBox();
    alpha.setWidth("4em");
    RootPanel.get("alpha").add(alpha);

    m_calc = new SchroedingerCalculator(new GWTParser());
    m_calc.recalculate();

    GChart.setCanvasFactory(new GWTCanvasBasedCanvasFactory());

    m_psiChart = new PsiChart();
    RootPanel.get("psiChart").add(m_psiChart);
    

    updatePlot();

    // Create a handler for the auto range check box
    class AutoRangeHandler implements ClickHandler {

      /**
       * Update plot.
       */
      private void runMe(boolean checked) {
        // First, we validate the input.
        if (checked) {
          m_calc.autoRange();
          xmin.setText(Double.toString(m_calc.xmin));
          xmax.setText(Double.toString(m_calc.xmax));
          Npoints.setText(Integer.toString(m_calc.N));
          xmin.setEnabled(false);
          xmax.setEnabled(false);
          Npoints.setEnabled(false);
        } else {
          xmin.setEnabled(true);
          xmax.setEnabled(true);
          Npoints.setEnabled(true);
        }
        updatePlot();
      }

      /**
       * Fired when the user clicks on the auto range.
       */
      public void onClick(ClickEvent event) {
          boolean checked = ((CheckBox) event.getSource()).getValue();
          runMe(checked);
      }
    }
    AutoRangeHandler autoRangeHandler = new AutoRangeHandler();
    autoRange.addClickHandler(autoRangeHandler);

    // Create a handler for the radial eq. check box
    class RadialHandler implements ClickHandler {

      /**
       * Update plot.
       */
      private void runMe(boolean checked) {
        // First, we validate the input.
        if (checked) {
          m_calc.setLogGrid(true);
          xmin.setText(Double.toString(m_calc.xmin));
          xmax.setText(Double.toString(m_calc.xmax));
          alpha.setText(Double.toString(m_calc.m_alpha));
          pot.setText("-1/x+0.25/(2*x*x)"); // set potential to the Hydrogen potential
        } else {
          m_calc.setLogGrid(false);
          xmin.setText(Double.toString(m_calc.xmin));
          xmax.setText(Double.toString(m_calc.xmax));
          alpha.setText(Double.toString(m_calc.m_alpha));
          pot.setText(m_calc.def_pot); // set potential to the harmonic oscillator
        }
      }

      /**
       * Fired when the user clicks on the radial eq.
       */
      public void onClick(ClickEvent event) {
          boolean checked = ((CheckBox) event.getSource()).getValue();
          runMe(checked);
      }
    }
    RadialHandler radialHandler = new RadialHandler();
    logGrid.addClickHandler(radialHandler);

    // Create a handler for the button
    class UpdateHandler implements ClickHandler, KeyUpHandler {

      /**
       * Fired when the user clicks on the button.
       */
      public void onClick(ClickEvent event) {
          runMe();
      }

      /**
       * Fired when the user types in the nameField.
       */
      public void onKeyUp(KeyUpEvent event) {
        if (event.getNativeKeyCode() == KeyCodes.KEY_ENTER) {
          runMe();
        }
      }

      /**
       * Update plot.
       */
      private void runMe() {
        // First, we validate the input.
        m_calc.setN(Integer.parseInt(n.getText()));
        m_calc.setX(Double.parseDouble(xmin.getText()), Double.parseDouble(xmax.getText()), Integer.parseInt(Npoints.getText()));
        m_calc.setIterStep(Integer.parseInt(iter.getText()), Double.parseDouble(alpha.getText()));
        m_calc.recalculatePotential(pot.getText());
        m_calc.recalculate();
        updatePlot();
      }
    }

    // Add a handler to send the name to the server
    UpdateHandler handler = new UpdateHandler();
    go.addClickHandler(handler);
    pot.addKeyUpHandler(handler);
  }
}

