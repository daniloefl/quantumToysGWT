package QuantumToys.client;
import com.googlecode.gchart.client.GChartCanvasLite;
import com.googlecode.gchart.client.GChartCanvasFactory;

public class GWTCanvasBasedCanvasFactory 
   implements GChartCanvasFactory {

   public GChartCanvasLite create() {
      GChartCanvasLite result = new GWTCanvasBasedCanvasLite();
      return result;
   }
   
}

