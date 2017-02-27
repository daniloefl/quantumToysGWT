package QuantumToys.client;

import QuantumToys.shared.Parser;

/**
 * Class that implements math parsing using JavaScript.
 * @author Danilo Ferreira de Lima <daniloefl@gmail.com>
 */

public class GWTParser extends Parser {


    public GWTParser() { }

    /**
     * Evaluates mathematical expression using JavaScript.
     * @param op (required) String to evaluate with only numbers.
     * @return The double result.
     */
    public native double calculateMath(String op) /*-{ return eval(op); }-*/;

    /**
     * Evaluates mathematical expression.
     * @param s (required) Expression to evaluate.
     * @param sx (required) Variable.
     * @param vx (required) Value of x.
     * @return Expression value for x = vx.
     */
    public double evaluate(String s, String sx, Double vx) {
        // split string in sums or differences
        String dx = Double.toString(vx);
        String op = s.replace(sx, dx);
        return calculateMath(op);
    }

}

