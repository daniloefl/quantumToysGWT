package QuantumToys.jws;

import QuantumToys.shared.Parser;
import net.objecthunter.exp4j.*;

public class JWSParser extends Parser {

    Expression exp_V;

    public JWSParser() { }

    /**
     * Evaluates mathematical expression.
     * @param s (required) Expression to evaluate.
     * @param sx (required) Variable.
     * @param vx (required) Value of x.
     * @return Expression value for x = vx.
     */
    public double evaluate(String s, String sx, Double vx) {
        exp_V = new ExpressionBuilder(s)
	                                .variables(sx)
                                        .build();
        return exp_V.setVariable(sx, vx).evaluate();
    }

}

