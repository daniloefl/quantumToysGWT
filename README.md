# quantumToysGWT

This solves the Schroedinger equation in 1D using the Numerov method.
It has two interfaces: one is built using the Google Web Toolkit and compiles
Java into JavaScript to be shown in a Web browser. The other generates a Java
Web Start application, which can be loaded if the user has Java installed in his/her computer.

The Java Web Start application can be directly accessed from the documentation.
You can see this in action by opening this link in a browser:
<https://daniloefl.github.io/quantumToysGWT/QuantumToys/SchroedingerSolver.jnlp>

# Installing the Google Web Toolkit implementation

To compile, please unzip the Google Web Toolkit, update this line of build.xml to
point to your installation of the Google Web Toolkit:

```
  <property name="gwt.sdk" location="/home/daniloefl/workspace/gwt-2.8.0" />
```

and compile it as follows:

```
ant gwtc
```

You can then open war/SchroedingerWebApp.html in a browser.

# Installing the Java Web Start implementation

This assumes that your Java Web Start libraries are available in /usr/share/icedtea-web/netx.jar.
If that is not the case, please change the following line in build.xml:

```
    <pathelement location="/usr/share/icedtea-web/netx.jar"/>
```

To compile, simply type:

```
ant jws
```

The compiled result will be in the directory `jws_bin`. To run it, one can either open the file
`jws_bin/SchroedingerSolver.jnlp` in a browser or type:

```
javaws jws_bin/SchroedingerSolver.jnlp
```

# More information

Please read the Javadoc documentation for of the class SchroedingerCalculator
to understand how the calculation is performed. This documentation is
in <https://daniloefl.github.io/quantumToysGWT/QuantumToys/shared/SchroedingerCalculator.html>.

For more information, please contact the author at:
Danilo Ferreira de Lima <daniloefl@gmail.com>

