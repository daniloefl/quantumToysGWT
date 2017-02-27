# quantumToysGWT

This is a Web Application built using Google
Web Toolkit 2.8.0 to transform Java into JavaScript.

It solves the Schroedinger equation in 1D using the Numerov method
and shows the result on the web page.

To compile, please unzip the Google Web Toolkit, update this line of build.xml to
point to your installation of the Google Web Toolkit:

```
  <property name="gwt.sdk" location="/home/daniloefl/workspace/gwt-2.8.0" />
```

and compile it as follows:

```
ant build
```

You can then open war/SchroedingerWebApp.html in a browser.

For more information, please contact the author at:
Danilo Ferreira de Lima <daniloefl@gmail.com>

