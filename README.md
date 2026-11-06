# Splines
A collection of interpolating splines.
# Intro
This repository contains a few interpolating splines, including S9, my own G2-continuous, low curvature variation spline. Each solver has cyclic and acyclic variants and produces either cubic or quadratic Bezier curves. Implemented splines are these -
* [CubicC2](Splines/SplineCubicC2.h). C2-continuous cubic spline.
* [CubicG2](Splines/SplineCubicG2.h). G2-continuous cubic spline based on [this paper](https://spline.technology/paper1.pdf).
* [Hobby spline](Splines/SplineHobby.h).
* [Adobe κ-Curves](Splines/SplineKappa.h). G2-continuous spline used in Adobe Illustrator.
* [S9](Splines/SplineS9.h). G2-continuous spline based on non-polynomial curves with low curvature variation.

