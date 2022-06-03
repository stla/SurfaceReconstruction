# SurfaceReconstruction

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/SurfaceReconstruction/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/SurfaceReconstruction/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-valgrind](https://github.com/stla/SurfaceReconstruction/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/SurfaceReconstruction/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->

Wrapping the C++ library **CGAL** in R. Convex hull, Delaunay tessellation, surface reconstruction.

## Some examples of Poisson reconstruction

*Toroidal helix:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/ToroidalHelix.png)

*Spider cage:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/SpiderCage.png)

*Solid Möbius strip:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/SolidMobiusStrip.png)

*Hopf torus:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/HopfTorus.png)

The Hopf torus is not very smooth. We can make it a bit smoother by reducing 
the `spacing` parameter of the `PoissonReconstruction` function:

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/HopfTorusMesh_spacing02.png)

Here is a series of three images which show the effect of this `spacing` 
parameter (0.05, 0.02, 0.005):

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/SolidMobiusStrip_spacings.png)

*Dupin cyclide:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/cyclide.png)

*Clifford torus:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/CliffordTorus.gif)

*Orthocircle:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/Orthocircle.png)

*ICN5D's eight-like surface:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/ICN5D_eight.png)

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/StanfordBunny.png)

*Stanford dragon:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/StanfordDragon.png)


## Advanced front surface reconstruction

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/Bunny.png)

*Stanford dragon:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/StanfordDragon.png)

*Dummy head:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/DummyHead.png)

*Skull:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/Skull.png)


## License

This package is provided under the GPL-3 license. If you wish to use CGAL for 
commercial purposes, you must obtain a license from the 
[GeometryFactory](https://geometryfactory.com).


## Blog post

I wrote a [blog post](https://laustep.github.io/stlahblog/posts/SurfaceReconstruction.html) 
devoted to surface reconstruction (using **RCGAL**, the ancestor of **SurfaceReconstruction**).
