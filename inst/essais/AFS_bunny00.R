library(SurfaceReconstruction, lib.loc = "C:/SL/Rloclib")
library(rgl, lib.loc = "C:/SL/Rloclib")

data(bunny, package = "onion")
mesh <- AFSreconstruction(bunny)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
shade3d(mesh, color = "darkred")

pts <- uniformly::runif_on_sphere(500, 3)
mesh <- AFSreconstruction(pts)
shade3d(mesh, color = "red")

tp <- tessellation::teapot()
mesh <- AFSreconstruction(tp)
shade3d(mesh, color = "red")

data(dummyhead, package = "Rvcg")
dummyhead <- t(dummyhead.mesh$vb[-4L, ])
mesh <- AFSreconstruction(dummyhead)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
shade3d(mesh, color = "darksalmon")

# setwd("C:/SL/MyPackages/RCGAL/inst/essais")
#
# tp <- readLines("teapot.obj", n=3644)
# writeLines(tp, "teapotVertices.txt")
# dat <- read.table("teapotVertices.txt")
# pts <- as.matrix(dat[, c(2, 3, 4)])
# mesh <- AFSreconstruction(pts)
# shade3d(mesh, color = "red")

#tp <- as.matrix(read.table("teapot.txt"))
library(rgl)
afs <- RCGAL:::AFSreconstruction_perimeter_cpp(teapot, 3)
mesh <- tmesh3d(t(teapot), afs$triangles)
mesh <- Rvcg::vcgClean(mesh, c(0, 7))
shade3d(mesh, color = "red")


psr_mesh <- PoissonReconstruction(SolidMobiusStrip)
library(rgl)
wire3d(psr_mesh, color = "black")

library(rgl)
x <- Rvcg::vcgUpdateNormals(SolidMobiusStrip)
teapot_normals <- t(x$normals[-4L, ])
poiss <- RCGAL:::Poisson_reconstruction_cpp(SolidMobiusStrip, teapot_normals)
mesh <- tmesh3d(t(poiss$vertices), t(poiss$facets+1L))
mesh <- addNormals(mesh)#Rvcg::vcgClean(mesh, c(0, 7))
wire3d(mesh, color = "black")


data(bunny, package = "onion")
psr_mesh <- PoissonReconstruction(bunny, spacing = 0.0005)
library(rgl)
shade3d(psr_mesh, color = "green")
