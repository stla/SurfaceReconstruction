#ifndef _SURFACERECONSTRUCTIONHEADER_
#include "SurfaceReconstruction.h"
#endif

// [[Rcpp::export]]
Rcpp::NumericMatrix jet_normals_cpp(const Rcpp::NumericMatrix pts,
                                    const unsigned nb_neighbors) {
  const size_t npoints = pts.ncol();
  std::vector<P3wn> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i); 
    points[i] = std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)),
                               Vector3(0.0, 0.0, 0.0));
  }

  CGAL::jet_estimate_normals<Concurrency_tag>(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  CGAL::mst_orient_normals(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  Rcpp::NumericMatrix normals(3, npoints);
  for(size_t i = 0; i < npoints; i++) {
    Rcpp::NumericVector normal_i(3);
    const Vector3 normal = points[i].second;
    normal_i(0) = normal.x();
    normal_i(1) = normal.y();
    normal_i(2) = normal.z();
    normals(Rcpp::_, i) = normal_i;
  }

  return normals;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pca_normals_cpp(Rcpp::NumericMatrix pts,
                                    unsigned nb_neighbors) {
  const size_t npoints = pts.ncol();
  std::vector<P3wn> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i); 
    points[i] = std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)),
                               Vector3(0.0, 0.0, 0.0));
  }

  CGAL::pca_estimate_normals<Concurrency_tag>(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  CGAL::mst_orient_normals(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  Rcpp::NumericMatrix normals(3, npoints);
  for(size_t i = 0; i < npoints; i++) {
    Rcpp::NumericVector normal_i(3);
    const Vector3 normal = points[i].second;
    normal_i(0) = normal.x();
    normal_i(1) = normal.y();
    normal_i(2) = normal.z();
    normals(Rcpp::_, i) = normal_i;
  }

  return normals;
}
