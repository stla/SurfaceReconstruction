#ifndef _SURFACERECONSTRUCTIONHEADER_
#define _SURFACERECONSTRUCTIONHEADER_
#endif

#define CGAL_EIGEN3_ENABLED 1

// #include <CGAL/assertions.h>
// #undef CGAL_error
// #define CGAL_error
// #undef CGAL_error_msg
// #define CGAL_error_msg(msg)

// #include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>

#include <CGAL/Kernel/global_functions.h>

#include <CGAL/Vector_3.h>

//#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/squared_distance_2.h>

//#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/utility.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/tuple.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/poisson_surface_reconstruction.h>

#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>

//#include <CGAL/Projection_traits_xy_3.h>

//#include <CGAL/Triangulation_face_base_with_info_2.h>

//#include <CGAL/Cartesian.h>


#include <RcppEigen.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point3;

typedef CGAL::Advancing_front_surface_reconstruction<> AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3 AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;
typedef K::Vector_3 Vector3;

typedef CGAL::Simple_cartesian<double> KSC;

typedef std::pair<Point3, Vector3> P3wn;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

typedef Eigen::
    Matrix<unsigned, Eigen::Dynamic, 3, Eigen::RowMajor | Eigen::AutoAlign>
        Imatrix;
typedef Eigen::Matrix<unsigned, 1, 3, Eigen::RowMajor | Eigen::AutoAlign>
    Ivector;

// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
