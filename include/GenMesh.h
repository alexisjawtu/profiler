#pragma once

#include <string>
#include <list>
#include <vector>
#include <map>
#include <valarray>
#include <algorithm>
#include <fstream>
#include <math.h>

#include "the_map.h"
#include "Point.h"
#include "Point3D.h"
#include "CustomizedSet.h"


using namespace std;


struct Prism {
    int front_up;
    int front_dn;
    int radial_up;
    int radial_dn;
    bool diagonal_orientation;
    vector<int> back;
};

typedef tuple<vector<int>, vector<int>, int> Brick;

// TODO: usar esto para measure()
typedef tuple<Point3D, Point3D, Point3D, Point3D> Tetrahedron;


class GenMesh 
{
    private:

    const int wtype;
    const Point pA, pD;
    const double probe_theta;
    const float minimum_float = 1.17549e-38;

    bool debug;
    bool add_tetrahedra_within;
    string output_dir;

    Point pAD;
    Point Orifla_start, Orifla_mid, Orifla_end, Orifla_2nd_start, Orifla_2nd_mid, Orifla_2nd_end;
    double h_max, h_min, R;
    double innerR;
    double inner_max_h;
    double dist_probe_to_bound_max;
    double dist_probe_to_bound_min;
    double dist_AD;
    double L1, L2;
    double angle_difference;
    double angle_orifla_start, angle_orifla_end, angle_orifla_mid, half_orifla_angle;
    double angle_2nd_orifla_start, angle_2nd_orifla_end, angle_2nd_orifla_mid, half_2nd_orifla_angle;
    double angle_notch_start, angle_notch_end;
    double Owidth, Foffset;
    double arg_pAD;

    int number_of_layers;
    // TODO: rename the following two as "base case" "general case"
    int number_of_outer_walls = 0;
    int number_of_bricks_in_wall = 0;
    double height_of_layer;

    Point direc_AD; //From A to D. Unit
    Point pCenter;

    // TODO: estos 3 pueden sacarse y usar directo el map
    valarray<double> upper_z; 
    valarray<double> lower_z;
    valarray<double> delta_over_plane_xy;

    vector<bool> next_profile_diagonals;
    vector<int> radial_sorted_upper_bdr_indices;
    vector<int> radial_sorted_lower_bdr_indices;
    vector<Point3D> wall_of_Point3D;

    map<int, vector<int>> indices_for_outer_facets;
    map<int, vector<int>> new_indices_for_outer_facets;
    map<int, vector<Point3D>> physical_facets;
    map<int, Prism> prisms_map;
    map<pair<bool, bool>, pair<int, int>> case_translator;
    map<pair<bool, bool>, pair<bool, bool>> outer_diagonals_translator;

    Point get_right_orthogonal(Point vector2D);
    void radial_sort_boundary_points();
    
    void orient_profile_diagonals();
    void add_tetrahedra_within_prisms();
    void construct_front_nodes_of_brick(vector<Point3D>& back_wall, int iteration);
    void split_prism(int prisma_kind, vector<int> up, vector<int> dn);
    void reset_profile_objects();
    double measure(vector<int>& tetra);

    ofstream _test_;

    public:
    
    GenMesh(
        int n_layers,
        double h_layer,
        bool debug_flag,
        string prefix,
        map<string, valarray<double>>& profile_params
    );

    bool measure_flag = false;
    int cuadrilaterals = 0;
    
    vector<Point> pointlist;
    vector<Point> bdr_pointlist;
    the_map sorted_boundary_points;
    vector<Point3D> all_wafer_Point3D;
    vector<bool> profile_diagonals;
    vector<int> degenerate_elements;
    vector<int> upper_bdr_points_global_indices;
    map<int, list<int>> vertices_by_elements;
    map<int, vector<int>> elements_by_vertices;

    ~GenMesh();

    Point project2D(Point3D p);

    vector<Point> get_boundary_points() const;
    void stream_diagonals_out();
    void stream_elements_out();
    void stream_nodes_out();
    void stream_boundary_nodes_out();
    void print_vector(vector<int>& v);
    void print_vector(vector<std::size_t>& v);
    void print_vector(vector<double>& v);
    void print_vector(vector<bool>& v);
    void make_3D_points();
    // void new_make_3D_points(vector<size_t>& in_vector);
    void build_profile_mesh(int number_of_points_in_the_input);

};
