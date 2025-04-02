#ifndef _PROFILE_H
#define _PROFILE_H

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
#include "toolbox.h"


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

// We may use this to measure elements and work the matrix condition
typedef tuple<Point3D, Point3D, Point3D, Point3D> Tetrahedron;


class Profile {

    private:

        bool debug;
        bool add_tetrahedra_within;
        string output_dir;

        int number_of_layers;
        // TODO: rename the following two as "base case" "general case"
        // in order to develop the recursion
        int number_of_outer_walls = 0;
        int number_of_bricks_in_wall = 0;
        double height_of_layer;

        // TODO: remove the following three and use the map directly
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
    
        Profile(
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
        vector<Point3D> all_profile_Point3D;
        vector<bool> profile_diagonals;
        vector<int> degenerate_elements;
        vector<int> upper_bdr_points_global_indices;
        map<int, vector<int>> cylinder_verts_by_elems;
        map<int, vector<int>> profile_elems_by_verts;

        ~Profile();

        Point project2D(Point3D p);

        void stream_diagonals_out();
        void stream_elements_out();
        void stream_nodes_out();
        void stream_boundary_nodes_out();
        void print_vector(vector<int>& v);
        void print_vector(vector<std::size_t>& v);
        void print_vector(vector<double>& v);
        void print_vector(vector<bool>& v);
        void find_global_coordinates_for_boundary();
        void make_3D_points();
        void build_profile_mesh(int number_of_points_in_the_input);

};

#endif
