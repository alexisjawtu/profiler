#include "Profile.h"
#include "SetParameter.h"


using namespace std;


int main(int argc, char *argv[]) {

    /*
    print("\nProfiler v1.2.0 August 2023." 
          "\vReminder: The boundary points in the input are assumed\n"
          "to be the 2D points as stored in "
          "vector<size_t> Cylinder::bdr_pointlist.");
    */


    // begin{ all this -> SetParameter
    bool debug_flag {false};
    // number of horizontal mesh layers
    int tetrahedral_layers {1};
    /**
     * TODO the task here is to eliminate the absurd definition of
     * layers as being 2 (2D layers) when there is only one 3D layer
     *
     * Eliminate variable layer and
     * keep only variable tetrahedral_layers.
    **/ 
    int layer = tetrahedral_layers + 1;
    // end{ all this -> SetParameter
    
    SetParameter parameters = SetParameter(argc, argv);

    string input_folder {parameters.cylinder_folder};
    ifstream boundary_data (input_folder + filenames::bdr_vertices_2D);
    ifstream nodes_data    (input_folder + filenames::cylinder_verts);

    if (!boundary_data || !nodes_data)
    {
        print("\vSomething is wrong with the input files.\n");
        exit(1);
    }
    /** TODO 
     * 1- if we keep this design of class Profile,
     *    then the constructor simply takes an InputCollection, a.k.a. SetParameter,
     *   as in
     *  
     *             Profile profile = Profile(parameters);
     *             bla bla bla
    **/ 

    Profile profile = Profile(
                        layer,
                        parameters.user_thickness_of_inner_wafer,
                        debug_flag,
                        parameters.output_folder,
                        parameters.profile_parameters
                      );

    int current_node {0};
    int comma;                                  // Stores the places of the commas.
    double coordx;                              // 2D coordinates of boundary points.
    double coordy;                              
    int row;                                    // Stores the index of a tetrahedron.
    int number_of_bdr_points {0};               // Stores the number_of_bdr_points of a 2D point.
    vector<int> element;                        // Stores a tetrahedron.

    // TODO de-hardcode this 100
    // int length_of_one_line_in_the_input = 100;  // Three floats of 22 characters long. No more.
    char* raw_point = new char[100];              
    boundary_data.getline(raw_point, 100);
    string point (raw_point);

    while (!boundary_data.eof())
    {
        comma  = point.find_first_of(",");
        stringstream(point.substr(0, comma)) >> coordx;
        point  = point.substr(comma + 1);

        comma  = point.find_first_of(",");
        stringstream(point.substr(0, comma)) >> coordy;
        point  = point.substr(comma + 1);

        profile.bdr_pointlist.push_back(Point(coordx, coordy));
        profile.pointlist.push_back(Point(coordx, coordy));

        boundary_data.getline(raw_point, 100);
        point = string(raw_point);

        number_of_bdr_points ++;
    }

    delete[] raw_point;
    boundary_data.close();

    while (false)
    {
        /**
         *  TODOs: all this cycle body to receive the 3D nodes missing.
         *
         *  Look for the indices in the boundary and erase it
         *  
         *  TASK Research how to use <valarray> to traverse and search a double
         *  TASK important transform the code in this 'main' into methods of Profile.
                 leave a minimal 'int main () {}'
         *  TASK replace the prints and couts with template print<T>
        **/
    }

    // profile.find_global_coordinates_for_boundary();
    profile.make_3D_points();
    profile.build_profile_mesh(number_of_bdr_points);
    profile.stream_elements_out();
    profile.stream_nodes_out();
    
    return 0;
}
