#include "Profile.h"


using namespace std;


int main(int argc, char *argv[]) {

    SetParameter parameters {argc, argv};

    string input_folder {parameters.cylinder_folder};
    ifstream boundary_data {input_folder + filenames::bdr_vertices_2D};
    ifstream nodes_data {input_folder + filenames::cylinder_verts};

    if (!boundary_data || !nodes_data)
    {
        print("\vSomething is wrong with the input files.\n");
        exit(1);
    }

//    auto profile = Profile(parameters);
    Profile profile {parameters};
    
    int current_node {0};
    int comma;                                  // Stores the places of the commas.
    double coordx;                              // 2D coordinates of boundary points.
    double coordy;                              
    int row;                                    // Stores the index of a tetrahedron.
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

        boundary_data.getline(raw_point, 100);
        point = string(raw_point);

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
    profile.build_mesh();
    profile.stream_elements_out();
    profile.stream_nodes_out();
    
    return 0;
}
