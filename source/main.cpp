#include <string>
#include <iostream>

#include "GenMesh.h"
#include "SetParameter.h"


using namespace std;


/** 
    cppreference.com:
* 
    std::size_t is the unsigned integer type of the result of the sizeof operator.
*   It depends on the architecture.

*/


template <class T> void print(T data)
{
    cout << T << "\n";
}


template <class T> void print_vector(T& v)
{
    for (auto& d: v)
    {
        cout << boolalpha << d << ", ";
    }
    cout << "\n";
}


int main(int argc, char *argv[])
{
    print("\nProfiler v1.2.0 August 2023." 
          "\vReminder: The boundary points in the input are assumed\n"
          "to be the 2D points as stored in "
          "vector<size_t> GenMesh::bdr_pointlist.");

    // begin{ all this -> SetParameter
    string prefix(".");
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

    SetParameter parameters = SetParameter(argc, argv);
    double thickness_of_inner_wafer = parameters.user_thickness_of_inner_wafer;
    double dt = thickness_of_inner_wafer/tetrahedral_layers;
    // end{ all this -> SetParameter

    // TODO genmesh has to be allocated regularly, not in dynamic memory
    /** 
     * TODO if we keep this design of class GenMesh,
     * then the constructor simply takes a SetParameter, as in
     *  
     *             GenMesh gm = GenMesh(parameters);
     *             bla bla bla
    **/ 
    GenMesh* genmesh = new GenMesh(
                        layer,
                        parameters.user_thickness_of_inner_wafer,
                        debug_flag,
                        parameters.output_folder,
                        parameters.profile_parameters
                    );

    ifstream boundary_data ("bdr_nodes.dat");

    if (!boundary_data)
    {
        cout << "Something is wrong with the input file." << endl;
        exit(1);
    }

    int length_of_one_line_in_the_input = 100;  // Three floats of 22 characters long. No more.
    int comma;                                  // Stores the places of the commas.
    double coordx;
    double coordy;
    int index = 0;

    boundary_data.seekg(0, boundary_data.end);  // set position at the end
    int length = boundary_data.tellg();         // tell which position is the end
    boundary_data.seekg(0, boundary_data.beg);  // set position back to beginning

    char* raw_point = new char[100];              
    boundary_data.getline(raw_point, 100);
    string point (raw_point);

    while (!boundary_data.eof())
    {
        comma  = point.find_first_of(",");
        coordx = stod(point.substr(0, comma));
        point  = point.substr(comma + 1);

        comma  = point.find_first_of(",");
        coordy = stod(point.substr(0, comma));
        point  = point.substr(comma + 1);

        genmesh->bdr_pointlist.push_back(Point(coordx, coordy));
        genmesh->pointlist.push_back(Point(coordx, coordy));

        boundary_data.getline(raw_point, 100);
        point = string(raw_point);

        index ++;
    }

    /**
     * 
     * Get the program to read the profile diagonals
     * 
     *   Cases of diagonals:
     *   
     *   Viewing the profile walls from <<OUTSIDE>> the wafer, those are:
     *  
     *              o------o                  o------o
     *     case 0:  |   /  |         case 1:  |  \   |
     *              |  /   |                  |   \  |
     *              o------o                  o------o
     *  
     * Please make sure that the number of diagonals is the same as the number
     * of boundary 2D points.
     * 
     */
    ifstream diagonals_data ("bdr_diagonals.dat");

    if (!diagonals_data)
    {
        cout << "Something is wrong with the diagonals input file." << endl;
        exit(1);
    }

    for (int i = 0; i < index; ++i)
    {
        genmesh -> profile_diagonals.push_back(diagonals_data.get() == '1');
    }

    genmesh -> make_3D_points();
    genmesh -> build_profile_mesh(index);
    genmesh -> stream_elements_out();
    genmesh -> stream_nodes_out();
    
    return 0;
}
