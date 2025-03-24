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

    string prefix(".");
    bool debug_flag {false};
    // number of horizontal 2D layers
    int layer {2};

    SetParameter parameters = SetParameter(argc, argv);
    double thickness_of_inner_wafer = parameters.user_thickness_of_inner_wafer;
    double dt = thickness_of_inner_wafer/(layer - 1);

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
    int number_of_bdr_points = 0;               // Exactly the number of vertical walls.

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

        number_of_bdr_points ++;
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

    for (int i = 0; i < number_of_bdr_points; ++i)
    {
        genmesh -> profile_diagonals.push_back(diagonals_data.get() == '1');
    }

    genmesh -> make_3D_points();
    genmesh -> build_profile_mesh(number_of_bdr_points);
    genmesh -> stream_elements_out();
    genmesh -> stream_nodes_out();
    
    return 0;
}
