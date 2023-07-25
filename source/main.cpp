#include <iostream>
// to_delete #include <fstream>
// to_delete #include <iomanip>
#include "GenMesh.h"
#include "SetParameter.h"
#include <string>

// to_delete #include <cstdio>
// to_delete #include <cstdlib>
// to_delete #include <stdio.h>
// to_delete #include <stdlib.h>


using namespace std;


/** 
    cppreference.com:
* 
    std::size_t is the unsigned integer type of the result of the sizeof operator.
*   It depends on the architecture.

*/


int main(int argc, char *argv[])
{
    cout << endl << "Profiler v1.1.0 July 2023." << endl;
    cout << endl << "Reminder: The boundary points in the input are assumed\n"
                    "to be the 2D points as stored in "
                    "vector<size_t> GenMesh::bdr_pointlist." << endl << endl;

    bool debug_flag = false;
    // to_delete int wafer_type   = 1;
    int layer        = 2;
    // to_delete Point A          = Point(-14.15, 0);
    // to_delete Point D          = Point(-14.15, 0.3);
    // to_delete double rescaled_diameter = 15;  // WAS R
    // to_delete double h                 = 0.01*0.75;   
    // to_delete double H                 = 0.75;
    // to_delete double L1                = 0.0;
    // to_delete double L2                = 0.0;
    // to_delete double scale             = 2.0;
    // to_delete double theta_p = 0;
    // to_delete double pitch = 1;
    // to_delete double pitch_error_AB = 0;
    // to_delete double pitch_error_BC = 0;
    // to_delete double pitch_error_CD = 0;
    // to_delete double angle_diff = 90;

    string prefix(".");

    SetParameter parameters = SetParameter(argc, argv);
    double thickness_of_inner_wafer = parameters.user_thickness_of_inner_wafer;
    double dt = thickness_of_inner_wafer/(layer - 1);
    
    GenMesh* genmesh = new GenMesh(
                        // to_delete wafer_type,
                        // to_delete A, 
                        // to_delete D, 
                        // to_delete H, 
                        // to_delete h, 
                        // to_delete rescaled_diameter, 
                        // to_delete theta_p/180.0 * PI, 
                        // to_delete L1, 
                        // to_delete L2, 
                        // to_delete angle_diff/180.0 * PI,
                        layer,
                        dt,
                        debug_flag,
                        prefix,
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
        genmesh->profile_diagonals.push_back(diagonals_data.get() == '1');
    }

    genmesh->make_3D_points();
    genmesh->build_profile_mesh(index);
    genmesh->stream_elements_out();
    genmesh->stream_nodes_out();
    
    return 0;
}
