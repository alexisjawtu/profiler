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

    ifstream boundary_data ("bdr_nodes.dat");
    ifstream elements_data ("elements.dat");
    ifstream nodes_data    ("nodes.dat");

    if (!boundary_data || !nodes_data || !elements_data)
    {
        print("\vSomething is wrong with the input files.\n");
        exit(1);
    }

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

    /** TODO 
     * 1- genmesh has to be allocated regularly, not in dynamic memory
     * 2- if we keep this design of class GenMesh,
     *    then the constructor simply takes a SetParameter, as in
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

    int current_node {0};
    int comma;                                  // Stores the places of the commas.
    double coordx;                              // 2D coordinates of boundary points.
    double coordy;                              
    int row;                                    // Stores the index of a tetrahedron.
    int index {0};                              // Stores the index of a 2D point.
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

        genmesh->bdr_pointlist.push_back(Point(coordx, coordy));
        genmesh->pointlist.push_back(Point(coordx, coordy));

        boundary_data.getline(raw_point, 100);
        point = string(raw_point);

        index ++;
    }

    delete[] raw_point;
    boundary_data.close();

    comma = 0;
    elements_data.seekg(0, elements_data.end);  // set position at the end
    int length {elements_data.tellg()};         // tell which position is the end
    elements_data.seekg(0, elements_data.beg);  // set position back to beginning

    char* raw_tetrahedron = new char[length];              
    elements_data.getline(raw_tetrahedron, length);
    string tetrahedron(raw_tetrahedron);

    row = 0;
    while (!elements_data.eof())
    {
        comma = tetrahedron.find_first_of(",");
        stringstream(tetrahedron.substr(0, comma)) >> current_node;
        element.push_back(current_node);
        tetrahedron = tetrahedron.substr(comma + 1);

        comma = tetrahedron.find_first_of(",");
        stringstream(tetrahedron.substr(0, comma)) >> current_node;
        element.push_back(current_node);
        tetrahedron = tetrahedron.substr(comma + 1);

        comma = tetrahedron.find_first_of(",");
        stringstream(tetrahedron.substr(0, comma)) >> current_node;
        element.push_back(current_node);
        tetrahedron = tetrahedron.substr(comma + 1);

        comma = tetrahedron.find_first_of(",");
        stringstream(tetrahedron.substr(0, comma)) >> current_node;
        element.push_back(current_node);

        genmesh->elements_by_vertices.emplace(row, element);
        elements_data.getline(raw_tetrahedron, length);
        tetrahedron = string(raw_tetrahedron);        
        element     = {};
        row        += 1;
    }

    delete[] raw_tetrahedron;
    elements_data.close();

    while (false)
    {
        /**
         *  TODOs CONTINUE HERE: all this cycle body to receive the 3D nodes missing.
         *
         *  Look for the indices in the boundary and erase it
         *  
         *  TASK Research how to use <valarray> to traverse and search a double
         *  TASK important transform the code in this 'main' into methods of GenMesh.
                 leave a minimal 'int main () {}'
         *  TASK replace the prints and couts with template print<T>
        **/
    }

    genmesh->find_global_coordinates_for_boundary();
    genmesh->make_3D_points();
    genmesh->build_profile_mesh(index);
    genmesh->stream_elements_out();
    genmesh->stream_nodes_out();
    
    return 0;
}
