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

    bool debug_flag {false};
    int layer {2};

    SetParameter parameters = SetParameter(argc, argv);
    
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

    while ()
    {
        CONTINUE HERE: falta todo el parrafo para tomar los nodes en R3.
            
            * buscar los indices del borde y borrarlo
            * ATENCION pensar en como usar un VALARRAY para la busqueda
                       de un double.
            * TAREA importante: pasar a funciones en GenMesh.cpp todo lo
                                de este archivo.
            * cambiar todos los print o cout por la template print<T>
    }

    genmesh->find_global_coordinates_for_boundary();
    genmesh->make_3D_points();
    genmesh->build_profile_mesh(index);
    genmesh->stream_elements_out();
    genmesh->stream_nodes_out();
    
    return 0;
}
