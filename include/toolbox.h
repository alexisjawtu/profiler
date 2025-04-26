#ifndef _TOOLS_H
#define _TOOLS_H

#ifdef DEBUG  // #define DEBUG
#define __info(x) std::cout << #x << " " << x << '\n';
#else
#define __info(x)
#endif


#include <charconv>  // to_chars converter
#include <iomanip>  // manipulators of input and output
#include <iostream>
#include <cassert>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>
 
#include "Point3D.h"

typedef tuple<vector<int>, vector<int>, int> Brick;
// We may use this to measure elements and work the matrix condition
typedef tuple<Point3D, Point3D, Point3D, Point3D> Tetrahedron;


namespace constants
{
    const bool debug {false};
    constexpr int horizontal_mesh_layers {1};
    constexpr double pi {3.141592653589793};
    const double two_pi {6.283185307179586};
    const double radian {0.017453292519943295};
}


namespace filenames
{
    const char output_dir[] {"./"};
    const char scalar_parameters[] {"input.txt"};
    const char sorted_3D_bdr_vertices[] {"sorted_boundary_vertices_on_space.csv"};
    const char bdr_vertices_2D[] {"2D_bdr_vertices.dat"};
    const char prof_elems[] {"profile_elements.dat"};
    const char prof_verts[] {"profile_vertices.dat"};
    const char cylinder_elems[] {"cylinder_elements.dat"};
    const char cylinder_verts[] {"cylinder_vertices.dat"};
    const char cylinder_verts_by_elems[] {"cylinder_verts_by_elems.dat"};
}


template <class T> void print(T data)
{
    std::cout << std::boolalpha << data << "\n";
}


template <class T> void print_vector(T& v)
{
    for (auto& d: v)
    {
        std::cout << std::boolalpha << d << ", ";
    }
    std::cout << "\n";
}


std::string num_to_str (auto numeric) {
    const size_t buf_size = 10;
    char buffer[buf_size] {};
    std::to_chars_result result = std::to_chars(buffer, buffer + buf_size, numeric);

    if (result.ec != std::errc())
        std::cout << std::make_error_code(result.ec).message() << '\n';
    
    std::string r(std::string_view(buffer, result.ptr - buffer));
    return r;
}


void num_to_chars (auto numeric) {
    const size_t buf_size = 10;
    char buffer[buf_size] {};
    std::to_chars_result result = std::to_chars(buffer, buffer + buf_size, numeric);

    if (result.ec != std::errc())
        std::cout << std::make_error_code(result.ec).message() << '\n';
    else {
        std::string_view str(buffer, result.ptr - buffer);
        std::cout << std::quoted(str) << '\n';
    }
}


#endif
