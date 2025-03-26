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
 

namespace constants {
    const double pi = 3.141592653589793;
    const double two_pi = 6.283185307179586;
    const double radian = 0.017453292519943295;
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
