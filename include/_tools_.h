#ifndef _TOOLS_H
#define _TOOLS_H

#ifdef DEBUG  // #define DEBUG
#define __info(x) std::cout << #x << " " << x << '\n';
#else
#define __info(x)
#endif

#define pi 3.141592653589793
#define two_pi 6.283185307179586
#define radian 0.017453292519943295


#include <iostream>
#include <cassert>


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


#endif
