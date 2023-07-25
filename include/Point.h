#ifndef _POINT_H_
#define _POINT_H_

#include <iostream>
#include <math.h>


using namespace std;


class Point
{

public:

    double m_x;
    double m_y;

    Point(/* args */);
    Point(double x, double y);
    Point(const Point &point);
    ~Point();
    
    Point operator+(const Point& p) const;
    Point operator-(const Point& p) const;
    Point operator*(const double scale);
    double operator*(const Point& p);  // inner product of two vector.
    Point operator/(const double scale);

    // To sort the points around Z axis.
    bool operator<(const Point& p);
    friend ostream& operator<<(ostream& o, const Point& p);
    double norm();
};

#endif
