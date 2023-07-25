#ifndef _POINT_3D_
#define _POINT_3D_

#include <string>

#include <Point.h>

class Point3D
{
    public:
    
    double x;
    double y;
    double z;

    Point3D(/* args */);
    Point3D(const Point& point2D);
    Point3D(const Point& point2D, double c);
    Point3D(double a, double b, double c);
    Point3D(const Point3D& point);
    ~Point3D();

    friend ostream& operator<<(ostream& o, const Point3D& p);
    Point3D operator+(const Point3D& p) const;
    Point3D operator-(const Point3D& p) const;
    Point3D operator*(const double scale);
    Point3D operator/(const double scale);
    Point3D product(Point3D& v);  // vector product
    double operator*(const Point3D& p);  // inner product
    bool operator<(const Point3D& p);
    double norm();
    string split(char separator);
};

#endif
