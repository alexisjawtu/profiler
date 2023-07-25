#include "Point.h"
#include "Point3D.h"
#include <math.h>


using namespace std;


Point3D::Point3D(/* args */): x(0), y(0), z(0) {}

Point3D::Point3D(const Point& point2D): x(point2D.m_x), y(point2D.m_y), z(0) {}

Point3D::Point3D(const Point& point2D, double c): x(point2D.m_x), y(point2D.m_y), z(c) {}

Point3D::Point3D(double a, double b, double c): x(a), y(b), z(c) {}

Point3D::Point3D(const Point3D& point) {
    x = point.x;
    y = point.y;
    z = point.z;
}

Point3D::~Point3D() {}

Point3D Point3D::operator+(const Point3D& p) const {
    return Point3D(x + p.x, y + p.y, z + p.z);
}

Point3D Point3D::operator-(const Point3D& p) const {
    return Point3D(x - p.x, y - p.y, z - p.z);
}

Point3D Point3D::operator*(const double scale) {
    return Point3D(x * scale, y * scale, z * scale);
}

Point3D Point3D::product(Point3D& v) {
    return Point3D(y*v.z - v.y*z, z*v.x - x*v.z, x*v.y - v.x*y);
}

double Point3D::operator*(const Point3D& p) {
    return x * p.x + y * p.y + z * p.z;
}

Point3D Point3D::operator/(const double scale) {
    return Point3D(x/scale, y/scale, z/scale);
}

bool Point3D::operator<(const Point3D& p) {
    // This is still a radial comparison
    return atan2(y, x) < atan2(p.y, p.x);
}

ostream& operator<<(ostream& o, const Point3D& p) {
    o << p.x << " " << p.y << " " << p.z;
    return o;
}

double Point3D::norm() {
    return sqrt(x * x + y * y + z * z);
}

string Point3D::split(char separator = ' ') {
    return to_string(x) + separator +
           to_string(y) + separator +
           to_string(z);
}