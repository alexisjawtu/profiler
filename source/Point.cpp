#include "Point.h"
#include "Point3D.h"
#include <math.h>


using namespace std;


Point::Point(/* args */) : m_x(0), m_y(0) {}

Point::Point(double x, double y) : m_x(x), m_y(y) {}

Point::Point(const Point& point)
{
    m_x = point.m_x;
    m_y = point.m_y;
}

/*
Point::Point(const Point3D& point3D)
{
    m_x = point3D.x;
    m_y = point3D.y;
}
*/

Point::~Point() {}

Point Point::operator+(const Point& p) const {
    Point temp;
    temp.m_x = m_x + p.m_x;
    temp.m_y = m_y + p.m_y;
    return temp;
}

Point Point::operator-(const Point& p) const {
    Point temp;
    temp.m_x = m_x - p.m_x;
    temp.m_y = m_y - p.m_y;
    return temp;
}

Point Point::operator*(const double scale) {
    Point temp;
    temp.m_x = m_x*scale;
    temp.m_y = m_y*scale;
    return temp;
}

double Point::operator*(const Point& _p) {
    return m_x*_p.m_x + m_y*_p.m_y;
}

Point Point::operator/(const double scale) {
    Point temp;
    temp.m_x = m_x/scale;
    temp.m_y = m_y/scale;
    return temp;
}

bool Point::operator<(const Point& p) {
    return atan2(m_y, m_x) < atan2(p.m_y, p.m_x);
}

ostream& operator<<(ostream& o, const Point& p) {
    o << p.m_x << ", " << p.m_y;
    return o;
}

double Point::norm(){
    return sqrt(m_x*m_x+m_y*m_y);
}
