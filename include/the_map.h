#ifndef _THE_MAP_H_
#define _THE_MAP_H_

#include <map>
#include <iostream>

#include "Point.h"


using namespace std; 

class the_map: public map<int, Point>
{
	public:
	
	friend ostream& operator<<(ostream& o, const the_map& tm);
};

#endif
