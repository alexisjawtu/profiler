#include "the_map.h"

ostream& operator<<(ostream& o, const the_map& tm) 
{
	for (auto& i: tm) {
		o << i.second.m_x << ", " << i.second.m_y << endl;
	}

	return o;
}
