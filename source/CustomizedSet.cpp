#include "CustomizedSet.h"

ostream& operator<<(ostream& o, const CustomizedSet& cset)
{
	for (auto& x: cset) {
		o << x << ", ";
	}
	o << endl;

	return o;
}
