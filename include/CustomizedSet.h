#ifndef _CUSTOMIZED_SET_
#define _CUSTOMIZED_SET_

#include <unordered_set>
#include <iostream>


using namespace std;

class CustomizedSet: public unordered_set<int> {

	public:

	friend ostream& operator<<(ostream& o, const CustomizedSet& cset);
};

#endif
