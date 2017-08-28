
#include "Poly.h"
#include <vector>

using namespace std;

Poly::Poly()
{
	mypoly = vector<pair<ld, int>>();
}



Poly::Poly(const vector<ld>& t)
{
	for (int i = 0; i < (int)t.size(); i++)
	{
		if (t[i] != 0)
			mypoly.emplace_back(t[i], i);
	}
}
