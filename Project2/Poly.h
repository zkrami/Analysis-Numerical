#pragma once
#include <vector>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

typedef long double ld;

const ld error = 1e-12;

class Poly {
	vector<pair<ld, int>> mypoly;
public:

	Poly();
	Poly(ld n) { mypoly.push_back(make_pair(n, 0)); }

	Poly(const vector<ld>&);

	inline bool iszero() {
		return !(fabsl(mypoly[0].first) > error || mypoly[0].second);
	}

	void add(ld a, int b) {
		mypoly.push_back(make_pair(a, b));
		this->shorten();
	}

	void shorten() {
		for (int i = 0; i < (int)mypoly.size(); i++)
			for (int j = i + 1; j < (int)mypoly.size(); j++)
				if (mypoly[i].second == mypoly[j].second) {
					mypoly[i].first += mypoly[j].first;
					mypoly.erase(mypoly.begin() + j--);
				}
				for (int i = 0; i < (int)mypoly.size(); i++)
					if (fabsl(mypoly[i].first) <= error)
						mypoly.erase(mypoly.begin() + i--);

				if (mypoly.empty())
					mypoly.push_back(make_pair(0.0, 0));
				sort(mypoly.begin(), mypoly.end(), [](pair<ld, int> a, pair<ld, int> b) {return a.second > b.second; });
	}

	Poly operator + (Poly temp) {
		Poly res = *this;
		res.mypoly.insert(res.mypoly.end(), temp.mypoly.begin(), temp.mypoly.end());
		res.shorten();
		return res;
	}
	void operator += (Poly temp) { *this = *this + temp; }
	Poly operator + (ld temp) {
		Poly res = *this;
		res.mypoly.push_back(make_pair(temp, 0));
		res.shorten();
		return res;
	}
	Poly operator - (Poly temp) {
		Poly res = *this;
		for (int i = 0; i < (int)temp.mypoly.size(); i++)
			res.mypoly.push_back(make_pair(-temp.mypoly[i].first, temp.mypoly[i].second));
		res.shorten();
		return res;
	}
	void operator -= (Poly temp) { *this = *this - temp; }

	Poly operator * (Poly temp) {
		Poly res;
		for (int i = 0; i < (int)mypoly.size(); i++)
			for (int j = 0; j < temp.mypoly.size(); j++)
				res.mypoly.push_back(make_pair(mypoly[i].first * temp.mypoly[j].first, mypoly[i].second + temp.mypoly[j].second));
		res.shorten();
		return res;
	}
	void operator *= (Poly temp) { *this = *this * temp; }

	Poly operator * (ld temp) {
		Poly res = *this;
		for (int i = 0; i < mypoly.size(); i++)
			res.mypoly[i].first *= temp;
		return res;
	}
	void operator *= (ld temp) { *this = *this * temp; }

	Poly operator / (Poly divr) {
		Poly res_temp, res, div = *this;
		if (div.iszero())
			return res;
		if (divr.iszero()) {
			printf("Cannot divide by 0!\n");
			int q = 0; q /= q;
		}
		int n = 30;
		while (n && !div.iszero()) { //div.mypoly[0].second >= divr.mypoly[0].second
			if (div.mypoly[0].second <= 0)
				n--;
			res_temp.mypoly.clear();
			res_temp.mypoly.push_back(make_pair(div.mypoly[0].first / divr.mypoly[0].first, div.mypoly[0].second - divr.mypoly[0].second));
			res += res_temp;
			div -= res_temp * divr;
			div.shorten();
		}
		res.shorten();
		return res;
	}
	void operator /= (Poly temp) { *this = *this / temp; }

	Poly operator / (ld temp) { return *this * (1 / temp); }
	void operator /= (ld temp) { *this = *this / temp; }

	Poly operator ^ (int n) {
		Poly res = *this;
		while (n--) {
			for (int i = 0; i < res.mypoly.size(); i++) {
				res.mypoly[i].first *= (res.mypoly[i].second);
				res.mypoly[i].second--;
			}
			res.shorten();
		}
		return res;
	}
	void operator ^= (int n) { *this = *this ^ n; }
	void operator = (Poly temp) { mypoly = temp.mypoly; }
	void operator = (ld temp) { *this = Poly(temp); }
	bool operator == (Poly temp) { return (mypoly == temp.mypoly); }
	bool operator != (Poly temp) { return (!(mypoly == temp.mypoly)); }
	bool operator == (ld temp) { if (temp == 0) return (iszero()); }

	ld value(ld x) {
		ld res = 0.0;
		for (int i = 0; i < (int)mypoly.size(); i++)
			res += mypoly[i].first * powl(x, mypoly[i].second);
		return res;
	}
	string to_string() {
		shorten();
		ostringstream oss;
		oss << "P(x) = ";
		if (iszero())
			oss << "0";
		else {
			if (abs(mypoly[0].first - 1.0) <= error) {
				if (mypoly[0].second == 0) 
					oss << "1";
				else if (mypoly[0].second == 1) 
					oss << "x";
				else 
					oss << "x^" << mypoly[0].second;
			} else if (abs(mypoly[0].first + 1.0) <= error) {
				if (mypoly[0].second == 0) 
					oss << "-1";
				else if (mypoly[0].second == 1) 
					oss << "-x";
				else 
					oss << "-x^" << mypoly[0].second;
			} else {
				if (mypoly[0].second == 0) oss << mypoly[0].first;
				else if (mypoly[0].second == 1)  oss << mypoly[0].first << " x";
				else oss << mypoly[0].first << " x^" << mypoly[0].second;
			}
			for (int i = 1; i < (int)mypoly.size(); i++) {
				if (abs(mypoly[i].first - 1.0) <= error) {
					if (mypoly[i].second == 0)
						oss << " + 1";
					else if (mypoly[i].second == 1)
						oss << " + x";
					else
						oss << " + x^" << mypoly[i].second;
				} else if (abs(mypoly[i].first + 1.0) <= error) {
					if (mypoly[i].second == 0)
						oss << " - 1";
					else if (mypoly[i].second == 1)
						oss << " - x";
					else
						oss << " - x^" << mypoly[i].second;
				} else {
					if (mypoly[i].second == 0)
						mypoly[i].first < 0 ? oss << " - " << abs(mypoly[i].first) : oss << " + " << mypoly[i].first;
					else if (mypoly[i].second == 1) {
						mypoly[i].first < 0 ? oss << " - " << abs(mypoly[i].first) : oss << " + " << mypoly[i].first;
						oss << " x";
					} else {
						mypoly[i].first < 0 ? oss << " - " << abs(mypoly[i].first) : oss << " + " << mypoly[i].first;
						oss << " x^" << mypoly[i].second;
					}
				}
			}
		}
		return oss.str();
	}
	vector <pair<ld, ld> > Draws(ld x , ld y , ld factorx , ld factory )
	{
		vector <pair<ld, ld> > res;
		for (ld i = -1000; i <= 1000; i+=0.5)
		{
			ld tx = i; 
			ld ty = value(i);
			if (tx*factorx> x / 2 || tx*factorx < -x / 2) continue; 
			if (ty*factory > y / 2 || ty*factory < -y / 2) continue;
			res.emplace_back(tx*factorx, ty*factory);
		}
		return res;
	}
};

struct spline_solution {
	spline_solution(int n) {
		v.resize(n);
	}

	struct spline_equation {
		Poly p;
		pair<ld, ld> interval;
		string to_string() {
			ostringstream oss;
			string s = p.to_string();
			s.erase(s.begin(), s.begin() + 7);
			oss << s << setw(50) << "  ; " << interval.first << " <= x <= " << interval.second << endl;
			return oss.str();
		}
	};

	vector<spline_equation> v;
	string spline_print() {
		ostringstream oss;
		for (int i = 0; i < (int)v.size(); i++) {
			oss << "S" << i << "(x) = " << v[i].to_string();
		}
		return oss.str();
	}
	spline_equation& operator[](unsigned int i) {
		return v[i];
	}


	ld evaluate(ld x) {
		for (int i = 0; i < (int)v.size(); i++) {
			if (x >= v[i].interval.first && x <= v[i].interval.second)
				return v[i].p.value(x);
		}
		throw  exception("Value out of range!\n");
	}
};

