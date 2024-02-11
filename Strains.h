#pragma once
#include "Header.h"

class Strains
{
private:
	vector<double> strains;
public:
	Strains();
	Strains(vector<double> v);
	Strains(double a, double b, double c);
	Strains& operator=(Strains b);
	vector<double> GetStrains();
	friend ostream& operator<<(ostream& cout, const Strains& b);
};
