#pragma once
#include "Header.h"

class Stresses
{
private:
	vector<double> stresses;
public:
	Stresses();
	Stresses(vector<double> v);
	Stresses(double a, double b, double c);
	Stresses& operator=(Stresses b);
	vector<double> GetStrains();
	friend ostream& operator<<(ostream& cout, const Stresses& b);
};