#include "Stresses.h"

Stresses::Stresses() {}
Stresses::Stresses(vector<double> v)
{
	stresses.resize(3);
	for (int i = 0; i < 3; i++) stresses[i] = v[i];
}
Stresses::Stresses(double a, double b, double c)
{
	stresses.resize(3);
	stresses[0] = a; stresses[1] = b; stresses[2] = c;
}
Stresses& Stresses::operator=(Stresses b)
{
	stresses = b.stresses;
	return *this;
}
vector<double> Stresses::GetStrains() { return stresses; }

ostream& operator<<(ostream& cout, const Stresses& b)
{
	for (int i = 0; i < 2; i++)
	{
		cout << b.stresses[i] << " ";
	}
	cout << b.stresses[2] << endl;
	return cout;
}
