#include "Strains.h"

Strains::Strains() {}
Strains::Strains(vector<double> v)
{
	strains.resize(3);
	for (int i = 0; i < 3; i++) strains[i] = v[i];
}
Strains::Strains(double a, double b, double c)
{
	strains.resize(3);
	strains[0] = a; strains[1] = b; strains[2] = c;
}
Strains& Strains::operator=(Strains b)
{
	strains = b.strains;
	return *this;
}
vector<double> Strains::GetStrains() { return strains; }

ostream& operator<<(ostream& cout, const Strains& b)
{
	for (int i = 0; i < 3; i++)
	{
		cout << b.strains[i] << " ";
	}
	cout << endl;
	return cout;
}