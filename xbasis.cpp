#include "myutils.h"
#include "abcxyz.h"

#define TRACE	0

// Calculate rotation matrix to transform basis
// unit vectors to (x, y, z)
void GetBasisR(const vector<vector<double> > &Basis,
  vector<vector<double> > &R)
	{
	AssertUnitBasis(Basis);

	InvertMx(Basis, R);

#if DEBUG
	{
	vector<vector<double> > RotatedBasis;
	MulMx(R, Basis, RotatedBasis);
	AssertCanonicalUnitBasis(RotatedBasis);
	}
#endif

#if TRACE
	Log("___________________________\n");
	Log("GetBasisR");
	LogMx("Basis", Basis);
	LogMx("R", R);
#endif
	}

static void Test(
  double x1, double y1, double z1,
  double x2, double y2, double z2,
  double x3, double y3, double z3)
	{
	vector<vector<double> > Basis(3);
	Basis[0].resize(3);
	Basis[1].resize(3);
	Basis[2].resize(3);

	Basis[0][X] = x1;
	Basis[0][Y] = y1;
	Basis[0][Z] = z1;

	Basis[1][X] = x2;
	Basis[1][Y] = y2;
	Basis[1][Z] = z2;

	Basis[2][X] = x3;
	Basis[2][Y] = y3;
	Basis[2][Z] = z3;

	vector<vector<double> > R;
	GetBasisR(Basis, R);
	}

static void Test(const vector<vector<double> > &Basis)
	{
	vector<vector<double> > R;
	GetBasisR(Basis, R);
	}

static double GetRandomAngle()
	{
	uint r = randu32();
	double d = double(r);
	double Theta = remainder(r, 2*PI);
	return Theta;
	}

static void GetRandomUnitBasis(vector<vector<double> > &Basis)
	{
	Basis.resize(3);
	for (uint i = 0; i < 3; ++i)
		Basis[i].resize(3);

	double Alpha = GetRandomAngle();
	double Beta = GetRandomAngle();
	double Gamma = GetRandomAngle();

	double sA = sin(Alpha);
	double cA = cos(Alpha);

	double sB = sin(Beta);
	double cB = cos(Beta);

	double sG = sin(Gamma);
	double cG = cos(Gamma);

	Basis[0][0] = cA*cB;
	Basis[0][1] = cA*sB*sG - sA*cG;
	Basis[0][2] = cA*sB*cG + sA*sG;

	Basis[1][0] = sA*cB;
	Basis[1][1] = sA*sB*sG + cA*cG;
	Basis[1][2] = sA*sB*cG - cA*sG;

	Basis[2][0] = -sB;
	Basis[2][1] = cB*sG;
	Basis[2][2] = cB*cG;

	LogMx("basis", Basis);

	AssertUnitBasis(Basis);
	}

void cmd_test_xbasis()
	{
	const double S2 = sqrt(2.0)/2.0;
	if (1)
		Test(
		  1, 0, 0,
		  0, 1, 0,
		  0, 0, 1);

	if (1)
		Test(
		  S2, S2, 0,
		  S2, -S2, 0,
		  0, 0, 1);

	if (1)
		Test(
		  0, 0, 1,
		  S2, S2, 0,
		  S2, -S2, 0);

	vector<vector<double> > Basis;
	for (uint Iter = 0; Iter < 100; ++Iter)
		{
		GetRandomUnitBasis(Basis);
		Test(Basis);
		}
	}
