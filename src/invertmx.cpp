#include "myutils.h"
#include "abcxyz.h"

void MulMx(
  const vector<vector<double> > &A,
  const vector<vector<double> > &B,
  vector<vector<double> > &Prod)
	{
	Prod.resize(3);
	for (uint i = 0; i < 3; ++i)
		{
		Prod[i].resize(3);
		for (uint j = 0; j < 3; ++j)
			{
			double Sum = 0;
			for (uint k = 0; k < 3; ++k)
				Sum += A[i][k]*B[k][j];
			Prod[i][j] = Sum;
			}
		}
	}

double GetMxDeterminant(const vector<vector<double> > &Mx)
	{
	AssertMx3D(Mx);
	double Det = 0;
	for(uint i = 0; i < 3; ++i)
		Det += Mx[0][i]*(Mx[1][(i+1)%3]*Mx[2][(i+2)%3] - Mx[1][(i+2)%3]*Mx[2][(i+1)%3]);
	asserta(abs(Det) > 1e-6);
	return Det;
	}

void InvertMx(const vector<vector<double> > &Mx,
  vector<vector<double> > &InvMx)
	{
	AssertMx3D(Mx);
	InvMx.resize(3);

	//double Det = 0;
	//for(uint i = 0; i < 3; ++i)
	//	Det += Mx[0][i]*(Mx[1][(i+1)%3]*Mx[2][(i+2)%3] - Mx[1][(i+2)%3]*Mx[2][(i+1)%3]);
	//asserta(abs(Det) > 1e-6);
	double Det = GetMxDeterminant(Mx);

	InvMx.resize(3);
	for (uint i = 0; i < 3; ++i)
		{
		InvMx[i].resize(3);

		for(uint j = 0; j < 3; ++j)
			{
			double Top = (Mx[(j+1)%3][(i+1)%3]*Mx[(j+2)%3][(i+2)%3]) - (Mx[(j+1)%3][(i+2)%3]*Mx[(j+2)%3][(i+1)%3]);
			InvMx[i][j] = Top/Det;
			}
		}

#if DEBUG
	vector<vector<double> > Prod;
	MulMx(Mx, InvMx, Prod);

	for (uint i = 0; i < 3; ++i)
		{
		for (uint j = 0; j < 3; ++j)
			{
			if (i == j)
				assert(feq(Prod[i][j], 1));
			else
				assert(feq(Prod[i][j], 0));
			}
		}
#endif
	}

static void Test(
  double x1, double y1, double z1,
  double x2, double y2, double z2,
  double x3, double y3, double z3)
	{
	vector<vector<double> > Mx(3);
	Mx[0].resize(3);
	Mx[1].resize(3);
	Mx[2].resize(3);

	Mx[0][X] = x1;
	Mx[0][Y] = y1;
	Mx[0][Z] = z1;

	Mx[1][X] = x2;
	Mx[1][Y] = y2;
	Mx[1][Z] = z2;

	Mx[2][X] = x3;
	Mx[2][Y] = y3;
	Mx[2][Z] = z3;

	vector<vector<double> > Inv(3);
	InvertMx(Mx, Inv);
	}

void cmd_test_inv()
	{
	for (uint Iter = 0; Iter < 100; ++Iter)
		{
#define r	((randu32()%1000000)/1e6)
		Test(
			r, r, r,
			r, r, r,
			r, r, r
		  );
		}
	}
