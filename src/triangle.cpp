#include "myutils.h"
#include "abcxyz.h"

#define TRACE	0

void GetTriangleCentroid(const vector<vector<double> > &MotifCoords,
  vector<double> &CentroidCoords)
	{
	CentroidCoords.resize(3);

	CentroidCoords[X] = (MotifCoords[0][X] + MotifCoords[1][X] + MotifCoords[2][X])/3.0;
	CentroidCoords[Y] = (MotifCoords[0][Y] + MotifCoords[1][Y] + MotifCoords[2][Y])/3.0;
	CentroidCoords[Z] = (MotifCoords[0][Z] + MotifCoords[1][Z] + MotifCoords[2][Z])/3.0;
	}

/***
        A               Y
      . .  .            |
     .  .	 .          |
    .   .      .        |
   S    .   N    . B     ----- X
        .      .            
        .	 .
        .  .
        C

  N = centroid
  S = point on N->B line perpendicular to A
        N->B line is Basis[X]
  X axis in direction (N) ---> B
  Y axis in plane ABC, right angles to x, in direction of A
***/
void GetTriangleBasis(const vector<vector<double> > &MotifCoords,
  vector<double> &CentroidCoords, vector<vector<double> > &Basis)
	{
#if TRACE
	Log("\n");
	Log("____________________________________\n");
	LogMx("MotifCoords", MotifCoords);
#endif
	Resize3x3(Basis);

	GetTriangleCentroid(MotifCoords, CentroidCoords);
#if TRACE
	LogVec("CentroidCoords", CentroidCoords);
#endif

	const vector<double> &PtA = MotifCoords[A];
	const vector<double> &PtB = MotifCoords[B];

	Sub_Vecs(PtB, CentroidCoords, Basis[0]);
	NormalizeVec(Basis[0]);
#if DEBUG
	double ModBasis0 = GetMod_Vec(Basis[0]);
	assert(feq(ModBasis0, 1));
#endif
#if 0 // TRACE
	LogVec("Basis[0]", Basis[0]);
#endif

	vector<double> VecAB(3);
	Sub_Vecs(MotifCoords[B], MotifCoords[A], VecAB);
	double ModAB = GetMod_Vec(VecAB);
	double ThetaNB_AB = GetTheta_Vecs(Basis[0], VecAB);
	double ModSB = ModAB*cos(ThetaNB_AB);

	vector<double> PtS(3);
	PtS[X] = PtB[X] - ModSB*Basis[0][X];
	PtS[Y] = PtB[Y] - ModSB*Basis[0][Y];
	PtS[Z] = PtB[Z] - ModSB*Basis[0][Z];

	vector<double> VecSA;
	Sub_Vecs(PtA, PtS, Basis[1]);
	NormalizeVec(Basis[1]);
#if 0 // TRACE
	LogVec("Basis[1]", Basis[1]);
#endif
#if DEBUG
	double ModBasis1 = GetMod_Vec(Basis[1]);
	assert(feq(ModBasis1, 1));
	double Theta01 = GetTheta_Vecs(Basis[0], Basis[1]);
	if (!feq(Theta01, PI/2))
		{
		Log("\n");
		LogVec("Basis X", Basis[0]);
		LogVec("Basis Y", Basis[1]);
		Log("Theta %.1f degrees (should be 90)\n",
		  degrees(Theta01));
		Die("Theta01");
		}
#endif

	CrossProduct(Basis[0], Basis[1], Basis[2]);
#if TRACE
	double ModBasis2 = GetMod_Vec(Basis[2]);
//	LogVec("Basis[2]", Basis[2]);
	LogMx("Basis", Basis);
#endif
	AssertUnitBasis(Basis);
	}
