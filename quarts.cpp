#include "myutils.h"
#include "sort.h"
#include "quarts.h"

void GetQuarts(const vector<unsigned> &v, Quarts &Q)
	{
	const unsigned N = SIZE(v);
	Q.Min = 0;
	Q.LoQ = 0;
	Q.Med = 0;
	Q.HiQ = 0;
	Q.Max = 0;
	Q.Total = 0;
	Q.Avg = 0.0;
	if (N == 0)
		return;

	vector<unsigned> v2 = v;
	unsigned *vs = v2.data();
	QuickSortInPlace(vs, N);

	for (unsigned i = 0; i < N; ++i)
		Q.Total += vs[i];

	Q.Min = vs[0];
	Q.LoQ = vs[N/4];
	Q.Med = vs[N/2];
	Q.HiQ = vs[(3*N)/4];
	Q.Max = vs[N-1];
	Q.Avg = double(Q.Total)/N;
	}

void GetQuartsDouble(const vector<double> &v, QuartsDouble &Q)
	{
	const unsigned N = SIZE(v);
	Q.Min = 0.0f;
	Q.LoQ = 0.0f;
	Q.Med = 0.0f;
	Q.HiQ = 0.0f;
	Q.Max = 0.0f;
	Q.Total = 0.0f;
	Q.Avg = 0.0f;
	if (N == 0)
		return;

	vector<double> v2 = v;
	double *vs = v2.data();
	QuickSortInPlace(vs, N);

	for (unsigned i = 0; i < N; ++i)
		Q.Total += vs[i];

	double Mean = double(Q.Total)/N;
	double Sumd = 0.0f;
	for (unsigned i = 0; i < N; ++i)
		{
		double x = vs[i];
		double d = (x - Mean)*(x - Mean);
		Sumd += d;
		}
	double StdDev = (double) sqrt(Sumd/N);

	Q.Min = vs[0];
	Q.LoQ = vs[N/4];
	Q.Med = vs[N/2];
	Q.HiQ = vs[(3*N)/4];
	Q.Max = vs[N-1];
	Q.Avg = Mean;
	Q.StdDev = StdDev;
	}
