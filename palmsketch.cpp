#include "myutils.h"

static uint Geti(char c)
	{
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 3;
		}
	asserta(false);
	return 0;
	}

static char GetSymbol(char c1, char c2, char c3)
	{
	uint i1 = Geti(c1);
	uint i2 = Geti(c2);
	uint i3 = Geti(c3);

	vector<uint> Counts(4);

	Counts[i1] += 1;
	Counts[i2] += 1;
	Counts[i3] += 1;

	uint Maxi = 0;
	uint Maxn = 0;
	for (uint i = 0; i < 4; ++i)
		{
		if (Counts[i] > Maxn)
			{
			Maxn = Counts[Maxi];
			Maxi = i;
			}
		}
	return "hst~"[Maxi];
	}

void GetPalmSketch(const string &ss, uint PSL, string &Sketch)
	{
	Sketch.clear();
	const uint L = SIZE(ss);
	for (uint i = 0; i < PSL; ++i)
		{
		uint k = 2 + (i*(L-4))/PSL;
		asserta(k >= 0 && k + 2 < L);

		char c1 = ss[k];
		char c2 = ss[k+1];
		char c3 = ss[k+2];
		char c = GetSymbol(c1, c2, c3);
		Sketch += c;
		}
	asserta(SIZE(Sketch) == PSL);
	}

void GetPalmSubSketch(const string &ss, uint Lo, uint n, 
  uint PSL, string &Sketch)
	{
	Sketch.clear();
	const uint L = SIZE(ss);
	for (uint j = 0; j < PSL; ++j)
		{
		uint k = Lo + (j*n)/PSL;
		asserta(k >= 0 && k + 2 < L);

		char c1 = ss[k];
		char c2 = ss[k + 1];
		char c3 = ss[k + 2];
		char c = GetSymbol(c1, c2, c3);
		Sketch += c;
		}
	asserta(SIZE(Sketch) == PSL);
	}
