#include "myutils.h"
#include "diagbox.h"

#define TEST	0

/***
	i = 0..LA-1
	j = 0..LB-1
	d = LA - i + j = 1 .. LA+LB-1
	j = d - LA + i
	i = LA - d + j
***/

void GetDiagRange(unsigned LA, unsigned LB, unsigned d,
  unsigned &mini, unsigned &minj, unsigned &maxi, unsigned &maxj)
	{
	if (d >= LA)
		{
		mini = 0;
		maxi = min(LA+LB-1-d, LA-1);
		minj = d - LA;
		maxj = min(LB-1, d-1);
		}
	else
		{
		mini = LA-d;
		maxi = min(LA+LB-1-d, LA-1);
		minj = 0;
		maxj = min(LB-1, d-1);
		}
	}

void GetDiagBox(unsigned LA, unsigned LB, unsigned DiagLo, unsigned DiagHi, DiagBox &Box)
	{
	asserta(DiagLo <= DiagHi);
	asserta(DiagLo >= 1);
	asserta(DiagHi <= LA + LB - 1);

	Box.LA = LA;
	Box.LB = LB;

	Box.dlo = DiagLo;
	Box.dhi = DiagHi;

	GetDiagRange(LA, LB, DiagLo, Box.dlo_mini, Box.dlo_minj, Box.dlo_maxi, Box.dlo_maxj);
	GetDiagRange(LA, LB, DiagHi, Box.dhi_mini, Box.dhi_minj, Box.dhi_maxi, Box.dhi_maxj);
	}

void GetDiagLoHi(unsigned LA, unsigned LB, const char *Path,
  unsigned &dlo, unsigned &dhi)
	{
	dlo = UINT_MAX;
	dhi = UINT_MAX;

	unsigned i = 0;
	unsigned j = 0;
	for (unsigned k = 0; ; ++k)
		{
		char c = Path[k];
		if (c == 0)
			break;
		if (c == 'M')
			{
			unsigned d = LA - i + j;
			if (dlo == UINT_MAX)
				{
				dlo = d;
				dhi = d;
				}
			else
				{
				if (d < dlo)
					dlo = d;
				if (d > dhi)
					dhi = d;
				}
			}
		if (c == 'M' || c == 'D')
			++i;
		if (c == 'M' || c == 'I')
			++j;
		}
	}

#if TEST
static void Test2(unsigned LA, unsigned LB, unsigned DiagLo, unsigned DiagHi)
	{
	DiagBox Box;
	GetDiagBox(LA, LB, DiagLo, DiagHi, Box);
	Box.LogMe();
	Box.Validate();
	}

static void Test1(unsigned LA, unsigned LB, unsigned d,
  unsigned i, unsigned j, unsigned I, unsigned J)
	{
	unsigned mini, maxi, minj, maxj;
	GetDiagRange(LA, LB, d, mini, minj, maxi, maxj);
	Log("LA=%u LB=%u d=%u (%u,%u) (%u,%u) expected (%u,%u) (%u,%u)\n",
	  LA, LB, d, mini, minj, maxi, maxj, i, j, I, J);

	asserta(mini == i);
	asserta(maxi == I);
	asserta(minj == j);
	asserta(maxj == J);
	}

void TestDiagBox()
	{
	Test2(16, 19, 17, 37);

	Test1(5, 3, 1, 4, 0, 4, 0);
	Test1(5, 3, 2, 3, 0, 4, 1);
	Test1(5, 3, 3, 2, 0, 4, 2);
	Test1(5, 3, 4, 1, 0, 3, 2);
	Test1(5, 3, 5, 0, 0, 2, 2);
	Test1(5, 3, 6, 0, 1, 1, 2);
	Test1(5, 3, 7, 0, 2, 0, 2);

	Test1(3, 5, 1, 2, 0, 2, 0);
	Test1(3, 5, 2, 1, 0, 2, 1);
	Test1(3, 5, 3, 0, 0, 2, 2);
	Test1(3, 5, 4, 0, 1, 2, 3);
	Test1(3, 5, 5, 0, 2, 2, 4);
	Test1(3, 5, 6, 0, 3, 1, 4);
	Test1(3, 5, 7, 0, 4, 0, 4);

	Test1(5, 5, 1, 4, 0, 4, 0);
	Test1(5, 5, 2, 3, 0, 4, 1);
	Test1(5, 5, 3, 2, 0, 4, 2);
	Test1(5, 5, 4, 1, 0, 4, 3);
	Test1(5, 5, 5, 0, 0, 4, 4);
	Test1(5, 5, 6, 0, 1, 3, 4);
	Test1(5, 5, 7, 0, 2, 2, 4);
	Test1(5, 5, 8, 0, 3, 1, 4);
	Test1(5, 5, 9, 0, 4, 0, 4);

	for (unsigned LA = 2; LA <= 5; ++LA)
		for (unsigned LB = 2; LB <= 5; ++LB)
			for (unsigned dlo = 1; dlo <= LA+LB-1; ++dlo)
				for (unsigned dhi = dlo; dhi <= LA+LB-1; ++dhi)
					Test2(LA, LB, dlo, dhi);

	Log("\n");
	Log("ALL OK\n");
	}
#endif // TEST
