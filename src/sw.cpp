#include "myutils.h"
#include "mx.h"
#include "best.h"
#include "seqdb.h"

static const float LOG_ZERO = -2e20f;
static const float TransMD = -5.0f;
static const float TransMI = -5.0f;
static const float TransDD = -0.5f;
static const float TransII = -0.5f;
static const float BIAS = -0.5f;

//#define		m(x)	static Mx<float> g_Fwd##x; static Mx<char> g_TB##x;
//	m(M)
//	m(D)
//	m(I)
//#undef m

float SW(const Mx<float> &SMx,
  Mx<float> &a_FwdM,
  Mx<float> &a_FwdD,
  Mx<float> &a_FwdI,
  Mx<char> &a_TBM,
  Mx<char> &a_TBD,
  Mx<char> &a_TBI,
  uint &Starti, uint &Startj, string &Path)
	{
	Path.clear();
	const float * const *SimMx = SMx.GetData();

	const uint LA = SMx.m_RowCount;
	const uint LB = SMx.m_ColCount;

#define		m(x)	a_Fwd##x.Alloc("SAff_Fwd"#x, LA+1, LB+1);	\
					a_TB##x.Alloc("SWAff_TB"#x, LA+1, LB+1);	\
					float **Fwd##x = a_Fwd##x.GetData();											\
					char **TB##x = a_TB##x.GetData();
	m(M)
	m(D)
	m(I)
#undef m

	for (uint i = 0; i <= LA; ++i)
		{
		FwdM[i][0] = 0;
		FwdD[i][0] = LOG_ZERO;
		FwdI[i][0] = LOG_ZERO;
		TBM[i][0] = 'S';
		TBD[i][0] = '?';
		TBI[i][0] = '?';
		}

	for (uint j = 0; j <= LB; ++j)
		{
		FwdM[0][j] = 0;
		FwdD[0][j] = LOG_ZERO;
		FwdI[0][j] = LOG_ZERO;
		TBM[0][j] = 'S';
		TBD[0][j] = '?';
		TBI[0][j] = '?';
		}

// Main loop
	float BestScore = LOG_ZERO;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;
	for (uint i = 0; i < LA; ++i)
		{
		const float *SimMxRow = SimMx[i];
		for (uint j = 0; j < LB; ++j)
			{
		// xM
			{
			float Match = SimMxRow[j] + BIAS;
			float MM = FwdM[i][j] + Match;
			float DM = FwdD[i][j] + Match;
			float IM = FwdI[i][j] + Match;
			float SM = Match;
			float Score;
			Best4(MM, DM, IM, SM, 'M', 'D', 'I', 'S', Score, TBM[i+1][j+1]);
			FwdM[i+1][j+1] = Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				Besti = i+1;
				Bestj = j+1;
				}
			}
			
		// xD
			{
			float MD = FwdM[i][j+1] + TransMD;
			float DD = FwdD[i][j+1] + TransDD;
			Best2(MD, DD, 'M', 'D', FwdD[i+1][j+1], TBD[i+1][j+1]);
			}
			
		// xI
			{
			float MI = FwdM[i+1][j] + TransMI;
			float II = FwdI[i+1][j] + TransII;
			Best2(MI, II, 'M', 'I', FwdI[i+1][j+1], TBI[i+1][j+1]);
			}
			}
		}

	if (opt_trace)
		{
		SMx.LogMe();
		a_FwdM.LogMe();
		a_FwdD.LogMe();
		a_FwdI.LogMe();
		a_TBM.LogMe();
		a_TBD.LogMe();
		a_TBI.LogMe();
		}

	if (Besti == UINT_MAX)
		return LOG_ZERO;

	uint i = Besti;
	uint j = Bestj;
	if (opt_trace)
		{
		Log("BestScore=%g\n", BestScore);
		Log("Besti,j=%u,%u\n", Besti, Bestj);
		Log("    i      j  State\n");
		Log("-----  -----  -----\n");
		}

	Path.clear();
	char State = 'M';
	for (;;)
		{
		if (opt_trace)
			Log("%5u  %5u  %5c\n", i, j, State);

		Path.push_back(State);
		switch (State)
			{
		case 'M':
			{
			State = TBM[i][j];
			--i;
			--j;
			break;
			}
		case 'D':
			{
			State = TBD[i][j];
			--i;
			break;
			}
		case 'I':
			{
			State = TBI[i][j];
			--j;
			break;
			}
		default:
			asserta(false);
			}

		if (i == 0 || j == 0 || State == 'S')
			break;
		}

	Starti = i;
	Startj = j;
	reverse(Path.begin(), Path.end());
	return BestScore;
	}
