#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "abcxyz.h"
#include <map>

/***
DaliLite v5
comparemodules.f, line 1436
===========================
        enveloperadius=20.0
        x=1/(enveloperadius*enveloperadius)
        do i=0,100
                wght(i)=exp(-x*i*i)
        end do
***/
static double Weight(double y)
	{
	const double D = 20.0;
	const double x = 1.0/(D*D);
	double w = exp(-x*y*y);
	return w;
	}

/***
comparemodules.f, line 1397
  a, b are integer distances in units of 1/10 Angstrom,
  so multiply by 10 to get Angstroms.
===========================
        function dpscorefun(a,b) result(s)
        implicit none
        include 'parsizes.for'
        real s
        integer*2 a,b
c
        real x,y,d0
        logical lela
        parameter(lela=.true.)
        parameter(d0=0.20)
c !!!   elastic uses weights !!!
        x=float(abs(a-b))/10
        if(lela) then
                y=float(a+b)/20
                if(y.gt.100) then
                        s=0.0
                else
                        if(y.gt.0) then
                          s=wght(nint(y))*(d0-x/y)
                        else
                          s=wght(nint(y))*d0
                        end if
                end if
        end if
***/
static double dpscorefun(double a, double b)
	{
	double Score = 0;
	const double d0 = 0.20;
	double x = fabs(a - b);
	double y = (a + b)/2;
	if (y > 100)
		Score = 0;
	else
		{
		if (y > 0)
			Score = Weight(y)*(d0 - x/y);
		else
			Score = Weight(y)*d0;
		}
	return Score;
	}

static double GetDALIScore_Col_Band(
  const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Col, uint Radius)
	{
	const int Lali = (int) SIZE(PosQs);
	asserta(SIZE(PosTs) == (uint) Lali);

	const double Theta = 0.2;
	double Sum = 0;
	uint PosQi = PosQs[Col];
	uint PosTi = PosTs[Col];
	int LoCol = int(Col) - int(Radius);
	int HiCol = int(Col) + int(Radius);
	if (LoCol < 0)
		LoCol = 0;
	if (HiCol >= Lali)
		HiCol = Lali - 1;
	for (int Col2 = LoCol; Col2 <= HiCol; ++Col2)
		{
		if ((int) Col == Col2)
			Sum += Theta;
		else
			{
			uint PosQj = PosQs[Col2];
			uint PosTj = PosTs[Col2];

			double dij_Q = Q.GetDist(PosQi, PosQj);
			double dij_T = T.GetDist(PosTi, PosTj);
			double x = dpscorefun(dij_Q, dij_T);
			Sum += x;
			}
		}
	int n = (HiCol - LoCol + 1);
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

static double GetDALIScore_Col(
  const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Col)
	{
	const uint Lali = SIZE(PosQs);
	asserta(SIZE(PosTs) == Lali);

	const double Theta = 0.2;
	double Sum = 0;
	uint PosQi = PosQs[Col];
	uint PosTi = PosTs[Col];
	for (uint Col2 = 0; Col2 < Lali; ++Col2)
		{
		if (Col == Col2)
			Sum += Theta;
		else
			{
			uint PosQj = PosQs[Col2];
			uint PosTj = PosTs[Col2];

			double dij_Q = Q.GetDist(PosQi, PosQj);
			double dij_T = T.GetDist(PosTi, PosTj);
			double x = dpscorefun(dij_Q, dij_T);
			Sum += x;
			}
		}
	return Sum;
	}

double GetDALIScore_Cols_Band(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Radius, vector<double> &ColScores)
	{
	ColScores.clear();
	const uint Lali = SIZE(PosQs);
	asserta(SIZE(PosTs) == Lali);

	double Sum = 0;
	for (uint Col = 0; Col < Lali; ++Col)
		{
		double ColScore = GetDALIScore_Col_Band(Q, T, PosQs, PosTs, Col, Radius);
		Sum += ColScore;
		ColScores.push_back(ColScore);
		}
	return Sum;
	}

double GetDALIScore_LSB_Pair(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Lo1, uint Hi1, uint Lo2, uint Hi2)
	{
	double Sum = 0;
	uint n = 0;
	asserta(SIZE(PosQs) == SIZE(PosTs));
	asserta(Hi1 < SIZE(PosQs));
	asserta(Hi2 < SIZE(PosQs));
	for (uint MatchCol1 = Lo1; MatchCol1 <= Hi1; ++MatchCol1)
		{
		uint PosQ1 = PosQs[MatchCol1];
		uint PosT1 = PosTs[MatchCol1];
		for (uint MatchCol2 = Lo2; MatchCol2 <= Hi2; ++MatchCol2)
			{
			uint PosQ2 = PosQs[MatchCol2];
			uint PosT2 = PosTs[MatchCol2];
			double dij_Q = Q.GetDist(PosQ1, PosQ2);
			double dij_T = T.GetDist(PosT1, PosT2);
			double x = dpscorefun(dij_Q, dij_T);
			Sum += x;
			++n;
			}
		}
	return Sum/n;
	}

double GetDALIScore_Cols(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  vector<double> &ColScores)
	{
	ColScores.clear();
	const uint Lali = SIZE(PosQs);
	asserta(SIZE(PosTs) == Lali);

	double Sum = 0;
	for (uint i = 0; i < Lali; ++i)
		{
		double ColScore = GetDALIScore_Col(Q, T, PosQs, PosTs, i);
		Sum += ColScore;
		ColScores.push_back(ColScore);
		}
	return Sum;
	}

static double GetDALIScore(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint Lali = SIZE(PosQs);
	asserta(SIZE(PosTs) == Lali);

	const double Theta = 0.2;
	double Sum = 0;
	for (uint i = 0; i < Lali; ++i)
		{
		uint PosQi = PosQs[i];
		uint PosTi = PosTs[i];
		for (uint j = 0; j < Lali; ++j)
			{
			if (i == j)
				Sum += Theta;
			else
				{
				uint PosQj = PosQs[j];
				uint PosTj = PosTs[j];

				double dij_Q = Q.GetDist(PosQi, PosQj);
				double dij_T = T.GetDist(PosTi, PosTj);
				double x = dpscorefun(dij_Q, dij_T);
				Sum += x;
				}
			}
		}
#if DEBUG
	{
	vector<double> ColScores;
	double Sum2 = GetDALIScore_Cols(Q, T, PosQs, PosTs, ColScores);
	asserta(feq(Sum2, Sum));
	}
#endif
	return Sum;
	}

/***
"./src/comparemodules.f" line 1473
===========================
        function zscore(l1,l2,score) result(z)
        implicit none
        real z,score
        integer l1,l2
c
        real n12,mean,sigma,x
c
        n12=sqrt(float(l1*l2))
        x=min(n12,400.0)
        mean=7.9494+0.70852*x+2.5895e-4*x*x-1.9156e-6*x*x*x
        if(n12.gt.400.0) mean=mean+(n12-400.0)*1.0              ! hack !
        sigma=0.50*mean
        z=(score-mean)/max(1.0,sigma)

        return
        end function zscore

***/

static double GetDALIZFromScoreAndLengths(double DALIScore, uint QL, uint TL)
	{
	double n12 = sqrt(QL*TL);
	double x = min(n12, 400.0);
	double mean = 7.9494 + 0.70852*x + 2.5895e-4*x*x - 1.9156e-6*x*x*x;
	if (n12 > 400)
		mean += n12 - 400.0;
	double sigma = 0.5*mean;
	double z = (DALIScore - mean)/max(1.0, sigma);
	return z;
	}

void GetPosVecs(const string &QRow, const string &TRow,
  vector<uint> &MatchColToPosQs, vector<uint> &MatchColToPosTs,
  vector<uint> &MatchColToCol, vector<uint> &ColToMatchCol)
	{
	MatchColToPosQs.clear();
	MatchColToPosTs.clear();
	MatchColToCol.clear();
	ColToMatchCol.clear();
	uint ColCount = SIZE(QRow);
	asserta(SIZE(TRow) == ColCount);
	uint PosQ = 0;
	uint PosT = 0;
	uint MatchColCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = QRow[Col];
		char t = TRow[Col];

		if (q != '-' && t != '-' &&
		  isupper(q) && isupper(t))
			{
			MatchColToPosQs.push_back(PosQ);
			MatchColToPosTs.push_back(PosT);
			MatchColToCol.push_back(Col);
			ColToMatchCol.push_back(MatchColCount);
			++MatchColCount;
			}
		else
			ColToMatchCol.push_back(UINT_MAX);
		if (q != '-')
			++PosQ;
		if (t != '-')
			++PosT;
		}
	}

double GetDALIZ_PosVecs(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint QL = Q.GetSeqLength();
	const uint TL = T.GetSeqLength();

	double DALI = GetDALIScore(Q, T, PosQs, PosTs);
	double Z = GetDALIZFromScoreAndLengths(DALI, QL, TL);
	return Z;
	}

double GetDALIZ(const PDBChain &Q, const PDBChain &T,
  const string &QRow, const string &TRow)
	{
	const uint QL = Q.GetSeqLength();
	const uint TL = T.GetSeqLength();

	vector<uint> PosQs;
	vector<uint> PosTs;
	vector<uint> MatchColToCol;
	vector<uint> ColToMatchCol;
	GetPosVecs(QRow, TRow, PosQs, PosTs, MatchColToCol, ColToMatchCol);
	double DALI = GetDALIScore(Q, T, PosQs, PosTs);
	double Z = GetDALIZFromScoreAndLengths(DALI, QL, TL);
	return Z;
	}

/***
palmscan2
  -daliz d:/a/res/daliz/data/palmprints.cal \
  -ref d:/a/res/daliz/data/palmprints.tsv \
  -fieldnrs 0,1,7,8 \
  -tsv daliz.tsv

palmprints.tsv from dali2tsv.py:
   0       1       2      3   4   5  6  7 .. 8
   Q       T       Z                    QRow .. TRow
1hhs    4ieg    14.5    1.8 124 128 19  VATDVSDHDTFWPGW...
***/

void cmd_daliz()
	{
	const string &QFN = opt_daliz;
	const string &TsvFN = opt_ref;

	vector<PDBChain *> Chains;
	ReadChains(QFN, Chains);

	const uint NQ = SIZE(Chains);
	map<string, uint> LabelToIndex;
	for (uint i = 0; i < NQ; ++i)
		{
		const string &Label = Chains[i]->m_Label;
		LabelToIndex[Label] = i;
		}

	vector<string> Fields;
	Split(opt_fieldnrs, Fields, ',');
	assert(SIZE(Fields) == 4);

	uint Qfn = StrToUint(Fields[0]);
	uint Tfn = StrToUint(Fields[1]);
	uint QRowfn = StrToUint(Fields[2]);
	uint TRowfn = StrToUint(Fields[3]);

	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');

		const string &Query = Fields[Qfn];
		const string &Target = Fields[Tfn];
		map<string, uint>::const_iterator pQ = LabelToIndex.find(Query);
		map<string, uint>::const_iterator pT = LabelToIndex.find(Target);
		if (pQ == LabelToIndex.end() || pT == LabelToIndex.end())
			continue;

		uint iQ = pQ->second;
		uint iT = pT->second;

		const PDBChain &Q = *Chains[iQ];
		const PDBChain &T = *Chains[iT];
		asserta(Q.m_Label == Query && T.m_Label == Target);

		const string &QRow = Fields[QRowfn];
		const string &TRow = Fields[TRowfn];

		double Z = GetDALIZ(Q, T, QRow, TRow);
		if (g_ftsv)
			{
			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			fprintf(g_ftsv, "\t%s", T.m_Label.c_str());
			fprintf(g_ftsv, "\t%.3g", Z);
			fprintf(g_ftsv, "\n");
			}
		}
	CloseStdioFile(f);
	}
