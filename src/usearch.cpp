#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "viterbi.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "usearch.h"
#include "alpha.h"

static XDPMem *g_Mem;
static ObjMgr *g_OM;
static SeqInfo *g_QSI;
static SeqInfo *g_TSI;
static PathInfo *g_PI;
static AlnParams g_AP;
static AlnHeuristics g_AH;
static bool g_Nucleo;
static const byte *g_CharToLetter;
static unsigned g_AlphaSize;

static void Init()
	{
	if (g_OM != 0)
		return;

	g_OM = new ObjMgr;
	g_Mem = new XDPMem;
	g_QSI = g_OM->GetSeqInfo();
	g_TSI = g_OM->GetSeqInfo();
	g_PI = g_OM->GetPathInfo();

	g_AP.InitFromCmdLine(true);
	g_AH.InitFromCmdLine(g_AP);

	g_Nucleo = g_AP.GetIsNucleo();
	if (g_Nucleo)
		{
		g_AlphaSize = 4;
		g_CharToLetter = g_CharToLetterNucleo;
		}
	else
		{
		g_AlphaSize = 20;
		g_CharToLetter = g_CharToLetterAmino;
		}
	}

static void GetAlnStats(const byte *QSeq, unsigned QL, const byte *TSeq, unsigned TL,
  const char *Path, unsigned &ColCount, unsigned &IdCount, unsigned &GapCount)
	{
	IdCount = 0;
	GapCount = 0;
	const unsigned N = (unsigned) strlen(Path);
	unsigned QPos = 0;
	unsigned TPos = 0;
	unsigned FirstM = UINT_MAX;
	unsigned LastM = UINT_MAX;
	for (unsigned Col = 0; Col < N; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = Col;
			LastM = Col;
			}
		}

	asserta(FirstM != UINT_MAX);
	ColCount = LastM - FirstM + 1;
	for (unsigned Col = 0; Col < N; ++Col)
		{
		char c = Path[Col];

		if (c == 'M')
			{
			byte q = QSeq[QPos];
			byte t = TSeq[TPos];
			byte qlet = g_CharToLetter[q];
			byte tlet = g_CharToLetter[t];
			if (qlet == tlet && qlet < g_AlphaSize)
				++IdCount;
			}
		else
			{
			asserta(c == 'D' || c == 'I');
			++GapCount;
			}
		
		if (c == 'M' || c == 'D')
			++QPos;
		if (c == 'M' || c == 'I')
			++TPos;
		}
	asserta(QPos == QL);
	asserta(TPos == TL);
	}

void Usearch(const char *QLabel, const byte *QSeq, unsigned QL,
  SeqDB &DB, vector<UsearchHit> &Hits)
	{
	Hits.clear();

	Init();

	g_QSI->m_Label = QLabel;
	g_QSI->m_Seq = QSeq;
	g_QSI->m_L = QL;

	vector<unsigned> WordCounts;
	vector<unsigned> Order;
	USort(*g_QSI, DB, WordCounts, Order);
	unsigned BestIdCount = UINT_MAX;
	const unsigned N = SIZE(WordCounts);
	UsearchHit Hit;
	Hit.QLabel = QLabel;
	Hit.QSeq = QSeq;
	Hit.QL = QL;
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];

		const char *TLabel = DB.GetLabel(i).c_str();
		const byte *TSeq = (const byte *) DB.GetSeq(i).c_str();
		unsigned TL = DB.GetSeqLength(i);

		ViterbiFastMainDiagMem(*g_Mem, QSeq, QL, TSeq, TL, g_AH.BandRadius, g_AP, *g_PI);

		const char *Path = g_PI->m_Path;
		unsigned ColCount;
		unsigned IdCount;
		unsigned GapCount;
		GetAlnStats(QSeq, QL, TSeq, TL, Path, ColCount, IdCount, GapCount);

		if (k == 0)
			BestIdCount = IdCount;
		else
			{
			if (IdCount > BestIdCount)
				BestIdCount = IdCount;
			else if (IdCount + opt_iddrop < BestIdCount)
				return;
			}

		Hit.TLabel = TLabel;
		Hit.TSeq = TSeq;
		Hit.TL = TL;
		
		Hit.Path = string(Path);

		Hit.ColCount = ColCount;
		Hit.IdCount = IdCount;
		Hit.GapCount = GapCount;

		Hits.push_back(Hit);
		}
	}
