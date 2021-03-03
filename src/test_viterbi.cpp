#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "viterbi.h"
#include "objmgr.h"
#include "seqinfo.h"

void Test_Viterbi()
	{
	const string &FileName = string(opt_input);

	SeqDB Input;
	Input.FromFasta(FileName);
	asserta(Input.GetSeqCount() >= 1);

	XDPMem Mem;

	AlnParams AP;
	AlnHeuristics AH;

	AP.InitFromCmdLine(true);
	AH.InitFromCmdLine(AP);

	ObjMgr *OM = new ObjMgr;
	PathInfo *PI = OM->GetPathInfo();

	const byte *Q = (const byte *) Input.GetSeq(0).c_str();
	unsigned QL = Input.GetSeqLength(0);

	SeqDB DB;
	DB.FromFasta(opt_db);

	SeqInfo *QSI = OM->GetSeqInfo();
	SeqInfo *TSI = OM->GetSeqInfo();

	QSI->m_Seq = Q;
	QSI->m_L = QL;
	QSI->m_Label = Input.GetLabel(0).c_str();

	vector<unsigned> WordCounts;
	vector<unsigned> Order;
	USort(*QSI, DB, WordCounts, Order);

	Log("\n");
	Log("  Words  Label\n");
	Log("-------  -----\n");
	const unsigned N = SIZE(WordCounts);
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		Log("%7u  %s\n", WordCounts[i], DB.GetLabel(i).c_str());
		}

	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];

		const byte *T = (const byte *) DB.GetSeq(i).c_str();
		unsigned TL = DB.GetSeqLength(i);
		const char *TLabel = DB.GetLabel(i).c_str();

		ViterbiFastMainDiagMem(Mem, Q, QL, T, TL, AH.BandRadius, AP, *PI);

		Log("\n");
		Log(">%s\n", TLabel);
		unsigned ColCount = (unsigned) strlen(PI->m_Path);
		LogAlnPretty(Q, T, PI->m_Path, true);
		}

	OM->Down(PI);
	}
