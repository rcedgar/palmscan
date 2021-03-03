#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "rdrpmodel.h"
#include "seqdb.h"
#include "pssmsearch.h"
#include "pathinfo.h"
#include "objmgr.h"

static void Search1(RdRpModel &Mod, const string &Label,
  const string &Seq)
	{
	ObjMgr OM;
	const uint PSSMCount = SIZE(Mod.m_PSSMs);
	for (uint PSSMIndex = 0; PSSMIndex < PSSMCount; ++PSSMIndex)
		{
		const string &Name = Mod.m_PSSMGroupNames[PSSMIndex];
		const char MotifLetter = Mod.m_PSSMMotifLetters[PSSMIndex];

		RPHit HitU;
		RPHit HitG;
		Mod.SearchAA1_Ungapped(PSSMIndex, Seq, HitU);
		Mod.SearchAA1_Gapped(PSSMIndex, Seq, HitG);

		Log(">%s", Label.c_str());
		Log("  %s", Name.c_str());
		Log(".%c", MotifLetter);
		Log("  %.4f", HitU.m_Score);
		Log("  %.4f", HitG.m_Score);
		Log("  %u", HitU.m_QPos);
		Log("  %u", HitG.m_QPos);
		Log("\n");

		Mod.LogAlnViterbi(PSSMIndex, Seq, HitG);
		}
	}

static void Search2(RdRpModel &Mod, const string &QueryLabel,
  const string &QuerySeq)
	{
	Mod.SearchAA(QueryLabel, QuerySeq);
	Mod.LogHitTable();
	}

void Test_Viterbi_PSSM()
	{
	const string &QueryFileName = opt_test;
	const string &ModelFileName = opt_model;

	RdRpModel Mod;
	Mod.FromModelFile(ModelFileName);

	SeqDB Input;
	Input.FromFasta(QueryFileName);

	const uint SeqCount = Input.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Searching");
		const string Label = string(Input.GetLabel(SeqIndex));
		const string &Seq = Input.GetSeq(SeqIndex);
		Search1(Mod, Label, Seq);
		}
	}
