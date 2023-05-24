#include "myutils.h"
#include "rdrpsearcher.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "omplock.h"
#include "abcxyz.h"
#include "outputfiles.h"

extern vector<string> g_ExcludeNames;
void SetExcludes();

static uint g_QueryCount;
static uint g_FoundCount;

void RdRpSearcher::MergeCFilterHits(const vector<RPHit> &Hits)
	{
	const uint N = SIZE(Hits);
	for (uint i = 0; i < N; ++i)
		{
		const RPHit &NewHit = Hits[i];
		for (uint j = 0; j < SIZE(m_CFilterHits); ++j)
			{
			RPHit &OldHit = m_CFilterHits[j];
			if (NewHit.m_QPos == OldHit.m_QPos)
				{
				if (NewHit.m_Score > OldHit.m_Score)
					OldHit = NewHit;
				return;
				}
			}
		m_CFilterHits.push_back(NewHit);
		}
	}

void RdRpSearcher::WriteCFilterHit(const RPHit &Hit) const
	{
	const uint QL = SIZE(m_QuerySeq);
	const uint QPos = Hit.m_QPos;
	asserta(QPos + CL <= QL);
	int iLo = int(QPos) - int(m_LeftFlankCFilter);
	if (iLo < 0)
		iLo = 0;
	uint Lo = uint(iLo);
	uint Hi = Hit.m_QPos + CL + m_LeftFlankCFilter;
	if (Hi >= QL)
		Hi = QL - 1;
	asserta(Hi > Lo);
	uint Len = Hi - Lo + 1;

	string CSeq;
	for (uint i = 0; i < CL; ++i)
		CSeq += m_QuerySeq[QPos + i];
	const char *GroupName = Hit.m_PSSM->m_GroupName.c_str();
	const char *ConsSeq = Hit.m_PSSM->m_ConsSeq.c_str();

	if (g_ftsv != 0)
		{
		fprintf(g_ftsv, "%s", m_QueryLabel.c_str());
		fprintf(g_ftsv, "\t%.1f", Hit.m_Score);
		fprintf(g_ftsv, "\t%u", Lo + 1);
		fprintf(g_ftsv, "\t%u", Hit.m_QPos + 1);
		fprintf(g_ftsv, "\t%u", Hi + 1);
		fprintf(g_ftsv, "\t%s", CSeq.c_str());
		fprintf(g_ftsv, "\t%s", ConsSeq);
		fprintf(g_ftsv, "\t%s", GroupName);
		fprintf(g_ftsv, "\n");
		}

	string OutLabel;
	Ps(OutLabel, "%s [%u-%u]", m_QueryLabel, Lo + 1, Hi + 1);
	SeqToFasta(g_ffasta, OutLabel.c_str(), m_QuerySeq.c_str() + 1, Len);
	}

bool RdRpSearcher::WriteCFilterOutput() const
	{
	if (m_CFilterHits.empty())
		return false;
	bool Any = false;
#pragma omp critical
	{
	for (uint i = 0; i < SIZE(m_CFilterHits); ++i)
		{
		const RPHit &Hit = m_CFilterHits[i];
		if (Hit.m_Score >= m_MinScore_CFiler)
			{
			Any = true;
			WriteCFilterHit(Hit);
			}
		}
	}
	return Any;
	}

void cmd_cfilter()
	{
	const string &QueryFileName = opt_cfilter;

	asserta(optset_model);
	const string &ModelFileName = opt_model;
	SetExcludes();

	if (!opt_notrunclabels)
		opt_trunclabels = true;

	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	const uint ThreadCount = GetRequestedThreadCount();
	vector<ObjMgr *> OMs;
	vector<RdRpSearcher *> Mods;
	for (int i = 0; i < int(ThreadCount); ++i)
		{
		ObjMgr *OM = new ObjMgr;
		OMs.push_back(OM);
		RdRpSearcher *RS = new RdRpSearcher;
		RS->Init(Model);
		Mods.push_back(RS);
		}

	FASTASeqSource *SS = new FASTASeqSource;
	SS->Open(QueryFileName);

	//RdRpSearcher::InitOutput();
	uint LastElapsedSecs = 0;
	uint CurrElapsedSecs = 0;
	ProgressStep(0, 1000, "Searching");
#pragma omp parallel num_threads(ThreadCount)
	{
	int ThreadIndex = GetThreadIndex();
	RdRpSearcher &RS = *Mods[ThreadIndex];
	ObjMgr *OM = OMs[ThreadIndex];
	for (;;)
		{
		SeqInfo *QSI = OM->GetSeqInfo();

		if (ThreadIndex == 0)
			{
			if (g_QueryCount%100 == 0)
				{
				CurrElapsedSecs = GetElapsedSecs();
				if (CurrElapsedSecs > LastElapsedSecs)
					{
					uint Pct10 = SS->GetPctDoneX10();
					double HitPct = GetPct(g_FoundCount, g_QueryCount);
					ProgressStep(Pct10, 1000, "Searching %u / %u hits (%.1f%%)",
					  g_FoundCount, g_QueryCount, HitPct);
					LastElapsedSecs = CurrElapsedSecs;
					}
				}
			}

		bool Ok = SS->GetNext(QSI);
		if (!Ok)
			break;

		const string Label = string(QSI->m_Label);
		string Seq;
		for (uint i = 0; i < QSI->m_L; ++i)
			Seq += char(toupper(QSI->m_Seq[i]));

		RS.CFilter(Label, Seq);
		bool HasHit = RS.WriteCFilterOutput();
		OM->Down(QSI);

		Lock();
		++g_QueryCount;
		if (HasHit)
			++g_FoundCount;
		Unlock();
		}
	}

	double HitPct = GetPct(g_FoundCount, g_QueryCount);
	ProgressStep(999, 1000, "Searching %u/%u hits (%.1f%%)",
	  g_FoundCount, g_QueryCount, HitPct);
	}
