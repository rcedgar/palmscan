#pragma once

#include "pssm.h"
#include "pssmsearch.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"
#include "rphit.h"
#include "rpresult.h"
#include "heuristics.h"

class RdRpModel
	{
public:
	int m_Thread;
	string m_QueryLabel;
	string m_QuerySeq;

	vector<string> m_GroupNames;
	vector<char> m_MotifLetters;

	vector<PSSM> m_PSSMs;
	vector<string> m_PSSMGroupNames;
	vector<char> m_PSSMMotifLetters;
	vector<RPHit> m_HitsU;
	vector<RPHit> m_HitsG;

	vector<vector<uint> > m_PIV;

	bool m_Gapped;
	ObjMgr m_OM;
	AlnParams m_AP;
	XDPMem m_DPMem;

	uint m_MaxX;

	RPResult m_Result;

public:
	RdRpModel()
		{
		m_Gapped = opt_gapped;
		m_MaxX = 10;
		if (optset_maxx)
			m_MaxX = opt_maxx;
		float GapOpen = -10;
		float GapExt = -1;
		float TermGapOpen = 0;
		float TermGapExt = 0;
		m_AP.Init4(0, GapOpen, GapExt, TermGapOpen, TermGapExt);
		m_AP.LOpenA = -999;
		m_AP.ROpenA = -999;
		}

	void Clear()
		{
		m_PSSMs.clear();
		m_PSSMGroupNames.clear();
		m_PSSMMotifLetters.clear();
		}

	const PSSM &GetGetPSSM(uint Index) const
		{
		asserta(Index < SIZE(m_PSSMs));
		return m_PSSMs[Index];
		}

	uint GetMotifCount() const
		{
		return SIZE(m_MotifLetters);
		}

	uint GetGroupCount() const
		{
		return SIZE(m_GroupNames);
		}

	uint GetPSSMCount() const
		{
		return SIZE(m_PSSMs);
		}

	const PSSM &GetPSSM(uint PSSMIndex) const
		{
		asserta(PSSMIndex < SIZE(m_PSSMs));
		return m_PSSMs[PSSMIndex];
		}

	uint GetPSSMIndex(uint MotifIndex, uint GroupIndex) const
		{
		asserta(MotifIndex < SIZE(m_PIV));
		asserta(GroupIndex < SIZE(m_PIV[0]));
		uint PSSMIndex = m_PIV[MotifIndex][GroupIndex];
		return PSSMIndex;
		}

	const PSSM *GetPSSM(uint MotifIndex, uint GroupIndex) const
		{
		uint PSSMIndex = GetPSSMIndex(MotifIndex, GroupIndex);
		if (PSSMIndex == UINT_MAX)
			return 0;
		asserta(PSSMIndex < SIZE(m_PSSMs));
		return &m_PSSMs[PSSMIndex];
		}

	void GetGroupName(uint GroupIndex, string &GroupName) const
		{
		if (GroupIndex >= SIZE(m_GroupNames))
			{
			GroupName = "*";
			return;
			}
		GroupName = m_GroupNames[GroupIndex];
		}

	uint GetPSSMLength(char MotifLetter, uint GroupIndex) const
		{
		return GetPSSMLength(GetMotifIndex(MotifLetter), GroupIndex);
		}

	const RPHit *GetHit(char MotifLetter, uint GroupIndex) const
		{
		return GetHit(GetMotifIndex(MotifLetter), GroupIndex);
		}

	uint GetTopGroup() const
		{
		uint TopIx;
		float Score;
		GetTopGroup(TopIx, Score);
		return TopIx;
		}

	const RPHit *GetHit(char MotifLetter) const
		{
		uint GroupIndex = GetTopGroup();
		if (GroupIndex == UINT_MAX)
			return 0;
		return GetHit(MotifLetter, GroupIndex);
		}

	uint GetPSSMLength(char MotifLetter) const
		{
		uint GroupIndex = GetTopGroup();
		if (GroupIndex == UINT_MAX)
			return 0;
		return GetPSSMLength(MotifLetter, GroupIndex);
		}

public:
	void FromSpecFile(const string &FileName);
	void FromModelFile(const string &FileName);
	void FromStrings(const vector<string> &Lines);
	void ToModelFile(const string &FileName) const;

	void SearchAA(const string &QueryLabel, const string &QuerySeq);
	void SearchAA_Seeded(const string &QueryLabel, const string &QuerySeq);

public:
	void SetPIV();
	uint GetMotifIndex(char MotifLetter) const;
	uint GetGroupIndex(const string &GroupName) const;
	uint GetPSSMIndex(char MotifLetter, const string &GroupName) const;
	uint GetPSSMLength(uint MotifIndex, uint GroupIndex) const;
	void SearchAA1(uint PSSMIndex, const string &Seq, RPHit &Hit);
	void SearchAA1_Ungapped(uint PSSMIndex, const string &Seq, RPHit &Hit);
	void SearchAA1_Gapped(uint PSSMIndex, const string &Seq, RPHit &Hit);
	void LogAlnViterbi(uint PSSMIndex, const string &Seq, const RPHit &Hit) const;
	void GetAln(uint MotifIndex, uint GroupIndex,
	  string &Q, string &P, string &Annot,
	  uint &QLo, uint &QHi) const;
	void GetXSeq(char MotifLetter, string &XSeq) const;
	int GetDist(const RPHit *Hit1, const RPHit *Hit2) const;
	int GetDist(char MotifLetter1, char MotifLetter2) const;
	uint GetSpan(const RPHit *Hit1, const RPHit *Hit2) const;

	void LogHitTable() const;
	float GetTotalScore(uint GroupIndex) const;
	float GetMinScore() const;
	uint GetTotalGaps() const;
	void GetTopGroup(uint &TopIx, float &TopScore) const;
	const RPHit *GetHit(uint MotifIndex, uint GroupIndex) const;
	void GetABCOrder(string &XXX) const;
	void GetOrderedHits(string &XXX, const RPHit *&Hit1,
	  const RPHit *&Hit2, const RPHit *&Hit3) const;
	bool GetABCRange(uint &Start, uint &End) const;
	void GetAlnRows(vector<string> &Rows) const;
	void GetFullAln(vector<string> &Rows) const;
	void GetMotifsSeq(string &Seq) const;
	void GetMotifsSeq2(string &Seq) const;
	void GetTrimmedSeq(string &Seq) const;
	float GetCScore() const;
	void GetSuperMotif(const string &MotifsSeq, string &s) const;

public:
	static uint GetNextSeedPos(const string &QuerySeq, uint SeedPos);
	static void SetResult_NoNtHit(const string &QueryLabel, uint QLnt,
	  RPResult &Result);
	static void SetResult_NoAaHit(const string &QueryLabel, uint QLnt,
	  RPResult &Result);

private:
	void SetResult();
	void CheckThread() const { asserta(omp_get_thread_num() == m_Thread); }
	};
