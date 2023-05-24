#pragma once

#include "rdrpmodel.h"
#include "palmhit.h"

#define MotLet "ABC"[MotifIndex]

class RdRpSearcher
	{
public:
	const RdRpModel *m_Model = 0;

	string m_QueryLabel;
	string m_QuerySeq;

	PalmHit m_TopPalmHit;
	PalmHit m_SecondPalmHit;

	uint m_MaxX = 10;
	float m_MinCScore = 3.0;
	float m_MinPSSMScore = 2.0;

	uint m_MinPPLength = 90;
	uint m_MaxPPLength = 180;

	uint m_MinInsert = 10;

	bool m_Trace = false;

	double m_MinScore_Core = 20.0;
	uint m_LeftFlank_Core = 149;	// see also MPCluster
	uint m_RightFlank_Core = 150;	// see also MPCluster
	uint m_MinFlankLen = 0;

	double m_MinScore_CFiler = 5.0;
	uint m_LeftFlankCFilter = 150;
	uint m_RightFlankCFilter = 300;

	vector<RPHit> m_CFilterHits;

public:
	RdRpSearcher()
		{
		m_MaxX = 10;
		if (optset_maxx)
			m_MaxX = opt_maxx;
		if (optset_mincscore)
			m_MinCScore = (float) opt_mincscore;
		if (optset_minflanklen)
			m_MinFlankLen = opt_minflanklen;
		asserta(m_MinCScore > 0);
		}

	void Clear()
		{
		m_QueryLabel.clear();
		m_QuerySeq.clear();
		m_CFilterHits.clear();
		m_TopPalmHit.Clear();
		m_SecondPalmHit.Clear();
		}

	void ClearSearch()
		{
		m_QueryLabel.clear();
		m_QuerySeq.clear();
		m_CFilterHits.clear();
		m_TopPalmHit.Clear();
		m_SecondPalmHit.Clear();
		}

	uint GetGroupCount() const
		{
		asserta(m_Model != 0);
		return m_Model->GetGroupCount();
		}

	void GetGroupName(uint GroupIndex, string &Name) const
		{
		asserta(m_Model != 0);
		return m_Model->GetGroupName(GroupIndex, Name);
		}

public:
	bool IsHit() const;
	void Init(const RdRpModel &Model);
	void Search(const string &QueryLabel, const string &QuerySeq);
	void SearchGroup(uint GroupIndex);
	void SearchMotif(uint GroupIndex, uint MotifIndex,
	  uint QLo, uint QHi, float MinScore, vector<RPHit> &Hits);
	void SearchMotif_TopHitOnly(uint GroupIndex, uint MotifIndex,
	  uint QLo, uint QHi, float MinScore, RPHit &Hit);
	void SearchAB(uint GroupIndex, const RPHit &Hit_C);
	bool SearchAB_NotPermuted(uint GroupIndex, const RPHit &Hit_C);
	bool SearchAB_Permuted(uint GroupIndex, const RPHit &Hit_C);
	void OnPalmHit(uint GroupIndex, const RPHit &Hit_A, const RPHit &Hit_B,
	  const RPHit &Hit_C, bool Permuted);
	void GetAlnRows(vector<string> &Rows) const;
	void GetMotifsSeq(const string &Sep, string &Seq) const;
	uint GetMotifPos(uint MotifIndex) const;
	void GetTrimmedSeq(string &Seq) const;
	void GetMotifPositions(uint &APos, uint &BPos, uint &CPos) const;
	void GetAln(const RPHit &Hit, string &Q, string &P, string &Annot,
	  uint &QLo, uint &QHi) const;
	uint GetPSSMLength(uint GroupIndex, uint MotifIndex) const;
	void GetSpan(uint &QLo, uint &QHi) const;
	void GetMotifSeq(uint MotifIndex, string &Seq) const;
	void WriteOutput() const;
	void WriteReport(FILE *f) const;
	void WriteTsv(FILE *f) const;
	void WriteFev(FILE *f) const;
	void WriteMotifs(FILE *fABC, FILE *fCAB) const;
	void WritePalmprintFasta(FILE *fABC, FILE *fCAB) const;
	void WriteFlankFasta(FILE *fABC, FILE *fCAB, bool Left) const;
	bool WriteCore(FILE *fABC, FILE *fCAB) const;
	void GetCoreSeq(string &Seq) const;

	void CFilter(const string &QueryLabel, const string &QuerySeq);
	void CFilterGroup(uint GroupIndex);
	void MergeCFilterHits(const vector<RPHit> &Hits);
	bool WriteCFilterOutput() const;
	void WriteCFilterHit(const RPHit &Hit) const;

	const PSSM &GetPSSM(uint GroupIndex, uint MotifIndex) const;

public:
	//static void InitOutput();
	};

typedef void (*fn_OnPalmHit)(const RdRpSearcher &RS);
void SearchPSSMs(const string &QueryFileName, fn_OnPalmHit OnHit);
