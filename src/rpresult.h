#pragma once

class RPResult
	{
public:
	string m_QueryLabel;
	uint m_QL;
	uint m_QLnt;
	string m_Gene;
	bool m_HiConf;
	string m_Comments;
	string m_XXX;
	int m_Frame;
	float m_PSSMTotalScore;
	float m_PSSMMinScore;
	uint m_StartPos;
	uint m_SegLength;
	uint m_StartPosNt;
	uint m_SegLengthNt;
	int m_Dist12;
	int m_Dist23;
	double m_FinalScore;
	string m_TopGroupName;
	string m_TrimmedSeq;
	string m_TrimmedSeqNt;
	string m_MotifsSeq;
	string m_MotifsSeq2;
	string m_SuperMotif;
	vector<string> m_Aln;
	vector<string> m_FullAln;

public:
	void Clear()
		{
		m_QueryLabel.clear();
		m_QL = UINT_MAX;
		m_QLnt = UINT_MAX;
		m_Gene.clear();
		m_HiConf = false;
		m_Comments.clear();
		m_XXX.clear();

		m_StartPos = UINT_MAX;
		m_StartPosNt = UINT_MAX;
		m_SegLength = UINT_MAX;
		m_SegLengthNt = UINT_MAX;
		m_Frame = -999;

		m_PSSMTotalScore = -999;
		m_PSSMMinScore = -999;
		m_Dist12 = -999;
		m_Dist23 = -999;
		m_SegLength = UINT_MAX;
		m_TopGroupName.clear();
		m_TrimmedSeq.clear();
		m_TrimmedSeqNt.clear();
		m_MotifsSeq.clear();
		m_MotifsSeq2.clear();
		m_SuperMotif.clear();
		m_Aln.clear();
		m_FullAln.clear();
		m_FinalScore = -999;
		}

	void ToFEVStr(string &s) const;
	void GetCat(string &Cat) const
		{
		if (m_Gene == "unclassified")
			Cat = "unclassified";
		else
			{
			if (m_HiConf)
				Cat = "high-confidence-" + m_Gene;
			else
				Cat = "low-confidence-" + m_Gene;
			}		
		}
	};
