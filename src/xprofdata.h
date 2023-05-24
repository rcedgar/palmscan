#pragma once

#include "pdbchain.h"

class XProfData : public PDBChain
	{
public:
	//string m_Label;
	//string m_Seq;
	//vector<double> m_Xs;
	//vector<double> m_Ys;
	//vector<double> m_Zs;
	vector<vector<double> > m_PosToFeatureVec;
	vector<vector<uint> > m_PosToIntFeatureVec;

public:
	void Clear()
		{
		PDBChain::Clear();
		//m_Label.clear();
		//m_Seq.clear();
		//m_Xs.clear();
		//m_Ys.clear();
		//m_Zs.clear();
		m_PosToFeatureVec.clear();
		m_PosToIntFeatureVec.clear();
		}

	bool FromCfv(FILE *f);
	//uint GetSeqLength() const { return SIZE(m_Seq); }
	uint GetIntFeature(uint FeatureIndex, uint Pos) const
		{
		assert(Pos < SIZE(m_PosToIntFeatureVec));
		assert(FeatureIndex < SIZE(m_PosToIntFeatureVec[Pos]));
		return m_PosToIntFeatureVec[Pos][FeatureIndex];
		}
	void Shuffle();
	};
