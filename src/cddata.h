#pragma once

#include "cdinfo.h"

class CDData
	{
public:
	const CDInfo *m_Info;
	vector<vector<double> > m_Data;

public:
	void Clear()
		{
		m_Info = 0;
		m_Data.clear();
		}

	void ToTsv(const string &Name, FILE *f) const;
	void FromTsv(const CDInfo &Info, FILE *f);
	void Init(const CDInfo &Info);
	double Get(uint MotifIndex1, uint i1,
	  uint MotifIndex2, uint i2) const;
	void Set(uint MotifIndex1, uint i1,
	  uint MotifIndex2, uint i2, double Value);
	double GetByIx(uint Ix1, uint Ix2) const;
	void SetByIx(uint Ix1, uint Ix2, double Value);
	};
