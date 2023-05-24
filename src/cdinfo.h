#pragma once

// Contact Distance info -- shape of contact
//   map distance data
class CDInfo
	{
public:
	uint m_MotifCount = 0;
	vector<string> m_MotifNames;
	vector<uint> m_MotifLengths;
	vector<uint> m_Offsets;
	vector<uint> m_MotifIndexes;
	vector<uint> m_MotifCoords;

public:
	void Clear()
		{
		m_MotifCount = 0;
		m_MotifNames.clear();
		m_MotifLengths.clear();
		m_Offsets.clear();
		m_MotifIndexes.clear();
		m_MotifCoords.clear();
		}

	void Init(const vector<string> &MotifNames,
	  vector<uint> &MotifLengths);
	void LogMe() const;
	uint GetIx(uint MotifIndex, uint i) const;
	void GetCoord(uint Ix, uint &MotifIndex, uint &i) const;
	void ToTsv(FILE *f) const;
	void FromTsv(FILE *f);
	uint GetSize() const;
	uint GetMotifCount() const { return m_MotifCount; }
	uint GetMotifLength(uint MotifIndex) const;
	const char *GetMotifName(uint MotifIndex) const;
	};
