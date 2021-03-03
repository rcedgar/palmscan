#ifndef seqdb_h
#define seqdb_h

class SeqDB
	{
public:
	bool m_IsAligned;
	bool m_IsNucleo;
	bool m_IsNucleoSet;
	unsigned m_ColCount;
	vector<string> m_Labels;
	vector<string> m_Seqs;

public:
	SeqDB()
		{
		m_IsAligned = false;
		m_IsNucleo = false;
		m_IsNucleoSet = false;
		m_ColCount = UINT_MAX;
		}

	unsigned AddSeq(const string &Label, const string &Seq);
	const string &GetSeq(unsigned SeqIndex) const;
	const string &GetLabel(unsigned SeqIndex) const;
	unsigned GetSeqLength(unsigned SeqIndex) const;
	bool IsAligned() const;
	unsigned GetColCount() const;
	bool GetIsNucleo();
	unsigned GetSeqCount() const { return SIZE(m_Seqs); }
	void FromFasta(const string &FileName);
	void WritePretty(FILE *f) const;
	void WriteMSAPretty(FILE *f) const;
	void LogMe() const;

private:
	void SetIsNucleo();
	};

#endif // seqdb_h
