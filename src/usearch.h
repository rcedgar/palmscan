#ifndef usearch_h
#define usearch_h

class UsearchHit
	{
public:
	const char *QLabel;
	const byte *QSeq;
	unsigned QL;

	const char *TLabel;
	const byte *TSeq;
	unsigned TL;

	string Path;
	unsigned ColCount;
	unsigned IdCount;
	unsigned GapCount;
	};

void Usearch(const char *QLabel, const byte *QSeq, unsigned QL,
  SeqDB &DB, vector<UsearchHit> &Hits);

#endif // usearch_h
