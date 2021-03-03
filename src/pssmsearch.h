#ifndef pssmsearch_h
#define pssmsearch_h

#include "usearch.h"
#include "pssms.h"

class PSSM;
class SFasta;

struct PSSMHitAA
	{
	const PSSM *P;
	const string *Label;
	const char *AASeq;
	unsigned L;
	unsigned AAPos;
	float Score;
	};

struct PSSMHitNt
	{
	const PSSM *P;
	const string *Label;
	const char *NtSeq;
	unsigned Lnt;

// AAPos is 0-based position from start of translated sequence
// NtPos is 0-based position from start of seq or rev-comp'd seq.
	unsigned AAPos;
	unsigned NtPos;
	int Frame;
	float Score;
	};

struct PSSMHitPairG
	{
	const PSSM *P1;
	const PSSM *P2;
	const string *Label;
	const char *NtSeq;
	unsigned Lnt;
	unsigned AAPos1;
	unsigned AAPos2;
	float Score1;
	float Score2;
	int Frame;
	char *AASeg1;
	char *AASeg2;
	char *AASeg;
	unsigned AASegLen;
	unsigned NtLo;
	unsigned NtHi;
	string NtSeg;
	string NtSeg2;
	};

struct PSSMHitPairNt
	{
	const PSSM *P1;
	const PSSM *P2;
	const string *Label;
	const char *NtSeq;
	unsigned Lnt;
	unsigned AAPos1;
	unsigned AAPos2;
	float Score1;
	float Score2;
	int Frame1;
	int Frame2;
	const char *AASeg1;
	const char *AASeg2;
	const char *AASeg;
	unsigned AASegLen;
	};

void GetTopHitX(const PSSM &P, float MinScore,
  const string &Label, const char *NtSeq, unsigned Lnt,
  PSSMHitNt &Hit);

void GetTopHitAA(const PSSM &P, float MinScore, 
  const string &Label, const char *AASeq, unsigned Laa,
  PSSMHitAA &Hit);

void GetTopHitNrPair(const PSSM &P1, float MinScore1,
  const PSSM &P2, float MinScore2,
  const string &Label, const char *NtSeq, unsigned Lnt,
  PSSMHitPairNt &Hit);

const char *GetSeg(const PSSMHitNt &Hit1, const PSSMHitNt &Hit2,
  int &Frame, unsigned &Lo, unsigned &Hi);

extern FILE *g_fOut;
extern FILE *g_fRep;

#endif // pssmsearch_h
