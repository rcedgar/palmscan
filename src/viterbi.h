#ifndef viterbi_h
#define viterbi_h

#include "xdpmem.h"
#include "alnparams.h"
#include "pathinfo.h"
#include "alnheuristics.h"
#include "seqinfo.h"

void LogAln(const byte *A, const byte *B, const char *Path);
void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps);

float ViterbiFastMainDiagMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, unsigned BandRadius, const AlnParams &AP,
  PathInfo &PI);

void USort(const SeqInfo &Query, const SeqDB &DB, vector<unsigned> &WordCounts, 
  vector<unsigned> &Order);

#endif // viterbi_h
