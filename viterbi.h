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

float ViterbiFastMem(XDPMem &Mem, float **ScoreMx,
  uint LA, uint LB, string &Path);

#endif // viterbi_h
