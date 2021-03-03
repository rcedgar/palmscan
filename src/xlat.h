#ifndef xlat_h
#define xlat_h

extern unsigned g_MaxLaa;
extern byte *g_AASeqs[7];

unsigned XlatX(const byte *nt, unsigned Lnt, int Frame, byte *aa);
void AllocXlat(unsigned Lnt);
unsigned AAPosToNtPos(unsigned AAPos, unsigned Lnt, int Frame);
const char *XlatStr(const char *nt, unsigned Lnt, unsigned NtPos, 
  unsigned CodonCount, int Frame, string &s);
void PrintFastaSegNt(FILE *f, const char *nt, unsigned Lnt, unsigned NtPos, 
  unsigned CodonCount, int Frame);
void GetNtSeg(const char *nt, unsigned Lnt, unsigned NtPos, 
  unsigned CodonCount, int Frame, string &s);

#endif // xlat_h
