#include "myutils.h"
#include "pssm.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "alnparams.h"
#include "objmgr.h"
#include "seqdb.h"

#if 0
float Viterbi_PSSM(XDPMem &Mem, const PSSM &PA, const PSSM &PB,
  const AlnParams &AP, PathInfo &PI);

void cmd_align_msas()
	{
	const string &FileNameA = opt_align_msas;
	const string &FileNameB = opt_msa2;

	AlnParams AP;
	AP.InitFromCmdLine(true);

	SeqDB MSAA;
	SeqDB MSAB;
	MSAA.FromFasta(FileNameA);
	MSAB.FromFasta(FileNameB);

	PSSM PA;
	PSSM PB;
	PA.FromFasta(FileNameA);
	PB.FromFasta(FileNameB);

	ObjMgr *OM = new ObjMgr;
	PathInfo *PI = OM->GetPathInfo();

	XDPMem Mem;
	float Score = Viterbi_PSSM(Mem, PA, PB, AP, *PI);

	const char *Path = PI->GetPath();

	ProgressLog("X=%s;Y=%s;Score=%.3g;Path=%s\n",
	  FileNameA.c_str(), FileNameB.c_str(), Score, Path);
	}
#endif // 0
