#include "myutils.h"

static const uint g_SegCount = 0
#define S(x)	+ 1
#include "segs.h"
	;

#define S(x)	FILE *g_f##x = 0;
#include "segs.h"

void OpenSegFiles()
	{
	string Prefix = opt_seg_fasta_prefix;
	if (Prefix.empty())
		return;

#define S(x)	g_f##x = CreateStdioFile(Prefix + #x);
#include "segs.h"
	}

void CloseSegFiles()
	{
#define S(x)	CloseStdioFile(g_f##x);
#include "segs.h"
	}
