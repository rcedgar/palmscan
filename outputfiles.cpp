#include "myutils.h"
#include "omplock.h"

#define F(x)	FILE *g_f##x;
#include "ofiles.h"

void OpenOutputFiles()
	{
#define F(x)	g_f##x = CreateStdioFile(opt_##x);
#include "ofiles.h"
	}

void CloseOutputFiles()
	{
#define F(x)	CloseStdioFile(g_f##x);
#include "ofiles.h"
	}

void LockOutput()
	{
	Lock("Output");
	}

void UnlockOutput()
	{
	Unlock("Output");
	}
