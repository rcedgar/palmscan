#include "myutils.h"
#include "sort.h"
#include "omplock.h"

#undef myalloc
#undef myfree

#include "allocid.h"

static unsigned g_AllocCounts[AllocIdCount];
static unsigned g_FreeCounts[AllocIdCount];
static float g_NetBytes[AllocIdCount];
static bool g_Trace;

void myalloc_trace(bool On)
	{
	g_Trace = On;
	}

static double g_TotalAlloc;
static double g_TotalFree;

static const char *AllocIdToStr(unsigned Id)
	{
	switch (Id)
		{
#define A(x)	case ALLOCID_##x: return #x;
#include "allocids.h"
		}
	return "?";
	}

void LogAllocs()
	{
	unsigned Order[AllocIdCount];

//	SortDescending(g_NetBytes, AllocIdCount, Order);
	QuickSortOrderDesc<float>(g_NetBytes, AllocIdCount, Order);

	float TotalBytes = 0.0f;
	for (unsigned k = 0; k < AllocIdCount; ++k)
		TotalBytes += g_NetBytes[k];

	Log("\n");
	Log("              Id      Allocs       Frees     Pct         Bytes\n");
	Log("----------------  ----------  ----------  ------  ------------\n");
	unsigned TotalAllocs = 0;
	unsigned TotalFrees = 0;
	double TotalPct = 0.0f;
	for (unsigned k = 0; k < AllocIdCount; ++k)
		{
		unsigned Id = Order[k];
		float Bytes = g_NetBytes[Id];
		if (Bytes == 0.0f)
			break;

		double pct = GetPct(Bytes, TotalBytes);
		TotalAllocs += g_AllocCounts[Id];
		TotalFrees+= g_FreeCounts[Id];
		TotalPct += pct;
		Log("%16.16s", AllocIdToStr(Id));
		Log("  %10u", g_AllocCounts[Id]);
		Log("  %10u", g_FreeCounts[Id]);
		Log("  %5.1f%%", pct);
		Log("  %12.12s", MemBytesToStr(Bytes));
		Log("\n");
		}
	Log("----------------  ----------  ----------  ------  ------------\n");
	Log("           Total  %10u  %10u  %5.1f%%  %12.12s\n",
	  TotalAllocs,
	  TotalFrees,
	  TotalPct,
	  MemBytesToStr(TotalBytes));

	double Curr = GetMemUseBytes();
	double Peak = GetPeakMemUseBytes();
	Log("\n");
	Log("%12.12s  Curr mem\n", MemBytesToStr(Curr));
	Log("%12.12s  Peak mem\n", MemBytesToStr(Peak));
	Log("%12.12s  Total alloc\n", MemBytesToStr(g_TotalAlloc));
	Log("%12.12s  Total free\n", MemBytesToStr(g_TotalFree));
	Log("%12.12s  Net\n", MemBytesToStr(g_TotalAlloc - g_TotalFree));

	double Excess = Curr - TotalBytes;
	Log("%12.12s  %cExcess\n", MemBytesToStr(fabs(Excess)), Excess > 0 ? '+' : '-');
	}

void *myalloc_track(uint32 Bytes, uint32 Id, const char *FileName, int LineNr)
	{
	LOCK();
	if (g_Trace && g_fLog != 0)
		Log("myalloc(%u, %s) %s(%d)\n",
		  Bytes, AllocIdToStr(Id), FileName, LineNr);

	g_TotalAlloc += (double) Bytes;

	uint32 *pu = (uint32 *) mymalloc(Bytes + 8);
	*pu++ = Id;
	*pu++ = Bytes;

	++g_AllocCounts[Id];
	g_NetBytes[Id] += Bytes;
	UNLOCK();
	return pu;
	}

void myfree_track(void *p, const char *FileName, int LineNr)
	{
	if (p == 0)
		return;

	LOCK();
	uint32 *pu = (uint32 *) p;

	unsigned Bytes = pu[-1];
	unsigned Id = pu[-2];
	if (g_Trace && g_fLog != 0)
		Log("myfree(%u, %s) %s(%d)\n",
		  Bytes, AllocIdToStr(Id), FileName, LineNr);

	g_TotalFree += (double) Bytes;

	asserta(Id < AllocIdCount);
	++g_FreeCounts[Id];
	g_NetBytes[Id] -= Bytes;

#if	RCE_MALLOC
	rce_free(pu - 2, __FILE__, __LINE__);
#else
	free(pu - 2);
#endif
	UNLOCK();
	}
