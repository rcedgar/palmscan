#include "myutils.h"
#include "timing.h"
#include <time.h>
#include <algorithm>

#ifdef _MSC_VER
#include <Windows.h>	// for IsDebuggerPresent()
#endif

#if	TIMING

const double MIN_TIMER_PCT = 0.1;

TIMER g_CurrTimer = TIMER_None;
TIMER g_LastTimer;
TICKS g_LastTicks;

#if	TIMING_2D
TICKS g_PrefixTicks2D[TimerCount][TimerCount];
unsigned g_PrefixCounts2D[TimerCount][TimerCount];
#else
#define T(x)							\
TICKS g_PrefixTicks##x;
#include "timers.h"
#endif

unsigned g_Counters[CounterCount];

#define T(x)							\
	TICKS g_TotalTicks##x;				\
	unsigned g_TimerStartCount##x;
#include "timers.h"

#define T(x)	unsigned g_STimerTicks##x;
#include "stimers.h"

static TICKS g_StartTicks;
static time_t g_StartTime;

static double GetTimerOverheadTicks()
	{
	const unsigned N = 1000*1000;

	TICKS t1 = GetClockTicks();
	for (unsigned i = 0; i < N; ++i)
		{
		StartTimer(GetTimerOverhead);
		EndTimer(GetTimerOverhead);
		}
	TICKS t2 = GetClockTicks();

	double Ticks = double(t2 - t1);
	return Ticks/N;
	}

void InitTiming()
	{
	g_StartTime = time(0);
	TICKS Now = GetClockTicks();
	g_StartTicks = Now;
	g_LastTicks = Now;
	}

static const char *CounterToStr(COUNTER c)
	{
	switch (c)
		{
#define C(x)	case COUNTER_##x: return #x;
#include "counters.h"
#undef C
		}
	return "?";
	}

static void LogCounters()
	{
	Log("\n");
	Log("                        Counter            N\n");
	Log("--------------------------------  ----------\n");
	for (unsigned i = 0; i < CounterCount; ++i)
		{
		double N = g_Counters[i];
		if (N == 0)
			continue;
		Log("%32.32s  %10.10s\n", CounterToStr((COUNTER) i), FloatToStr(N));
		}
	}

const char *TimerToStr(TIMER t)
	{
	switch (t)
		{
#define T(x)	case TIMER_##x: return #x;
#include "timers.h"
#undef T
		}
	return "?";
	}

static void LogTimers(const GlobalTimingData &GTD,
  const vector<TimerData> &TDs)
	{
	const unsigned TDCount = SIZE(TDs);

	double TotalTimerCalls = 0.0;
	double TotalTimerTicks = 0.0;
	for (unsigned i = 0; i < TDCount; ++i)
		{
		const TimerData &TD = TDs[i];
		TotalTimerTicks += TD.Ticks;
		if (!TD.Is1D && !TD.Is2D)
			TotalTimerCalls += TD.Calls;
		}

	double TicksPerSec = GTD.ElapsedTicks/GTD.ElapsedSecs;
	double TotalTimerSecs = TotalTimerTicks/TicksPerSec;
	double TotalTimerOverheadTicks = TotalTimerCalls*GTD.TimerOverheadTicks;
	double TotalTimerOverheadSecs = TotalTimerOverheadTicks/TicksPerSec;

	double TotalPct = 0.0;
	double TotalPctShown = 0.0;
	double TotalCPct = 0.0;
	double TotalSecs = 0.0;

	Log("\n");
	Log("    Pct   TotPct        Ticks        Secs       Calls  Ticks/Call  Timer\n");
	Log("-------  -------  -----------  ----------  ----------  ----------  -----\n");

	for (unsigned i = 0; i < TDCount; ++i)
		{
		const TimerData &TD = TDs[i];

		double pct = Pct(TD.Ticks, TotalTimerTicks);
		TotalPct += pct;

		double Secs = TD.Ticks/TicksPerSec;
		TotalSecs += Secs;

		double CTicks = TD.Ticks;
		if (!TD.Is1D && !TD.Is2D)
			CTicks -= TD.Calls*GTD.TimerOverheadTicks;
		double TicksPerCall = TD.Calls == 0.0 ? 0.0 : CTicks/TD.Calls;

		if (pct >= MIN_TIMER_PCT)
			{
			TotalPctShown += pct;

			Log("%6.1f%%", pct);
			Log("  %6.1f%%", TotalPct);
			Log("  %11.4e", TD.Ticks);
			Log("  %10.10s", FloatToStr(Secs));
			Log("  %10.10s", FloatToStr(TD.Calls));
			Log("  %10.10s", TD.Is1D ? "" : FloatToStr(TicksPerCall));
			Log("  %s", TD.Name.c_str());
			Log("\n");
			}
		}

	Log("=======  =======  ===========  ==========  ==========  ==========  =====\n");

	Log("%6.1f%%", TotalPctShown);
	Log("  %6.1f%%", TotalPct);
	Log("  %11.11s", FloatToStr(TotalTimerTicks));
	Log("  %10.10s", FloatToStr(TotalSecs));
	Log("  %10.10s", FloatToStr(TotalTimerCalls));
	Log("  %10.10s", "");
	Log("  TOTAL\n");
	Log("\n");

	Log("\n");
	Log("%11.4g  Ticks elapsed (%.1f secs)\n",
	  GTD.ElapsedTicks, GTD.ElapsedTicks/TicksPerSec);
	Log("%11.4g  Ticks total timers (%.1f%%)\n",
	  TotalTimerTicks,
	  Pct(TotalTimerTicks, GTD.ElapsedTicks));
	Log("%11.11s  Secs elapsed\n", SecsToStr(GTD.ElapsedSecs));
	Log("%11.11s  Timer overhead ticks\n", FloatToStr(GTD.TimerOverheadTicks));

	double TimerOverheadSecs = TotalTimerCalls*GTD.TimerOverheadTicks/TicksPerSec;
	double EstdSecs = GTD.ElapsedSecs - TimerOverheadSecs;
	Log("%11.11s  Secs timer overhead (%.1f%%)\n",
	  SecsToStr(TimerOverheadSecs),
	  Pct(TimerOverheadSecs, GTD.ElapsedSecs));
	Log("%11.11s  Secs estd. without timers\n",
	  SecsToStr(EstdSecs));
	}

void LogTiming()
	{
	TICKS ExitTicks  = GetClockTicks();
	time_t ExitTime = time(0);

#if	TIMING_2D
	g_PrefixTicks2D[g_LastTimer][TIMER_ExitTiming] = ExitTicks - g_LastTicks;
#else
	g_TotalTicksExitTiming = ExitTicks - g_LastTicks;
#endif

	GlobalTimingData GTD;
	double ElapsedSecs = double(ExitTime - g_StartTime);

	if (ElapsedSecs == 0.0)
		ElapsedSecs = 1.0;

	double ElapsedTicks = double(ExitTicks - g_StartTicks);

	vector<TimerData> TDs;
#define T(x)											\
	{													\
	double Ticks = (double) g_TotalTicks##x;			\
	if (Ticks > 0.0)									\
		{												\
		TimerData TD;									\
		TD.Name = #x;									\
		TD.Ticks = Ticks;								\
		TD.Calls = g_TimerStartCount##x;				\
		TD.Is1D = false;								\
		TD.Is2D = false;								\
		TDs.push_back(TD);								\
		}												\
	}
#include "timers.h"

#if	TIMING_2D

#define T(x)											\
	{													\
	for (unsigned i = 0; i < TimerCount; ++i)			\
		{												\
		double Ticks = (double) g_PrefixTicks2D[TIMER_##x][i];\
		if (Ticks > 0.0)								\
			{											\
			TimerData TD;								\
			TD.Name = string(#x) + string(" - ")		\
			  + string(TimerToStr((TIMER) i));			\
			TD.Ticks = Ticks;							\
			TD.Calls = g_PrefixCounts2D[i][TIMER_##x];	\
			TD.Is1D = false;							\
			TD.Is2D = true;								\
			TDs.push_back(TD);							\
			}											\
		}												\
	}
#include "timers.h"

#else

#define T(x)											\
	{													\
	double Ticks = (double) g_PrefixTicks##x;			\
	if (Ticks > 0.0)									\
		{												\
		TimerData TD;									\
		TD.Name = "> " #x;								\
		TD.Ticks = Ticks;								\
		TD.Calls = g_TimerStartCount##x;				\
		TD.Is1D = true;									\
		TD.Is2D = false;								\
		TDs.push_back(TD);								\
		}												\
	}
#include "timers.h"

#endif

	sort(TDs.begin(), TDs.end());

	Log("\n");
#ifdef _MSC_VER
	Log("Debugger %s\n", IsDebuggerPresent() ? "Yes " : "No");
#endif

	GTD.ElapsedSecs = ElapsedSecs;
	GTD.ElapsedTicks = ElapsedTicks;
	GTD.TimerOverheadTicks = GetTimerOverheadTicks();

	LogCounters();
	LogTimers(GTD, TDs);

	double TicksPerSec = GTD.ElapsedTicks/GTD.ElapsedSecs;

	bool Any = false;
#define T(x)	\
	if (g_STimerTicks##x > 0)	\
		{ \
		if (!Any) \
			{ \
			Any = true; \
			Log("\n"); \
			Log("                          STimer        Ticks  Secs\n"); \
			Log("--------------------------------  -----------  ----\n"); \
			} \
		Log("%32.32s  %11.4g  %s\n", \
		  #x, \
		  (double) g_STimerTicks##x, \
		  SecsToStr(g_STimerTicks##x/TicksPerSec)); \
		}
#include "stimers.h"
	}

void ZeroTimers()
	{
#define T(x)	g_TimerStartCount##x = 0;
#include "timers.h"
	}

#else

void InitTiming()
	{
	}

void LogTiming()
	{
	}

#endif // TIMING
