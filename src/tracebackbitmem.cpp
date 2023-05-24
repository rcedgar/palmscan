#include "myutils.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "pathinfo.h"

void TraceBackBitMem(XDPMem &Mem, unsigned LA, unsigned LB,
  char State, string &Path)
	{
	Path.clear();

	byte **TB = Mem.GetTBBit();

	size_t i = LA;
	size_t j = LB;
	for (;;)
		{
		if (i == 0 && j == 0)
			break;

		Path += State;

		byte t;
		switch (State)
			{
		case 'M':
			asserta(i > 0 && j > 0);
			t = TB[i-1][j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else
				State = 'M';
			--i;
			--j;
			break;
		case 'D':
			asserta(i > 0);
			t = TB[i-1][j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			t = TB[i][j-1];
			if (t & TRACEBITS_MI)
				State = 'M';
			else
				State = 'I';
			--j;
			break;

		default:
			Die("TraceBackBit, invalid state %c", State);
			}
		}
	reverse(Path.begin(), Path.end());
	}
