#include "myutils.h"

const double A = 0.1;
const double DIAG = 1.0 - 3.0*A;

// Jukes-Cantor
double g_NtMx_JP_Rows[4][4] =
	{
//	A,		C,		G,		T
	{ DIAG,	A,		A,		A		},	// A
	{ A,	DIAG,	A,		A		},	// C
	{ A,	A,		DIAG,	A		},	// G
	{ A,	A,		A,		DIAG	},	// T
	};
