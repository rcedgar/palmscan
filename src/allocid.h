#ifndef allocid_h
#define allocid_h

enum ALLOCID
	{
#define	A(x)	ALLOCID_##x,
#include "allocids.h"
	};

const unsigned AllocIdCount = 0
#define A(x)	+ 1
#include "allocids.h"
	;

#endif // allocid_h
