#pragma once

enum PF
	{
#define F(x)	PF_##x,
#include "pppfeatures.h"
	};
