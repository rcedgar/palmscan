#pragma once

#define S(x)	extern FILE *g_f##x;
#include "segs.h"

void OpenSegFiles();
void CloseSegFiles();
