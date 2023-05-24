#pragma once

#define F(x)	extern FILE *g_f##x;
#include "ofiles.h"

void OpenOutputFiles();
void CloseOutputFiles();
void LockOutput();
void UnlockOutput();
