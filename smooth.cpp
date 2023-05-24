#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"

void cmd_smooth()
	{
	const string &FN = opt_smooth;

	ChainReader CR;
	CR.Open(FN);

	const uint ShapeLength = 64;
	const uint w = 2;

	PDBChain Chain;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		if (!Ok)
			break;

		Log("\n");
		Log(">%s\n", Chain.m_Label.c_str());
		for (uint i = 0; i < ShapeLength; ++i)
			{
			double Coord_X = 
			  Chain.GetSmoothedCoord(X, i, ShapeLength, w);

			double Coord_Y = 
			  Chain.GetSmoothedCoord(Y, i, ShapeLength, w);

			double Coord_Z = 
			  Chain.GetSmoothedCoord(Z, i, ShapeLength, w);

			Log("[%5u]", i);
			Log("  %8.2f", Coord_X);
			Log("  %8.2f", Coord_Y);
			Log("  %8.2f", Coord_Z);
			Log("\n");
			}
		}
	}
