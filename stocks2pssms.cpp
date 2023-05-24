#include "myutils.h"
#include "stock.h"
#include "pssm.h"
#include "rdrpmodel.h"

static void GetNameFromPath(const string &Path, string &Name)
	{
	vector<string> Fields;
	Split(Path, Fields, '/');
	const string s1 = Fields[SIZE(Fields) - 1];

	Split(s1, Fields, '\\');
	const string s2 = Fields[SIZE(Fields) - 1];

	Split(s2, Fields, '.');
	Name = Fields[0];
	}

void cmd_stocks2pssms()
	{
	const string &InputFN = opt_stocks2pssms;
	RdRpModel Model;
	
	vector<string> StockFileNames;
	ReadLinesFromFile(InputFN, StockFileNames);
	const uint N = SIZE(StockFileNames);
	vector<Stock *> Stocks;
	vector<string> GroupNames;
	vector<PSSM> PAs;
	vector<PSSM> PBs;
	vector<PSSM> PCs;
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = StockFileNames[i];
		string GroupName;
		GetNameFromPath(FN, GroupName);

		ProgressStep(i, N, "%s", FN.c_str());
		Stock * S = new Stock;
		S->FromFile(FN);
		Stocks.push_back(S);

		vector<string> As, Bs, Cs;
		S->GetMotifsStrings(As, Bs, Cs);

		PSSM PA, PB, PC;
		PA.FromSeqs(As);
		PB.FromSeqs(Bs);
		PC.FromSeqs(Cs);

		GroupNames.push_back(GroupName);
		PA.m_GroupName = GroupName;
		PB.m_GroupName = GroupName;
		PC.m_GroupName = GroupName;

		string ConsA, ConsB, ConsC;
		PA.CalcConsSeq(ConsA);
		PB.CalcConsSeq(ConsB);
		PC.CalcConsSeq(ConsC);
		Log("%s  %s  %s  %s\n",
		  ConsA.c_str(),  ConsA.c_str(),  ConsA.c_str(),
		  GroupName.c_str());

		PAs.push_back(PA);
		PBs.push_back(PB);
		PCs.push_back(PC);

		vector<double> Probs;
		for (uint i = 0; i < 8; ++i)
			{
			PA.CalcProbs(i, Probs);
			double Sum = 0;
			for (uint Letter = 0; Letter < 20; ++Letter)
				{
				char c = g_LetterToCharAmino[Letter];
				double P = Probs[Letter];
				Sum += P;
				Log(" %c=%.4f", c, P);
				}
			Log("  Sum=%.4f\n", Sum);
			}
		}

	Model.FromPSSMs(GroupNames, PAs, PBs, PCs);
	Model.ToModelFile(opt_model);
	}
