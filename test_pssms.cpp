#include "myutils.h"
#include "pssm.h"
#include "rdrpsearcher.h"
#include "seqdb.h"

void cmd_test_pssms()
	{
	const string &QueryFileName = opt_test_pssms;

	RdRpModel Model;
	Model.FromModelFile(opt_model);

	RdRpSearcher RS;
	RS.Init(Model);

	SeqDB Query;
	Query.FromFasta(QueryFileName);

	float MinScore = 1.0f;
	const uint SeqCount = Query.GetSeqCount();
	const uint GroupCount = RS.GetGroupCount();
	vector<RPHit> Hits;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = Query.GetLabel(SeqIndex);
		const string &Q = Query.GetSeq(SeqIndex);
		const uint QL = SIZE(Q);
		if (QL < 8)
			continue;

		Log("\n");
		Log("__________________________________________________\n");
		Log(">%s\n", Label.c_str());

		RS.m_QueryLabel = Label;
		RS.m_QuerySeq = Q;
		RS.m_QueryLabel = Label;

		for (uint MotifIndex = 0; MotifIndex < 3; ++MotifIndex)
			{
			for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
				{
				RS.SearchMotif(GroupIndex, MotifIndex, 0, QL-1, MinScore, Hits);
				const uint HitCount = SIZE(Hits);
				for (uint HitIndex = 0; HitIndex < HitCount; ++HitIndex)
					{
					const RPHit &Hit = Hits[HitIndex];
					uint QLo, QHi;
					string QRow, PRow, ARow;
					RS.GetAln(Hit, QRow, PRow, ARow, QLo, QHi);

					char MotifLetter = "ABC"[MotifIndex];
					string GroupName;
					RS.GetGroupName(GroupIndex, GroupName);

					Log("\n");
					Log("%s.%c", GroupName.c_str(), MotifLetter);
					Log(" (%.1f)", Hit.m_Score);
					Log(" %u-%u", QLo, QHi);
					Log("\n");

					Log("  %s\n", QRow.c_str());
					Log("  %s\n", ARow.c_str());
					Log("  %s\n", PRow.c_str());
					}
				}
			}
		}
	}
