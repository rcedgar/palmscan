#include "myutils.h"
#include "mpcluster.h"
#include "pssm.h"
#include "rdrpmodel.h"

void MPCluster::ReadSeqsVec(const string &TsvFileName,
  vector<vector<string> > &SeqsVec)
	{
	SeqsVec.clear();

	FILE *f = OpenStdioFile(TsvFileName);
	string Line;
	vector<string> Fields;
/***
C       0       11477   ICTDFSKFDQHFSGSGGTNADETLAHCLGDDGVL	ABC
M       0       5.15    VCTDFSSFDQHFSGSGGTNCDETIAHCLGDDGIL	ABC
M       0       4.91    VCTDFSRFEPHFSGSGGTNADETLVHCLGDDGII	ABC
M       0       4.82    VCTDFSRFDQHFSGSGETNADESLAHCNGDDGIL	ABC
***/
	vector<string> MemberSeqs;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 5);
		const string &Type = Fields[0];
		const string &Seq = Fields[3];
		if (Type == "C")
			{
			if (!MemberSeqs.empty())
				{
				SeqsVec.push_back(MemberSeqs);
				MemberSeqs.clear();
				MemberSeqs.push_back(Seq);
				}
			}
		else if (Type == "M")
			MemberSeqs.push_back(Seq);
		else
			asserta(false);
		}
	asserta(!MemberSeqs.empty());
	SeqsVec.push_back(MemberSeqs);

	CloseStdioFile(f);
	}

void MPCluster::ReadMPs(const string &TsvFileName,
  vector<MotifProfile *> &MPs)
	{
	MPs.clear();

	vector<vector<string> > SeqsVec;
	ReadSeqsVec(TsvFileName, SeqsVec);

	const uint ClusterCount = SIZE(SeqsVec);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount;
	  ++ClusterIndex)
		{
		MotifProfile *MP = new MotifProfile;
		const vector<string> &MemberSeqs = SeqsVec[ClusterIndex];
		MP->FromSeqs(MemberSeqs);
		MPs.push_back(MP);
		}
	}

static void SeqToMotifSeqs(const string &Seq,
  string &SeqA, string &SeqB, string &SeqC)
	{
	asserta(SIZE(Seq) == AL + BL + CL);
	SeqA = Seq.substr(0, AL);
	SeqB = Seq.substr(AL, BL);
	SeqC = Seq.substr(AL+BL, CL);
	asserta(SIZE(SeqA) == AL);
	asserta(SIZE(SeqB) == BL);
	asserta(SIZE(SeqC) == CL);
	}

static void SeqsToMotifSeqs(const vector<string> &Seqs, 
  vector<string> &SeqsA,
  vector<string> &SeqsB,
  vector<string> &SeqsC)
	{
	SeqsA.clear();
	SeqsB.clear();
	SeqsC.clear();

	const uint N = SIZE(Seqs);
	for (uint i = 0; i < N; ++i)
		{
		string SeqA;
		string SeqB;
		string SeqC;
		SeqToMotifSeqs(Seqs[i], SeqA, SeqB, SeqC);
		SeqsA.push_back(SeqA);
		SeqsB.push_back(SeqB);
		SeqsC.push_back(SeqC);
		}
	}

void cmd_clusters2model()
	{
	const string &InputFileName = opt_clusters2model;
	vector<vector<string> > SeqsVec;
	MPCluster::ReadSeqsVec(InputFileName, SeqsVec);
	const uint ClusterCount = SIZE(SeqsVec);
	vector<string> SeqsA;
	vector<string> SeqsB;
	vector<string> SeqsC;
	uint N = ClusterCount;
	if (optset_topn)
		N = min(ClusterCount, opt_topn);

	RdRpModel RM;
	for (uint ClusterIndex = 0; ClusterIndex < N; ++ClusterIndex)
		{
		const vector<string> &Seqs = SeqsVec[ClusterIndex];
		SeqsToMotifSeqs(Seqs, SeqsA, SeqsB, SeqsC);

		PSSM *PA = new PSSM;
		PSSM *PB = new PSSM;
		PSSM *PC = new PSSM;

		PA->FromSeqs(SeqsA);
		PB->FromSeqs(SeqsB);
		PC->FromSeqs(SeqsC);

		string ConsSeqA;
		string ConsSeqB;
		string ConsSeqC;
		PA->CalcConsSeq(ConsSeqA);
		PB->CalcConsSeq(ConsSeqB);
		PC->CalcConsSeq(ConsSeqC);
		Log("[%5u]  %s  %s  %s\n",
		  ClusterIndex, ConsSeqA.c_str(), ConsSeqB.c_str(), ConsSeqC.c_str());

		string GroupName;
		Ps(GroupName, "Cluster_%u", ClusterIndex+1);
		RM.m_GroupNames.push_back(GroupName);
		RM.m_PSSMs.push_back(*PA);
		RM.m_PSSMs.push_back(*PB);
		RM.m_PSSMs.push_back(*PC);
		}
	RM.ToModelFile(opt_model);
	}
