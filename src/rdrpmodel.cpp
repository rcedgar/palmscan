#include "myutils.h"
#include "rdrpmodel.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "alnparams.h"
#include "pathinfo.h"
#include "sort.h"
#include <set>

void SeqToFasta(FILE *f, const char *Label, const char *Seq, unsigned L);

void RdRpModel::FromPSSMs(const vector<string> &GroupNames,
  const vector<PSSM> &PAs,
  const vector<PSSM> &PBs,
  const vector<PSSM> &PCs)
	{
	Clear();
	const uint GroupCount = SIZE(GroupNames);
	asserta(SIZE(PAs) == GroupCount);
	asserta(SIZE(PBs) == GroupCount);
	asserta(SIZE(PCs) == GroupCount);

	m_GroupNames = GroupNames;
	m_PSSMs.resize(GroupCount*3);
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		m_PSSMs[3*GroupIndex + 0] = PAs[GroupIndex];
		m_PSSMs[3*GroupIndex + 1] = PBs[GroupIndex];
		m_PSSMs[3*GroupIndex + 2] = PCs[GroupIndex];
		}
	}

// PSSMs must be in current directory, filenames GroupName.X where X=A,B,C.
void RdRpModel::FromPSSMDir(const string &PSSMDir,
  const vector<string> &GroupNames)
	{
	Clear();

	m_GroupNames = GroupNames;
	const uint GroupCount = SIZE(m_GroupNames);
	m_PSSMs.resize(GroupCount*3);
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		string FileNameA = PSSMDir + m_GroupNames[GroupIndex] + ".A";
		string FileNameB = PSSMDir + m_GroupNames[GroupIndex] + ".B";
		string FileNameC = PSSMDir + m_GroupNames[GroupIndex] + ".C";

		m_PSSMs[3*GroupIndex + 0].FromFile(FileNameA);
		m_PSSMs[3*GroupIndex + 1].FromFile(FileNameB);
		m_PSSMs[3*GroupIndex + 2].FromFile(FileNameC);
		}
	}

uint RdRpModel::GetPSSMLength(uint GroupIndex, uint MotifIndex) const
	{
	const PSSM &P = GetPSSM(GroupIndex, MotifIndex);
	uint L = P.GetColCount();
	return L;
	}

void RdRpModel::GetMotifProfile(uint GroupIndex, MotifProfile &MP) const
	{
	MP.Clear();
	asserta(GroupIndex < SIZE(m_GroupNames));
	const PSSM &PA = GetPSSM(GroupIndex, 0);
	const PSSM &PB = GetPSSM(GroupIndex, 1);
	const PSSM &PC = GetPSSM(GroupIndex, 2);
	MP.FromPSSMs(PA, PB, PC);
	}
