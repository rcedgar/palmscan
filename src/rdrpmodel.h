#pragma once

#include "pssm.h"
#include "pssmsearch.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"
#include "rphit.h"
#include "heuristics.h"
#include "motifprofile.h"

class RdRpModel
	{
public:
	vector<string> m_GroupNames;
	vector<PSSM> m_PSSMs;
	uint m_LineIndex = 0;
	const vector<string> *m_Lines = 0;

public:

	void Clear()
		{
		m_PSSMs.clear();
		m_GroupNames.clear();
		m_Lines = 0;
		m_LineIndex = 0;
		}

	uint GetGroupCount() const
		{
		uint GroupCount = SIZE(m_GroupNames);
		return SIZE(m_GroupNames);
		}

	uint GetPSSMCount() const
		{
		return SIZE(m_PSSMs);
		}

	uint GetPSSMIndex(uint GroupIndex, uint MotifIndex) const
		{
		asserta(MotifIndex < 3);
		asserta(GroupIndex < SIZE(m_GroupNames));
		uint PSSMIndex = GroupIndex*3 + MotifIndex;
		asserta(PSSMIndex < SIZE(m_PSSMs));
		return PSSMIndex;
		}

	const PSSM &GetPSSM(uint GroupIndex, uint MotifIndex) const
		{
		uint PSSMIndex = GetPSSMIndex(GroupIndex, MotifIndex);
		return m_PSSMs[PSSMIndex];
		}

	void GetGroupName(uint GroupIndex, string &GroupName) const
		{
		asserta(GroupIndex < SIZE(m_GroupNames));
		GroupName = m_GroupNames[GroupIndex];
		}

	void FromPSSMs(const vector<string> &GroupNames,
	  const vector<PSSM> &PAs,
	  const vector<PSSM> &PBs,
	  const vector<PSSM> &PCs);
	void FromPSSMDir(const string &PSSMDir,
	  const vector<string> &GroupNames);
	void FromLines(const vector<string> &Lines);
	void FromModelFile(const string &FileName);
	bool GetNextLine(string &Line);
	void ToModelFile(const string &FileName) const;

	uint GetPSSMLength(uint GroupIndex, uint MotifIndex) const;
	void GetMotifProfile(uint GroupIndex, MotifProfile &MP) const;

public:
	void GetDefaultModelLines(vector<string> &Lines);
	};

extern vector<string> g_ModelStrings;
extern vector<string> g_ModelStringsRdRp;
