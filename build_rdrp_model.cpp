#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "rdrpmodel.h"

bool RdRpModel::GetNextLine(string &Line)
	{
	asserta(m_Lines != 0);
	if (m_LineIndex >= SIZE(*m_Lines))
		return false;
	Line = (*m_Lines)[m_LineIndex++];
	return true;
	}

void RdRpModel::FromLines(const vector<string> &Lines)
	{
	Clear();
	m_Lines = &Lines;
	m_LineIndex = 0;

	string Line;
	vector<string> Fields;
	GetNextLine(Line);
	asserta(Line == "palmprint_model");

	GetNextLine(Line);
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2 || Fields[0] != "groups")
		Die("Bad model (groups) '%s'", Line.c_str());
	uint GroupCount = (uint) atoi(Fields[1].c_str());
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		GetNextLine(Line);
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3 || Fields[0] != "group"
		  || (uint) atoi(Fields[1].c_str()) != GroupIndex)
			Die("Bad model expected group %u got '%s'",
			  GroupIndex, Line.c_str());

		const string &Name = Fields[2];
		m_GroupNames.push_back(Name);
		}

	uint PSSMCount = 3*GroupCount;
	m_PSSMs.resize(PSSMCount);
	for (uint PSSMIndex = 0; PSSMIndex < PSSMCount; ++PSSMIndex)
		{
		GetNextLine(Line);
		Split(Line, Fields, '\t');

		uint GroupIndex = PSSMIndex/3;
		string GroupName = m_GroupNames[GroupIndex];

		uint MotifIndex = PSSMIndex%3;
		char MotifLetter = "ABC"[MotifIndex];
		string MotifLetterStr;
		MotifLetterStr.push_back(MotifLetter);

		if (SIZE(Fields) != 3
		  || Fields[0] != "pssm"
		  || Fields[1] != GroupName
		  || Fields[2] != MotifLetterStr)
			Die("Bad PSSM prefix line '%s'", Line.c_str());
	
		vector<string> PSSMStrings;
		for (;;)
			{
			bool Ok = GetNextLine(Line);
			if (!Ok)
				Die("Missing // at end of PSSM strings");
			PSSMStrings.push_back(Line);
			if (Line == "//")
				break;
			}

		m_PSSMs[PSSMIndex].FromStrings(GroupName, PSSMStrings);
		m_PSSMs[PSSMIndex].m_GroupName = GroupName;
		}
	GetNextLine(Line);
	if (Line != ".")
		Die("Expected . in last line got '%s'", Line.c_str());
	}

void RdRpModel::FromModelFile(const string &FileName)
	{
	asserta(FileName != "");

	Clear();

	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(FileName);
	ReadLineStdioFile(f, Line);
	asserta(Line == "palmprint_model");

	ReadLineStdioFile(f, Line);
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2 || Fields[0] != "groups")
		Die("Bad model (groups) '%s'", Line.c_str());
	uint GroupCount = (uint) atoi(Fields[1].c_str());
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		ReadLineStdioFile(f, Line);
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3 || Fields[0] != "group"
		  || (uint) atoi(Fields[1].c_str()) != GroupIndex)
			Die("Bad model expected group %u got '%s'",
			  GroupIndex, Line.c_str());

		const string &Name = Fields[2];
		m_GroupNames.push_back(Name);
		}

	uint PSSMCount = 3*GroupCount;
	m_PSSMs.resize(PSSMCount);
	for (uint PSSMIndex = 0; PSSMIndex < PSSMCount; ++PSSMIndex)
		{
		ReadLineStdioFile(f, Line);
		Split(Line, Fields, '\t');

		uint GroupIndex = PSSMIndex/3;
		string GroupName = m_GroupNames[GroupIndex];

		uint MotifIndex = PSSMIndex%3;
		char MotifLetter = "ABC"[MotifIndex];
		string MotifLetterStr;
		MotifLetterStr.push_back(MotifLetter);

		if (SIZE(Fields) != 3
		  || Fields[0] != "pssm"
		  || Fields[1] != GroupName
		  || Fields[2] != MotifLetterStr)
			Die("Bad PSSM prefix line '%s'", Line.c_str());

		m_PSSMs[PSSMIndex].FromFile(f);
		m_PSSMs[PSSMIndex].m_GroupName = GroupName;
		}
	ReadLineStdioFile(f, Line);
	if (Line != ".")
		Die("Expected . in last line got '%s'", Line.c_str());

	CloseStdioFile(f);
	}

void RdRpModel::ToModelFile(const string &FileName) const
	{
	if (FileName == "")
		return;

	uint GroupCount = GetGroupCount();

	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "palmprint_model\n");
	fprintf(f, "groups	%u\n", GroupCount);
	for (uint i = 0; i < GroupCount; ++i)
		fprintf(f, "group\t%u\t%s\n", i, m_GroupNames[i].c_str());

	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		for (uint MotifIndex = 0; MotifIndex < 3; ++MotifIndex)
			{
			char MotifLetter = "ABC"[MotifIndex];
			string GroupName;
			GetGroupName(GroupIndex, GroupName);

			fprintf(f, "pssm\t%s\t%c\n",GroupName.c_str(), MotifLetter);
			m_PSSMs[GroupIndex*3 + MotifIndex].ToFile(f);
			}
		}

	fprintf(f, ".\n");
	CloseStdioFile(f);
	}

void cmd_build_rdrp_model()
	{
	const string &NamesFileName = opt_build_rdrp_model;
	string PSSMDir = "./";
	if (optset_pssmdir)
		PSSMDir = opt_pssmdir;
	uint n = SIZE(PSSMDir);
	if (PSSMDir[n-1] != '/')
		PSSMDir += "/";

	vector<string> Names;
	string Name;
	FILE *f = OpenStdioFile(NamesFileName);
	while (ReadLineStdioFile(f, Name))
		Names.push_back(Name);
	CloseStdioFile(f);

	RdRpModel Mod;
	Mod.FromPSSMDir(PSSMDir, Names);
	Mod.ToModelFile(opt_output);
	}
