#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "rdrpmodel.h"

void RdRpModel::FromStrings(const vector<string> &Lines)
	{
	Clear();

	string Line;
	uint LineNr = 0;
#define NEXT_LINE()	{ asserta(LineNr < SIZE(Lines)); Line = Lines[LineNr++]; }

	vector<string> Fields;
	NEXT_LINE();
	asserta(Line == "rdrp_model");

	NEXT_LINE();
	Split(Line, Fields, '\t');
	if (SIZE(Fields) < 2 || Fields[0] != "motifs")
		Die("Bad motif line (fields) '%s'", Line.c_str());
	uint MotifCount = (uint) atoi(Fields[1].c_str());
	if (MotifCount == 0 || SIZE(Fields) != MotifCount+2)
		Die("Bad motif line (MotifCount) '%s'", Line.c_str());
	for (uint i = 0; i < MotifCount; ++i)
		{
		const string &s = Fields[i+2];
		if (SIZE(s) != 1)
			Die("Bad motif line (motif len) '%s'", Line.c_str());
		m_MotifLetters.push_back(s[0]);
		}

	NEXT_LINE();
	Split(Line, Fields, '\t');
	if (SIZE(Fields) < 2 || Fields[0] != "groups")
		Die("Bad groups line (fields) '%s'", Line.c_str());
	uint GroupCount = (uint) atoi(Fields[1].c_str());
	if (GroupCount == 0 || SIZE(Fields) != GroupCount+2)
		Die("Bad groups line (GroupCount) '%s'", Line.c_str());
	for (uint i = 0; i < GroupCount; ++i)
		{
		const string &s = Fields[i+2];
		m_GroupNames.push_back(s);
		}

	NEXT_LINE();
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2 || Fields[0] != "PSSMs")
		Die("Bad PSSMs line '%s'", Line.c_str());
	uint PSSMCount = (uint) atoi(Fields[1].c_str());
	m_PSSMs.resize(PSSMCount);
	for (uint i = 0; i < PSSMCount; ++i)
		{
		NEXT_LINE();
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 4 || Fields[0] != "PSSM" || SIZE(Fields[2]) != 1)
			Die("Bad PSSM prefix line '%s'", Line.c_str());

		char MotifLetter = Fields[2][0];
		const string &GroupName = Fields[3];
		m_PSSMMotifLetters.push_back(MotifLetter);
		m_PSSMGroupNames.push_back(GroupName);

		vector<string> PSSMStrings;
		for (;;)
			{
			NEXT_LINE();
			PSSMStrings.push_back(Line);
			if (Line == "//")
				break;
			}
		// m_PSSMs[i].FromFile(f);
		string PSSMName = GroupName + "." + MotifLetter;
		m_PSSMs[i].FromStrings(PSSMName, PSSMStrings);
		}
	NEXT_LINE();
	if (Line != ".")
		Die("Expected . in last line got '%s'", Line.c_str());

	SetPIV();
	}

void RdRpModel::FromModelFile(const string &FileName)
	{
	asserta(FileName != "");

	Clear();

	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(FileName);
	ReadLineStdioFile(f, Line);
	asserta(Line == "rdrp_model");

	ReadLineStdioFile(f, Line);
	Split(Line, Fields, '\t');
	if (SIZE(Fields) < 2 || Fields[0] != "motifs")
		Die("Bad motif line (fields) '%s'", Line.c_str());
	uint MotifCount = (uint) atoi(Fields[1].c_str());
	if (MotifCount == 0 || SIZE(Fields) != MotifCount+2)
		Die("Bad motif line (MotifCount) '%s'", Line.c_str());
	for (uint i = 0; i < MotifCount; ++i)
		{
		const string &s = Fields[i+2];
		if (SIZE(s) != 1)
			Die("Bad motif line (motif len) '%s'", Line.c_str());
		m_MotifLetters.push_back(s[0]);
		}

	ReadLineStdioFile(f, Line);
	Split(Line, Fields, '\t');
	if (SIZE(Fields) < 2 || Fields[0] != "groups")
		Die("Bad groups line (fields) '%s'", Line.c_str());
	uint GroupCount = (uint) atoi(Fields[1].c_str());
	if (GroupCount == 0 || SIZE(Fields) != GroupCount+2)
		Die("Bad groups line (GroupCount) '%s'", Line.c_str());
	for (uint i = 0; i < GroupCount; ++i)
		{
		const string &s = Fields[i+2];
		m_GroupNames.push_back(s);
		}

	ReadLineStdioFile(f, Line);
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2 || Fields[0] != "PSSMs")
		Die("Bad PSSMs line '%s'", Line.c_str());
	uint PSSMCount = (uint) atoi(Fields[1].c_str());
	m_PSSMs.resize(PSSMCount);
	for (uint i = 0; i < PSSMCount; ++i)
		{
		//fprintf(f, "PSSM	%u	%c	%s\n",
		//  i+1, m_PSSMMotifLetters[i], m_PSSMGroupNames[i].c_str());
		ReadLineStdioFile(f, Line);
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 4 || Fields[0] != "PSSM" || SIZE(Fields[2]) != 1)
			Die("Bad PSSM prefix line '%s'", Line.c_str());

		char MotifLetter = Fields[2][0];
		const string &GroupName = Fields[3];
		m_PSSMMotifLetters.push_back(MotifLetter);
		m_PSSMGroupNames.push_back(GroupName);

		m_PSSMs[i].FromFile(f);
		}
	ReadLineStdioFile(f, Line);
	if (Line != ".")
		Die("Expected . in last line got '%s'", Line.c_str());

	CloseStdioFile(f);

	SetPIV();
	}

void RdRpModel::ToModelFile(const string &FileName) const
	{
	asserta(FileName != "");

	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "rdrp_model\n");
	fprintf(f, "motifs	%u", GetMotifCount());
	for (uint i = 0; i < GetMotifCount(); ++i)
		fprintf(f, "\t%c", m_MotifLetters[i]);
	fprintf(f, "\n");

	fprintf(f, "groups	%u", GetGroupCount());
	for (uint i = 0; i < GetGroupCount(); ++i)
		fprintf(f, "\t%s", m_GroupNames[i].c_str());
	fprintf(f, "\n");

	fprintf(f, "PSSMs	%u\n", SIZE(m_PSSMs));
	for (uint i = 0; i < SIZE(m_PSSMs); ++i)
		{
		fprintf(f, "PSSM	%u	%c	%s\n",
		  i+1, m_PSSMMotifLetters[i], m_PSSMGroupNames[i].c_str());
		m_PSSMs[i].ToFile(f);
		}

	fprintf(f, ".\n");
	CloseStdioFile(f);
	}

void BuildRdRpModel()
	{
	const string &SpecFileName = opt_build_rdrp_model;

	RdRpModel Mod;
	Mod.FromSpecFile(SpecFileName);
	Mod.ToModelFile(opt_output);
	}
