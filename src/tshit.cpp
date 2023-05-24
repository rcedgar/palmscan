#include "myutils.h"
#include "tshit.h"
#include "abcxyz.h"

float GetBlosum62Score(char a, char b);
void GetPalmSketch(const string& ss, uint PSL, string& Sketch);
void GetPalmSubSketch(const string& ss, uint Lo, uint n, 
  uint PSL, string& Sketch);

const uint V_SKETCH_LENGTH = 20;

static void MakeAnnotSketch(const string &QSketch, const string &RSketch,
  string &AnnotSketch)
	{
	AnnotSketch.clear();
	const uint L = SIZE(QSketch);
	asserta(SIZE(RSketch) == L);
	for (uint i = 0; i < L; ++i)
		{
		char q = QSketch[i];
		char r = RSketch[i];
		if (q == r  && q != ' ' && r != ' ')
			AnnotSketch += '|';
		else
			AnnotSketch += ' ';
		}
	}

// D:\src\py\conf_score_test.py
double GetConfidenceScore(double MotifRMSD)
	{
	const double T = 2.0;	// RMSD where results start to degrade
	double Ratio = MotifRMSD*MotifRMSD/(T*T);
	double ex = exp(Ratio);
	double Conf = (1.0 + ex)/(2.0*ex);
	return Conf;
	}

static char GetAnnotChar(char a, char b)
	{
	if (!isalpha(a) || !isalpha(b))
		return ' ';
	if (a == b)
		return '|';

	float Score = GetBlosum62Score(a, b);

	if (Score >= 0.5)
		return '+';
	else if (Score > 0)
		return '.';
	return ' ';
	}

static void MakeAnnot(const string &A, const string &B, 
  string &Annot)
	{
	const uint L = SIZE(A);
	assert(SIZE(B) == L);
	Annot.clear();
	for (uint i = 0; i < L; ++i)
		{
		char c = GetAnnotChar(A[i], B[i]);
		Annot += c;
		}
	}

void TSHit::WriteAln(FILE *f) const
	{
	asserta(m_Query != 0 && m_Ref != 0);
	
	const string &QLabel = m_Query->m_Label;
	const string &RLabel = m_Ref->m_Label;

	uint QPosA = m_QPosA;
	uint QPosB = m_QPosB;
	uint QPosC = m_QPosC;

	string QASeq, QBSeq, QCSeq;
	m_Query->GetSubSeq(QPosA, AL, true, QASeq);
	m_Query->GetSubSeq(QPosB, BL, true, QBSeq);
	m_Query->GetSubSeq(QPosC, CL, true, QCSeq);

	asserta(m_Ref->m_MotifPosVec.size() == 3);
	uint RPosA = m_Ref->m_MotifPosVec[A];
	uint RPosB = m_Ref->m_MotifPosVec[B];
	uint RPosC = m_Ref->m_MotifPosVec[C];

	string RASeq, RBSeq, RCSeq;
	m_Ref->GetSubSeq(RPosA, AL, true, RASeq);
	m_Ref->GetSubSeq(RPosB, BL, true, RBSeq);
	m_Ref->GetSubSeq(RPosC, CL, true, RCSeq);

	string AnnotA, AnnotB, AnnotC;
	MakeAnnot(QASeq, RASeq, AnnotA);
	MakeAnnot(QBSeq, RBSeq, AnnotB);
	MakeAnnot(QCSeq, RCSeq, AnnotC);

	string QRow = QASeq + "   " + QBSeq + "   " + QCSeq;
	string RRow = RASeq + "   " + RBSeq + "   " + RCSeq;
	string AnnotRow = AnnotA + "   " + AnnotB + "   " + AnnotC;

	QRow += "  " + QLabel;
	RRow += "  " + RLabel;

	string CoordStr;
	string CoordsA;
	string CoordsB;
	string CoordsC;

	Ps(CoordsA, "A:%u", QPosA+1);
	Ps(CoordsB, "B:%u", QPosB+1);
	Ps(CoordsC, "C:%u", QPosC+1);

	uint PPL = QPosC + CL - QPosA;
	Ps(CoordStr, "%-12.12s   %-14.14s   %-8.8s  Palmprint(%u-%u, %u aa)",
	  CoordsA.c_str(), CoordsB.c_str(), CoordsC.c_str(),
	  QPosA+1, QPosC+CL, PPL);

	fprintf(f, "\n");
	fprintf(f, "%s\n", CoordStr.c_str());
	fprintf(f, "%s\n", QRow.c_str());
	fprintf(f, "%s\n", AnnotRow.c_str());
	fprintf(f, "%s\n", RRow.c_str());
	}

void TSHit::WriteTsv(FILE *f) const
	{
	if (f == 0)
		return;

	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "Query");
		fprintf(f, "\tRef");
		fprintf(f, "\tRMSDm");
		fprintf(f, "\tPosA");
		fprintf(f, "\tPosB");
		fprintf(f, "\tPosC");
		fprintf(f, "\tMotifA");
		fprintf(f, "\tMotifB");
		fprintf(f, "\tMotifC");
		fprintf(f, "\n");
		HdrDone = true;
		}

	asserta(m_Query != 0);
	asserta(m_Ref != 0);

	const string &QLabel = m_Query->m_Label;
	const string &RLabel = m_Ref->m_Label;

	double MotifRMSD = sqrt(m_MotifRMSD2);
	double Score = m_Score;

	string QASeq, QBSeq, QCSeq;
	m_Query->GetSubSeq(m_QPosA, AL, true, QASeq);
	m_Query->GetSubSeq(m_QPosB, BL, true, QBSeq);
	m_Query->GetSubSeq(m_QPosC, CL, true, QCSeq);

	fprintf(f, "%s", QLabel.c_str());
	fprintf(f, "\t%s", RLabel.c_str());
	fprintf(f, "\t%.3g", MotifRMSD);
	fprintf(f, "\t%u", m_QPosA + 1);
	fprintf(f, "\t%u", m_QPosB + 1);
	fprintf(f, "\t%u", m_QPosC + 1);
	fprintf(f, "\t%s", QASeq.c_str());
	fprintf(f, "\t%s", QBSeq.c_str());
	fprintf(f, "\t%s", QCSeq.c_str());
	fprintf(f, "\n");
	}

void TSHit::WritePalmprintFasta(FILE *f) const
	{
	if (f == 0)
		return;
	asserta(m_Query != 0);
	const string &Q = m_Query->m_Seq;
	const uint QL = SIZE(Q);
	uint PPLo = m_QPosA;
	uint PPHi = m_QPosC + CL - 1;
	uint PPL = PPHi - PPLo + 1;
	asserta(PPHi > PPLo && PPHi < QL);
	string Label;
	const string &QLabel = m_Query->m_Label;
	SeqToFasta(f, Label.c_str(), Q.c_str() + PPLo, PPL);
	}

void TSHit::SetSketch()
	{
	//asserta(SIZE(m_Query->m_SS) == SIZE(m_Query->m_Seq));
	m_Query->m_MotifPosVec.clear();
	m_Query->m_MotifPosVec.push_back(m_QPosA);
	m_Query->m_MotifPosVec.push_back(m_QPosB);
	m_Query->m_MotifPosVec.push_back(m_QPosC);

	if (m_Query->m_SS.empty())
		m_Query->SetSS();

	assert(m_Query->CheckMotifCoords());
	assert(m_Ref->CheckPPCMotifCoords());

	m_SketchPct = 0;
	m_QSketch.clear();
	m_RSketch.clear();

	uint RPosA = m_Ref->m_MotifPosVec[A];
	uint RPosB = m_Ref->m_MotifPosVec[B];
	uint RPosC = m_Ref->m_MotifPosVec[C];

	const string& QSS = m_Query->m_SS;
	const string& RSS = m_Ref->m_SS;

	uint QL = m_Query->GetSeqLength();
	uint RL = m_Ref->GetSeqLength();

	string QASk = QSS.substr(m_QPosA + 2, AL - 2);
	string RASk = RSS.substr(RPosA + 2, AL - 2);

	string QBSk = QSS.substr(m_QPosB, BL);
	string RBSk = RSS.substr(RPosB, BL);

	string QCSk = QSS.substr(m_QPosC, CL - 2);
	string RCSk = RSS.substr(RPosC, CL - 2);

	uint QV1L = m_QPosB - (m_QPosA + AL);
	uint RV1L = RPosB - (RPosA + AL);

	uint QV2L = m_QPosC - (m_QPosB + BL);
	uint RV2L = RPosC - (RPosB + BL);

	string QV1Sk;
	string RV1Sk;
	GetPalmSubSketch(QSS, m_QPosA + 1, QV1L, V_SKETCH_LENGTH, QV1Sk);
	GetPalmSubSketch(RSS, RPosA + 1, RV1L, V_SKETCH_LENGTH, RV1Sk);

	string QV2Sk;
	string RV2Sk;
	GetPalmSubSketch(QSS, m_QPosB + 1, QV2L, V_SKETCH_LENGTH, QV2Sk);
	GetPalmSubSketch(RSS, RPosB + 1, RV2L, V_SKETCH_LENGTH, RV2Sk);

	m_QSketch = QASk + "  " + QV1Sk + "  " + QBSk + "  " + QV2Sk + "  " + QCSk;
	m_RSketch = RASk + "  " + RV1Sk + "  " + RBSk + "  " + RV2Sk + "  " + RCSk;
	m_AnnotSketch;
	MakeAnnotSketch(m_QSketch, m_RSketch, m_AnnotSketch);

	const uint ColCount = SIZE(m_QSketch);
	asserta(SIZE(m_RSketch) == ColCount);
	uint N = 0;
	uint n = 0;
	for (uint i = 0; i < ColCount; ++i)
		{
		char q = m_QSketch[i];
		char r = m_RSketch[i];
		if (q == ' ' || r == ' ')
			continue;
		N += 1;
		if (q == r)
			n += 1;
		}
	m_SketchPct = GetPct(n, N);
	}

void TSHit::WriteSketch(FILE *f) const
	{
	if (f == 0)
		return;
	asserta(!m_QSketch.empty());

	fprintf(f, "\n");
	fprintf(f,
"_____A____  _________V1_________  ______B_______  _________V2_________  ___C__\n");

	fprintf(f, "%s  %s\n", m_QSketch.c_str(), m_Query->m_Label.c_str());
	fprintf(f, "%s\n", m_AnnotSketch.c_str());
	fprintf(f, "%s  %s\n", m_RSketch.c_str(), m_Ref->m_Label.c_str());
	}

void TSHit::WriteReport(FILE *f) const
	{
	if (f == 0)
		return;

	const string &QLabel = m_Query->m_Label;
	uint QL = m_Query->GetSeqLength();

	fprintf(f, "\n");
	fprintf(f, "_______________________________________________\n");
	fprintf(f, "Query >%s (%u aa)\n", QLabel.c_str(), QL);

	WriteAln(f);
	WriteSketch(f);
	double RMSD = sqrt(m_MotifRMSD2);

	fprintf(f, "\n");
	fprintf(f, "Score %.3f (%.2f, %.1f%%)",
	  m_Score, RMSD, m_SketchPct);

	if (m_Score >= 0.75)
		fprintf(f, " high confidence");
	else if (m_Score >= 0.50)
		fprintf(f, " low confidence");
	else
		fprintf(f, " likely false positive");
	fprintf(f, "\n");
	}
