#include "myutils.h"
#include "seqdb.h"
#include "abcxyz.h"
#include "mx.h"
#include "alpha.h"
#include "rdrpsearcher.h"
#include "mpcluster.h"
#include "outputfiles.h"
#include "sort.h"
#include "segfiles.h"

void SetBLOSUM62();

extern float **g_SubstMx;
uint MPCluster::m_LeftFlank_Core = 149;		// see also RdRpSearcher
uint MPCluster::m_RightFlank_Core = 150;	// see also RdRpSearcher

void MPCluster::TrimLeft(string &Left)
	{
	uint n = SIZE(Left);
	if (n > m_LeftFlank_Core)
		Left = Left.substr(n-m_LeftFlank_Core, m_LeftFlank_Core);
	}

void MPCluster::TrimRight(string &Right)
	{
	uint n = SIZE(Right);
	if (n > m_RightFlank_Core)
		Right.resize(m_RightFlank_Core);
	}

char GetAnnotChar(uint Letter1, uint Letter2)
	{
	if (Letter1 == Letter2)
		return '|';
	float Score = g_SubstMx[Letter1][Letter2];
	if (Score >= 1)
		return '+';
	if (Score > 0)
		return '.';
	return ' ';
	}

void MPCluster::LogLogos(const vector<MotifProfile *> &MPs)
	{
	for (uint i = 0; i < SIZE(MPs); ++i)
		{
		string Logo;
		MPs[i]->GetLogo(Logo);
		Log("[%5u]  %s\n", i, Logo.c_str());
		}
	}

float MPCluster::GetScore_DotProduct(const MotifProfile &MP1,
  const MotifProfile &MP2) const
	{
#if DEBUG
	MP1.ValidateFreqs();
	MP2.ValidateFreqs();
#endif
	float Sum = 0;
	for (uint i = 0; i < MPL; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			float f1 = MP1.m_FreqVec[i][j];
			float f2 = MP2.m_FreqVec[i][j];
			Sum += f1*f2;
			}
		}
	return Sum;
	}

float MPCluster::GetScore(const MotifProfile &MP1,
  const MotifProfile &MP2) const
	{
	if (opt_dotproduct)
		return GetScore_DotProduct(MP1, MP2);

#if DEBUG
	MP1.ValidateFreqs();
	MP2.ValidateFreqs();
#endif
	float Sum = 0;
	for (uint i = 0; i < MPL; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			float f1 = MP1.m_FreqVec[i][j];
			if (f1 < 1e-6)
				continue;
			byte c1 = g_LetterToCharAmino[j];
			for (uint k = 0; k < 20; ++k)
				{
				byte c2 = g_LetterToCharAmino[k];
				float f2 = MP2.m_FreqVec[i][k];
				float Score = f1*f2*g_SubstMx[c1][c2];
				Sum += Score;
				}
			}
		}
	float Score = Sum/MPL;
	return Score;
	}

void MPCluster::LogPair(const MotifProfile &MP1,
  const MotifProfile &MP2) const
	{
	string Logo1;
	string Logo2;
	MP1.GetLogo(Logo1);
	MP2.GetLogo(Logo2);

	float Score = GetScore(MP1, MP2);

	asserta(Logo1.size() == MPL);
	asserta(Logo2.size() == MPL);
	string Annot;
	for (uint i = 0; i < MPL; ++i)
		{
		char a = GetAnnotChar(Logo1[i], Logo2[i]);
		Annot.push_back(a);
		}
	Log("\n");
	Log("%s\n", Logo1.c_str());
	Log("%s\n", Annot.c_str());
	Log("%s\n", Logo2.c_str());
	Log("Score = %.3f\n", Score);
	}

void MotifProfile::GetLettersFromSeq(const string &Seq,
  vector<uint> &Letters)
	{
	Letters.clear();
	asserta(SIZE(Seq) == AL+BL+CL);

	for (uint i = 0; i < MPL; ++i)
		{
		byte c = (byte) Seq[i];
		assert(c != 0);
		uint Letter = g_CharToLetterAmino[c];
		Letters.push_back(Letter);
		}
	}

void MotifProfile::FromSeq(uint SeqIndex, uint PosA, uint PosB, uint PosC,
  const string &Seq)
	{
	Clear();
	m_InputSeqIndex = SeqIndex;
	m_PosA = PosA;
	m_PosB = PosB;
	m_PosC = PosC;

	vector<uint> Letters;
	GetLettersFromSeq(Seq, Letters);

	asserta(SIZE(Letters) == MPL);
	for (uint i = 0; i < MPL; ++i)
		{
		byte c = (byte) Seq[i];
		assert(c != 0);
		uint Letter = g_CharToLetterAmino[c];
		if (Letter >= 20)
			{
			for (uint j = 0; j < 20; ++j)
				m_FreqVec[i][j] = 1.0f/20.0f;
			continue;
			}
		m_FreqVec[i][Letter] = 1;
		}
	}

void MotifProfile::FromPSSMs(const PSSM &PA, const PSSM &PB, const PSSM &PC)
	{
	Clear();
	vector<float> Probs;
	for (uint i = 0; i < 12; ++i)
		{
		PA.CalcProbs(i, Probs);
		for (uint j = 0; j < 20; ++j)
			m_FreqVec[i][j] = Probs[j];
		}
	for (uint i = 0; i < 14; ++i)
		{
		PB.CalcProbs(i, Probs);
		for (uint j = 0; j < 20; ++j)
			m_FreqVec[12 + i][j] = Probs[j];
		}
	for (uint i = 0; i < 8; ++i)
		{
		PC.CalcProbs(i, Probs);
		for (uint j = 0; j < 20; ++j)
			m_FreqVec[12 + 14 + i][j] = Probs[j];
		}
	}

void MotifProfile::FromSeqs(const vector<string> &InputSeqs)
	{
	Clear();
	m_InputSeqIndex = UINT_MAX;

	const uint InputSeqCount = SIZE(InputSeqs);
	set<string> UniqueSeqs;
	for (uint SeqIndex = 0; SeqIndex < InputSeqCount; ++SeqIndex)
		{
		const string &Seq = InputSeqs[SeqIndex];
		asserta(SIZE(Seq) == MPL);
		UniqueSeqs.insert(Seq);
		}

	vector<string> Seqs;
	for (set<string>::const_iterator p = UniqueSeqs.begin();
	  p != UniqueSeqs.end(); ++p)
		Seqs.push_back(*p);

	const uint SeqCount = SIZE(Seqs);

	asserta(SeqCount > 0);
	vector<uint> Letters;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Seq = Seqs[SeqIndex];
		GetLettersFromSeq(Seq, Letters);
		asserta(SIZE(Letters) == MPL);
		for (uint i = 0; i < MPL; ++i)
			{
			byte c = (byte) Seq[i];
			assert(c != 0);
			uint Letter = g_CharToLetterAmino[c];
			if (Letter >= 20)
				{
				for (uint j = 0; j < 20; ++j)
					m_FreqVec[i][j] += 1.0f/(20.0f*SeqCount);
				continue;
				}
			m_FreqVec[i][Letter] += 1.0f/SeqCount;
			}
		}
	ValidateFreqs();
	}

void MotifProfile::GetMaxLetter(uint i, uint &MaxLetter, float &MaxFreq) const
	{
	asserta(i < SIZE(m_FreqVec));
	const vector<float> &fs = m_FreqVec[i];
	MaxFreq = 0;
	MaxLetter = 0;
	for (uint j = 0; j < 20; ++j)
		{
		if (fs[j] > MaxFreq)
			{
			MaxFreq = fs[j];
			MaxLetter = j;
			}
		}
	}

void MotifProfile::GetConsSeq(string &ConsSeq) const
	{
	ConsSeq.clear();
	double Sum = 0;
	for (uint i = 0; i < MPL; ++i)
		{
		uint Letter;
		float Freq;
		GetMaxLetter(i, Letter, Freq);
		char c = (char) g_LetterToCharAmino[Letter];
		if (Freq < 0.1)
			c = 'X';
		else if (Freq < 0.25)
			c = tolower(c);
		ConsSeq.push_back(c);
		}
	}

void MotifProfile::GetLogo(string &Logo) const
	{
	Logo.clear();
	double Sum = 0;
	for (uint i = 0; i < MPL; ++i)
		{
		if (i == AL || i == AL+BL)
			Logo += "  ";
		uint Letter;
		float Freq;
		GetMaxLetter(i, Letter, Freq);
		char c = (char) g_LetterToCharAmino[Letter];
		if (Freq < 0.1)
			c = 'X';
		else if (Freq < 0.25)
			c = '.';
		else if (Freq < 0.5)
			c = tolower(c);
		Logo.push_back(c);
		}
	}

void MotifProfile::LogMe() const
	{
#if 0
	vector<uint> Order(20);
	for (uint i = 0; i < MPL; ++i)
		{
		const vector<float> &fs = m_FreqVec[i];
		QuickSortOrderDesc(fs.data(), 20, Order.data());
		Log("%c[%02u]  ", PosToMotif(i), PosToOffset(i));
		for (uint k = 0; k < 3; ++k)
			{
			uint Letter = Order[k];
			float f = fs[Letter];
			if (f < 1e-6)
				continue;
			if (k > 0)
				Log("  ");
			Log("%c/%.4f", g_LetterToCharAmino[Letter], fs[Letter]);
			}
		Log("\n");
		}
#endif
	string Logo;
	GetLogo(Logo);
	Log("%s\n", Logo.c_str());
	}

void MPCluster::GreedyCluster(const vector<MotifProfile *> &Input,
  float MinScore)
	{
	Clear();
	m_Input = &Input;
	m_MinScore = MinScore;
	const uint InputCount = SIZE(Input);
	for (uint i = 0; i < InputCount; ++i)
		m_PendingIndexes.insert(i);

	uint DoneCount = 0;
	vector<uint> ClusterSizes;
	uint Counter = 0;
	while (m_PendingIndexes.size() > 0)
		{
		ProgressStep(DoneCount, InputCount + 1, "Clustering");
		uint SizeStart = SIZE(m_PendingIndexes);
//		uint CentroidIndex = *m_PendingIndexes.begin();
		uint CentroidIndex = GetNextGreedyCentroid();
		++DoneCount;
		m_PendingIndexes.erase(CentroidIndex);
		m_CentroidIndexes.push_back(CentroidIndex);

		MotifProfile &CP = GetProfile(CentroidIndex);
		vector<uint> MemberIndexes;
		MemberIndexes.push_back(CentroidIndex);
		for (set<uint>::const_iterator p = m_PendingIndexes.begin();
		  p != m_PendingIndexes.end(); ++p)
			{
			++Counter;
			if (Counter%1000 == 0)
				ProgressStep(DoneCount, InputCount + 1, "Clustering");

			uint Index = *p;
			MotifProfile &P = GetProfile(Index);
			float Score = GetScore(CP, P);
			if (Score >= m_MinScore)
				{
				++DoneCount;
				MemberIndexes.push_back(Index);
				}
			}

		uint Size = SIZE(MemberIndexes);
		ClusterSizes.push_back(Size);
		m_CentroidIndexToMemberIndexes.push_back(MemberIndexes);
		for (uint i = 0; i < SIZE(MemberIndexes); ++i)
			{
			uint Index = MemberIndexes[i];
			m_PendingIndexes.erase(Index);
			}
		uint SizeEnd = SIZE(m_PendingIndexes);
		asserta(SizeEnd < SizeStart);
		uint ClusterCount = SIZE(m_CentroidIndexes);
		asserta(SIZE(m_CentroidIndexToMemberIndexes) == ClusterCount);
		if (ClusterCount == m_TopN)
			break;
		asserta(m_TopN == UINT_MAX || ClusterCount < m_TopN);
		}
	ProgressStep(InputCount, InputCount + 1, "Clustering");

	uint ClusterCount = SIZE(ClusterSizes);
	m_ClusterSizeOrder.resize(ClusterCount);
	QuickSortOrderDesc(ClusterSizes.data(), ClusterCount,
	  m_ClusterSizeOrder.data());
	}

void MPCluster::LogCluster(uint ClusterIndex) const
	{
	asserta(ClusterIndex < SIZE(m_CentroidIndexes));
	asserta(ClusterIndex < SIZE(m_CentroidIndexToMemberIndexes));

	uint CentroidIndex = m_CentroidIndexes[ClusterIndex];
	const vector<uint> &MemberIndexes = m_CentroidIndexToMemberIndexes[ClusterIndex];
	const uint Size = SIZE(MemberIndexes);

	Log("\n");
	Log("Cluster %u, size %u\n", ClusterIndex, Size);

	const MotifProfile &CP = GetProfile(CentroidIndex);
	string Logo;
	CP.GetLogo(Logo);
	//Log("[Centroid]  %s\n", Logo.c_str());
	vector<float> Scores;
	for (uint i = 0; i < Size; ++i)
		{
		uint Index = MemberIndexes[i];
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		Scores.push_back(Score);
		}
	vector<uint> Order(Size);
	QuickSortOrderDesc(Scores.data(), Size, Order.data());

	for (uint k = 0; k < Size; ++k)
		{
		uint i = Order[k];
		uint Index = MemberIndexes[i];
		if (Index == CentroidIndex)
			continue;
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		P.GetLogo(Logo);
//		Log("[%8.4f]  %s\n", Score, Logo.c_str());
		}
	}

void MPCluster::LogClusters() const
	{
	const uint ClusterCount = SIZE(m_CentroidIndexes);
	asserta(SIZE(m_CentroidIndexToMemberIndexes) == ClusterCount);
	Log("\n");
	Log("%u clusters\n", ClusterCount);
	for (uint i = 0; i < ClusterCount; ++i)
		LogCluster(i);
	}

void MPCluster::WriteFasta(const MotifProfile &MP) const
	{
	if (m_FLSeqDB.GetSeqCount() == 0)
		return;

	uint SeqIndex = MP.m_InputSeqIndex;
	asserta(SeqIndex != UINT_MAX);

	const string &Label = m_FLSeqDB.GetLabel(SeqIndex);
	const string &Seq = m_FLSeqDB.GetSeq(SeqIndex);
	uint L = SIZE(Seq);

	uint PosA = MP.m_PosA;
	uint PosB = MP.m_PosB;
	uint PosC = MP.m_PosC;

	asserta(PosA + AL <= L);
	asserta(PosB + BL <= L);
	asserta(PosC + CL <= L);

	string SeqA = Seq.substr(PosA, AL);
	string SeqB = Seq.substr(PosB, BL);
	string SeqC = Seq.substr(PosC, CL);

	string OutLabel = Label;
	Psa(OutLabel, " A:%u:%s", PosA + 1, SeqA.c_str());
	Psa(OutLabel, " B:%u:%s", PosB + 1, SeqB.c_str());
	Psa(OutLabel, " C:%u:%s", PosC + 1, SeqC.c_str());

	string Sep;
	if (optset_sep)
		Sep = string(opt_sep);

	bool Perm = (PosC < PosA);
	string Left, V1, V2, Right;
	string OutSeq;
	if (Perm)
		{
		if (PosC > 0)
			Left = Seq.substr(0, PosC);

		uint V1Lo = PosC + CL;
		uint V1Hi = PosA - 1;
		if (V1Lo < V1Hi)
			{
			uint V1L = V1Hi - V1Lo + 1;
			V1 = Seq.substr(V1Lo, V1L);
			}

		uint V2Lo = PosA + AL;
		uint V2Hi = PosB - 1;
		if (V2Lo < V2Hi)
			{
			uint V2L = V2Hi - V2Lo + 1;
			V2 = Seq.substr(V2Lo, V2L);
			}

		if (PosB+BL < L)
			Right = Seq.substr(PosB+BL, string::npos);

		string CheckSeq = Left + SeqC + V1 + SeqA + V2 + SeqB + Right;
		if (CheckSeq != Seq)
			{
			Log("\n");
			Log("%s\n", CheckSeq.c_str());
			Log("%s\n", Seq.c_str());
			Die("Seq2!=Seq");
			}

		TrimLeft(Left);
		TrimRight(Right);
		if (Sep == "1")
			{
			OutSeq = Left;

			OutSeq += "5";
			OutSeq += SeqC;
			OutSeq += "6";

			OutSeq += V1;

			OutSeq += "1";
			OutSeq += SeqA;
			OutSeq += "2";

			OutSeq += V2;

			OutSeq += "3";
			OutSeq += SeqB;
			OutSeq += "4";

			OutSeq += Right;
			}
		else
			{
			OutSeq = Left;

			OutSeq += Sep;
			OutSeq += SeqC;
			OutSeq += Sep;

			OutSeq += V1;

			OutSeq += Sep;
			OutSeq += SeqA;
			OutSeq += Sep;

			OutSeq += V2;

			OutSeq += Sep;
			OutSeq += SeqB;
			OutSeq += Sep;

			OutSeq += Right;
			}
		SeqToFasta(g_fcluster_fasta, OutLabel.c_str(), OutSeq.c_str(), SIZE(OutSeq));
		SeqToFasta(g_fcluster_fasta_cab, OutLabel.c_str(), OutSeq.c_str(), SIZE(OutSeq));

		string PP = SeqC + V1 + SeqA + V2 + SeqB;
		SeqToFasta(g_fCAB_PP, Label.c_str(), PP.c_str());

		SeqToFasta(g_fCAB_A, Label.c_str(), SeqA.c_str());
		SeqToFasta(g_fCAB_B, Label.c_str(), SeqB.c_str());
		SeqToFasta(g_fCAB_C, Label.c_str(), SeqC.c_str());

		if (!V1.empty())
			SeqToFasta(g_fCAB_V1, Label.c_str(), V1.c_str());
		if (!V2.empty())
			SeqToFasta(g_fCAB_V2, Label.c_str(), V2.c_str());

		if (!Left.empty())
			SeqToFasta(g_fCAB_Left, Label.c_str(), Left.c_str());
		if (!Right.empty())
			SeqToFasta(g_fCAB_Right, Label.c_str(), Right.c_str());
		}
	else
		{
		if (PosA > 0)
			Left = Seq.substr(0, PosA);

		uint V1Lo = PosA + AL;
		uint V1Hi = PosB - 1;
		if (V1Lo < V1Hi)
			{
			uint V1L = V1Hi - V1Lo + 1;
			V1 = Seq.substr(V1Lo, V1L);
			}

		uint V2Lo = PosB + BL;
		uint V2Hi = PosC - 1;
		if (V2Lo < V2Hi)
			{
			uint V2L = V2Hi - V2Lo + 1;
			V2 = Seq.substr(V2Lo, V2L);
			}

		if (PosC+CL < L)
			Right = Seq.substr(PosC+CL, string::npos);

		string CheckSeq = Left + SeqA + V1 + SeqB + V2 + SeqC + Right;
		if (CheckSeq != Seq)
			{
			Log("\n");
			Log("%s\n", CheckSeq.c_str());
			Log("%s\n", Seq.c_str());
			Die("CheckSeq!=Seq");
			}

		TrimLeft(Left);
		TrimRight(Right);
		if (Sep == "1")
			{
			OutSeq = Left;

			OutSeq += "1";
			OutSeq += SeqA;
			OutSeq += "2";

			OutSeq += V1;

			OutSeq += "3";
			OutSeq += SeqB;
			OutSeq += "4";

			OutSeq += V2;

			OutSeq += "5";
			OutSeq += SeqC;
			OutSeq += "6";

			OutSeq += Right;
			}
		else
			{
			OutSeq = Left;

			OutSeq += Sep;
			OutSeq += SeqA;
			OutSeq += Sep;

			OutSeq += V1;

			OutSeq += Sep;
			OutSeq += SeqB;
			OutSeq += Sep;

			OutSeq += V2;

			OutSeq += Sep;
			OutSeq += SeqC;
			OutSeq += Sep;

			OutSeq += Right;
			}
		SeqToFasta(g_fcluster_fasta, OutLabel.c_str(), OutSeq.c_str(), SIZE(OutSeq));
		SeqToFasta(g_fcluster_fasta_abc, OutLabel.c_str(), OutSeq.c_str(), SIZE(OutSeq));

		string PP = SeqA + V1 + SeqB + V2 + SeqC;
		SeqToFasta(g_fABC_PP, Label.c_str(), PP.c_str());

		SeqToFasta(g_fABC_A, Label.c_str(), SeqA.c_str());
		SeqToFasta(g_fABC_B, Label.c_str(), SeqB.c_str());
		SeqToFasta(g_fABC_C, Label.c_str(), SeqC.c_str());

		if (!V1.empty())
			SeqToFasta(g_fABC_V1, Label.c_str(), V1.c_str());
		if (!V2.empty())
			SeqToFasta(g_fABC_V2, Label.c_str(), V2.c_str());

		if (!Left.empty())
			SeqToFasta(g_fABC_Left, Label.c_str(), Left.c_str());
		if (!Right.empty())
			SeqToFasta(g_fABC_Right, Label.c_str(), Right.c_str());
		}
	}

void MPCluster::WriteCluster(uint OrderIndex) const
	{
	asserta(OrderIndex < SIZE(m_ClusterSizeOrder));
	uint ClusterIndex = m_ClusterSizeOrder[OrderIndex];
	asserta(ClusterIndex < SIZE(m_CentroidIndexes));
	asserta(ClusterIndex < SIZE(m_CentroidIndexToMemberIndexes));

	uint CentroidIndex = m_CentroidIndexes[ClusterIndex];
	const vector<uint> &MemberIndexes = m_CentroidIndexToMemberIndexes[ClusterIndex];
	const uint Size = SIZE(MemberIndexes);

	const MotifProfile &CP = GetProfile(CentroidIndex);
	string CentroidSeq;
	CP.GetConsSeq(CentroidSeq);
	const char *ABC = (CP.m_PosA > CP.m_PosC ? "CAB" : "ABC");
		
	Pf(g_fcluster_tsv, "C\t%u\t%u\t%s\t%s\n", 
	  OrderIndex, Size, CentroidSeq.c_str(), ABC);
	WriteFasta(CP);

	vector<float> Scores;
	for (uint i = 0; i < Size; ++i)
		{
		uint Index = MemberIndexes[i];
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		Scores.push_back(Score);
		}
	vector<uint> Order(Size);
	QuickSortOrderDesc(Scores.data(), Size, Order.data());

	for (uint k = 0; k < Size; ++k)
		{
		uint i = Order[k];
		uint Index = MemberIndexes[i];
		if (Index == CentroidIndex)
			continue;
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		string ConsSeq;
		P.GetConsSeq(ConsSeq);
		const char *ABC = (P.m_PosA > P.m_PosC ? "CAB" : "ABC");
		Pf(g_fcluster_tsv, "M\t%u\t%.2f\t%s\t%s\n",
		  OrderIndex, Score, ConsSeq.c_str(), ABC);
		WriteFasta(P);
		}
	}

void MPCluster::WriteOutput() const
	{
	const uint ClusterCount = SIZE(m_CentroidIndexes);
	asserta(SIZE(m_CentroidIndexToMemberIndexes) == ClusterCount);
	uint SumSize = 0;
	for (uint OrderIndex = 0; OrderIndex < ClusterCount; ++OrderIndex)
		{
		asserta(OrderIndex < SIZE(m_ClusterSizeOrder));
		uint ClusterIndex = m_ClusterSizeOrder[OrderIndex];
		const vector<uint> &MemberIndexes = m_CentroidIndexToMemberIndexes[ClusterIndex];
		const uint Size = SIZE(MemberIndexes);
		SumSize += Size;
		Log("Cluster %5u  total  %7u\n", OrderIndex, SumSize);

		WriteCluster(OrderIndex);
		}
	}

void MPCluster::FindNN(uint &Index1, uint &Index2) const
	{
	asserta(SIZE(m_PendingIndexes) >= 2);
	Index1 = UINT_MAX;
	Index2 = UINT_MAX;
	float MaxScore = -9999;
	for (set<uint>::const_iterator p = m_PendingIndexes.begin();
	  p != m_PendingIndexes.end(); ++p)
		{
		uint i1 = *p;
		set<uint>::const_iterator q = p;
		for (;;)
			{
			++q;
			if (q == m_PendingIndexes.end())
				break;
			uint i2 = *q;
			float Score = GetScoreNNPair(i1, i2);
			if (Score > MaxScore)
				{
				MaxScore = Score;
				Index1 = i1;
				Index2 = i2;
				}
			}
		}
	}

void MPCluster::LogNN() const
	{
	const uint N = SIZE(m_MPs);
	Log("%u input, %u tree MPs\n", m_InputMPCount, N);
	const uint M = SIZE(m_Parents);
	Log("%u internal nodes\n", M);
	asserta(SIZE(m_Lefts) == M);
	asserta(SIZE(m_Rights) == M);
	for (uint i = 0; i < M; ++i)
		{
		Log("[%4u]", i);
		Log("  P=%4u", m_Parents[i]);
		Log("  L=%4u", m_Lefts[i]);
		Log("  R=%4u", m_Rights[i]);
		Log("\n");
		}
	}

void MPCluster::NNNodeToNewickFile(FILE *f, uint NodeIndex) const
	{
	if (NodeIndex >= m_InputMPCount)
		{
	// Internal node
		uint Left = m_Lefts[NodeIndex - m_InputMPCount];
		uint Right = m_Rights[NodeIndex - m_InputMPCount];
		fprintf(f, "(\n");
		NNNodeToNewickFile(f, Left);
		fprintf(f, ",\n");
		NNNodeToNewickFile(f, Right);
		fprintf(f, ")\n");
		}
	else
		{
	// Leaf node
		const MotifProfile &MP = *m_MPs[NodeIndex];
		fprintf(f, "%s\n", MP.m_Name.c_str());
		}
	}

void MPCluster::NNToNewickFile(const string &FileName) const
	{
	LogNN();
	if (FileName == "")
		return;

	FILE *f = CreateStdioFile(FileName);
	asserta(SIZE(m_PendingIndexes) == 1);
	uint RootNode = *m_PendingIndexes.begin();
	NNNodeToNewickFile(f, RootNode);
	fprintf(f, ";\n");
	CloseStdioFile(f);
	}

void MPCluster::Join(uint Index1, uint Index2)
	{
	MotifProfile &P = CreateProfileNN(Index1, Index2);
	uint Parent = SIZE(m_MPs);
	uint LeftSize = m_Sizes[Index1];
	uint RightSize = m_Sizes[Index2];
	uint ParentSize = LeftSize + RightSize;
	m_MPs.push_back(&P);
	m_Parents.push_back(Parent);
	m_Lefts.push_back(Index1);
	m_Rights.push_back(Index2);
	m_PendingIndexes.erase(Index1);
	m_PendingIndexes.erase(Index2);
	m_PendingIndexes.insert(Parent);
	m_Sizes.push_back(ParentSize);
	}

float MPCluster::GetScoreNNPair(uint i1, uint i2) const
	{
	asserta(i1 < SIZE(m_MPs));
	asserta(i2 < SIZE(m_MPs));

	const MotifProfile &MP1 = *m_MPs[i1];
	const MotifProfile &MP2 = *m_MPs[i2];
	float Score = GetScore(MP1, MP2);
	return Score;
	}

MotifProfile &MPCluster::CreateProfileNN(uint i1, uint i2) const
	{
	asserta(i1 < SIZE(m_MPs));
	asserta(i2 < SIZE(m_MPs));

	const MotifProfile &MP1 = *m_MPs[i1];
	const MotifProfile &MP2 = *m_MPs[i2];

	uint Size1 = m_Sizes[i1];
	uint Size2 = m_Sizes[i2];
	float w1 = (float) sqrt(Size1);
	float w2 = (float) sqrt(Size2);
	float w12 = w1 + w2;
	w1 /= w12;
	w2 /= w12;
	asserta(feq(w1 + w2, 1.0));

	MotifProfile &P = *new MotifProfile;
	for (uint i = 0; i < MPL; ++i)
		{
		const vector<float> &v1 = MP1.m_FreqVec[i];
		const vector<float> &v2 = MP2.m_FreqVec[i];
		vector<float> &v = P.m_FreqVec[i];

		float Sum = 0;
		for (uint j = 0; j < 20; ++j)
			{
			float f1 = v1[j];
			float f2 = v2[j];
			float f = w1*f1 + w2*f2;
			v[j] = f;
			Sum += f;
			}
		asserta(feq(Sum, 1.0f));
		}

	string Logo1;
	string Logo2;
	string LogoP;
	MP1.GetLogo(Logo1);
	MP2.GetLogo(Logo2);
	P.GetLogo(LogoP);
	
	Log("  L[%5u, %6.4f]  %s\n", Size1, w1, Logo1.c_str());
	Log("  R[%5u, %6.4f]  %s\n", Size2, w2, Logo2.c_str());

	return P;
	}

void MPCluster::NNCluster(const vector<MotifProfile *> &Input,
  float MinScore)
	{
	Clear();
	m_Input = &Input;
	const uint N = SIZE(Input);
	m_InputMPCount = N;
	for (uint i = 0; i < N; ++i)
		{
		m_MPs.push_back(Input[i]);
		m_PendingIndexes.insert(i);
		m_Sizes.push_back(1);
		}

	uint JoinCount = N - 1;
	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		ProgressStep(JoinIndex, JoinCount, "NN cluster");
		uint Index1, Index2;
		FindNN(Index1, Index2);
		Log("\nJoin %u\n", JoinIndex);
		Join(Index1, Index2);
		}
	}
void MPCluster::GetMembers(uint Centroid, vector<uint> &Members) const
	{
	Members.clear();
	MotifProfile &CP = GetProfile(Centroid);
	for (set<uint>::const_iterator p = m_PendingIndexes.begin();
		p != m_PendingIndexes.end(); ++p)
		{
		uint Index = *p;
		MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		if (Score >= m_MinScore)
			Members.push_back(Index);
		}
	}

void MPCluster::GetRandomPending(uint n, vector<uint> &v) const
	{
	uint M = SIZE(m_PendingIndexes);
	uint k = M/n;
	if (k <= 3)
		{
		for (set<uint>::const_iterator p = m_PendingIndexes.begin();
		  p != m_PendingIndexes.end(); ++p)
			{
			v.push_back(*p);
			if (SIZE(v) == n)
				break;
			}
		return;
		}

	for (set<uint>::const_iterator p = m_PendingIndexes.begin();
		p != m_PendingIndexes.end(); ++p)
		{
		if (randu32()%k == 0)
			v.push_back(*p);
		if (SIZE(v) == n)
			break;
		}
	if (v.empty())
		v.push_back(*m_PendingIndexes.begin());
	}

uint MPCluster::GetBestCentroid(const vector<uint> &v) const
	{
	asserta(!v.empty());
	uint BestSize = 0;
	uint BestCentroid = v[0];
	vector<uint> Members;
	for (uint i = 0; i < SIZE(v); ++i)
		{
		GetMembers(v[i], Members);
		uint Size = SIZE(Members);
		if (Size > BestSize)
			{
			BestSize = Size;
			BestCentroid = v[i];
			}
		}
	return BestCentroid;
	}

uint MPCluster::GetNextGreedyCentroid() const
	{
	asserta(!m_PendingIndexes.empty());
	vector<uint> v;
	GetRandomPending(m_Sample1, v);
	asserta(!v.empty());
	uint Centroid = GetBestCentroid(v);
	vector<uint> Members;
	GetMembers(Centroid, Members);
	random_shuffle(Members.begin(), Members.end());
	Members.resize(m_Sample2);
	uint BestCentroid = GetBestCentroid(Members);
	return BestCentroid;
	}

SeqDB MPCluster::m_FLSeqDB;
SeqDB MPCluster::m_MotifSeqDB;
vector<uint> MPCluster::m_PosAs;
vector<uint> MPCluster::m_PosBs;
vector<uint> MPCluster::m_PosCs;

static void OnPalmHit(const RdRpSearcher &RS)
	{
	float Score = RS.m_TopPalmHit.m_Score;
	asserta(Score >= opt_min_palm_score);

	uint PosA = RS.GetMotifPos(0);
	uint PosB = RS.GetMotifPos(1);
	uint PosC = RS.GetMotifPos(2);
	asserta(PosA != UINT_MAX && PosB != UINT_MAX && PosC != UINT_MAX);

	string SeqA;
	string SeqB;
	string SeqC;
	RS.GetMotifSeq(0, SeqA);
	RS.GetMotifSeq(1, SeqB);
	RS.GetMotifSeq(2, SeqC);

	string Motifs = SeqA + SeqB + SeqC;
	asserta(Motifs.size() == 34);

	MPCluster::m_FLSeqDB.AddSeq(RS.m_QueryLabel, RS.m_QuerySeq);
	MPCluster::m_MotifSeqDB.AddSeq(RS.m_QueryLabel, Motifs);
	MPCluster::m_PosAs.push_back(PosA);
	MPCluster::m_PosBs.push_back(PosB);
	MPCluster::m_PosCs.push_back(PosC);
	}

static void ClusterMotifs(const string &InputFileName,
  const string &Strategy)
	{
	OpenSegFiles();

	SetBLOSUM62();

	asserta(optset_motif_cluster_minscore);
	const float MinScore = (float) opt_motif_cluster_minscore;

	SearchPSSMs(InputFileName, OnPalmHit);

	const uint InputCount = MPCluster::m_MotifSeqDB.GetSeqCount();
	asserta(MPCluster::m_FLSeqDB.GetSeqCount() == InputCount);
	asserta(SIZE(MPCluster::m_PosAs) == InputCount);
	asserta(SIZE(MPCluster::m_PosBs) == InputCount);
	asserta(SIZE(MPCluster::m_PosCs) == InputCount);
	ProgressLog("%u palm hits\n", InputCount);

	vector<MotifProfile *> MPs;
	for (uint SeqIndex = 0; SeqIndex < InputCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, InputCount, "Build motif db");

		const string &Motifs = MPCluster::m_MotifSeqDB.GetSeq(SeqIndex);
		uint PosA = MPCluster::m_PosAs[SeqIndex];
		uint PosB = MPCluster::m_PosBs[SeqIndex];
		uint PosC = MPCluster::m_PosCs[SeqIndex];

		MotifProfile *MP = new MotifProfile;
		MP->FromSeq(SeqIndex, PosA, PosB, PosC, Motifs);
		MPs.push_back(MP);
		}

	MPCluster MC;
	if (Strategy == "greedy")
		{
		if (optset_topn)
			MC.m_TopN = opt_topn;
		else
			MC.m_TopN = UINT_MAX;
		MC.GreedyCluster(MPs, MinScore);
		MC.LogClusters();
		ProgressLog("%u motifs, %u clusters\n",
		  InputCount, SIZE(MC.m_CentroidIndexes));
		MC.WriteOutput();
		}
	else if (Strategy == "nn")
		MC.NNCluster(MPs, MinScore);
	else
		Die("Invalid strategy '%s'", Strategy.c_str());

	CloseSegFiles();
	}

// Input is scanned by PSSMs
void cmd_cluster_motifs_greedy()
	{
	const string &InputFileName = opt_cluster_motifs_greedy;
	ClusterMotifs(InputFileName, "greedy");
	}

// Input is scanned by PSSMs
void cmd_cluster_motifs_nn()
	{
	const string &InputFileName = opt_cluster_motifs_nn;
	ClusterMotifs(InputFileName, "nn");
	}

// Input is 34aa ABC sequences
void cmd_cluster_motifs_greedy3()
	{
	const string &InputFileName = opt_cluster_motifs_greedy3;
	OpenSegFiles();

	SetBLOSUM62();

	asserta(optset_motif_cluster_minscore);
	const float MinScore = (float) opt_motif_cluster_minscore;

	MPCluster::m_MotifSeqDB.FromFasta(InputFileName);

	const uint InputCount = MPCluster::m_MotifSeqDB.GetSeqCount();
	ProgressLog("%u input motifs\n", InputCount);

	vector<MotifProfile *> MPs;
	for (uint SeqIndex = 0; SeqIndex < InputCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, InputCount, "Build motif db");

		const string &Motifs = MPCluster::m_MotifSeqDB.GetSeq(SeqIndex);
		uint L = SIZE(Motifs);
		if (L != AL+BL+CL)
			Die("Sequence length %u != 32 >%s",
			  L, MPCluster::m_MotifSeqDB.GetLabel(SeqIndex).c_str());

		uint PosA = 0;
		uint PosB = AL;
		uint PosC = AL+BL;

		MotifProfile *MP = new MotifProfile;
		MP->FromSeq(SeqIndex, PosA, PosB, PosC, Motifs);
		MPs.push_back(MP);
		}

	MPCluster MC;
	if (optset_topn)
		MC.m_TopN = opt_topn;
	else
		MC.m_TopN = UINT_MAX;
	MC.GreedyCluster(MPs, MinScore);
	MC.LogClusters();
	ProgressLog("%u motifs, %u clusters\n",
		InputCount, SIZE(MC.m_CentroidIndexes));
	MC.WriteOutput();

	CloseSegFiles();
	}
