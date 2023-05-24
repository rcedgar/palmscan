#include "myutils.h"
#include "mx.h"

float SW(const Mx<float> &SMx,
  Mx<float> &a_FwdM,
  Mx<float> &a_FwdD,
  Mx<float> &a_FwdI,
  Mx<char> &a_TBM,
  Mx<char> &a_TBD,
  Mx<char> &a_TBI,
  uint &Starti, uint &Startj, string &Path);

float GetBlosum62Score(char a, char b);

void cmd_test_sw()
	{
	const string &A = "SMITHWATERMAN";
	const string &B = "SMITHATERMAN";

	const uint LA = SIZE(A);
	const uint LB = SIZE(B);

	Mx<float> SMx;
	SMx.Alloc("SMx", LA, LB);
	float **Sim = SMx.GetData();

	Mx<float> a_FwdM;
	Mx<float> a_FwdD;
	Mx<float> a_FwdI;
	Mx<char> a_TBM;
	Mx<char> a_TBD;
	Mx<char> a_TBI;

	for (unsigned i = 0; i < LA; ++i)
		{
		char a = A[i];
		float *SimRow = Sim[i];
		for (unsigned j = 0; j < LB; ++j)
			{
			char b = B[j];
			SimRow[j] = GetBlosum62Score(a, b);
			}
		}

	uint Starti, Startj;
	string Path;
	float Score = SW(SMx, a_FwdM, a_FwdD, a_FwdI,
	  a_TBM, a_TBD, a_TBI,
	  Starti, Startj, Path);

	Log("Score = %.3g, Path %s\n", Score, Path.c_str());

	const uint ColCount = SIZE(Path);
	uint i = Starti;
	uint j = Startj;
	string ARow;
	string BRow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			{
			ARow += A[i];
			++i;
			}
		else
			ARow += '-';

		if (c == 'M' || c == 'I')
			{
			BRow += B[j];
			++j;
			}
		else
			BRow += '-';
		}
	Log("A  %s\n", ARow.c_str());
	Log("B  %s\n", BRow.c_str());
	}
