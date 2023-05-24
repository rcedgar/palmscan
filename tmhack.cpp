#include "myutils.h"
#include "tma.h"

static void XFormPt_TM(
  const vector<double> &Pt,
  const vector<double> &t,
  const vector<vector<double> > &u,
  vector<double> &XPt)
	{
	XPt.resize(3);
	double x = Pt[0];
	double y = Pt[1];
	double z = Pt[2];

	XPt[0] = t[0] + u[0][0]*x + u[0][1]*y + u[0][2]*z;
	XPt[1] = t[1] + u[1][0]*x + u[1][1]*y + u[1][2]*z;
	XPt[2] = t[2] + u[2][0]*x + u[2][1]*y + u[2][2]*z;
	}

void cmd_tmhack()
	{
	PDBChain ChainA;
	PDBChain ChainB;

	ChainA.FromCal("2ijd-poli.cal");
	ChainB.FromCal("6u5o_L-mono.cal");

/***
==> Code for rotating Structure A from (x,y,z) to (X,Y,Z):
for(i=0; i<L; i++)
{
   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];
   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];
   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];
}

***/
	vector<double> t;
	vector<vector<double> > u(3);

#if 0
	{
/***
m               t[m]        u[m][0]        u[m][1]        u[m][2]
0      -0.7258910570   0.9766654312  -0.0901997765  -0.1949067361
1      -0.4388934162  -0.1703278250   0.2274721337  -0.9587725801
2      -0.9047675209   0.1308169236   0.9695980759   0.2068006378
***/
	t.push_back(-0.7258910570);
	t.push_back(-0.4388934162);
	t.push_back(-0.9047675209);

	u[0].push_back( 0.9766654312); u[0].push_back(-0.0901997765); u[0].push_back(-0.1949067361);
	u[1].push_back(-0.1703278250); u[1].push_back( 0.2274721337); u[1].push_back(-0.9587725801);
	u[2].push_back( 0.1308169236); u[2].push_back( 0.9695980759); u[2].push_back( 0.2068006378);
	}
#endif

#if 1
	{
/***
m               t[m]        u[m][0]        u[m][1]        u[m][2]
0       0.7949676501   0.9764519473  -0.1707121042   0.1319051629
1       0.8819205160  -0.0904220553   0.2312715589   0.9686781292
2      -0.3991706881  -0.1958709944  -0.9577947815   0.2103894247
***/
	t.push_back( 0.7949676501);
	t.push_back( 0.8819205160);
	t.push_back(-0.3991706881);

	u[0].push_back( 0.9764519473); u[0].push_back(-0.1707121042); u[0].push_back(0.1319051629);
	u[1].push_back(-0.0904220553); u[1].push_back( 0.2312715589); u[1].push_back(0.9686781292);
	u[2].push_back(-0.1958709944); u[2].push_back(-0.9577947815); u[2].push_back(0.2103894247);
	}
#endif

	const uint LA = ChainA.GetSeqLength();
	for (uint Pos = 0; Pos < LA; ++Pos)
		{
		vector<double> Pt;
		ChainA.GetPt(Pos, Pt);

		vector<double> XPt;
		XFormPt_TM(Pt, t, u, XPt);

		ChainA.m_Xs[Pos] = XPt[0];
		ChainA.m_Ys[Pos] = XPt[1];
		ChainA.m_Zs[Pos] = XPt[2];
		}

	ChainA.LogMe(true);
	}
