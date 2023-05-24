#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "sort.h"

void GetCavityLatLong(double x, double y, double z,
  double &r, double &lat_deg, double &long_deg);

const double LAT_D = 8;
const double LONG_D = 16;

static bool Hides(double Lat_j, double Long_j,
  double Lat_i, double Long_i)
	{
	double MaxLat = max(Lat_i, Lat_j);
	double MinLat = min(Lat_i, Lat_j);
	double dLat = MaxLat - MinLat;
	if (dLat > 90)
		dLat = 180 - dLat;

	double MaxLong = max(Long_i, Long_j);
	double MinLong = min(Long_i, Long_j);
	double dLong = MaxLong - MinLong;
	if (dLong > 180)
		dLong = 360 - dLong;

	bool h = dLat <= LAT_D && dLong <= LONG_D;
	if (h)
		Log("");//@@
	return h;
	}

void RemoveHidden(const PDBChain &Chain, 
  const vector<uint> &PosVec, 
  vector<uint> &PosVec2)
	{
	const uint N = SIZE(PosVec);
	vector<double> Rs;
	vector<double> Lats;
	vector<double> Longs;
	vector<bool> Hiddens;
	for (uint i = 0; i < N; ++i)
		{
		uint Pos = PosVec[i];
		vector<double> Pt;
		Chain.GetPt(Pos, Pt);

		double R, Lat, Long;
		GetCavityLatLong(Pt[X], Pt[Y], Pt[Z], R, Lat, Long);
		Rs.push_back(R);
		Lats.push_back(Lat);
		Longs.push_back(Long);
		Hiddens.push_back(false);
		}

	vector<uint> Order(N);
	QuickSortOrderDesc<double>(Rs.data(), N, Order.data());

	double LastR = DBL_MAX;
	uint HideCount = 0;
	for (uint k = 0; k < N; ++k)
		{
		uint i = Order[k];
		asserta(!Hiddens[i]);

		double R_i = Rs[i];
		asserta(R_i <= LastR);
		LastR = R_i;

		double Lat_i = Lats[i];
		double Long_i = Longs[i];

		for (uint j = 0; j < N; ++j)
			{
			if (j == i)
				continue;

			double R_j = Rs[j];
			double Lat_j = Lats[j];
			double Long_j = Longs[j];
			if (R_j >= R_i)
				continue;

			if (Hides(Lat_j, Long_j, Lat_i, Long_i))
				{
				Hiddens[i] = true;
				++HideCount;
				continue;
				}
			}
		}

	ProgressLog("Hidden %u\n", HideCount);
	for (uint i = 0; i < N; ++i)
		{
		if (!Hiddens[i])
			PosVec2.push_back(PosVec[i]);
		}
	}
