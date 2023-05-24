#include "myutils.h"
#include "quarts.h"
#include "outputfiles.h"
#include <map>
#include <list>
#include <set>

static map<string, size_t> LabelToIndex;
static vector<string> Labels;
static vector<vector<double> > DistMx;

static vector<size_t> IndexToCluster;
static vector<size_t> ClusterSizes;
static vector<set<size_t> > Clusters;
static set<size_t> UnassignedIndexes;
static set<size_t> AssignedIndexes;

static bool CmpClusterSize(size_t C1, size_t C2)
	{
	size_t Size1 = ClusterSizes[C1];
	size_t Size2 = ClusterSizes[C2];
	return Size1 > Size2;
	}

static size_t GetIndex(const string &Label)
	{
	auto p = LabelToIndex.find(Label);
	if (p == LabelToIndex.end())
		{
		size_t Index = LabelToIndex.size();
		Labels.push_back(Label);
		LabelToIndex[Label] = Index;
		return Index;
		}
	size_t Index = p->second;
	return Index;
	}

static size_t GetCentroid(size_t ClusterIndex, QuartsDouble &DistQuarts)
	{
	DistQuarts.Clear();

	asserta(ClusterIndex < Clusters.size());
	const set<size_t> &Cluster = Clusters[ClusterIndex];

	size_t BestIndex = *Cluster.begin();
	if (Cluster.size() <= 2)
		return BestIndex;

	double BestSum = 0;
	for (auto Index1 : Cluster)
		{
		double Sum = 0;
		for (auto Index2 : Cluster)
			{
			double d12 = DistMx[Index1][Index2];
			if (d12 != DBL_MAX)
				Sum += d12;
			}
		if (Sum > BestSum)
			{
			BestIndex = Index1;
			BestSum = Sum;
			}
		}

	const size_t Centroid = BestIndex;

	vector<double> Dists;
	for (auto Index : Cluster)
		{
		if (Index != Centroid)
			{
			double d = DistMx[Centroid][Index];
			if (d != DBL_MAX)
				Dists.push_back(d);
			}
		}

	GetQuartsDouble(Dists, DistQuarts);
	return Centroid;
	}

void cmd_cluster_cl()
	{
	const string &DistMxFileName = opt_cluster_cl;
	asserta(optset_threshold);
	double Threshold = opt_threshold;

	string Line;
	FILE *f = OpenStdioFile(DistMxFileName);
	vector<string> Fields;
	uint PairCount = 0;
	while (ReadLineStdioFile(f, Line))
		{
		if (PairCount%10000 == 0)
			Progress("Pass 1 (%s, %u)\r",
			  IntToStr(PairCount), (uint) Labels.size());
		++PairCount;

		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		double Dist = StrToFloat(Fields[2].c_str());

		size_t Index1 = GetIndex(Label1);
		size_t Index2 = GetIndex(Label2);
		}
	Progress("Pass 1 (%u, %u)\n", PairCount, (uint) Labels.size());

	SetStdioFilePos64(f, 0);

	size_t LabelCount = Labels.size();
	DistMx.resize(LabelCount);
	for (uint i = 0; i < LabelCount; ++i)
		{
		ProgressStep(i, LabelCount, "Init mem");
		DistMx[i].resize(0);
		DistMx[i].resize(LabelCount, DBL_MAX);
		}

	uint PairCount2 = 0;
	while (ReadLineStdioFile(f, Line))
		{
		if (PairCount2%10000 == 0)
			Progress("Pass 2 (%s)\r", IntToStr(PairCount2));
		++PairCount2;

		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		double Dist = StrToFloat(Fields[2].c_str());

		size_t Index1 = GetIndex(Label1);
		size_t Index2 = GetIndex(Label2);

		DistMx[Index1][Index2] = Dist;
		}
	Progress("Pass 2 (%u, %u)\n", PairCount2, (uint) Labels.size());
	asserta(PairCount2 == PairCount);

	for (size_t i = 0; i < LabelCount; ++i)
		{
		ProgressStep(i, LabelCount, "Fix up matrix");
		DistMx[i][i] = DBL_MAX;
		for (size_t j = 0; j < i; ++j)
			{
			double dij = DistMx[i][j];
			double dji = DistMx[j][i];
			if (dij == DBL_MAX && dji == DBL_MAX)
				continue;
			else if (dij == DBL_MAX)
				DistMx[i][j] = dji;
			else if (dji == DBL_MAX)
				DistMx[j][i] = dij;
			else
				{
				double d = (dij + dji)/2.0;
				DistMx[j][i] = d;
				DistMx[i][j] = d;
				}
			}
		}

	for (size_t i = 0; i < LabelCount; ++i)
		{
		ProgressStep(i, LabelCount, "Check symmetry");
		DistMx[i][i] = DBL_MAX;
		for (size_t j = 0; j < i; ++j)
			asserta(DistMx[i][j] == DistMx[j][i]);
		}

	Progress("Init more mem...");
	IndexToCluster.clear();
	IndexToCluster.resize(LabelCount, UINT_MAX);

	UnassignedIndexes.clear();
	for (size_t i = 0; i < LabelCount; ++i)
		UnassignedIndexes.insert(i);
	Progress(" done.\n");

	size_t SumClusterSizes = 0;
	ClusterSizes.clear();
	uint ClusterCount10 = 0;
	uint ClusterCount1 = 0;
	while (!UnassignedIndexes.empty())
		{
		size_t ClusterIndex = Clusters.size();
		size_t UnassignedCount = SIZE(UnassignedIndexes);

		Progress("%u clusters, %u unassigned\r", ClusterIndex, UnassignedCount);

		Clusters.resize(ClusterIndex+1);
		set<size_t> &Cluster = Clusters.back();
		size_t Seed = *UnassignedIndexes.begin();
		UnassignedIndexes.erase(Seed);

		set<size_t> PendingIndexesThisCluster;
		PendingIndexesThisCluster.insert(Seed);
		while (!PendingIndexesThisCluster.empty())
			{
			size_t Index = *PendingIndexesThisCluster.begin();
			PendingIndexesThisCluster.erase(Index);
			Cluster.insert(Index);
			asserta(AssignedIndexes.find(Index) == AssignedIndexes.end());

			AssignedIndexes.insert(Index);
			IndexToCluster[Index] = ClusterIndex;
			UnassignedIndexes.erase(Index);

			for (auto p : UnassignedIndexes)
				{
				size_t Index2 = p;
				double d = DistMx[Index][Index2];
				if (d != DBL_MAX && d >= Threshold)
					PendingIndexesThisCluster.insert(Index2);
				}

			for (auto p: PendingIndexesThisCluster)
				{
				size_t Index2 = p;

				Cluster.insert(Index2);
				IndexToCluster[Index2] = ClusterIndex;
				UnassignedIndexes.erase(Index2);
				}
			}

		size_t ClusterSize = Cluster.size();
		ClusterSizes.push_back(ClusterSize);

		SumClusterSizes += ClusterSize;
		if (ClusterSize >= 10)
			++ClusterCount10;
		if (ClusterSize == 1)
			++ClusterCount1;
		}
	size_t ClusterCount = Clusters.size();

	vector<size_t> Order;
	for (size_t i = 0; i < ClusterCount; ++i)
		Order.push_back(i);
	std::sort(Order.begin(), Order.end(), CmpClusterSize);

	for (auto i : Order)
		Log("Cluster %5u   size %5u\n", i, ClusterSizes[i]);

	asserta(SumClusterSizes == LabelCount);
	ProgressLog("Threshold %.1f: %u clusters, %u >ten, %u singletons\n",
	  Threshold, ClusterCount, ClusterCount10, ClusterCount1);

	uint ErrCount = 0;
	for (size_t Index1 = 0; Index1 < LabelCount; ++Index1)
		{
		ProgressStep(Index1, LabelCount, "Validating");
		const string &Label1 = Labels[Index1];
		size_t C1 = IndexToCluster[Index1];

		for (size_t Index2 = 0; Index2 < LabelCount; ++Index2)
			{
			size_t C2 = IndexToCluster[Index2];
			if (C1 != C2)
				{
				double d12 = DistMx[Index1][Index2];
				if (d12 == DBL_MAX)
					continue;
				asserta(d12 < Threshold);
				}
			}

		bool Found = false;
		const set<size_t> &Cluster = Clusters[C1];
		size_t ClusterSize = Cluster.size();
		for (auto Index2 : Cluster)
			{
			const string &Label2 = Labels[Index2];
			double d12 = DistMx[Index1][Index2];
			if (d12 == DBL_MAX)
				continue;
			if (d12 >= Threshold)
				{
				Found = true;
				break;
				}
			}
		asserta(Found || ClusterSize == 1);
		}

	for (auto ClusterIndex : Order)
//	for (size_t ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		if (g_ftsv != 0)
			{
			size_t ClusterSize = ClusterSizes[ClusterIndex];

			QuartsDouble DistQuarts;
			size_t Centroid = GetCentroid(ClusterIndex, DistQuarts);
			const string &CentroidLabel = Labels[Centroid];

			fprintf(g_ftsv, "C");
			fprintf(g_ftsv, "\t%u", (uint) ClusterIndex);
			fprintf(g_ftsv, "\t%u", (uint) ClusterSize);
			fprintf(g_ftsv, "\t%u", (uint) Centroid);
			fprintf(g_ftsv, "\t%s", CentroidLabel.c_str());

			fprintf(g_ftsv, "\t%.4g", DistQuarts.Min);
			fprintf(g_ftsv, "\t%.4g", DistQuarts.Avg);
			fprintf(g_ftsv, "\t%.4g", DistQuarts.Med);
			fprintf(g_ftsv, "\t%.4g", DistQuarts.Max);

			fprintf(g_ftsv, "\n");
			}
		}

	for (size_t LabelIndex = 0; LabelIndex < LabelCount; ++LabelIndex)
		{
		size_t ClusterIndex = IndexToCluster[LabelIndex];
		asserta(IndexToCluster[LabelIndex] < ClusterCount);
		const string &Label = Labels[LabelIndex];
		if (g_ftsv != 0)
			{
			if (ClusterIndex == UINT_MAX)
				fprintf(g_ftsv, "M	%s	*\n", Label.c_str());
			else
				fprintf(g_ftsv, "M	%s	%u\n", Label.c_str(), (uint) ClusterIndex);
			}
		}

	CloseStdioFile(g_ftsv);
	}
