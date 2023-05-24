#include "myutils.h"
#include "rdrpmodel.h"
#include "mpcluster.h"

void SetBLOSUM62();

void cmd_cluster_ppm()
	{
	const string &ModelFileName = opt_cluster_ppm;

	SetBLOSUM62();

	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	const uint GroupCount = Model.GetGroupCount();
	vector<MotifProfile *> MPs;
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		MotifProfile *MP = new MotifProfile;
		Model.GetMotifProfile(GroupIndex, *MP);
		Model.GetGroupName(GroupIndex, MP->m_Name);
		MP->ValidateFreqs();
		MPs.push_back(MP);
		}

	Log("%u groups\n", GroupCount);
	MPCluster MPC;
	MPC.NNCluster(MPs, -999);
	MPC.NNToNewickFile(opt_tree);
	}
