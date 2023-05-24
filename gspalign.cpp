#include "myutils.h"
#include "gspaligner.h"

void cmd_gspalign()
	{
	const string &QueryFileName = opt_gspalign;
	const string &TargetFileName = opt_db;

	GSProf Q;
	GSProf T;

	Q.FromTsv(QueryFileName);
	T.FromTsv(TargetFileName);

	GSPAligner GA;
	GA.Align(Q, T);
	}
