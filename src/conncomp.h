#ifndef conncomp_h
#define conncomp_h

#include <map>
#include <set>

void GetConnectedComponents(const vector<vector<unsigned> > &AdjMx_, vector<vector<unsigned> > &CCs,
  bool ShowProgress);

void GetConnectedComponents2(const vector<unsigned> &FromVec, const vector<unsigned> &ToVec, 
  vector<vector<unsigned> > &CCs,  bool ShowProgress);

void GetConnectedComponents3(const map<unsigned, set<unsigned> > &Edges,
  vector<vector<unsigned> > &CCs,  bool ShowProgress);

#endif // conncomp_h
