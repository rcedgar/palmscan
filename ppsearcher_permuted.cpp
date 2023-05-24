#include "myutils.h"
#include "cmpsearcher.h"

#define TRACE 0

//// Manually extended to allow outliers
//static const uint min_aadist_AdBg = 10;
//static const uint max_aadist_AdBg = 150;
//
//static const uint min_aadist_BgCd = 10;
//static const uint max_aadist_BgCd = 150;

// Manual guesses, lack of training data
static const uint min_aadist_CdAd_permuted = 10;
static const uint max_aadist_CdAd_permuted = 150;

static const uint min_aadist_AdBg_permuted = 10;
static const uint max_aadist_AdBg_permuted = 150;

void CMPSearcher::GetAdLoHi_Permuted(uint Cd, uint &AdLo, uint &AdHi) const
	{
	AdLo = Cd + min_aadist_CdAd_permuted;
	AdHi = Cd + max_aadist_CdAd_permuted;
	}

void CMPSearcher::GetBgLoHi_Permuted(uint Ad, uint &BgLo, uint &BgHi) const
	{
	BgLo = Ad + min_aadist_AdBg_permuted;
	BgHi = Ad + max_aadist_AdBg_permuted;
	}

void CMPSearcher::Search_CAB(const PDBChain &Query)
	{
	const uint QL = SIZE(Query.m_Seq);

	vector<uint> AdVec;
	vector<uint> BgVec;
	vector<uint> CdVec;

	SearchCd(0, QL, CdVec);
	const uint ncd = SIZE(CdVec);
	for (uint icd = 0; icd < ncd; ++icd)
		{
		uint Cd = CdVec[icd];
		uint AdLo, AdHi;
		GetAdLoHi_Permuted(Cd, AdLo, AdHi);

		SearchAd(AdLo, AdHi, AdVec);
		const uint nad = SIZE(AdVec);
		for (uint iad = 0; iad < nad; ++iad)
			{
			uint Ad = AdVec[iad];
			uint BgLo, BgHi;
			GetBgLoHi_Permuted(Ad, BgLo, BgHi);

			SearchBg(BgLo, BgHi, BgVec);
			const uint nbg = SIZE(BgVec);
			for (uint ibg = 0; ibg < nbg; ++ibg)
				{
				uint Bg = BgVec[ibg];
				double Score = GetScore(Ad, Bg, Cd);
#if TRACE
				Log("[Score3] %.4f\n", Score);
#endif
				CheckHit(Ad, Bg, Cd, Score);
				}
			}
		}
	}
