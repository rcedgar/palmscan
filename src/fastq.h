#ifndef fastq_h
#define fastq_h

#include "alpha.h"

// Hard-code Illumina variant

static const unsigned MAX_FASTQ_LEN = 16*1024;

enum FASTQ_FILTER
	{
	FF_Good = 0,
	FF_Short = 1,
	FF_HighErr = 2,
	};

class FastQ
	{
public:
	static double m_Wildcard_P_e;
	static byte m_ASCII_Offset;
	static byte m_IntQual_Min;
	static byte m_IntQual_Max;
	static byte m_IntQualOut_Max;
	static double m_CharToProb[256];
	static double m_IntQualToProb[256];
	static byte **m_PairMatch;
	static byte **m_PairMismatch;

public:
	static void InitFromCmdLine();
	static void Init(byte Q_Min, byte Q_Max, byte ASCII_Offset);
	static void SetDefaultVariant();
	static void SetVariant(const string &Name);
	static void InitTables();
	static void InitMerge();

	static inline byte CharToIntQual(byte Ch)
		{
		int IntQual = (int) Ch - (int) m_ASCII_Offset;
		if (IntQual < (int) m_IntQual_Min || IntQual > (int) m_IntQual_Max)
			Die("Phred score %u out of range %u..%u", IntQual, m_IntQual_Min, m_IntQual_Max);
		return (byte) IntQual;
		}

	static inline byte IntQualToCharOut(byte IntQual)
		{
		if (IntQual < m_IntQual_Min || IntQual > m_IntQualOut_Max)
			Die("Phred score %u out of range %u..%u", IntQual, m_IntQual_Min, m_IntQual_Max);
		return IntQual + m_ASCII_Offset;
		}

	static inline byte IntQualToChar(byte IntQual)
		{
		if (IntQual < m_IntQual_Min || IntQual > m_IntQual_Max)
			Die("Phred score %u out of range %u..%u", IntQual, m_IntQual_Min, m_IntQual_Max);
		return IntQual + m_ASCII_Offset;
		}

	static inline double IntQualToProb(byte IntQual)
		{
		return m_IntQualToProb[IntQual];
		}

	static inline double CharToProb(byte Ch)
		{
		return m_CharToProb[Ch];
		}

	static inline double ProbToFloatQual(double P)
		{
		asserta(P >= 0.0 && P <= 1.0);
		if (P == 0.0)
			return m_IntQualOut_Max;
		double FloatQual = -10.0*log10(P);
		if (FloatQual > (double) m_IntQualOut_Max)
			FloatQual = (double) m_IntQualOut_Max;
		return FloatQual;
		}
	};

#endif // fastq_h
