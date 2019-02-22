// VerticalSHResponse.h: CVerticalSHResponse �N���X�̃C���^�[�t�F�C�X
//
// ���w�\���ɐ������˂���g�ɑ΂�����g���������v�Z����N���X�ł��B
// �w�̋��E�ɒu���āA�ψʂƉ��̘͂A�����E�������������Ƃɂ���ċ��߂Ă��܂��B
//
// �R���X�g���N�^�Ŕ}���\����^���ăI�u�W�F�N�g���\�z���A
// CalcResponse()�Ń�t�A�T���v������^���Ď��g���������v�Z���܂��B
// ���̌�A�C�ӂ̑w�ԍ����w�肵�āAGetSpectrum()��ʂ���std::vector< std::complex<double> >��
// �������擾���Ă��������B
// ��x�\�z�������g�������͔j������܂ŗL���ł��B
//
// �}���\���͔g�̈ʑ����x�A���x�A�w���ŋK�肳��܂��B�����l�̒P�ʌn�͑����Ă��������B
//
// �Q�l�����F
// �����וv�A1978�A�e���g���_�A��g���X�A118--120
//
// Ver. 1.02  2002 S.Senna(MSS)
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VERTICALSHRESPONSE_H__79A9A070_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_)
#define AFX_VERTICALSHRESPONSE_H__79A9A070_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_

#if _MSC_VER > 1000
//#pragma once
#endif // _MSC_VER > 1000

//#pragma warning(disable:4786)

#include <vector>
#include <complex>


class CVerticalSHResponse
{
public:
	CVerticalSHResponse( const std::vector<double>& velocity, const std::vector<double>& density, const std::vector<double>& thickness);
	virtual ~CVerticalSHResponse();

	void CalcResponse( unsigned int numofSampling, double delta_time );
	void GetSpectrum( unsigned int layerNo, std::vector< std::complex<double> >& spectrum, double delta_time );
	void GetWave( unsigned int layerNo, const std::vector<double>& inWave, std::vector<double>& outWave );
    
	static inline std::complex<double> exp( std::complex<double> x ){ static const double base_e = 2.718281828459045; return pow( base_e, x ); }
	//static inline std::complex<double> exp( const std::complex<double> *x ){ static const double base_e = 2.718281828459045; return std::pow( base_e, x ); }

	double DeltaT(){ return m_dt; }
	unsigned int NumofSample(){ return m_a[0].size(); }

private:
	const std::vector<double>& m_velocity;
	const std::vector<double>& m_density;
	const std::vector<double>& m_thickness;

	std::vector< std::vector< std::complex<double> > > m_a; //����������t�[���G�W���@m_a[�w�ԍ�][����or���g��]
	std::vector< std::vector< std::complex<double> > > m_b; //�����������t�[���G�W���@m_b[�w�ԍ�][����or���g��]

	double m_dt;
};

#endif // !defined(AFX_VERTICALSHRESPONSE_H__79A9A070_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_)
