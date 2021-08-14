#if !defined(AFX_STD_FORMAT__INCLUDED_)
#define AFX_STD_FORMAT__INCLUDED_

#include <stdio.h>
#include <string>

#ifndef __WIN32
#include <stdarg.h>
#endif

namespace std{

std::string format( const char* pszFormat, ... )
{
	class X{ // �֐����N���X
	public:
		char* _ptr;
		X( size_t size ){ _ptr = new char[ size ]; };
		~X(){ delete[] _ptr; };
	};
	 
	const int iInit = 80;//�o�b�t�@�̏����T�C�Y
	const int iSpan = 32;//�o�b�t�@�̊g��X�p��

	for( int iSize = iInit; ; iSize += iSpan ) {
		X tmp( iSize +1 );

		va_list args;
		va_start(args, pszFormat);
#ifdef _WIN32
		int nResult = _vsnprintf( tmp._ptr, iSize, pszFormat, args);
#else
		int nResult = vsprintf( tmp._ptr, pszFormat, args);
#endif
		va_end( args );

		if( nResult == -1 ){ //�o�b�t�@���蓖�Ď��s
			continue;
		}else if( nResult < 0 ){ //�o�̓G���[
			break; // ���͓K���ȗ�O�𑗏o�Ƃ��A�T�[�g����Ƃ��B
		}

		return tmp._ptr; // ����
	}

	return std::string();
}

}

#endif
