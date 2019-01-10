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
	class X{ // 関数内クラス
	public:
		char* _ptr;
		X( size_t size ){ _ptr = new char[ size ]; };
		~X(){ delete[] _ptr; };
	};
	 
	const int iInit = 80;//バッファの初期サイズ
	const int iSpan = 32;//バッファの拡大スパン

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

		if( nResult == -1 ){ //バッファ割り当て失敗
			continue;
		}else if( nResult < 0 ){ //出力エラー
			break; // 又は適当な例外を送出とかアサートするとか。
		}

		return tmp._ptr; // 成功
	}

	return '\0';
}

}

#endif
