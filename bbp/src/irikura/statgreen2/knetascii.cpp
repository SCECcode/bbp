// knetascii.cpp: knetascii クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////

#include "knetascii.h"
#include "format.h"

#include <algorithm>
#include <cmath>
//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

knetascii::knetascii()
{
    origin_time_			= "1996/06/03 00:00";
    epic_latitude_			= 0.0;
    epic_longitude_			= 0.0;
    epic_depth_				= 0.0;
    magnitude_				= 0.0;
    site_code_				= "";
    site_latitude_			= 0.0;
    site_longitude_			= 0.0;
    site_height_			= 0.0;
    record_time_			= "1996/06/03 00:00:00";
    sampling_frequency_		= 0;
    duration_				= 0.0;
    direction_				= "N-S";
    numerator_				= 2000.0;
    denominator_			= 8388608.0;
    max_value_				= 0.0;
    last_correction_		= "1996/06/03 00:00:00";
    memo_					= "";
	kind_of_data_			= KNET_ACC;
}

knetascii::~knetascii()
{

}

string_list knetascii::make_header()
{
	string	strtmp;

	string_list	slist;
	string	line01=string("Origin Time       ");
	string	line02=string("Lat.              ");
	string	line03=string("Long.             ");
	string	line04=string("Depth. (km)       ");
	string	line05=string("Mag.              ");
	string	line06=string("Station Code      ");
	string	line07=string("Station Lat.      ");
	string	line08=string("Station Long.     ");
	string	line09=string("Station Height(m) ");
	string	line10=string("Record Time       ");
	string	line11=string("Sampling Freq(Hz) ");
	string	line12=string("Duration Time(s)  ");
	string	line13=string("Dir.              ");
	string	line14=string("Scale Factor      ");
	string	line15;
	string	line16=string("Last Correction   ");
	string	line17=string("Memo.             ");

	switch(kind_of_data_){
	case KNET_ACC:
		line15=string("Max. Acc. (gal)   ");
		break;
	case KNET_VEL:
		line15=string("Max. Vel. (kine)  ");
		break;
	case KNET_DIS:
		line15=string("Max. Dis. (cm)    ");
		break;
	case KNET_OTHER:
		line15=string("Max. Value        ");
		break;
	default:
		break;
	}

	line01	+= get_origin_time();
	line01.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%5.2f",get_epic_latitude());
	line02	+= strtmp;
	line02.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%5.2f",get_epic_longitude());
	line03	+= strtmp;
	line03.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%3.0f",get_epic_depth());
	line04	+= strtmp;
	line04.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%3.1f",get_magnitude());
	line05	+= strtmp;
	line05.resize(LINE_LENGTH, ' ');

	line06	+= site_code_;
	line06.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%8.4f",get_site_latitude());
	line07	+= strtmp;
	line07.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%8.4f",get_site_longitude());
	line08	+= strtmp;
	line08.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%.0f",get_site_height());
	line09	+= strtmp;
	line09.resize(LINE_LENGTH, ' ');

	line10	+= record_time_;
	line10.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%dHz",get_sampling_frequency());
	line11	+= strtmp;
	line11.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%f",get_duration());
	line12	+= strtmp;
	line12.resize(LINE_LENGTH, ' ');

	line13	+= get_direction();
	line13.resize(LINE_LENGTH, ' ');

	switch(kind_of_data_){
	case KNET_ACC:
		strtmp	= std::format("%f(gal)/%f",get_denominator(), get_numerator());
		break;
	case KNET_VEL:
		strtmp	= std::format("%f(kine)/%f",get_denominator(), get_numerator());
		break;
	case KNET_DIS:
		strtmp	= std::format("%f(cm)/%f",get_denominator(), get_numerator());
		break;
	case KNET_OTHER:
		strtmp	= std::format("%f/%f",get_denominator(), get_numerator());
		break;
	default:
		break;
	}
	line14	+= strtmp;
	line14.resize(LINE_LENGTH, ' ');

	strtmp	= std::format("%f",get_max_value());
	line15	+= strtmp;
	line15.resize(LINE_LENGTH, ' ');

	line16	+= get_last_correction();
	line16.resize(LINE_LENGTH, ' ');

	line17	+= get_memo();
	line17.resize(LINE_LENGTH, ' ');

    slist.push_back(line01);
    slist.push_back(line02);
    slist.push_back(line03);
    slist.push_back(line04);
    slist.push_back(line05);
    slist.push_back(line06);
    slist.push_back(line07);
    slist.push_back(line08);
    slist.push_back(line09);
    slist.push_back(line10);
    slist.push_back(line11);
    slist.push_back(line12);
    slist.push_back(line13);
    slist.push_back(line14);
    slist.push_back(line15);
    slist.push_back(line16);
    slist.push_back(line17);

    return slist;
}

void knetascii::get_data(vector<int>& data)
{
	int	i, n;
	n	= data_.size();
	data.clear();
	data.resize(n);
	for(i=0;i<n;i++) data[i] = data_[i];
}

void knetascii::get_data(vector<double>& data)
{
	int	i, n;
	n	= data_.size();
	data.clear();
	data.resize(n);
	for(i=0;i<n;i++) data[i] = denominator_/numerator_*(double)data_[i];
}

void knetascii::set_data(vector<double> data)
{
	double min = fabs(*min_element(data.begin(), data.end()));
	double max = fabs(*max_element(data.begin(), data.end()));
	double ful = (min > max ? min : max);
	set_denominator(ful);
	int i;
	data_.clear();
	for(i=0;i<data.size();i++){
		int n = (int)(data[i]*numerator_/denominator_);
		data_.push_back(n);
	}
}

ostream& operator<<(ostream& strm, knetascii& knet)
{
	vector<int>	data;
	knet.get_data(data);

    string_list	slist	= knet.make_header();
    string_list::iterator i;
    for(i=slist.begin();i!=slist.end();++i) strm << *i << endl;
	int	j, n=data.size(), m=8-(n%8);
	m = (m == 8 ? 0 : m);
	const int	k = 0;
	for(j=0;j<n;j++){
		strm << std::format("%9d",data[j]);

		if((j+1)%8 == 0) strm << endl;
	}
	for(j=0;j<m;j++){
		strm << std::format("%9d",k);
	}
    return strm;
}

istream& operator>>(istream& strm, knetascii& knet)
{
	int	i, n;
	string	strLine, strItem, strData;
	char	buf[1024], *token;

	knet.clear_data();

	for(i=0;i<HEADER_LINE_NUMBER;i++){
		getline(strm, strLine);
		if(strm.eof()) return strm;
		strItem	= strLine.substr(0, 17);
		strData	= strLine.substr(18);
		strData	= strData.substr(0,strData.find_last_not_of(' ')+1);
		switch(i){
		case 0:
			knet.set_origin_time(strData);
			break;
		case 1:
			knet.set_epic_latitude(atof(strData.c_str()));
			break;
		case 2:
			knet.set_epic_longitude(atof(strData.c_str()));
			break;
		case 3:
			knet.set_epic_depth(atof(strData.c_str()));
			break;
		case 4:
			knet.set_magnitude(atof(strData.c_str()));
			break;
		case 5:
			knet.set_site_code(strData);
			break;
		case 6:
			knet.set_site_latitude(atof(strData.c_str()));
			break;
		case 7:
			knet.set_site_longitude(atof(strData.c_str()));
			break;
		case 8:
			knet.set_site_height(atof(strData.c_str()));
			break;
		case 9:
			knet.set_record_time(strData.c_str());
			break;
		case 10:
			knet.set_sampling_frequency(atoi(strData.c_str()));
			break;
		case 11:
			knet.set_duration(atof(strData.c_str()));
			break;
		case 12:
			knet.set_direction(strData);
			break;
		case 13:
			n	= strData.find('/');
			knet.set_denominator(atof((strData.substr(0,n-1)).c_str()));
			knet.set_numerator(atof((strData.substr(n+1)).c_str()));
			break;
		case 14:
			knet.set_max_value(atof(strData.c_str()));
			break;
		case 15:
			knet.set_last_correction(strData);
			break;
		case 16:
			knet.set_memo(strData);
			break;
		}
	}
	i=0;
	while(!strm.eof()){
		getline(strm, strData);
		strcpy(buf, strData.c_str());

		token	= strtok(buf, " \n");
		while( token != NULL ){
			knet.push_back(atoi(token));
			token	= strtok( NULL, " \n");
		}
	}
    return strm;
}
