#include "string_conversion.h"

#include <limits.h>

using namespace std;

size_t str_to_uint32_t(const string &m_str)
{
	size_t ret = 0;
	size_t power = 1;

	if( (m_str.size() > 2) && (m_str[0] == '0') && ( toupper(m_str[1]) == 'X') ){

		const string tmp = m_str.substr(2, m_str.size() - 2);

		// Read as a hexadecimal number
		for(string::const_reverse_iterator i = tmp.rbegin();i != tmp.rend();++i){
	
			switch(*i){
				case '0':
					break;
				case '1':
					ret += power;
					break;
				case '2':
					ret += power*2;
					break;
				case '3':
					ret += power*3;
					break;
				case '4':
					ret += power*4;
					break;
				case '5':
					ret += power*5;
					break;
				case '6':
					ret += power*6;
					break;
				case '7':
					ret += power*7;
					break;
				case '8':
					ret += power*8;
					break;
				case '9':
					ret += power*9;
					break;
				case 'A': case 'a':
					ret += power*10;
					break;
				case 'B': case 'b':
					ret += power*11;
					break;
				case 'C': case 'c':
					ret += power*12;
					break;
				case 'D': case 'd':
					ret += power*13;
					break;
				case 'E': case 'e':
					ret += power*14;
					break;
				case 'F': case 'f':
					ret += power*15;
					break;
				default:
					throw __FILE__ ":str_to_uint32_t: Illegal character in hex form";
			};

			power *= 16;
		}
	}
	else{ // Parse as a decimal number

		for(string::const_reverse_iterator i = m_str.rbegin();i != m_str.rend();++i){
	
			if( (*i < '0') || (*i > '9') ){
				throw __FILE__ ":str_to_uint32_t: Illegal character in decimal form";
			}
		
			ret += power*(*i - '0');
			power *= 10;		
		}
	}

	if(ret > UINT_MAX){
		throw __FILE__ ":str_to_uint32_t: Overflow!";
	}
	
	return uint32_t(ret);
}

size_t str_to_size_t(const string &m_str)
{
	size_t ret = 0;
	size_t power = 1;

	if( (m_str.size() > 2) && (m_str[0] == '0') && ( toupper(m_str[1]) == 'X') ){

		const string tmp = m_str.substr(2, m_str.size() - 2);

		// Read as a hexadecimal number
		for(string::const_reverse_iterator i = tmp.rbegin();i != tmp.rend();++i){
	
			switch(*i){
				case '0':
					break;
				case '1':
					ret += power;
					break;
				case '2':
					ret += power*2;
					break;
				case '3':
					ret += power*3;
					break;
				case '4':
					ret += power*4;
					break;
				case '5':
					ret += power*5;
					break;
				case '6':
					ret += power*6;
					break;
				case '7':
					ret += power*7;
					break;
				case '8':
					ret += power*8;
					break;
				case '9':
					ret += power*9;
					break;
				case 'A': case 'a':
					ret += power*10;
					break;
				case 'B': case 'b':
					ret += power*11;
					break;
				case 'C': case 'c':
					ret += power*12;
					break;
				case 'D': case 'd':
					ret += power*13;
					break;
				case 'E': case 'e':
					ret += power*14;
					break;
				case 'F': case 'f':
					ret += power*15;
					break;
				default:
					throw __FILE__ ":str_to_size_t: Illegal character in hex form";
			};

			power *= 16;
		}
	}
	else{ // Parse as a decimal number

		for(string::const_reverse_iterator i = m_str.rbegin();i != m_str.rend();++i){
		
			if( (*i < '0') || (*i > '9') ){
				throw __FILE__ ":str_to_size_t: Illegal character in decimal form";
			}
		
			ret += power*(*i - '0');
			power *= 10;		
		}
	}
	
	return ret;
}