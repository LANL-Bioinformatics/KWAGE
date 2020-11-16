#include "sra_accession.h"
#include <algorithm>

// DEBUG
#include <iostream>
using namespace std;

#define	LETTER		26	// Letters use 5 bits
#define	DIGIT		10	// Digits use 4 bits

#define	NUM_LETTERS	3

// Accession is a 64 bit unsigned integer, i.e.
// [63,62,61, ..., 1, 0]
// 4 bits [0-3] -> number of digits
// 60 bits [4 - 63] -> payload (binary accession)

// The payload is partitioned as: 
// [4 bit num digits][ 3 letters * 5 bits/letter][ (num digits)* 4 bits/digit]
//
// 4 bits for the number of digits
// 15 bits for the three letters
// 41 bits/ (4 bits/digit) = max of 10 digits
#define		MAX_NUM_DIGITS		10

// Pack an accession string 
SraAccession str_to_accession(const std::string &m_accession)
{	
	size_t num_letter = 0;
	size_t num_digit = 0;
	size_t data = 0;
	
	for(std::string::const_iterator i = m_accession.begin();i != m_accession.end();++i){
		
        switch( toupper(*i) ){
            case 'A': case 'B': case 'C': case 'D': case 'E': case 'F': case 'G': case 'H':
            case 'I': case 'J': case 'K': case 'L': case 'M': case 'N': case 'O': case 'P':
            case 'Q': case 'R': case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
            case 'Y': case 'Z':
                
                ++num_letter;
                data = data*LETTER + ( toupper(*i) - 'A' );
                break;
            case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
            case '8': case '9':
                
                ++num_digit;
                data = data*DIGIT + ( *i - '0' );
                break;
        };
	}
    
    if( (num_letter != NUM_LETTERS) || (num_digit == 0) || (num_digit > MAX_NUM_DIGITS) ){
        throw __FILE__ ":str_to_accession: Unable to parse accession string";
    }
    
	const SraAccession ret = (num_digit - 1) | (data << 4);

	if(ret == INVALID_ACCESSION){
		throw __FILE__ ":str_to_accession: Mapped input string to INVALID_ACCESSION";
	}

	return ret;
}

// Unpack an accession
std::string accession_to_str(const SraAccession &m_accession)
{
	std::string ret;

	const size_t num_letter = 3;
	const size_t num_digit = (m_accession & 0xF) + 1;
	size_t data = (m_accession >> 4) & 0x0FFFFFFFFFFFFFFF;
	
	for(size_t i = 0;i < num_digit;++i){
		
		const char c = (data % DIGIT) + '0';
		
		data /= DIGIT;
		
		ret.push_back(c);
	}
	
	for(size_t i = 0;i < num_letter;++i){
		
		const char c = (data % LETTER) + 'A';
		
		data /= LETTER;
		
		ret.push_back(c);
	}
	
	std::reverse( ret.begin(), ret.end() );
	
	return ret;
}
