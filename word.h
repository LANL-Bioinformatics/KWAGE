#ifndef __KMER_TOOLBOX
#define __KMER_TOOLBOX

#include <string>
#include <deque>
#include <vector>

// Storing words in 64-bits (at two bits per base) yeilds a max word
// size of 32.
#define		MAX_WORD_LEN	32
#define		MAX_WORD_VALUE	0xFFFFFFFFFFFFFFFF
#define		BITS_PER_BASE	2
#define		BASE_MASK	3
typedef		size_t	Word;

// A two bit binary base encoding. Order the bases lexographically to ensure that
// the definition of canonical k-mer matches the definition used by the bigsi python 
// implementation. 
enum {BASE_A, BASE_C, BASE_G, BASE_T};

// Fast bit to base conversion (make sure the order of the characters in the string
// matches the order of the enumeration).
inline char bits_to_base(unsigned char m_binary_value)
{
	return "ACGT"[m_binary_value];
}

////////////////////////////////////////////////////////////////////////////////////
// In word.cpp
std::string word_to_string(const Word &m_w, const size_t &m_k);
Word kmer_word_mask(const size_t &m_len);
void find_abundant_kmers(std::deque<Word> &m_out, std::vector< std::deque<Word> > &m_in,
	const size_t &m_threshold);

// Similar to std::unique(), but with the requirement that valid elements must occur at least
// m_threshold-times in the input array.
// - This funtion assumes that the input array has been sorted!
template<typename _ITERATOR>
_ITERATOR thresholded_unique(_ITERATOR m_begin, _ITERATOR m_end, const size_t &m_threshold)
{
	_ITERATOR head = m_begin;
	size_t count = 0;
	
	for(_ITERATOR i = m_begin;i != m_end;++i){
		
		if(*head == *i){
			++count;
		}
		else{ // *head != *i
			
			head += (count >= m_threshold);
			
			*head = *i;
			count = 1;
		}
	}
	
	// Handle the last element
	head += (count >= m_threshold);

	return head;
}

// Macros for digesting sequences into kmers
#define ForEachDuplexWord(__SEQ, __LEN)\
	{\
		const Word __comp_shift = BITS_PER_BASE*( (__LEN) - 1 );\
		const Word __mask = kmer_word_mask(__LEN);\
		Word __w = 0;\
		Word __comp_w = 0;\
		unsigned int __word_len = 0;\
		size_t __k = __LEN;\
		size_t __index = 0;\
		for(std::string::const_iterator __i = (__SEQ).begin();__i != (__SEQ).end();++__i,++__index){ \
			++__word_len;\
			switch(*__i){\
				case 'A': case 'a':\
					__w = (__w << BITS_PER_BASE) | BASE_A;\
					__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_T) << __comp_shift);\
					break;\
				case 'T': case 't':\
					__w = (__w << BITS_PER_BASE) | BASE_T;\
                                	__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_A) << __comp_shift);\
					break;\
				case 'G': case 'g':\
					__w = (__w << BITS_PER_BASE) | BASE_G;\
                                	__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_C) << __comp_shift);\
					break;\
				case 'C': case 'c':\
					__w = (__w << BITS_PER_BASE) | BASE_C;\
                                	__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_G) << __comp_shift);\
					break;\
				default:\
					__word_len = 0;\
					break;\
			};

#define ForEachSenseWord(__SEQ, __LEN)\
	{\
		const Word __mask = kmer_word_mask(__LEN);\
		Word __w = 0;\
		unsigned int __word_len = 0;\
		size_t __k = __LEN;\
		size_t __index = 0;\
		for(std::string::const_iterator __i = (__SEQ).begin();__i != (__SEQ).end();++__i,++__index){ \
			++__word_len;\
			switch(*__i){\
				case 'A': case 'a':\
					__w = (__w << BITS_PER_BASE) | BASE_A;\
					break;\
				case 'T': case 't':\
					__w = (__w << BITS_PER_BASE) | BASE_T;\
					break;\
				case 'G': case 'g':\
					__w = (__w << BITS_PER_BASE) | BASE_G;\
					break;\
				case 'C': case 'c':\
					__w = (__w << BITS_PER_BASE) | BASE_C;\
					break;\
				default:\
					__word_len = 0;\
					break;\
			};

#define ForEachAntisenseWord(__SEQ, __LEN)\
	{\
		const Word __comp_shift = BITS_PER_BASE*( (__LEN) - 1 );\
		const Word __mask = kmer_word_mask(__LEN);\
		Word __comp_w = 0;\
		unsigned int __word_len = 0;\
		size_t __k = __LEN;\
		size_t __index = 0;\
		for(std::string::const_iterator __i = (__SEQ).begin();__i != (__SEQ).end();++__i,++__index){ \
			++__word_len;\
			switch(*__i){\
				case 'A': case 'a':\
					__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_T) << __comp_shift);\
					break;\
				case 'T': case 't':\
                                	__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_A) << __comp_shift);\
					break;\
				case 'G': case 'g':\
                                	__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_C) << __comp_shift);\
					break;\
				case 'C': case 'c':\
                                	__comp_w = (__comp_w >> BITS_PER_BASE) | (Word(BASE_G) << __comp_shift);\
					break;\
				default:\
					__word_len = 0;\
					break;\
			};
				
#define	ErrorWord	(__word_len == 0)
#define	ValidWord	(__word_len >= __k)
#define	SenseWord	(__w & __mask)
#define	AntisenseWord	(__comp_w & __mask)
#define	CanonicalWord	std::min(SenseWord, AntisenseWord)
#define	CurrentBase	(*__i)
#define	Loc3		(__index)		// Location of last base in current kmer
#define	Loc5		( (__index + 1) - __k)	// Location of the first base in the current kmer

#define EndWord\
		}\
	}

#endif // __KMER_TOOLBOX
