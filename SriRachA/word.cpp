#include "word.h"
#include "sort.h"

#include <algorithm>
#include <iostream>

using namespace std;

Word kmer_word_mask(const size_t &m_len)
{
	Word ret = 0;
	
	// The following code is commented out because it does not
	// work when the m_kmer_len is equal to full word length
	// (we overflow the Word).
        //const Word ret = ( 1UL << (2*m_len) ) - 1;

	for(size_t i = 0;i < BITS_PER_BASE*m_len;++i){
		ret |= (1UL << i);
	}

	return ret;
}

string word_to_string(const Word &m_w, const size_t &m_k)
{
	Word w = m_w;
	deque<char> buffer;
	
	for(size_t i = 0;i < m_k;++i){
		
		switch(w & BASE_MASK){
			case BASE_A:
				buffer.push_back('A');
				break;
			case BASE_T:
				buffer.push_back('T');
				break;
			case BASE_G:
				buffer.push_back('G');
				break;
			case BASE_C:
				buffer.push_back('C');
				break;
			default:
				throw __FILE__ ":word_to_string: Unknown base!";
		};
		
		w = (w >> BITS_PER_BASE);
	}
	
	// Since we've added bases in 3' -> 5' order, we need to reverse the sequence
	reverse( buffer.begin(), buffer.end() );
	
	return string( buffer.begin(), buffer.end() );
}

void find_abundant_kmers(deque<Word> &m_out, vector< deque<Word> > &m_in, const size_t &m_threshold)
{
	// Sort all inputs in ascending order
	for(vector< deque<Word> >::iterator i = m_in.begin();i != m_in.end();++i){
		SORT( i->begin(), i->end() );
	}
	
	typedef deque<Word>::iterator I;
	const unsigned int len = m_in.size();
	
	vector<I> iter(len);
	vector<I> end(len);
	
	for(unsigned int i = 0;i < len;++i){
		
		iter[i] = m_in[i].begin();
		end[i] = m_in[i].end();
	}
	
	bool valid = true;
	
	while(valid){
		
		Word curr = MAX_WORD_VALUE;
		valid = false;
		
		for(unsigned int i = 0;i < len;++i){
		
			if(iter[i] != end[i]){
				
				curr = min(curr, *iter[i]);
				valid = true;
			}
		}
		
		size_t count = 0;
		
		for(unsigned int i = 0;i < len;++i){
			
			while( (iter[i] != end[i]) && (*iter[i] == curr) ){
				
				++count;
				++iter[i];
			}
		}
		
		if(count >= m_threshold){
			m_out.push_back(curr);
		}		
	}
		
	// Since the inputs were sorted, the output will also be sorted
}
