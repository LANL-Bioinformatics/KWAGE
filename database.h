#ifndef __BLOOM_DATABASE
#define __BLOOM_DATABASE

#include <deque>
#include <fstream>

#include "options.h"
#include "bloom.h"
#include "bigsi++.h"
#include "update.h"

void merge_bloom_filters(const std::deque<BloomFilter> &m_db, 
	const BloomParam &m_param,
	UpdateInfo &m_progress,
	const BuildOptions &m_opt);

#endif // __BLOOM_DATABASE
