#ifndef __ACCESSION
#define __ACCESSION

#include <string>

// A Compact form for storing NCBI SRA accession values of the form:
// [3] capital letters + [1-12] digits

typedef size_t SraAccession;
#define	INVALID_ACCESSION	0x0

SraAccession str_to_accession(const std::string &m_accession);
std::string accession_to_str(const SraAccession &m_accession);

#endif // __ACCESSION
