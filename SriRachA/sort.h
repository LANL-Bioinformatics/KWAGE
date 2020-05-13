#ifndef __SORT
#define __SORT

// While we can get OpenMP working on clang, I have not yet found a
// way to use the Gnu parallel C++ api ...
#if defined(_OPENMP) && !defined(__clang__)

	#include <parallel/algorithm>

	// The default parallel sort algorithm used in __gnu_parallel::sort (potentially
	// called via the SORT macro) is the "multiway merge sort" (MWMS) and
	// requires *twice* the memory of input array to sort in parallel! The other
	// sort algorithm implemented by __gnu_parallel::sort, i.e. quick sort
	// (and balanced quick sort), does not require any extra RAM, but offers a
	// much worse parallel speed up (two fold at best, unless nested parallelism
	// is enabled for OpenMP).
	// * Fast/High memory using MWMS: __gnu_parallel::multiway_mergesort_tag()
	// * Slow/Low memory using balanced QS: __gnu_parallel::balanced_quicksort_tag()

	// Enable OpenMP-based parallel sorting	
	#define	SORT	__gnu_parallel::sort
	//#define	SORT(X, Y)	__gnu_parallel::sort( (X), (Y), __gnu_parallel::balanced_quicksort_tag() )
#else
	#include <algorithm>
	
	// Use standard serial-based sorting
	#define	SORT	std::sort
#endif // _OPENMP

#endif // __SORT
