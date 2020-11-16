#ifndef __SRA_DATE
#define __SRA_DATE

#include <string>
#include <iostream>

// Needed to keep g++ happy (clang does not need these includes)
#include "binary_io.h"
#include "mpi_util.h"

class Date
{
	private:

		// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
		// ensure that structure variable are correctly serialized
		#define DATE_MEMBERS \
			VARIABLE(unsigned int, day) \
			VARIABLE(unsigned int, month) \
			VARIABLE(unsigned int, year)
	
		#define VARIABLE(A, B) A B;
			DATE_MEMBERS
		#undef VARIABLE

	public:
		
		Date()
		{
			day = month = year = 0;
		};

		Date(const std::string &m_str)
		{
			// Dates are formatted as YYYY-MM-DDThh:mm:ssZ
			// For example: "2010-03-24T03:10:22Z"

			// Make sure that we have YYYY-MM-DD
			if(m_str.size() < 10){
				std::cerr << "Malformed date: " <<  m_str << std::endl;
				throw __FILE__ ":Date: Malformed date string";
			}

			day = month = year = 0;

			///////////////////////////////////////////////
			// Parse the four digit year
			///////////////////////////////////////////////
			if( (m_str[0] < '0') || (m_str[0] > '9') ){
				throw __FILE__ ":Date: Malformed year[0]";
			}

			year += 1000*(m_str[0] - '0');

			if( (m_str[1] < '0') || (m_str[1] > '9') ){
				throw __FILE__ ":Date: Malformed year[1]";
			}

			year += 100*(m_str[1] - '0');

			if( (m_str[2] < '0') || (m_str[2] > '9') ){
				throw __FILE__ ":Date: Malformed year[2]";
			}

			year += 10*(m_str[2] - '0');

			if( (m_str[3] < '0') || (m_str[3] > '9') ){
				throw __FILE__ ":Date: Malformed year[3]";
			}

			year += (m_str[3] - '0');

			if(m_str[4] != '-'){
				throw __FILE__ ":Date: Malformed year-month separator";
			}

			///////////////////////////////////////////////
			// Parse the two digit month
			///////////////////////////////////////////////

			if( (m_str[5] < '0') || (m_str[5] > '9') ){
				throw __FILE__ ":Date: Malformed month[5]";
			}

			month += 10*(m_str[5] - '0');

			if( (m_str[6] < '0') || (m_str[6] > '9') ){
				throw __FILE__ ":Date: Malformed month[6]";
			}

			month += (m_str[6] - '0');

			if(m_str[7] != '-'){
				throw __FILE__ ":Date: Malformed month-day separator";
			}

			///////////////////////////////////////////////
			// Parse the two digit day
			///////////////////////////////////////////////

			if( (m_str[8] < '0') || (m_str[8] > '9') ){
				throw __FILE__ ":Date: Malformed day[8]";
			}

			day += 10*(m_str[8] - '0');

			if( (m_str[9] < '0') || (m_str[9] > '9') ){
				throw __FILE__ ":Date: Malformed month[9]";
			}

			day += (m_str[9] - '0');
		};

		inline bool operator==(const Date &m_rhs) const
		{
			return (day == m_rhs.day) &&
				(month == m_rhs.month) &&
				(year == m_rhs.year);
		};

		inline bool operator<(const Date &m_rhs) const
		{
			if(year == m_rhs.year){
				
				if(month == m_rhs.month){
					return (day < m_rhs.day);
				}

				return (month < m_rhs.month);
			}
			
			return (year < m_rhs.year);
		};

		inline bool operator<=(const Date &m_rhs) const
		{
			if(year == m_rhs.year){
				
				if(month == m_rhs.month){
					return (day <= m_rhs.day);
				}

				return (month < m_rhs.month);
			}
			
			return (year < m_rhs.year);
		};

		inline bool operator>(const Date &m_rhs) const
		{
			if(year == m_rhs.year){
				
				if(month == m_rhs.month){
					return (day > m_rhs.day);
				}

				return (month > m_rhs.month);
			}
			
			return (year > m_rhs.year);
		};

		inline bool operator>=(const Date &m_rhs) const
		{
			if(year == m_rhs.year){
				
				if(month == m_rhs.month){
					return (day >= m_rhs.day);
				}

				return (month > m_rhs.month);
			}
			
			return (year > m_rhs.year);
		};

		inline bool is_valid() const
		{
			return (year != 0) && (month != 0) && (day != 0);
		};

		inline unsigned int get_year() const
		{
			return year;
		};

		inline unsigned int get_month() const
		{
			return month;
		};

		inline unsigned int get_day() const
		{
			return day;
		};

		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
		
		template<class T> friend void binary_write(std::ostream &m_out,
			const T &m_obj);
		template<class T> friend void binary_read(std::istream &m_in, 
			T &m_obj);
};

template<> size_t mpi_size(const Date &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Date &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Date &m_obj);

template<> void binary_write(std::ostream &m_out, const Date &m_obj);
template<> void binary_read(std::istream &m_in, Date &m_obj);

std::ostream& operator <<(std::ostream &m_stream, const Date &m_date);

#endif // __SRA_DATE
