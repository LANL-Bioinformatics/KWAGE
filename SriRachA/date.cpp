#include "date.h"

using namespace std;

ostream& operator <<(ostream &m_stream, const Date &m_date)
{
	m_stream << m_date.get_year() << '-' << m_date.get_month() << '-' << m_date.get_day();

	return m_stream;
}