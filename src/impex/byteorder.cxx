#include <algorithm>
#include <stdexcept>
#include "byteorder.hxx"

namespace vigra
{
    const byteorder::host byteorder::m_host;

    byteorder::host::host()
    {
        // byteorder check: if the first byte is the least significant,
        // we have little endian byteorder.
        unsigned int testint = 0x01;
        unsigned char * testchar
            = reinterpret_cast< unsigned char * >(&testint);
        if ( testchar[0] == 0x01 )
            m_string = "little endian";
        else
            m_string = "big endian";
    }

    const std::string & byteorder::host::get() const
    {
        return m_string;
    }

    const std::string & byteorder::get_host_byteorder() const
    {
        return m_host.get();
    }

    byteorder::byteorder()
    {
        set( m_host.get() );
    }

    byteorder::byteorder( const std::string & s )
    {
        set(s);
    }

    const std::string & byteorder::get() const
    {
        return m_string;
    }

    void byteorder::set( const std::string & s )
    {
        m_string = s;
        native = get_host_byteorder() == s;
    }
}
