#ifndef CGP2D_MACROS_HXX
#define CGP2D_MACROS_HXX

#include <sstream>

#define CGP_ASSERT_OP(a,comp,b) \
    if(!  static_cast<bool>( a comp b )   ) { \
       std::stringstream s; \
       s << "OpenGM assertion " << #a <<#comp <<#b<< " failed:\n"; \
       s << #a "="<<a<<"\n"; \
       s << #b "="<<b<<"\n"; \
       s << "in file " << __FILE__ << ", line " << __LINE__ << "\n"; \
       throw std::runtime_error(s.str()); \
    }


#define VIGRA_ASSERT_OP(a,comp,b) \
    if(!  static_cast<bool>( a comp b )   ) { \
       std::stringstream s; \
       s << "OpenGM assertion " << #a <<#comp <<#b<< " failed:\n"; \
       s << #a "="<<a<<"\n"; \
       s << #b "="<<b<<"\n"; \
       s << "in file " << __FILE__ << ", line " << __LINE__ << "\n"; \
       throw std::runtime_error(s.str()); \
    }

#endif // CGP2D_MACROS_HXX