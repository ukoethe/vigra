#ifndef VIGRA_IMPEX_ERROR_HXX
#define VIGRA_IMPEX_ERROR_HXX

#include "vigra/error.hxx"

#define VIGRA_IMPEX_FINALIZED(p) { if (p) \
    vigra_precondition( false, "encoder settings were already finalized" ); }

#endif // VIGRA_IMPEX_ERROR_HXX
