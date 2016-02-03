#include <hdf5.h>

#if defined(CHECK_VERSION)
# if (H5_VERS_MAJOR < MIN_MAJOR) || \
     ((H5_VERS_MAJOR == MIN_MAJOR) && (H5_VERS_MINOR < MIN_MINOR))
#  error "insufficient HDF5 version"
# endif
#endif

int main()
{}
