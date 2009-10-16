#include <hdf5.h>

int main()
{
    if(H5_VERS_MAJOR > MIN_MAJOR)
        return 0;
   
    if(H5_VERS_MAJOR == MIN_MAJOR && H5_VERS_MINOR >= MIN_MINOR)
        return 0;
      
    return 1;
}
