#include <H5pubconf.h>

int main()
{
#if defined(H5_SOMETHING) && H5_SOMETHING == 1
    return 0;
#else   
    return 1;
#endif
}
