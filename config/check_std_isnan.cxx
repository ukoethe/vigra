#include <cmath>
int main()
{
    // Some old versions of libc++ do not provide isnan() for integers,
    //  in which case this will fail to compile.
    std::isnan(int(0));
    return 0;
}
