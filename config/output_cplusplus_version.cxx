#include <iostream>

int main() {
    if (__cplusplus == 1)
    {
        // This is a workaround for a gcc bug:
        // http://stackoverflow.com/questions/7530047/gnu-c-macro-cplusplus-standard-conform
#ifdef __GXX_EXPERIMENTAL_CXX0X__
        std::cout << 201103;
#else
        std::cout << 199711;
#endif
    }
    else
    {
        std::cout << __cplusplus;
    }
    return 0;
}
