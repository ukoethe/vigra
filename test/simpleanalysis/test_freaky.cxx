#include <iostream>
#include <cmath>

int main()
{
    int x,y;
    
    for(y=0; y<7; ++y)
    {
        for(x=0; x<7; ++x)
        {
            double dist1 = std::sqrt((2.0 - x)*(2.0 - x) +
                                     (2.0 - y)*(2.0 - y));
            double dist2 = std::sqrt((5.0 - x)*(5.0 - x) +
                                     (5.0 - y)*(5.0 - y));
            
            if(dist1 != dist2)
            {
                std::cerr << 1
                          << ((dist1 == dist2) ? "compiler broken\n" : "");
            }
        }
    }
}
