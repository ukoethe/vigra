template <int N>
struct Sum
{
    static const int sum = 1 + Sum<N-1>::sum;
};

template <>
struct Sum<0>
{
    static const int sum = 0;
};

int main()
{
    static const int sum = Sum<DEPTH>::sum;
}
