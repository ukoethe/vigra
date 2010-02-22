#include "../multi_array.hxx"
#include <set>
#include <vector>

namespace vigra
{

template<class T>
struct SampleRange
{
    SampleRange(int start,int end,int num_features)
    {
        this->start=start;
        this->end=end;
        this->min_boundaries.resize(num_features,-FLT_MAX);
        this->max_boundaries.resize(num_features,FLT_MAX);
    }
    int start;
    mutable int end;
    mutable std::vector<T> max_boundaries;
    mutable std::vector<T> min_boundaries;
    bool operator<(const SampleRange& o) const
    {
        return o.start<start;
    }
};

template<class T>
class OnlinePredictionSet
{
public:
    template<class U>
    OnlinePredictionSet(MultiArrayView<2,T,U>& features,int num_sets)
    {
        this->features=features;
        std::vector<int> init(features.shape(0));
        for(int i=0;i<init.size();++i)
            init[i]=i;
        indices.resize(num_sets,init);
        std::set<SampleRange<T> > set_init;
        set_init.insert(SampleRange<T>(0,init.size(),features.shape(1)));
        ranges.resize(num_sets,set_init);
	cumulativePredTime.resize(num_sets,0);
    }
    int get_worsed_tree()
    {
        int result=0;
	for(int i=0;i<cumulativePredTime.size();++i)
	{
	    if(cumulativePredTime[i]>cumulativePredTime[result])
	    {
	        result=i;
	    }
	}
	return result;
    }
    void reset_tree(int index)
    {
        index=index % ranges.size();
        std::set<SampleRange<T> > set_init;
        set_init.insert(SampleRange<T>(0,features.shape(0),features.shape(1)));
        ranges[index]=set_init;
	cumulativePredTime[index]=0;
    }
    std::vector<std::set<SampleRange<T> > > ranges;
    std::vector<std::vector<int> > indices;
    std::vector<int> cumulativePredTime;
    MultiArray<2,T> features;
};

}

