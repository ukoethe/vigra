#ifndef RN_FEATURES_HXX
#define RN_FEATURES_HXX





class FeatureBase
{
public:
    MultiArrayShape<2>::type shape() = 0;
    int shape(int index)
    {
        return shape()[index];
    }
    double & operator() (int i, int j) = 0;
};


class CompositFeatures : public FeatureBase
{
public:
    typedef typename MultiArrayShape<2>::type
        Shp;
    ArrayVector<Shp> 
        ext2int;
    ArrayVector<FeatureBase > sub_feats;
    void add(FeatureBase & feat)
    {
        if(feat.shape(0) != this->shape(0))
            throw std::runtime_error("Wrong Number Of samples");

        sub_feats.push_back(feat);
        for(int ii = 0; ii < feat.shape(1); ++ii)
            ext2int.push_back(Shp(sub_feats.size()-1,
                                  ii));

    }

    MultiArrayShape<2>::type shape()
    {
        return MultiArrayShape<2>::type(sub_feats[0].shape(0),
                                        ext2int.size());
    }

    double & operator() (int i, int j)
    {
        return sub_feats[ext2int[j][0]](i, ext2int[j][1]);
    }
};

template<int N, class T, class C>
class NeighborFeatures : public FeatureBase
{
public:
    typedef typename MultiArrayShape<N>::type Shp;
    MultiArrayView<N, T, C> raw_data;
    ArrayVector<Shp > 
                            feat_coos;
    MultiArrayShape<2>::type shape()
    {
        return MultiArrayShape<2>::type(raw_data.size(),
                                        feat_coos.size());
    }

    double & operator() (int i, int j)
    {
        return raw_data(raw_data.scanOrderIndexToCoordinate(i) + feat_coos(j));
    }
    NeighborFeatures(MultiArrayView<N, T, C> & in, MultiArrayView<2, int> const & coos)
        : raw_data(in)
    {
        for(int ii = 0; ii < coos.shape(0); ++ii)
            feat_coos.push_back(Shp());
            for(int jj = 0; jj < coos.shape(1); ++jj)
                feat_coos.back()[jj] = coos(ii, jj);
    }

};


class BindFeatureColumn : public FeatureBase
{
    typedef typename MultiArrayShape<N>::type Shp;
    int index; 
    FeatureBase & underlying;

    MultiArrayShape<2>::type shape()
    {
        return MultiArrayShape<2>::type(underlying.shape(0), 1);
    }

    double & operator() (int i, int j)
    {
        return  underlying(i, index);
    }
    double & operator[](int i)
    {
        return underlying(i, index); 
    }
    
    BindFeatureColumn(FeaetureBase & in, int index_)
        : index(index_), underlying(in) 
    {
		;
	}
};

FeatureBase columnVector(FeatureBase & in, int ii)
{
    return BindFeatureColumn(in, ii);
}

#endif
