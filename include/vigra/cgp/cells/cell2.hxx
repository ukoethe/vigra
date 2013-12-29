#ifndef CGP2D_CELL_2_HXX
#define CGP2D_CELL_2_HXX

/* std library */
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>

/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>

/* vigra cgp */
#include <vigra/cgp/cells/cell_base.hxx>
#include <vigra/cgp/tgrid.hxx>


namespace vigra{


template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cell<COORDINATE_TYPE,LABEL_TYPE,2> : public CellBase<COORDINATE_TYPE,LABEL_TYPE,2>{
public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
private:
};


}

#endif // CGP2D_CELL_2_HXX