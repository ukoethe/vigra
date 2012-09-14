#pragma once
#ifndef VIGRA_TWS_HXX
#define VIGRA_TWS_HXX

#include <stddef.h>

#include <vector>
#include <queue>
#include <algorithm> // max
#include <map>

#include <vigra/windows.h>
#include <vigra/multi_array.hxx>
#include <vigra/bucket_queue.hxx>

namespace vigra {

template<class T>
void twsc
(
	const MultiArrayView<3, UInt8, StridedArrayTag>&,
	MultiArrayView<3, T, StridedArrayTag>&,
	MultiArrayView<3, UInt8, StridedArrayTag>&
);

template<class T>
inline bool isAtSeedBorder
(
	const MultiArrayView<3, T, StridedArrayTag>& labeling,
	const MultiArrayIndex& index
);

template<class S1, class T, class S2>
void tws
(
	const MultiArrayView<3, UInt8, S1>& vol,
	MultiArrayView<3, T, S2>& labeling
)
{
    size_t numVoxels = vol.size();

	// define 256 queues, one for each gray level.
	std::vector<std::queue<MultiArrayIndex> > queues(256);

    std::cout << "uint8 version\n" << std::flush;

	// add each unlabeled pixels which is adjacent to a seed
	// to the queue corresponding to its gray level
	for(MultiArrayIndex j = 0; j < labeling.size(); ++j) {
        if(j % 1000 == 0) {
            std::cout << "\r  initializing queues " << (j/float(numVoxels)*100) << "%                    " << std::flush;
        }
		if(isAtSeedBorder<T>(labeling, j)) {
			queues[vol[j]].push(j);
		}
	}
    std::cout << std::endl;

	// flood
	UInt8 grayLevel = 0;

    size_t voxelsProcessed = 0;

	for(;;) {
		while(!queues[grayLevel].empty()) {
			// label pixel and remove from queue
			MultiArrayIndex j = queues[grayLevel].front();
			queues[grayLevel].pop();

            ++voxelsProcessed;

            if(voxelsProcessed % 1000 == 0) {
                std::cout << "\r  watersheds " << (voxelsProcessed/float(numVoxels)*100) << "%                    " << std::flush;
            }

			// add unlabeled neighbors to queues
			// left, upper, and front voxel
			typename MultiArrayView<3, UInt32>::difference_type coordinate = labeling.scanOrderIndexToCoordinate(j);
			for(unsigned short d = 0; d<3; ++d) {
				if(coordinate[d] != 0) {
					--coordinate[d];
					if(labeling[coordinate] == 0) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						UInt8 queueIndex = std::max(vol[coordinate], grayLevel);
						labeling[k] = labeling[j]; // label pixel
						queues[queueIndex].push(k);
					}
					++coordinate[d];
				}
			}
			// right, lower, and back voxel
			for(unsigned short d = 0; d<3; ++d) {
				if(coordinate[d] < labeling.shape(d)-1) {
					++coordinate[d];
					if(labeling[coordinate] == 0) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						UInt8 queueIndex = std::max(vol[coordinate], grayLevel);
						labeling[k] = labeling[j]; // label pixel
						queues[queueIndex].push(k);
					}
					--coordinate[d];
				}
			}
		}
		if(grayLevel == 255) {
			break;
		}
		else {
			queues[grayLevel] = std::queue<MultiArrayIndex>(); // free memory
			++grayLevel;
		}
	}
    std::cout << std::endl;
}

template<class S1, class T, class S2>
void tws
(
	const MultiArrayView<3, float, S1>& vol,
	MultiArrayView<3, T, S2>& labeling
)
{
    size_t numVoxels = vol.size();

	// define 256 queues, one for each gray level.
    typedef PriorityQueue<MultiArrayIndex, float, true> PQ;
	std::vector<PQ> queues(256);

    std::cout << "float version\n" << std::flush;

	// add each unlabeled pixels which is adjacent to a seed
	// to the queue corresponding to its gray level
	for(MultiArrayIndex j = 0; j < labeling.size(); ++j) {
        if(j % 1000 == 0) {
            std::cout << "\r  initializing queues " << (j/float(numVoxels)*100) << "%                    " << std::flush;
        }
		if(isAtSeedBorder<T>(labeling, j)) 
        {
			queues[(int)vol[j]].push(j, vol[j]);
		}
	}
    std::cout << std::endl;

	// flood
	UInt8 grayLevel = 0;

    size_t voxelsProcessed = 0;

	for(;;) {
		while(!queues[grayLevel].empty()) {
			// label pixel and remove from queue
			MultiArrayIndex j = queues[grayLevel].top();
            float p = queues[grayLevel].topPriority();
			queues[grayLevel].pop();

            ++voxelsProcessed;

            if(voxelsProcessed % 1000 == 0) {
                std::cout << "\r  watersheds " << (voxelsProcessed/float(numVoxels)*100) << "%                    " << std::flush;
            }

			// add unlabeled neighbors to queues
			// left, upper, and front voxel
			typename MultiArrayView<3, UInt32>::difference_type coordinate = labeling.scanOrderIndexToCoordinate(j);
			for(unsigned short d = 0; d<3; ++d) {
				if(coordinate[d] != 0) {
					--coordinate[d];
					if(labeling[coordinate] == 0) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						float priority = std::max(vol[coordinate], p);
						labeling[k] = labeling[j]; // label pixel
						queues[(int)priority].push(k, priority);
					}
					++coordinate[d];
				}
			}
			// right, lower, and back voxel
			for(unsigned short d = 0; d<3; ++d) {
				if(coordinate[d] < labeling.shape(d)-1) {
					++coordinate[d];
					if(labeling[coordinate] == 0) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						float priority = std::max(vol[coordinate], p);
						labeling[k] = labeling[j]; // label pixel
						queues[(int)priority].push(k, priority);
					}
					--coordinate[d];
				}
			}
		}
		if(grayLevel == 255) {
			break;
		}
		else {
			queues[grayLevel] = PQ(); // free memory
			++grayLevel;
		}
	}
    std::cout << std::endl;
}

template <class T>
struct TWS
{
    template<class S1, class U, class S2>
    static void exec(const MultiArrayView<3, T, S1>&,
                     MultiArrayView<3, U, S2>&)
    {
        vigra_precondition(false,
           "tws(): Turbo watershed needs input type UInt8 or float32.");
    }
};

template <>
struct TWS<UInt8>
{
    template<class S1, class U, class S2>
    static void exec(const MultiArrayView<3, UInt8, S1>& vol,
                     MultiArrayView<3, U, S2>& labeling)
    {
        tws(vol, labeling);
    }
};

template <>
struct TWS<float>
{
    template<class S1, class U, class S2>
    static void exec(const MultiArrayView<3, float, S1>& vol,
                     MultiArrayView<3, U, S2>& labeling)
    {
        tws(vol, labeling);
    }
};

template<class T>
void twsc
(
	const MultiArrayView<3, UInt8>& vol,
	MultiArrayView<3, T>& labeling, 
	MultiArrayView<3, UInt8>& directions,
	std::map<std::pair<T, T>, std::pair<MultiArrayIndex, MultiArrayIndex> >& adjacency
)
{
	// define 256 queues, one for each gray level.
	std::vector<std::queue<MultiArrayIndex> > queues(256);

	// add each unlabeled pixels which is adjacent to a seed
	// to the queue corresponding to its gray level
	for(MultiArrayIndex j = 0; j < labeling.size(); ++j) {
		if(isAtSeedBorder<T>(labeling, j)) {
			queues[vol[j]].push(j);
		}
	}

	// flood
	UInt8 grayLevel = 0;
	for(;;) {
		while(!queues[grayLevel].empty()) {
			// label pixel and remove from queue
			MultiArrayIndex j = queues[grayLevel].front();
			queues[grayLevel].pop();

			// add unlabeled neighbors to queues
			// left, upper, and front voxel
			typename MultiArrayView<3, UInt32>::difference_type coordinate = labeling.scanOrderIndexToCoordinate(j);
			for(unsigned short d = 0; d<3; ++d) {
				if(coordinate[d] != 0) {
					--coordinate[d];
					if(labeling[coordinate] == 0) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						UInt8 queueIndex = std::max(vol[coordinate], grayLevel);
						labeling[k] = labeling[j]; // label pixel
						directions[k] = d+1; // save direction
						queues[queueIndex].push(k);
					}
					else if(labeling[coordinate] != labeling[j]) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						if(labeling[j] < labeling[k]) {
							std::pair<T, T> p(labeling[j], labeling[k]);
							if(adjacency.count(p) == 0) {
								adjacency[p] = std::pair<MultiArrayIndex, MultiArrayIndex>(j, k);
							}
						}
						else {
							std::pair<T, T> p(labeling[k], labeling[j]);
							if(adjacency.count(p) == 0) {
								adjacency[p] = std::pair<MultiArrayIndex, MultiArrayIndex>(k, j);
							}
						}
					}
					++coordinate[d];
				}
			}
			// right, lower, and back voxel
			for(unsigned short d = 0; d<3; ++d) {
				if(coordinate[d] < labeling.shape(d)-1) {
					++coordinate[d];
					if(labeling[coordinate] == 0) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						UInt8 queueIndex = std::max(vol[coordinate], grayLevel);
						labeling[k] = labeling[j]; // label pixel
						directions[k] = d+4; // save direction
						queues[queueIndex].push(k);
					}
					else if(labeling[coordinate] != labeling[j]) {
						MultiArrayIndex k = labeling.coordinateToScanOrderIndex(coordinate);
						if(labeling[j] < labeling[k]) {
							std::pair<T, T> p(labeling[j], labeling[k]);
							if(adjacency.count(p) == 0) {
								adjacency[p] = std::pair<MultiArrayIndex, MultiArrayIndex>(j, k);
							}
						}
						else {
							std::pair<T, T> p(labeling[k], labeling[j]);
							if(adjacency.count(p) == 0) {
								adjacency[p] = std::pair<MultiArrayIndex, MultiArrayIndex>(k, j);
							}
						}
					}
					--coordinate[d];
				}
			}
		}
		if(grayLevel == 255) {
			break;
		}
		else {
			queues[grayLevel] = std::queue<MultiArrayIndex>(); // free memory
			++grayLevel;
		}
	}
}

template<class T>
inline bool isAtSeedBorder
(
	const MultiArrayView<3, T, StridedArrayTag>& labeling,
	const MultiArrayIndex& index
)
{
	if(labeling[index] == 0) {	
		return false; // not a seed voxel
	}
	else {
		typename MultiArrayView<3, UInt32>::difference_type coordinate
			= labeling.scanOrderIndexToCoordinate(index);
		// check left, upper, and front voxel for zero label
		for(unsigned short d = 0; d<3; ++d) {
			if(coordinate[d] != 0) {
				--coordinate[d];
				if(labeling[coordinate] == 0) {
					return true;
				}
				++coordinate[d];
			}
		}
		// check right, lower, and back voxel for zero label
		for(unsigned short d = 0; d<3; ++d) {
			if(coordinate[d] < labeling.shape(d)-1) {
				++coordinate[d];
				if(labeling[coordinate] == 0) {
					return true;
				}
				--coordinate[d];
			}
		}
		return false;
	}
}

} // namespace vigra

#endif // VIGRA_TWS_HXX
