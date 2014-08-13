
// Author(s): Thomas Walter
// $Date$
// $Rev$
// $URL$

#ifndef MORPHO_WATERSHED_HXX_
#define MORPHO_WATERSHED_HXX_


#include "morpho_utilities.hxx"
#include "morpho_geodesy.hxx"

namespace vigra{
namespace morpho{

template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE>
int morphoBifurcationPoints(Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
                            Iterator2 destUpperLeft, Accessor2 desta,
                            NBTYPE & nbOffset)
{
    int BIF_CONDITION = 2;

    int width  = srcLowerRight.x - srcUpperLeft.x;
    int height = srcLowerRight.y - srcUpperLeft.y;

    typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
    typedef typename NBTYPE::SIZETYPE SIZETYPE;

    Diff2D o0(0,0);

    for(o0.y = 0; o0.y < height; ++o0.y)
    {
      for(o0.x = 0; o0.x < width; ++o0.x)
      {
        if(srca(srcUpperLeft, o0) != 0)
        {
          int nb_neighbors=0;

          // count the non-zero neighbors
          for(ITERATORTYPE iter = nbOffset.begin();
              iter != nbOffset.end();
              ++iter)
          {
            Diff2D o1 = o0 + *iter;

            // if the neighbor is not outside the image
            if(!nbOffset.isOutsidePixel(o1))
            {
              if(srca(srcUpperLeft,o1) != 0)
              {
                nb_neighbors++;
              }
            }
          }

          // if the number of neighbors exceeds 2, the point is a bifurcation point
          if(nb_neighbors>BIF_CONDITION) {
            desta.set(srca(srcUpperLeft, o0), destUpperLeft, o0);
          }
        } // end if src has non-zero pixel
      } // end x-loop
    } // end y-loop
  return(0);
} // end of function

//////////////////////////
// morphoBifurcationPoints
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE>
inline
int morphoBifurcationPoints(vigra::triple<Iterator1, Iterator1, Accessor1> src,
                            vigra::pair<Iterator2, Accessor2> dest,
                            NBTYPE & neighborOffset)
{
  //clock_t startTime = clock();
  return(morphoBifurcationPoints(src.first, src.second, src.third,
                                 dest.first, dest.second,
                                 neighborOffset));
  //StopTime(startTime);
}

template<class Image1, class Image2, class NB>
inline
int morphoBifurcationPoints(const Image1 & imin, Image2 & imout, NB & nb)
{
  return(morphoBifurcationPoints(srcImageRange(imin), destImage(imout), nb));
}


template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE>
int morphoLabel(Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
        Iterator2 destUpperLeft, Accessor2 desta,
      NBTYPE & nbOffset)
{
  std::queue<Diff2D> Q;

  int label = 1;

  int width  = srcLowerRight.x - srcUpperLeft.x;
    int height = srcLowerRight.y - srcUpperLeft.y;

  typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
  typedef typename NBTYPE::SIZETYPE SIZETYPE;

  Diff2D o0(0,0);

  for(o0.y = 0; o0.y < height; ++o0.y)
  {
    for(o0.x = 0; o0.x < width; ++o0.x)
    {
      if( (desta(destUpperLeft,o0) == 0) && (srca(srcUpperLeft, o0) != 0))
      {
        desta.set(label, destUpperLeft, o0);
        Q.push(o0);
        while(!Q.empty())
        {
          // take the first pixel out of the queue.
          Diff2D o1 = Q.front(); Q.pop();

          // look to the neighborhood.
          for(ITERATORTYPE iter = nbOffset.begin();
            iter != nbOffset.end();
            ++iter)
          {
            Diff2D o2 = o1 + *iter;
            // if the neighbor is not outside the image
            if(!nbOffset.isOutsidePixel(o2))
            {
              if( (srca(srcUpperLeft,o2) != 0) && (desta(destUpperLeft,o2) == 0))
              {
                desta.set(label, destUpperLeft, o2);
                Q.push(o2);
              }
            }  // end if not outside pixel
          }  // end for (neighborhood)
        } // end while !Q.empty()

        label++;
      } // end if controlimage != DONE
    } // end x-loop
  }  // end y-loop

  return(label-1);
} // end of function morphoLabel

////////////////
// morphoLabel
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE>
inline
int morphoLabel(vigra::triple<Iterator1, Iterator1, Accessor1> src,
      vigra::pair<Iterator2, Accessor2> dest,
      NBTYPE & neighborOffset)
{
  clock_t startTime = clock();
  return(morphoLabel(src.first, src.second, src.third,
               dest.first, dest.second,
              neighborOffset));
  //StopTime(startTime);
}

template<class Image1, class Image2, class NB>
inline
int morphoLabel(const Image1 & imin, Image2 & imout, NB & nb)
{
  return(morphoLabel(srcImageRange(imin), destImage(imout), nb));
}


//////////////////
// morphoWatershed
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE,
     class MinmaxFunctor,
     class PriorityFunctor>
void morphoWatershed(Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
         Iterator2 destUpperLeft, Accessor2 desta,
         NBTYPE & nbOffset,
         MinmaxFunctor minmax,
         PriorityFunctor priority
         )
{
  clock_t startTime = clock();

  // for the work image (negative labels are allowed)
  const int WS_QUEUED = -1;
  const int WS_NOT_PROCESSED = 0;
  const int WS_WSLABEL = -2;
  // for the output image: negative labels are not allowed).
  const int OUT_WSLABEL = 0;

  unsigned long insertionOrder = 0;

  typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
  typedef typename NBTYPE::SIZETYPE SIZETYPE;
  typedef typename Accessor1::value_type VALUETYPE;

  typedef Pixel2D<VALUETYPE> PIX;


  int width  = srcLowerRight.x - srcUpperLeft.x;
  int height = srcLowerRight.y - srcUpperLeft.y;

  vigra::BasicImage<int> labelImage(width, height);
  vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft();
  vigra::BasicImage<int>::Accessor lab;
  typedef int LABTYPE;

  utilityMinMaxLabel(srcUpperLeft, srcLowerRight, srca,
               labUpperLeft, lab,
             minmax,
             nbOffset);

  std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);
  std::queue<Diff2D> Q;

  Diff2D o0(0,0);

  VALUETYPE maxval = srca(srcUpperLeft, o0);

  // initialization of the hierarchical queue
  for(o0.y = 0; o0.y < height; ++o0.y)
  {
    for(o0.x = 0; o0.x < width; ++o0.x)

    {

      maxval = std::max(srca(srcUpperLeft, o0), maxval);

      LABTYPE label = lab(labUpperLeft, o0);
      desta.set(label, destUpperLeft, o0);
      if(label > WS_NOT_PROCESSED)
      {
        // look to the neighborhood.
        for(ITERATORTYPE iter = nbOffset.begin();
          iter != nbOffset.end();
          ++iter)
        {
          Diff2D o1 = o0 + *iter;
          // if the neighbor is not outside the image
          // and if it has no label and if it is not in the queue
          if(    (!nbOffset.isOutsidePixel(o1))
            && (lab(labUpperLeft, o1) == WS_NOT_PROCESSED))
          {
            VALUETYPE priority = std::max(srca(srcUpperLeft, o1), srca(srcUpperLeft, o0));
            PQ.push(PIX(priority, o1, insertionOrder++));
            lab.set(WS_QUEUED, labUpperLeft, o1);
          }
        } // end for neighborhood
      } // end if label
    } // end x-loop
  } // end y-loop

  // until the hierarchical queue is empty ...
  while(!PQ.empty())
  {
    PIX px = PQ.top();
    PQ.pop();
    Diff2D o0 = px.offset;

    // normal flooding procedure
    int label1 = WS_NOT_PROCESSED;
    int label2 = WS_NOT_PROCESSED;

    VALUETYPE currentval = px.value;

    // look to the neighborhood to determine the label of pixel o0.
    for(ITERATORTYPE iter = nbOffset.begin();
      iter != nbOffset.end();
      ++iter)
    {
      Diff2D o1 = o0 + *iter;
      if(!nbOffset.isOutsidePixel(o1))
      {
        LABTYPE label_o1 = lab(labUpperLeft, o1);
        // first case: neighbor pixel has not been processed.
        if(label_o1 == WS_NOT_PROCESSED)
        {
          VALUETYPE priority = std::max(srca(srcUpperLeft, o1), currentval);
          PQ.push(PIX(priority, o1, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o1);
        }

        // second case: neighbor pixel is already in the queue
        // or part of the watershed line
        // nothing is to be done, then.

        // third case: the neighbor has a label
        // but a real label (not a WSL label)
        if(label_o1 > WS_NOT_PROCESSED)
        {
          label2 = label_o1;

          if(label1 == 0)
          {
            // in this case, the label is the first
            // which has been found in the neighborhood.
            label1 = label2;

            lab.set(label1, labUpperLeft, o0);
            desta.set(label1, destUpperLeft, o0);
          }
          else
          {
            // in this case, a label has already been assigned to o0.
            // o0 is part of the watershed line.
            if(label1 != label2)
            {
              lab.set(WS_WSLABEL, labUpperLeft, o0);
              desta.set(OUT_WSLABEL, destUpperLeft, o0);
            }

          }

        }
      }
    } // end for neighborhood

    // if there was no label assigned to the pixel
    // this happens if the all neighbors
    // are either part of the watershed line or did not
    // yet come out of the queue. There is at least one neighbor
    // which is part of the WSL in this case.
    if(label1 == WS_NOT_PROCESSED)
    {
      if(currentval < maxval) {
        // reinsertion of the point into the queue with incremented priority
          PQ.push(PIX(currentval+1, o0, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o0);
      }
      else {
        // insertion of the point into the final queue
        Q.push(o0);
      }
    }

  } // end of PRIORITY QUEUE

  while(!Q.empty()) {
    // The queue corresponds to the points that only have WSL neighbors.
    // They are also set to WSL
    Diff2D o0 = Q.front(); Q.pop();
    desta.set(OUT_WSLABEL, destUpperLeft, o0);
  }

} // end of function

/////////////////
// Watershed
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE>
void morphoWatershed(vigra::triple<Iterator1, Iterator1, Accessor1> src,
           vigra::pair<Iterator2, Accessor2> dest,
           NBTYPE & neighborOffset)
{

  typedef typename Accessor1::value_type val_type;
  typedef Pixel2D<val_type> PIX;

  morphoWatershed(src.first, src.second, src.third,
            dest.first, dest.second,
           neighborOffset,
            IsSmaller<typename Accessor1::value_type, typename Accessor1::value_type>(),
            PriorityBottomUp<val_type>());

}

template<class Image1, class Image2, class NBTYPE>
void morphoWatershed(const Image1 & imin, Image2 & imout, NBTYPE & nbOffset)
{
  morphoWatershed(srcImageRange(imin), destImage(imout), nbOffset);
}



//////////////////
// morphoWatershed
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE,
     class MinmaxFunctor,
     class PriorityFunctor>
void morphoSelectiveWatershed(Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
                              Iterator2 destUpperLeft, Accessor2 desta,
                              int dyn_thresh,
                              NBTYPE & nbOffset,
                              MinmaxFunctor minmax,
                              PriorityFunctor priority)
{
  clock_t startTime = clock();

  // for the work image (negative labels are allowed)
  const int WS_QUEUED = -1;
  const int WS_NOT_PROCESSED = 0;
  const int WS_WSLABEL = -2;
  // for the output image: negative labels are not allowed).
  const int OUT_WSLABEL = 0;

  unsigned long insertionOrder = 0;

  typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
  typedef typename NBTYPE::SIZETYPE SIZETYPE;
  typedef typename Accessor1::value_type VALUETYPE;
  typedef typename Accessor2::value_type OUTVALUETYPE;

  typedef Pixel2D<VALUETYPE> PIX;


  int width  = srcLowerRight.x - srcUpperLeft.x;
  int height = srcLowerRight.y - srcUpperLeft.y;

  vigra::BasicImage<int> labelImage(width, height);
  vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft();
  vigra::BasicImage<int>::Accessor lab;
  typedef int LABTYPE;

  int numberOfMinima = utilityMinMaxLabel(srcUpperLeft, srcLowerRight, srca,
                                          labUpperLeft, lab,
                                          minmax,
                                          nbOffset);

  // structures for selectivity

  // equivalence takes the label of the lake with which it has been fused.
  // at the moment, this is simply i.
  std::vector<int> equivalence(numberOfMinima + 1);
  for(std::vector<int>::size_type i = 0; i != equivalence.size(); ++i)
      equivalence[i] = i;

  // the vector containing the dynamics
  std::vector<VALUETYPE> dynamics(numberOfMinima + 1);
  for(int i = 0; i != numberOfMinima + 1; ++i)
      dynamics.push_back(0);
  VALUETYPE maxDyn = 0;
  LABTYPE labelOfMaxDyn = -1;

  // to take the values of the minma
  std::vector<VALUETYPE> valOfMin(numberOfMinima + 1);

  // Queue for watershed
  std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);
  std::queue<Diff2D> Q;

  Diff2D o0(0,0);

  VALUETYPE maxval = srca(srcUpperLeft, o0);

  // initialization of the hierarchical queue
  for(o0.y = 0; o0.y < height; ++o0.y)
  {
    for(o0.x = 0; o0.x < width; ++o0.x)

    {

      maxval = std::max(srca(srcUpperLeft, o0), maxval);

      LABTYPE label = lab(labUpperLeft, o0);
      //desta.set(label, destUpperLeft, o0);
      if(label > WS_NOT_PROCESSED)
      {

        // the value of the minimum is assigned.
        valOfMin[label] = srca(srcUpperLeft, o0);

        // look to the neighborhood.
        for(ITERATORTYPE iter = nbOffset.begin();
          iter != nbOffset.end();
          ++iter)
        {
          Diff2D o1 = o0 + *iter;
          // if the neighbor is not outside the image
          // and if it has no label and if it is not in the queue
          if(    (!nbOffset.isOutsidePixel(o1))
            && (lab(labUpperLeft, o1) == WS_NOT_PROCESSED))
          {
            VALUETYPE priority = std::max(srca(srcUpperLeft, o1), srca(srcUpperLeft, o0));
            PQ.push(PIX(priority, o1, insertionOrder++));
            lab.set(WS_QUEUED, labUpperLeft, o1);
          }
        } // end for neighborhood
      } // end if label
    } // end x-loop
  } // end y-loop

  // until the hierarchical queue is empty ...
  while(!PQ.empty())
  {
    PIX px = PQ.top();
    PQ.pop();
    Diff2D o0 = px.offset;
    VALUETYPE level = px.value;

    // normal flooding procedure
    int label1 = WS_NOT_PROCESSED;
    int label2 = WS_NOT_PROCESSED;
    //int losing_label = WS_NOT_PROCESSED;

    VALUETYPE currentval = px.value;

    // look to the neighborhood to determine the label of pixel o0.
    for(ITERATORTYPE iter = nbOffset.begin();
      iter != nbOffset.end();
      ++iter)
    {
      Diff2D o1 = o0 + *iter;
      if(!nbOffset.isOutsidePixel(o1))
      {
        LABTYPE label_o1 = lab(labUpperLeft, o1);
        // first case: pixel has not been processed.
        if(label_o1 == WS_NOT_PROCESSED)
        {
          VALUETYPE priority = std::max(srca(srcUpperLeft, o1), currentval);
          PQ.push(PIX(priority, o1, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o1);
        }

        // second case: neighbor pixel is already in the queue
        // or part of the watershed line
        // nothing is to be done then.

        // third case: the neighbor has a label
        // but a real label (not a WSL label)
        if(label_o1 > WS_NOT_PROCESSED)
        {
          label2 = label_o1;

          // find the original label to which the lake has been fused.
          while(label2 != equivalence[label2])
              label2 = equivalence[label2];

          if(label1 == WS_NOT_PROCESSED)
          {
            // in this case, the label is the first
            // which has been found in the neighborhood.
            label1 = label2;

            lab.set(label1, labUpperLeft, o0);
            //desta.set(label1, destUpperLeft, o0);
          }
          else
          {

            // in this case, a label has already been assigned to o0.
            // o0 is part of the watershed line.
            if(label1 != label2)
            {
              // in this case, we have a meeting point of two lakes.
              // we therefore have to fuse the two lakes or to assign the WSL label
              if(minmax(valOfMin[label1], valOfMin[label2]))
              {
                VALUETYPE curr_dyn = level - valOfMin[label2];

                // If the lower dynamic is below the thresh, no water shed line is built.
                // we build a watershed line only if the dynamic is higher than the threshold
                if(curr_dyn > dyn_thresh) {

                  // Watershed line.
                  lab.set(WS_WSLABEL, labUpperLeft, o0);
                  desta.set(OUT_WSLABEL, destUpperLeft, o0);
                }
                else {

                  // fusion: the pixel keeps its label.
                  dynamics[label2] = curr_dyn;
                  equivalence[label2] = label1;
                }
              }
              else
              {
                VALUETYPE curr_dyn = level - valOfMin[label1];

                // If the lower dynamic is below the thresh, no water shed line is built.
                // we build a watershed line only if the dynamic is higher than the threshold
                if(curr_dyn > dyn_thresh) {

                  // Watershed line.
                  lab.set(WS_WSLABEL, labUpperLeft, o0);
                  desta.set(OUT_WSLABEL, destUpperLeft, o0);
                }
                else {

                  // fusion: the pixel keeps its label.
                  dynamics[label1] = curr_dyn;
                  equivalence[label1] = label2;
                  label1 = label2;
                }
              }

            } // end of label1 != label2

          } // end of label1 != 0

        } // end of label_o1 > 0
      } // end of (px not outside)
    } // end for neighborhood

    // if there was no label assigned to the pixel
    // this happens if the all neighbors
    // are either part of the watershed line or did not
    // yet come out of the queue. There is at least one neighbor
    // which is part of the WSL in this case.
    if(label1 == WS_NOT_PROCESSED)
    {
      if(currentval < maxval) {
        // reinsertion of the point into the queue with incremented priority
          PQ.push(PIX(currentval+1, o0, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o0);
      }
      else {
        // insertion of the point into the final queue
        Q.push(o0);
      }
    }

  } // end of PRIORITY QUEUE

  // assignment of final labels (without jumps)
  // first assigning the final labels to the original labels
  // which have been maintained
  // note that equivalence[0], k=0, final_label[0] is not used and irrelevant.
  int k = 1;
  std::vector<int> final_label(numberOfMinima + 1);
  for(std::vector<int>::size_type i = 1; i != equivalence.size(); ++i) {
    if(i == equivalence[i]) {
      final_label[i] = k;
      k++;
    }
  }

  // assignment of final labels to the fused labels
  for(std::vector<int>::size_type i = 0; i != equivalence.size(); ++i) {
    int label = equivalence[i];

    // find the original label to which the lake has been fused.
    while(label != equivalence[label])
      label = equivalence[label];

    final_label[i] = final_label[label];
    //std::cout << i << " : " << final_label[i] << std::endl;
  }


  // assigning the final labels to pixels.
  for(o0.y = 0; o0.y < height; ++o0.y)
  {
    for(o0.x = 0; o0.x < width; ++o0.x)
    {
      LABTYPE label = lab(labUpperLeft, o0);
      if(label > WS_NOT_PROCESSED) {
        OUTVALUETYPE val = (OUTVALUETYPE)final_label[label];
        desta.set(val, destUpperLeft, o0);
      }
    }
  }

  // setting the watershed label to the points still in the queue
  // (the points that could not be assigned)
  while(!Q.empty()) {
    // The queue corresponds to the points that only have WSL neighbors.
    // They are also set to WSL
    Diff2D o0 = Q.front(); Q.pop();
    desta.set(OUT_WSLABEL, destUpperLeft, o0);
    //desta.set(23, destUpperLeft, o0);
  }

} // end of function

//////////////////////
// Selective Watershed
template<class Iterator1, class Accessor1,
         class Iterator2, class Accessor2,
         class NBTYPE>
void morphoSelectiveWatershed(vigra::triple<Iterator1, Iterator1, Accessor1> src,
                              vigra::pair<Iterator2, Accessor2> dest,
                              int dyn_thresh,
                              NBTYPE & neighborOffset)
{

  typedef typename Accessor1::value_type val_type;
  typedef Pixel2D<val_type> PIX;

  morphoSelectiveWatershed(src.first, src.second, src.third,
                           dest.first, dest.second,
                           dyn_thresh,
                           neighborOffset,
                           IsSmaller<typename Accessor1::value_type, typename Accessor1::value_type>(),
                           PriorityBottomUp<val_type>());

}

template<class Image1, class Image2, class NBTYPE>
void morphoSelectiveWatershed(const Image1 & imin, Image2 & imout,
                              int dyn_thresh,
                              NBTYPE & nbOffset)
{
  morphoSelectiveWatershed(srcImageRange(imin), destImage(imout),
                           dyn_thresh,
                           nbOffset);
}


/////////////////
// Thalweg
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class NBTYPE>
void morphoThalweg(vigra::triple<Iterator1, Iterator1, Accessor1> src,
         vigra::pair<Iterator2, Accessor2> dest,
         NBTYPE & neighborOffset)
{
  typedef typename Accessor1::value_type val_type;
  typedef Pixel2D<val_type> PIX;

  morphoWatershed(src.first, src.second, src.third,
            dest.first, dest.second,
            neighborOffset,
            IsGreater<val_type, val_type>(),
            PriorityTopDown<val_type>());
}

template<class Image1, class Image2, class NBTYPE>
void morphoThalweg(const Image1 & imin, Image2 & imout, NBTYPE & nbOffset)
{
  morphoThalweg(srcImageRange(imin), destImage(imout), nbOffset);
}

/////////////////////////
// ImConstrainedWatershed
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class Iterator3, class Accessor3,
     class NBTYPE,
     class PriorityFunctor>
void morphoConstrainedWatershed(Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
               Iterator2 markerUpperLeft, Iterator2 markerLowerRight, Accessor2 marka,
               Iterator3 destUpperLeft, Accessor3 desta,
               NBTYPE & nbOffset,
               PriorityFunctor priority
               )
{
  // for the work image (negative labels are allowed)
  const int WS_QUEUED = -1;
  const int WS_NOT_PROCESSED = 0;
  const int WS_WSLABEL = -2;
  // for the output image: negative labels are not allowed).
  const int OUT_WSLABEL = 0;

  unsigned long insertionOrder = 0;

  typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
  typedef typename NBTYPE::SIZETYPE SIZETYPE;
  typedef typename Accessor1::value_type VALUETYPE;

  typedef Pixel2D<VALUETYPE> PIX;

  int width  = srcLowerRight.x - srcUpperLeft.x;
  int height = srcLowerRight.y - srcUpperLeft.y;

  vigra::BasicImage<int> labelImage(width, height);
  vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft();
  vigra::BasicImage<int>::Accessor lab;

  morphoLabel(markerUpperLeft, markerLowerRight, marka,
              labUpperLeft, lab,
              nbOffset);

  std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);
  std::queue<Diff2D> Q;

  Diff2D o0(0,0);

  VALUETYPE maxval = srca(srcUpperLeft, o0);

  // initialization of the hierarchical queue
  for(o0.y = 0; o0.y < height; ++o0.y)
  {
    for(o0.x = 0; o0.x < width; ++o0.x)
    {
      desta.set(lab(labUpperLeft, o0), destUpperLeft, o0);
      if(lab(labUpperLeft, o0) > WS_NOT_PROCESSED)
      {
        // look to the neighborhood.
        for(SIZETYPE i = 0; i < nbOffset.numberOfPixels(); ++i)
        {
          // if the neighbor is not outside the image
          // and if it has no label and if it is not in the queue
          if(    (!nbOffset.isOutsidePixel(o0, nbOffset[i]))
            && (lab(labUpperLeft, o0 + nbOffset[i]) == WS_NOT_PROCESSED))
          {
            Diff2D o1 = o0 + nbOffset[i];
            VALUETYPE priority = std::max(srca(srcUpperLeft, o1), srca(srcUpperLeft, o0));
            PQ.push(PIX(priority, o1, insertionOrder++));
            lab.set(WS_QUEUED, labUpperLeft, o1);
          }
        } // end for neighborhood
      } // end if label
    } // end x-loop
  } // end y-loop

  // until the hierarchical queue is empty ...
  while(!PQ.empty())
  {
    PIX px = PQ.top();
    PQ.pop();
    Diff2D o0 = px.offset;

    // normal flooding procedure
    int label1 = WS_NOT_PROCESSED;
    int label2 = WS_NOT_PROCESSED;

    // the current flooding value is taken from the queue entry.
    // it is not necessarily the same value as in the original image,
    // because some lower regions might not have been flooded
    // (either because of a buttonhole or because of the constraint).
    VALUETYPE currentval = px.value;

    // look to the neighborhood to determine the label of pixel o0.
    for(SIZETYPE i = 0; i < nbOffset.numberOfPixels(); ++i)
    {
      if(!nbOffset.isOutsidePixel(o0, nbOffset[i]))
      {
        Diff2D o1 = o0 + nbOffset[i];
        // first case: pixel has not been processed.
        if(lab(labUpperLeft, o1) == WS_NOT_PROCESSED)
        {
          VALUETYPE priority = std::max(srca(srcUpperLeft, o1), currentval);
          PQ.push(PIX(priority, o1, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o1);
        }

        // second case: neighbor pixel is already in the queue:
        // nothing is to be done, then.

        // third case: the neighbor has a label
         if(lab(labUpperLeft, o1) > WS_NOT_PROCESSED)
        {
          label2 = lab(labUpperLeft, o1);

          if(label1 == 0)
          {
            // in this case, the label is the first
            // which has been found in the neighborhood.
            label1 = label2;

            lab.set(label1, labUpperLeft, o0);
            desta.set(label1, destUpperLeft, o0);
          }
          else
          {
            // in this case, a label has already been assigned to o0.
            // o0 is part of the watershed line.
            if(label1 != label2)
            {
              lab.set(WS_WSLABEL, labUpperLeft, o0);
              desta.set(OUT_WSLABEL, destUpperLeft, o0);
            }
          }

        }
      }
    } // end for neighborhood

    // if there was no label assigned to the pixel
    // (this can happen in some pathological but not uncommon situations)
    if(label1 == WS_NOT_PROCESSED)
    {
      if(currentval < maxval) {
        // in this case the pixel is pushed back to the queue with
        // increased priority level (so that it will be checked out again
        // at the next grey level).
          PQ.push(PIX(currentval+1, o0, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o0);
      }
      else {
        // if the maximum level is already reached, the pixel remains
        // outside the queue and is pushed to the final non-hierarchical queue.
          Q.push(o0);
      }
    }

  } // end of PRIORITY QUEUE

  // all points in the rest queue are given the watershedline-label.
  while(!Q.empty()) {
     Diff2D o0 = Q.front(); Q.pop();
     desta.set(OUT_WSLABEL, destUpperLeft, o0);
  }

} // end of function

/////////////////////////
// ImConstrainedWatershed
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class Iterator3, class Accessor3,
     class NBTYPE,
     class PriorityFunctor>
void morphoConstrainedWatershedOpt(Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
               Iterator2 markerUpperLeft, Iterator2 markerLowerRight, Accessor2 marka,
               Iterator3 destUpperLeft, Accessor3 desta,
               NBTYPE & nbOffset,
               PriorityFunctor priority
               )
{
  // for the work image (negative labels are allowed)
  const int WS_QUEUED = -1;
  const int WS_NOT_PROCESSED = 0;
  const int WS_WSLABEL = -2;
  // for the output image: negative labels are not allowed).
  const int OUT_WSLABEL = 0;

  unsigned long insertionOrder = 0;

  typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
  typedef typename NBTYPE::SIZETYPE SIZETYPE;
  typedef typename Accessor1::value_type VALUETYPE;

  typedef Pixel2D<VALUETYPE> PIX;

  int width  = srcLowerRight.x - srcUpperLeft.x;
    int height = srcLowerRight.y - srcUpperLeft.y;

  vigra::BasicImage<int> labelImage(width, height);
  vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft();
  vigra::BasicImage<int>::Accessor lab;
  typedef int LABTYPE;

  morphoLabel(markerUpperLeft, markerLowerRight, marka,
              labUpperLeft, lab,
              nbOffset);
  //globalDebEnv.DebugWriteImage(labelImage, "label");

  std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);

  Diff2D o0(0,0);

  // initialization of the hierarchical queue
  for(o0.y = 0; o0.y < height; ++o0.y)
  {
    for(o0.x = 0; o0.x < width; ++o0.x)
    {
      LABTYPE  label_o0 = lab(labUpperLeft, o0);
      desta.set(label_o0, destUpperLeft, o0);
      if(label_o0 > WS_NOT_PROCESSED)
      {
        // look to the neighborhood.
        for(ITERATORTYPE iter = nbOffset.begin();
          iter != nbOffset.end();
          ++iter)
        {
          Diff2D o1 = o0 + *iter;

          // if the neighbor is not outside the image
          // and if it has no label and if it is not in the queue
          if(    (!nbOffset.isOutsidePixel(o1))
            && (lab(labUpperLeft, o1) == WS_NOT_PROCESSED))
          {
            VALUETYPE priority = std::max(srca(srcUpperLeft, o1), srca(srcUpperLeft, o0));
            PQ.push(PIX(priority, o1, insertionOrder++));
            lab.set(WS_QUEUED, labUpperLeft, o1);
          }
        } // end for neighborhood
      } // end if label
    } // end x-loop
  } // end y-loop

  // until the hierarchical queue is empty ...
  while(!PQ.empty())
  {
    PIX px = PQ.top();
    PQ.pop();
    Diff2D o0 = px.offset;

    // normal flooding procedure
    int label1 = WS_NOT_PROCESSED;
    int label2 = WS_NOT_PROCESSED;

    // look to the neighborhood to determine the label of pixel o0.
    for(ITERATORTYPE iter = nbOffset.begin();
      iter != nbOffset.end();
      ++iter)
    {
      Diff2D o1 = o0 + *iter;
      if(!nbOffset.isOutsidePixel(o1))
      {
        label2 = lab(labUpperLeft, o1);
        // first case: pixel has not been processed.
        if(label2 == WS_NOT_PROCESSED)
        {
          // the priority is at least the current pixel value (value of o0).
          VALUETYPE priority = std::max(srca(srcUpperLeft, o1), srca(srcUpperLeft, o0));
          PQ.push(PIX(priority, o1, insertionOrder++));
          lab.set(WS_QUEUED, labUpperLeft, o1);
        }

        // second case: neighbor pixel is already in the queue:
        // nothing is to be done, then.

        // third case: the neighbor has a label
         if(label2 > WS_NOT_PROCESSED)
        {
          if(label1 == 0)
          {
            // in this case, the label is the first
            // which has been found in the neighborhood.
            label1 = label2;

            lab.set(label1, labUpperLeft, o0);
            desta.set(label1, destUpperLeft, o0);
          }
          else
          {
            // in this case, a label has already been assigned to o0.
            // o0 is part of the watershed line.
            if(label1 != label2)
            {
              lab.set(WS_WSLABEL, labUpperLeft, o0);
              desta.set(OUT_WSLABEL, destUpperLeft, o0);
            }
          }

        }
      }
    } // end for neighborhood

    // if the pixel has not been assigned a label
    if(label1 == WS_NOT_PROCESSED)
    {
      // we know that this is not correct, but we think
      // that the differences do not concern many pixels.
      lab.set(WS_WSLABEL, labUpperLeft, o0);
      desta.set(OUT_WSLABEL, destUpperLeft, o0);
    }

  } // end of PRIORITY QUEUE

} // end of function

/////////////////
// Watershed
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class Iterator3, class Accessor3,
     class NBTYPE>
void morphoConstrainedWatershed(vigra::triple<Iterator1, Iterator1, Accessor1> src,
              vigra::triple<Iterator2, Iterator2, Accessor2> marker,
                 vigra::pair<Iterator3, Accessor3> dest,
                NBTYPE & neighborOffset)
{
  typedef typename Accessor1::value_type val_type;

  morphoConstrainedWatershed(src.first, src.second, src.third,
               marker.first, marker.second, marker.third,
                   dest.first, dest.second,
                   neighborOffset,
                       PriorityBottomUp<val_type>());
}

template<class Image1, class Image2, class Image3, class NBTYPE>
void morphoConstrainedWatershed(const Image1 & imin, const Image2 & marker, Image3 & imout, NBTYPE & nbOffset)
{
  morphoConstrainedWatershed(srcImageRange(imin), srcImageRange(marker), destImage(imout), nbOffset);
}

};
};
#endif /*MORPHO_WATERSHED_HXX_*/




