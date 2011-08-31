#######################################################################
#                                                                      
#         Copyright 2009-2010 by Ullrich Koethe                        
#                                                                      
#    This file is part of the VIGRA computer vision library.           
#    The VIGRA Website is                                              
#        http://hci.iwr.uni-heidelberg.de/vigra/                       
#    Please direct questions, bug reports, and contributions to        
#        ullrich.koethe@iwr.uni-heidelberg.de    or                    
#        vigra@informatik.uni-hamburg.de                               
#                                                                      
#    Permission is hereby granted, free of charge, to any person       
#    obtaining a copy of this software and associated documentation    
#    files (the "Software"), to deal in the Software without           
#    restriction, including without limitation the rights to use,      
#    copy, modify, merge, publish, distribute, sublicense, and/or      
#    sell copies of the Software, and to permit persons to whom the    
#    Software is furnished to do so, subject to the following          
#    conditions:                                                       
#                                                                      
#    The above copyright notice and this permission notice shall be    
#    included in all copies or substantial portions of the             
#    Software.                                                         
#                                                                      
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     
#    OTHER DEALINGS IN THE SOFTWARE.                                   
#                                                                      
#######################################################################

import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

from nose.tools import assert_equal, raises
import numpy as np
import vigra.vigranumpycore # FIXME: why is this needed? (without, impex returns ndarray)
import vigra.impex as im
import vigra.arraytypes as at

#in the hope, that functions are tested in C++, we basicly test return types

image=at.RGBImage(np.random.rand(10,10,3)*255,dtype=np.float32, 
                  axistags=at.VigraArray.defaultAxistags(3, 'V'))
image2=at.RGBImage(np.random.rand(20,20,3)*255,dtype=np.uint8, 
                  axistags=at.VigraArray.defaultAxistags(3, 'V'))
scalar_image=at.ScalarImage(np.random.rand(10,10)*255,dtype=np.float32, 
                  axistags=at.AxisTags(at.AxisInfo.x, at.AxisInfo.y))
volume256=at.Volume(np.random.rand(8,9,10)*255,dtype=np.uint8, 
                  axistags=at.AxisTags(at.AxisInfo.x, at.AxisInfo.y, at.AxisInfo.z))
volumeFloat=at.Volume(np.random.rand(3,4,5,6)*100,dtype=np.float32, 
                  axistags=at.VigraArray.defaultAxistags(4, 'C'))

def checkEqualData(i1,i2):
    assert_equal(i1.shape, i2.shape)
    assert(np.all(i1==i2))

def checkUnequalData(i1,i2):
    assert_equal(i1.shape, i2.shape)
    assert(np.any(i1!=i2))

def test_multiImageTiff():
    if not 'TIFF' in im.listFormats():
        return

    filename = 'resimage.tif'

    # decompose the RGB image and write the three channels as individual images
    # to a multi-image TIFF
    for i in range(3):
        # the first image requires mode="w" in case the image already exists
        im.writeImage(image2[:,:,i], filename, mode="w" if i == 0 else "a")

    # test different dimensions and data types
    im.writeImage(image, filename, mode="a")
    im.writeImage(image2, filename, mode="a")
    im.writeImage(scalar_image, filename, mode="a")

    # check number of images contained in the file
    assert(im.numberImages(filename) == 6)

    # check for equal data
    for i in range(3):
        img_test = im.readImage(filename, index=i)
        checkEqualData(img_test.dropChannelAxis(), image2[:,:,i])
    checkEqualData(im.readImage(filename, index=3), image)
    checkEqualData(im.readImage(filename, index=4), image2)
    checkEqualData(im.readImage(filename, index=5).dropChannelAxis(), scalar_image)


def test_writeAndReadImageHDF5():
    try:
        import h5py
    except:
        print "Warning: 'import h5py' failed, not executing HDF5 import/export tests"
        return
    
    # positive tests
    # write and read image
    im.writeHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata")
    image_imp = im.readHDF5("hdf5test.hd5", "group/subgroup/imgdata")
    checkEqualData(image,image_imp)
    # write and read scalar image
    im.writeHDF5(scalar_image, "hdf5test.hd5", "group/subgroup/imgdata")
    scalar_image_imp = im.readHDF5("hdf5test.hd5", "group/subgroup/imgdata")
    scalar_image_imp = scalar_image_imp.dropChannelAxis()
    checkEqualData(scalar_image,scalar_image_imp)
    # write multiple sets and check if they are all there afterwards
    im.writeHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata")
    im.writeHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata2")
    image_imp1 = im.readHDF5("hdf5test.hd5", "group/subgroup/imgdata")
    image_imp2 = im.readHDF5("hdf5test.hd5", "group/subgroup/imgdata2")
    checkEqualData(image,image_imp1)
    checkEqualData(image,image_imp2)

    # negative tests
    # write and read image
    image_imp[1,1,1] = 100000
    checkUnequalData(image,image_imp)
    # write and read scalar image
    scalar_image_imp[1,1] = 100000
    checkUnequalData(scalar_image,scalar_image_imp)
    
def test_writeAndReadVolumeHDF5():
    try:
        import h5py
    except:
        return
    
    # positive tests
    # write and read volume
    im.writeHDF5(volume256, "hdf5test.hd5", "group/subgroup/voldata")
    volume256_imp = im.readHDF5("hdf5test.hd5", "group/subgroup/voldata")
    checkEqualData(volume256,volume256_imp)
    # write and read binary volume
    im.writeHDF5(volumeFloat, "hdf5test.hd5", "group/subgroup/voldata")
    volumeFloat_imp = im.readHDF5("hdf5test.hd5", "group/subgroup/voldata")
    checkEqualData(volumeFloat.transposeToDefaultOrder(), volumeFloat_imp)
    # write multiple sets and check if they are all there afterwards
    im.writeHDF5(volume256, "hdf5test.hd5", "group/subgroup/voldata")
    im.writeHDF5(volume256, "hdf5test.hd5", "group/subgroup/voldata2")
    volume256_imp1 = im.readHDF5("hdf5test.hd5", "group/subgroup/voldata")
    volume256_imp1 = volume256_imp1.dropChannelAxis()
    volume256_imp2 = im.readHDF5("hdf5test.hd5", "group/subgroup/voldata2")
    volume256_imp2 = volume256_imp2.dropChannelAxis()
    checkEqualData(volume256,volume256_imp1)
    checkEqualData(volume256,volume256_imp2)

    # negative tests
    # write and read volume
    volume256_imp[1,1,1] = 100000
    checkUnequalData(volume256,volume256_imp)
    # write and read binary volume
    volumeFloat_imp[1,1,1] = 100000
    checkUnequalData(volumeFloat.transposeToDefaultOrder(), volumeFloat_imp)
