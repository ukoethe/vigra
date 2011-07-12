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
import vigra.impex as im
import vigra.arraytypes as at

#in the hope, that functions are tested in C++, we basicly test return types

image=at.RGBImage(np.random.rand(10,10,3)*255,dtype=np.float32)
scalar_image=at.ScalarImage(np.random.rand(10,10)*255,dtype=np.float32)
volume256=at.Volume(np.random.rand(8,9,10)*255,dtype=np.uint8)
volumeFloat=at.Volume(np.random.rand(3,4,5,6)*100,dtype=np.float32)

def checkEqualData(i1,i2):
    assert(i1.shape==i2.shape)
    assert(np.all(i1==i2))

def checkUnequalData(i1,i2):
    assert(i1.shape==i2.shape)
    assert(np.any(i1!=i2))

def test_writeAndReadImageHDF5():
    if not hasattr(im, 'writeImageToHDF5'):
        return
    
    # positive tests
    # write and read image
    im.writeImageToHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata")
    image_imp = im.readImageFromHDF5("hdf5test.hd5", "group/subgroup/imgdata")
    checkEqualData(image,image_imp)
    # write and read scalar image
    im.writeImageToHDF5(scalar_image, "hdf5test.hd5", "group/subgroup/imgdata")
    scalar_image_imp = im.readImageFromHDF5("hdf5test.hd5", "group/subgroup/imgdata")
    checkEqualData(scalar_image,scalar_image_imp)
    # write multiple sets and check if they are all there afterwards
    im.writeImageToHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata")
    im.writeImageToHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata2")
    image_imp1 = im.readImageFromHDF5("hdf5test.hd5", "group/subgroup/imgdata")
    image_imp2 = im.readImageFromHDF5("hdf5test.hd5", "group/subgroup/imgdata2")
    checkEqualData(image,image_imp1)
    checkEqualData(image,image_imp2)

    # negative tests
    # write and read image
    image_imp[1,1,1] = 100000
    checkUnequalData(image,image_imp)
    # write and read scalar image
    scalar_image_imp[1,1] = 100000
    checkUnequalData(scalar_image,scalar_image_imp)
    
    # check that we do the same as h5py
    try:
        import h5py
    except:
        return
        
    h5py_file = h5py.File('hdf5test.hd5', 'w')
    h5py_file.create_dataset('imgdata', data=image.swapaxes(0, 1))
    h5py_file.close() 
    raise
    image_imp3 = im.readImageFromHDF5("hdf5test.hd5", "imgdata")
    checkEqualData(image,image_imp3)
        
    im.writeImageToHDF5(image, "hdf5test.hd5", "group/subgroup/imgdata")
    h5py_file = h5py.File('hdf5test.hd5', 'r')
    image_imp4 = h5py_file['/group/subgroup/imgdata']
    checkEqualData(image,image_imp4.value.swapaxes(0,1))

def test_writeAndReadVolumeHDF5():
    if not hasattr(im, 'writeVolumeToHDF5'):
        return
    
    # positive tests
    # write and read volume
    im.writeVolumeToHDF5(volume256, "hdf5test.hd5", "group/subgroup/voldata")
    volume256_imp = im.readVolumeFromHDF5("hdf5test.hd5", "group/subgroup/voldata")
    checkEqualData(volume256,volume256_imp)
    # write and read binary volume
    im.writeVolumeToHDF5(volumeFloat, "hdf5test.hd5", "group/subgroup/voldata")
    volumeFloat_imp = im.readVolumeFromHDF5("hdf5test.hd5", "group/subgroup/voldata")
    checkEqualData(volumeFloat,volumeFloat_imp)
    # write multiple sets and check if they are all there afterwards
    im.writeVolumeToHDF5(volume256, "hdf5test.hd5", "group/subgroup/voldata")
    im.writeVolumeToHDF5(volume256, "hdf5test.hd5", "group/subgroup/voldata2")
    volume256_imp1 = im.readVolumeFromHDF5("hdf5test.hd5", "group/subgroup/voldata")
    volume256_imp2 = im.readVolumeFromHDF5("hdf5test.hd5", "group/subgroup/voldata2")
    checkEqualData(volume256,volume256_imp1)
    checkEqualData(volume256,volume256_imp2)

    # negative tests
    # write and read volume
    volume256_imp[1,1,1] = 100000
    checkUnequalData(volume256,volume256_imp)
    # write and read binary volume
    volumeFloat_imp[1,1,1] = 100000
    checkUnequalData(volumeFloat,volumeFloat_imp)
    
    # check that we do the same as h5py
    try:
        import h5py
    except:
        return
        
    h5py_file = h5py.File('hdf5test.hd5', 'w')
    h5py_file.create_dataset('voldata', data=volumeFloat.swapaxes(0, 2))
    h5py_file.close() 
    volumeFloat_imp1 = im.readVolumeFromHDF5("hdf5test.hd5", "voldata")
    checkEqualData(volumeFloat, volumeFloat_imp1)
        
    im.writeVolumeToHDF5(volumeFloat, "hdf5test.hd5", "group/subgroup/voldata")
    h5py_file = h5py.File('hdf5test.hd5', 'r')
    volumeFloat_imp2 = h5py_file['/group/subgroup/voldata']
    checkEqualData(volumeFloat, volumeFloat_imp2.value.swapaxes(0,2))
