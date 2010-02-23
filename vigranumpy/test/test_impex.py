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

def test_writeAndReadVolumeHDF5():
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
