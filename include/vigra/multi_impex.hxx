/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_MULTI_IMPEX_HXX
#define VIGRA_MULTI_IMPEX_HXX

#include <memory>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <sys/types.h>
#include <libgen.h>
#include <dirent.h>
#include <errno.h>
#include <regex.h>
#include "vigra/basicimageview.hxx"
#include "vigra/impex.hxx"
#include "vigra/multi_array.hxx"

namespace vigra {

/** \addtogroup VolumeImpex Import/export of volume data.
*/

//@{

/********************************************************/
/*                                                      */
/*                    importVolume                      */
/*                                                      */
/********************************************************/

/** \brief Function for importing a 3D volume.

    The data are expected to be stored in a by-slice manner,
    where the slices are enumerated from <tt>name_base+"000"+name_ext</tt>
    (the number of zeros depends on the depth of the given volume).
    The <tt>volume</tt> will be reshaped to match the number and size of the images
    found. 

    <b>\#include</b>
    "<a href="multi_impex_8hxx-source.html">vigra/multi_impex.hxx</a>"

    Namespace: vigra
*/
template <class T, class Allocator>
void importVolume (MultiArray <3, T, Allocator> & volume,
                   const std::string &name_base,
                   const std::string &name_ext)
{
    // find out how many images we have
    char * name, * path, * base;
    name = std::strdup(name_base.c_str());
    base = basename(name);
    path = dirname(name);
    
    DIR * dir = opendir(path);
    if(!dir)
    {
        free(name);
        vigra_fail("importVolume(): Unable to open directory.");
    }
      
    regex_t    re;
    std::string pattern = std::string(base) + "[0-9]+" + name_ext;  
    if (regcomp(&re, pattern.c_str(), REG_EXTENDED|REG_NOSUB) != 0) 
    {
        closedir(dir);
        free(name);
        vigra_fail("importVolume(): Failed to compile regular expression.");
    }
    
    int count = 0;
    dirent * dp;
    errno = 0;
    while ((dp = readdir(dir)) != NULL) 
    {
        if(regexec(&re, dp->d_name, (size_t) 0, NULL, 0) == 0)
            ++count;
    }

    regfree(&re);
    closedir(dir);
    free(name);
    
    vigra_precondition(errno == 0,
          "importVolume(): Error while searching for images.");
    vigra_precondition(count > 0,
          "importVolume(): No image with given basename found.");
          
    int numlen = static_cast <int> (std::ceil (std::log10 (count)));
    for (unsigned int i = 0; i < count; ++i) 
    {
        // build the filename
        std::stringstream stream;
        stream << std::setfill ('0') << std::setw (numlen) << i;
        std::string name_num;
        stream >> name_num;
        std::string name = name_base + name_num + name_ext;

        // import the image
        ImageImportInfo info (name.c_str ());
        
        // reshape the array according to size of first image
        if(i == 0)
        {
            volume.reshape(typename 
              MultiArray <3, T>::difference_type(info.width(), info.height(), count));
        }

        // generate a basic image view to the current layer
        MultiArrayView <2, T> array_view (volume.bindOuter (i));
        BasicImageView <T> view = makeBasicImageView (array_view);
        vigra_precondition(view.size() == info.size(),
            "importVolume(): image size mismatch.");
            
        importImage (info, destImage(view));
    }
}


/********************************************************/
/*                                                      */
/*                    exportVolume                      */
/*                                                      */
/********************************************************/

/** \brief Function for exporting a 3D volume.

    The volume is exported in a by-slice manner, where the number of slices equals
    the depth of the volume. The file names will be enumerated like 
    <tt>name_base+"000"+name_ext</tt>, <tt>name_base+"001"+name_ext</tt> etc.
    (the actual number of zeros depends on the depth).

    <b>\#include</b>
    "<a href="multi_impex_8hxx-source.html">vigra/multi_impex.hxx</a>"

    Namespace: vigra
*/
template <class T, class Tag>
void exportVolume (MultiArrayView <3, T, Tag> const & volume,
                   const std::string &name_base,
                   const std::string &name_ext)
{
    
    const unsigned int depth = volume.shape (2);
    int numlen = static_cast <int> (std::ceil (std::log10 (depth)));
    for (unsigned int i = 0; i < depth; ++i) 
    {

        // build the filename
        std::stringstream stream;
        stream << std::setfill ('0') << std::setw (numlen) << i;
        std::string name_num;
        stream >> name_num;
        std::string name = name_base + name_num + name_ext;

        // generate a basic image view to the current layer
        MultiArrayView <2, T, Tag> array_view (volume.bindOuter (i));
        BasicImageView <T> view = makeBasicImageView (array_view);

        // import the image
        ImageExportInfo info (name.c_str ());
        exportImage (srcImageRange(view), info);
    }
}

//@}

} // namespace vigra

#endif // VIGRA_MULTI_IMPEX_HXX
