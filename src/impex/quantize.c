/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 

/*                                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%  This file contains source code adapted from ImageMagick                    %
%                                                                             %
%  ImageMagick is Copyright 1998 E. I. du Pont de Nemours and Company         %
%                                                                             %
%  Permission is hereby granted, free of charge, to any person obtaining a    %
%  copy of this software and associated documentation files ("ImageMagick"),  %
%  to deal in ImageMagick without restriction, including without limitation   %
%  the rights to use, copy, modify, merge, publish, distribute, sublicense,   %
%  and/or sell copies of ImageMagick, and to permit persons to whom the       %
%  ImageMagick is furnished to do so, subject to the following conditions:    %
%                                                                             %
%  The above copyright notice and this permission notice shall be included in %
%  all copies or substantial portions of ImageMagick.                         %
%                                                                             %
%  The software is provided "as is", without warranty of any kind, express or %
%  implied, including but not limited to the warranties of merchantability,   %
%  fitness for a particular purpose and noninfringement.  In no event shall   %
%  E. I. du Pont de Nemours and Company be liable for any claim, damages or   %
%  other liability, whether in an action of contract, tort or otherwise,      %
%  arising from, out of or in connection with ImageMagick or the use or other %
%  dealings in ImageMagick.                                                   %
%                                                                             %
%  Except as contained in this notice, the name of the E. I. du Pont de       %
%  Nemours and Company shall not be used in advertising or otherwise to       %
%  promote the sale, use or other dealings in ImageMagick without prior       %
%  written authorization from the E. I. du Pont de Nemours and Company.       %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
/*
%
%  Realism in computer graphics typically requires using 24 bits/pixel to
%  generate an image.  Yet many graphic display devices do not contain
%  the amount of memory necessary to match the spatial and color
%  resolution of the human eye.  The QUANTIZE program takes a 24 bit
%  image and reduces the number of colors so it can be displayed on
%  raster device with less bits per pixel.  In most instances, the
%  quantized image closely resembles the original reference image.
%
%  A reduction of colors in an image is also desirable for image
%  transmission and real-time animation.
%
%  Function Quantize takes a standard RGB or monochrome images and quantizes
%  them down to some fixed number of colors.
%
%  For purposes of color allocation, an image is a set of n pixels, where
%  each pixel is a point in RGB space.  RGB space is a 3-dimensional
%  vector space, and each pixel, pi,  is defined by an ordered triple of
%  red, green, and blue coordinates, (ri, gi, bi).
%
%  Each primary color component (red, green, or blue) represents an
%  intensity which varies linearly from 0 to a maximum value, cmax, which
%  corresponds to full saturation of that color.  Color allocation is
%  defined over a domain consisting of the cube in RGB space with
%  opposite vertices at (0,0,0) and (cmax,cmax,cmax).  QUANTIZE requires
%  cmax = 255.
%
%  The algorithm maps this domain onto a tree in which each node
%  represents a cube within that domain.  In the following discussion
%  these cubes are defined by the coordinate of two opposite vertices:
%  The vertex nearest the origin in RGB space and the vertex farthest
%  from the origin.
%
%  The tree's root node represents the the entire domain, (0,0,0) through
%  (cmax,cmax,cmax).  Each lower level in the tree is generated by
%  subdividing one node's cube into eight smaller cubes of equal size.
%  This corresponds to bisecting the parent cube with planes passing
%  through the midpoints of each edge.
%
%  The basic algorithm operates in three phases: Classification,
%  Reduction, and Assignment.  Classification builds a color
%  description tree for the image.  Reduction collapses the tree until
%  the number it represents, at most, the number of colors desired in the
%  output image.  Assignment defines the output image's color map and
%  sets each pixel's color by reclassification in the reduced tree.
%  Our goal is to minimize the numerical discrepancies between the original
%  colors and quantized colors (quantization error).
%
%  Classification begins by initializing a color description tree of
%  sufficient depth to represent each possible input color in a leaf.
%  However, it is impractical to generate a fully-formed color
%  description tree in the classification phase for realistic values of
%  cmax.  If colors components in the input image are quantized to k-bit
%  precision, so that cmax= 2k-1, the tree would need k levels below the
%  root node to allow representing each possible input color in a leaf.
%  This becomes prohibitive because the tree's total number of nodes is
%  1 + sum(i=1,k,8k).
%
%  A complete tree would require 19,173,961 nodes for k = 8, cmax = 255.
%  Therefore, to avoid building a fully populated tree, QUANTIZE: (1)
%  Initializes data structures for nodes only as they are needed;  (2)
%  Chooses a maximum depth for the tree as a function of the desired
%  number of colors in the output image (currently log2(colormap size)).
%
%  For each pixel in the input image, classification scans downward from
%  the root of the color description tree.  At each level of the tree it
%  identifies the single node which represents a cube in RGB space
%  containing the pixel's color.  It updates the following data for each
%  such node:
%
%    n1: Number of pixels whose color is contained in the RGB cube
%    which this node represents;
%
%    n2: Number of pixels whose color is not represented in a node at
%    lower depth in the tree;  initially,  n2 = 0 for all nodes except
%    leaves of the tree.
%
%    Sr, Sg, Sb: Sums of the red, green, and blue component values for
%    all pixels not classified at a lower depth. The combination of
%    these sums and n2  will ultimately characterize the mean color of a
%    set of pixels represented by this node.
%
%    E: The distance squared in RGB space between each pixel contained
%    within a node and the nodes' center.  This represents the quantization
%    error for a node.
%
%  Reduction repeatedly prunes the tree until the number of nodes with
%  n2 > 0 is less than or equal to the maximum number of colors allowed
%  in the output image.  On any given iteration over the tree, it selects
%  those nodes whose E count is minimal for pruning and merges their
%  color statistics upward. It uses a pruning threshold, Ep, to govern
%  node selection as follows:
%
%    Ep = 0
%    while number of nodes with (n2 > 0) > required maximum number of colors
%      prune all nodes such that E <= Ep
%      Set Ep to minimum E in remaining nodes
%
%  This has the effect of minimizing any quantization error when merging
%  two nodes together.
%
%  When a node to be pruned has offspring, the pruning procedure invokes
%  itself recursively in order to prune the tree from the leaves upward.
%  n2,  Sr, Sg,  and  Sb in a node being pruned are always added to the
%  corresponding data in that node's parent.  This retains the pruned
%  node's color characteristics for later averaging.
%
%  For each node, n2 pixels exist for which that node represents the
%  smallest volume in RGB space containing those pixel's colors.  When n2
%  > 0 the node will uniquely define a color in the output image. At the
%  beginning of reduction,  n2 = 0  for all nodes except a the leaves of
%  the tree which represent colors present in the input image.
%
%  The other pixel count, n1, indicates the total number of colors
%  within the cubic volume which the node represents.  This includes n1 -
%  n2  pixels whose colors should be defined by nodes at a lower level in
%  the tree.
%
%  Assignment generates the output image from the pruned tree.  The
%  output image consists of two parts: (1)  A color map, which is an
%  array of color descriptions (RGB triples) for each color present in
%  the output image;  (2)  A pixel array, which represents each pixel as
%  an index into the color map array.
%
%  First, the assignment phase makes one pass over the pruned color
%  description tree to establish the image's color map.  For each node
%  with n2  > 0, it divides Sr, Sg, and Sb by n2 .  This produces the
%  mean color of all pixels that classify no lower than this node.  Each
%  of these colors becomes an entry in the color map.
%
%  Finally,  the assignment phase reclassifies each pixel in the pruned
%  tree to identify the deepest node containing the pixel's color.  The
%  pixel's value in the pixel array becomes the index of this node's mean
%  color in the color map.
%
%  For efficiency, QUANTIZE requires that the reference image be in a
%  run-length encoded format.
%
%  With the permission of USC Information Sciences Institute, 4676 Admiralty
%  Way, Marina del Rey, California  90292, this code was adapted from module
%  ALCOLS written by Paul Raveling.
%
%  The names of ISI and USC are not used in advertising or publicity
%  pertaining to distribution of the software without prior specific
%  written permission from ISI.
%
%
*/

/*
  Include declarations.
*/
#include <stdio.h>
#include <stdlib.h>
#if defined(_MSC_VER)
#    include <direct.h>
#else
#    include <unistd.h>
#endif
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "vigra/impex.h"
#include "error.h"
#include "utility.h"

/*
  Define declarations.
*/
#define MaxNodes  266817
#define MaxSpan ((1 << MaxTreeDepth)-1)
#define MaxTreeDepth  8
#define NodesInAList  2048

/*
  Structures.
*/
typedef struct _NodeInfo
{
  unsigned char
    id,
    level,
    census;

  unsigned int
    color_number,
    number_unique;

  double
    quantization_error,
    total_red,
    total_green,
    total_blue;

  struct _NodeInfo
    *parent,
    *child[8];
} NodeInfo;

typedef struct _Nodes
{
  NodeInfo
    nodes[NodesInAList];

  struct _Nodes
    *next;
} Nodes;

typedef struct _CubeInfo
{
  NodeInfo
    *root;

  unsigned int
    depth;

  unsigned long
    colors;

  VigraImpexColorPacket
    color,
    *colormap;

  double
    distance,
    pruning_threshold,
    next_pruning_threshold;

  unsigned int
    *squares,
    nodes,
    free_nodes,
    color_number;

  NodeInfo
    *next_node;

  Nodes
    *node_queue;
} CubeInfo;

/*
  Function prototypes.
*/
static NodeInfo
  *vigraImpexInitializeNodeInfo(CubeInfo *,const unsigned int,const unsigned int,
    NodeInfo *);

static unsigned int
  DitherImage(CubeInfo *,VigraImpexImage *);

static void
  ClosestColor(CubeInfo *,const NodeInfo *),
  DefineColormap(CubeInfo *,NodeInfo *),
  PruneLevel(CubeInfo *,const NodeInfo *);

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   A s s i g n m e n t                                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure Assignment generates the output image from the pruned tree.  The
%  output image consists of two parts: (1)  A color map, which is an
%  array of color descriptions (RGB triples) for each color present in
%  the output image;  (2)  A pixel array, which represents each pixel as
%  an index into the color map array.
%
%  First, the assignment phase makes one pass over the pruned color
%  description tree to establish the image's color map.  For each node
%  with n2  > 0, it divides Sr, Sg, and Sb by n2 .  This produces the
%  mean color of all pixels that classify no lower than this node.  Each
%  of these colors becomes an entry in the color map.
%
%  Finally,  the assignment phase reclassifies each pixel in the pruned
%  tree to identify the deepest node containing the pixel's color.  The
%  pixel's value in the pixel array becomes the index of this node's mean
%  color in the color map.
%
%  The format of the Assignment routine is:
%
%      Assignment(color_cube,quantize_info,image)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o quantize_info: Specifies a pointer to an VigraImpexQuantizeInfo structure.
%
%    o image: Specifies a pointer to an VigraImpexImage structure;  returned from
%      ReadImage.
%
%
*/
static void Assignment(CubeInfo *color_cube,VigraImpexQuantizeInfo *quantize_info,
  VigraImpexImage *image)
{
#define AssignImageText  "  Assigning image colors...  "

  register int
    i;

  const NodeInfo
    *node_info;

  register VigraImpexRunlengthPacket
    *p;

  register unsigned short
    index;

  unsigned int
    id;

  /*
    Allocate image colormap.
  */
  if (image->colormap == (VigraImpexColorPacket *) NULL)
    image->colormap=(VigraImpexColorPacket *)
      malloc(color_cube->colors*sizeof(VigraImpexColorPacket));
  else
    image->colormap=(VigraImpexColorPacket *)
      realloc((char *) image->colormap,color_cube->colors*sizeof(VigraImpexColorPacket));
  if (image->colormap == (VigraImpexColorPacket *) NULL)
    vigraImpexMagickError(ResourceLimitError,"Unable to quantize image",
      "Memory allocation failed");
  color_cube->colormap=image->colormap;
  color_cube->colors=0;
  DefineColormap(color_cube,color_cube->root);
  if ((quantize_info->number_colors == 2) &&
      (quantize_info->colorspace == VigraImpexGRAYColorspace))
    {
      unsigned int
        polarity;

      /*
        Monochrome image.
      */
      polarity=Intensity(image->colormap[0]) > Intensity(image->colormap[1]);
      image->colormap[polarity].red=0;
      image->colormap[polarity].green=0;
      image->colormap[polarity].blue=0;
      image->colormap[!polarity].red=MaxRGB;
      image->colormap[!polarity].green=MaxRGB;
      image->colormap[!polarity].blue=MaxRGB;
    }
  if (quantize_info->colorspace != VigraImpexTransparentColorspace)
    {
      image->matte=False;
      image->c_class=VigraImpexPseudoClass;
    }
  image->colors=(unsigned int) color_cube->colors;
  /*
    Create a reduced color image.
  */
  if (quantize_info->dither)
    quantize_info->dither=!DitherImage(color_cube,image);
  p=image->pixels;
  if (!quantize_info->dither)
    for (i=0; i < image->packets; i++)
    {
      /*
        Identify the deepest node containing the pixel's color.
      */
      node_info=color_cube->root;
      for (index=MaxTreeDepth-1; (int) index > 0; index--)
      {
        id=((DownScale(p->red) >> index) & 0x01) << 2 |
           ((DownScale(p->green) >> index) & 0x01) << 1 |
           ((DownScale(p->blue) >> index) & 0x01);
        if ((node_info->census & (1 << id)) == 0)
          break;
        node_info=node_info->child[id];
      }
      /*
        Find closest color among siblings and their children.
      */
      color_cube->color.red=p->red;
      color_cube->color.green=p->green;
      color_cube->color.blue=p->blue;
      color_cube->distance=3.0*(MaxRGB+1)*(MaxRGB+1);
      ClosestColor(color_cube,node_info->parent);
      index=color_cube->color_number;
      if (image->c_class == VigraImpexPseudoClass)
        p->index=index;
      else
        {
          p->red=image->colormap[index].red;
          p->green=image->colormap[index].green;
          p->blue=image->colormap[index].blue;
        }
      p++;
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   C l a s s i f i c a t i o n                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure Classification begins by initializing a color description tree
%  of sufficient depth to represent each possible input color in a leaf.
%  However, it is impractical to generate a fully-formed color
%  description tree in the classification phase for realistic values of
%  cmax.  If colors components in the input image are quantized to k-bit
%  precision, so that cmax= 2k-1, the tree would need k levels below the
%  root node to allow representing each possible input color in a leaf.
%  This becomes prohibitive because the tree's total number of nodes is
%  1 + sum(i=1,k,8k).
%
%  A complete tree would require 19,173,961 nodes for k = 8, cmax = 255.
%  Therefore, to avoid building a fully populated tree, QUANTIZE: (1)
%  Initializes data structures for nodes only as they are needed;  (2)
%  Chooses a maximum depth for the tree as a function of the desired
%  number of colors in the output image (currently log2(colormap size)).
%
%  For each pixel in the input image, classification scans downward from
%  the root of the color description tree.  At each level of the tree it
%  identifies the single node which represents a cube in RGB space
%  containing It updates the following data for each such node:
%
%    n1 : Number of pixels whose color is contained in the RGB cube
%    which this node represents;
%
%    n2 : Number of pixels whose color is not represented in a node at
%    lower depth in the tree;  initially,  n2 = 0 for all nodes except
%    leaves of the tree.
%
%    Sr, Sg, Sb : Sums of the red, green, and blue component values for
%    all pixels not classified at a lower depth. The combination of
%    these sums and n2  will ultimately characterize the mean color of a
%    set of pixels represented by this node.
%
%    E: The distance squared in RGB space between each pixel contained
%    within a node and the nodes' center.  This represents the quantization
%    error for a node.
%
%  The format of the Classification routine is:
%
%      Classification(color_cube,image)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o image: Specifies a pointer to an VigraImpexImage structure;  returned from
%      ReadImage.
%
%
*/
static void Classification(CubeInfo *color_cube,VigraImpexImage *image)
{
#define ClassifyImageText  "  Classifying image colors...  "

  double
    distance_squared;

  float
    mid_red,
    mid_green,
    mid_blue;

  register float
    distance;

  register int
    i;

  NodeInfo
    *node_info;

  register VigraImpexRunlengthPacket
    *p;

  register unsigned int
    index;

  static float
    bisect[] = { 128.0, 64.0, 32.0, 16.0, 8.0, 4.0, 2.0, 1.0, 0.5 };

  unsigned int
    id,
    level;

  color_cube->root->quantization_error=
    3.0*(MaxRGB/2.0)*(MaxRGB/2.0)*image->columns*image->rows;
  p=image->pixels;
  for (i=0; i < image->packets; i++)
  {
    if (color_cube->nodes > MaxNodes)
      {
        /*
          Prune one level if the color tree is too large.
        */
        PruneLevel(color_cube,color_cube->root);
        color_cube->depth--;
      }
    /*
      Start at the root and descend the color cube tree.
    */
    node_info=color_cube->root;
    index=MaxTreeDepth-1;
    mid_red=MaxRGB/2.0;
    mid_green=MaxRGB/2.0;
    mid_blue=MaxRGB/2.0;
    for (level=1; level <= color_cube->depth; level++)
    {
      id=((DownScale(p->red) >> index) & 0x01) << 2 |
         ((DownScale(p->green) >> index) & 0x01) << 1 |
         ((DownScale(p->blue) >> index) & 0x01);
      mid_red+=id & 4 ? bisect[level] : -bisect[level];
      mid_green+=id & 2 ? bisect[level] : -bisect[level];
      mid_blue+=id & 1 ? bisect[level] : -bisect[level];
      if (node_info->child[id] == (NodeInfo *) NULL)
        {
          /*
            Set colors of new node to contain pixel.
          */
          node_info->census|=1 << id;
          node_info->child[id]=
            vigraImpexInitializeNodeInfo(color_cube,id,level,node_info);
          if (node_info->child[id] == (NodeInfo *) NULL)
            vigraImpexMagickError(ResourceLimitError,"Unable to quantize image",
              "Memory allocation failed");
          if (level == color_cube->depth)
            color_cube->colors++;
        }
      node_info=node_info->child[id];
      if (level != 8)
        {
          /*
            Approximate the quantization error represented by this node.
          */
          distance=AbsoluteValue(DownScale(p->red)-mid_red);
          distance_squared=color_cube->squares[(int) distance]+distance-0.25;
          distance=AbsoluteValue(DownScale(p->green)-mid_green);
          distance_squared+=color_cube->squares[(int) distance]+distance-0.25;
          distance=AbsoluteValue(DownScale(p->blue)-mid_blue);
          distance_squared+=color_cube->squares[(int) distance]+distance-0.25;
          node_info->quantization_error+=distance_squared*(p->length+1);
        }
      index--;
    }
    /*
      Sum RGB values for this leaf for later derivation of the mean cube color.
    */
    node_info->number_unique+=(p->length+1);
    node_info->total_red+=p->red*(p->length+1);
    node_info->total_green+=p->green*(p->length+1);
    node_info->total_blue+=p->blue*(p->length+1);
    p++;
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   C l o s e s t C o l o r                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure ClosestColor traverses the color cube tree at a particular node
%  and determines which colormap entry best represents the input color.
%
%  The format of the ClosestColor routine is:
%
%      ClosestColor(color_cube,node_info)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o node_info: The address of a structure of type NodeInfo which points to a
%      node in the color cube tree that is to be pruned.
%
%
*/
static void ClosestColor(CubeInfo *color_cube,const NodeInfo *node_info)
{
  register unsigned int
    id;

  /*
    Traverse any children.
  */
  if (node_info->census != 0)
    for (id=0; id < 8; id++)
      if (node_info->census & (1 << id))
        ClosestColor(color_cube,node_info->child[id]);
  if (node_info->number_unique != 0)
    {
      register float
        distance_squared;

      register int
        distance;

      register VigraImpexColorPacket
        *color;

      /*
        Determine if this color is "closest".
      */
      color=color_cube->colormap+node_info->color_number;
      distance=(int) color->red-(int) color_cube->color.red;
      distance_squared=color_cube->squares[distance];
      distance=(int) color->green-(int) color_cube->color.green;
      distance_squared+=color_cube->squares[distance];
      distance=(int) color->blue-(int) color_cube->color.blue;
      distance_squared+=color_cube->squares[distance];
      if (distance_squared < color_cube->distance)
        {
          color_cube->distance=distance_squared;
          color_cube->color_number=(unsigned short) node_info->color_number;
        }
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   D e f i n e C o l o r m a p                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure DefineColormap traverses the color cube tree and notes each
%  colormap entry.  A colormap entry is any node in the color cube tree where
%  the number of unique colors is not zero.
%
%  The format of the DefineColormap routine is:
%
%      DefineColormap(color_cube,node_info)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o node_info: The address of a structure of type NodeInfo which points to a
%      node in the color cube tree that is to be pruned.
%
%
*/
static void DefineColormap(CubeInfo *color_cube,NodeInfo *node_info)
{
  register double
    number_unique;

  register unsigned int
    id;

  /*
    Traverse any children.
  */
  if (node_info->census != 0)
    for (id=0; id < 8; id++)
      if (node_info->census & (1 << id))
        DefineColormap(color_cube,node_info->child[id]);
  if (node_info->number_unique != 0)
    {
      /*
        Colormap entry is defined by the mean color in this color_cube.
      */
      number_unique=node_info->number_unique;
      color_cube->colormap[color_cube->colors].red=(Quantum)
        ((node_info->total_red+number_unique/2)/number_unique);
      color_cube->colormap[color_cube->colors].green=(Quantum)
        ((node_info->total_green+number_unique/2)/number_unique);
      color_cube->colormap[color_cube->colors].blue=(Quantum)
        ((node_info->total_blue+number_unique/2)/number_unique);
      node_info->color_number=(unsigned int) color_cube->colors++;
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   D i t h e r I m a g e                                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure DitherImage uses the Floyd-Steinberg algorithm to dither the
%  image.  Refer to "An Adaptive Algorithm for Spatial GreySscale", Robert W.
%  Floyd and Louis Steinberg, Proceedings of the S.I.D., Volume 17(2), 1976.
%
%  First find the closest representation to the reference pixel color in the
%  colormap, the node pixel is assigned this color.  Next, the colormap color
%  is subtracted from the reference pixels color, this represents the
%  quantization error.  Various amounts of this error are added to the pixels
%  ahead and below the node pixel to correct for this error.  The error
%  proportions are:
%
%            P     7/16
%      3/16  5/16  1/16
%
%  The error is distributed left-to-right for even scanlines and right-to-left
%  for odd scanlines.
%
%  The format of the DitherImage routine is:
%
%      DitherImage(color_cube,image)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o image: Specifies a pointer to an VigraImpexImage structure;  returned from
%      ReadImage.
%
%
*/
static unsigned int DitherImage(CubeInfo *color_cube,VigraImpexImage *image)
{
#define CacheShift  (QuantumDepth-6)
#define DitherImageText  "  Dithering image...  "


  typedef struct _ErrorPacket
  {
    int
      red,
      green,
      blue;
  } ErrorPacket;

  ErrorPacket
    *error;

  int
    blue_error,
    green_error,
    red_error,
    step;

  NodeInfo
    *node_info;

  Quantum
    blue,
    green,
    red;

  register ErrorPacket
    *cs,
    *ns;

  register int
    *cache,
    *error_limit,
    i;

  register Quantum
    *range_limit;

  register VigraImpexRunlengthPacket
    *q;

  register unsigned short
    index;

  unsigned int
    id,
    quantum,
    x,
    y;

  /*
    VigraImpexImage must be uncompressed.
  */
  if (!vigraImpexUncondenseImage(image))
    return(True);
  /*
    Allocate memory.
  */
  cache=(int *) malloc((1 << 18)*sizeof(int));
  error=(ErrorPacket *) malloc(((image->columns+2) << 1)*sizeof(ErrorPacket));
  error_limit=(int *) malloc((MaxRGB*2+1)*sizeof(int));
  range_limit=(Quantum *) malloc(3*(MaxRGB+1)*sizeof(Quantum));
  if ((cache == (int *) NULL) || (error == (ErrorPacket *) NULL) ||
      (error_limit == (int *) NULL) || (range_limit == (Quantum *) NULL))
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to dither image",
        "Memory allocation failed");
      return(True);
    }
  /*
    Initialize color cache.
  */
  for (i=0; i < (1 << 18); i++)
    cache[i]=(-1);
  /*
    Initialize error tables.
  */
  for (i=0; i < ((image->columns+2) << 1); i++)
  {
    error[i].red=0;
    error[i].green=0;
    error[i].blue=0;
  }
  /*
    Initialize error limit (constrain error).
  */
  quantum=(unsigned int) Max(image->colors >> 4,1);
  if ((quantum > 1) && (QuantumDepth != 8))
    quantum>>=1;
  error_limit+=MaxRGB;
  step=0;
  for (i=0; i < ((MaxRGB+1)/quantum); i++)
  {
    error_limit[i]=step;
    error_limit[-i]=(-step);
    step++;
  }
  if (quantum > 3)
    for ( ; i < (3*(MaxRGB+1)/quantum); i++)
    {
      error_limit[i]=step;
      error_limit[-i]=(-step);
      step+=(i & 0x01) ? 0 : 1;
    }
  for ( ; i <= MaxRGB; i++)
  {
    error_limit[i]=step;
    error_limit[-i]=(-step);
  }
  /*
    Initialize range tables.
  */
  for (i=0; i <= MaxRGB; i++)
  {
    range_limit[i]=0;
    range_limit[i+(MaxRGB+1)]=(Quantum) i;
    range_limit[i+(MaxRGB+1)*2]=MaxRGB;
  }
  range_limit+=(MaxRGB+1);
  /*
    Dither image.
  */
  for (y=0; y < image->rows; y++)
  {
    q=image->pixels+image->columns*y;
    cs=error+1;
    ns=cs+image->columns+2;
    step=y & 0x01 ? -1 : 1;
    if (step < 0)
      {
        /*
          Distribute error right-to-left for odd scanlines.
        */
        q+=(image->columns-1);
        ns=error+image->columns;
        cs=ns+image->columns+2;
      }
    for (x=0; x < image->columns; x++)
    {
      red_error=error_limit[(cs->red+8)/16];
      green_error=error_limit[(cs->green+8)/16];
      blue_error=error_limit[(cs->blue+8)/16];
      red=range_limit[(int) q->red+red_error];
      green=range_limit[(int) q->green+green_error];
      blue=range_limit[(int) q->blue+blue_error];
      i=(blue >> CacheShift) << 12 | (green >> CacheShift) << 6 |
        (red >> CacheShift);
      if (cache[i] < 0)
        {
          /*
            Identify the deepest node containing the pixel's color.
          */
          node_info=color_cube->root;
          for (index=MaxTreeDepth-1; (int) index > 0; index--)
          {
            id=((DownScale(red) >> index) & 0x01) << 2 |
               ((DownScale(green) >> index) & 0x01) << 1 |
               ((DownScale(blue) >> index) & 0x01);
            if ((node_info->census & (1 << id)) == 0)
              break;
            node_info=node_info->child[id];
          }
          /*
            Find closest color among siblings and their children.
          */
          color_cube->color.red=red;
          color_cube->color.green=green;
          color_cube->color.blue=blue;
          color_cube->distance=3.0*(MaxRGB+1)*(MaxRGB+1);
          ClosestColor(color_cube,node_info->parent);
          cache[i]=color_cube->color_number;
        }
      index=(unsigned short) cache[i];
      if (image->c_class == VigraImpexPseudoClass)
        q->index=index;
      else
        {
          q->red=image->colormap[index].red;
          q->green=image->colormap[index].green;
          q->blue=image->colormap[index].blue;
        }
      q+=step;
      /*
        Propagate the error in these proportions:
                Q     7/16
          3/16  5/16  1/16
      */
      red_error=(int) red-(int) image->colormap[index].red;
      green_error=(int) green-(int) image->colormap[index].green;
      blue_error=(int) blue-(int) image->colormap[index].blue;
      cs->red=0;
      cs->green=0;
      cs->blue=0;
      cs+=step;
      cs->red+=7*red_error;
      cs->green+=7*green_error;
      cs->blue+=7*blue_error;
      ns-=step;
      ns->red+=3*red_error;
      ns->green+=3*green_error;
      ns->blue+=3*blue_error;
      ns+=step;
      ns->red+=5*red_error;
      ns->green+=5*green_error;
      ns->blue+=5*blue_error;
      ns+=step;
      ns->red+=red_error;
      ns->green+=green_error;
      ns->blue+=blue_error;
    }
  }
  /*
    Free up memory.
  */
  range_limit-=(MaxRGB+1);
  free((char *) range_limit);
  error_limit-=MaxRGB;
  free((char *) error_limit);
  free((char *) error);
  free((char *) cache);
  return(False);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t Q u a n t i z e I n f o                                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexGetQuantizeInfo initializes the VigraImpexQuantizeInfo structure.
%
%  The format of the vigraImpexGetQuantizeInfo routine is:
%
%      vigraImpexGetQuantizeInfo(quantize_info)
%
%  A description of each parameter follows:
%
%    o quantize_info: Specifies a pointer to a MontageInfo structure.
%
%
*/
Export void vigraImpexGetQuantizeInfo(VigraImpexQuantizeInfo *quantize_info)
{
  assert(quantize_info != (VigraImpexQuantizeInfo *) NULL);
  quantize_info->number_colors=0;
  quantize_info->tree_depth=8;
  quantize_info->dither=False;
  quantize_info->colorspace=VigraImpexRGBColorspace;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   I n i t i a l i z e C u b e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function InitializeCube initialize the Cube data structure.
%
%  The format of the InitializeCube routine is:
%
%      InitializeCube(color_cube,depth)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o depth: Normally, this integer value is zero or one.  A zero or
%      one tells Quantize to choose a optimal tree depth of Log4(number_colors).
%      A tree of this depth generally allows the best representation of the
%      reference image with the least amount of memory and the fastest
%      computational speed.  In some cases, such as an image with low color
%      dispersion (a few number of colors), a value other than
%      Log4(number_colors) is required.  To expand the color tree completely,
%      use a value of 8.
%
%
*/
static void InitializeCube(CubeInfo *color_cube,int depth)
{
  register int
    i;

  /*
    Initialize tree to describe color color_cube.
  */
  color_cube->node_queue=(Nodes *) NULL;
  color_cube->nodes=0;
  color_cube->free_nodes=0;
  if (depth > MaxTreeDepth)
    depth=MaxTreeDepth;
  if (depth < 2)
    depth=2;
  color_cube->depth=depth;
  /*
    Initialize root node.
  */
  color_cube->root=vigraImpexInitializeNodeInfo(color_cube,0,0,(NodeInfo *) NULL);
  color_cube->squares=(unsigned int *)
    malloc((MaxRGB+MaxRGB+1)*sizeof(unsigned int));
  if ((color_cube->root == (NodeInfo *) NULL) ||
      (color_cube->squares == (unsigned int *) NULL))
    vigraImpexMagickError(ResourceLimitError,"Unable to quantize image",
      "Memory allocation failed");
  color_cube->root->parent=color_cube->root;
  color_cube->root->quantization_error=0.0;
  color_cube->colors=0;
  color_cube->squares+=MaxRGB;
  for (i=(-MaxRGB); i <= MaxRGB; i++)
    color_cube->squares[i]=i*i;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   I n i t i a l i z e N o d e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexInitializeNodeInfo allocates memory for a new node in the color cube
%  tree and presets all fields to zero.
%
%  The format of the vigraImpexInitializeNodeInfo routine is:
%
%      node=vigraImpexInitializeNodeInfo(color_cube,node,id,level)
%
%  A description of each parameter follows.
%
%    o node: The vigraImpexInitializeNodeInfo routine returns this integer address.
%
%    o id: Specifies the child number of the node.
%
%    o level: Specifies the level in the classification the node resides.
%
%
*/
static NodeInfo *vigraImpexInitializeNodeInfo(CubeInfo *color_cube,const unsigned int id,
  const unsigned int level,NodeInfo *parent)
{
  register int
    i;

  NodeInfo
    *node_info;

  if (color_cube->free_nodes == 0)
    {
      Nodes
        *nodes;

      /*
        Allocate a new nodes of nodes.
      */
      nodes=(Nodes *) malloc(sizeof(Nodes));
      if (nodes == (Nodes *) NULL)
        return((NodeInfo *) NULL);
      nodes->next=color_cube->node_queue;
      color_cube->node_queue=nodes;
      color_cube->next_node=nodes->nodes;
      color_cube->free_nodes=NodesInAList;
    }
  color_cube->nodes++;
  color_cube->free_nodes--;
  node_info=color_cube->next_node++;
  node_info->parent=parent;
  for (i=0; i < 8; i++)
    node_info->child[i]=(NodeInfo *) NULL;
  node_info->id=id;
  node_info->level=level;
  node_info->census=0;
  node_info->number_unique=0;
  node_info->quantization_error=0.0;
  node_info->total_red=0.0;
  node_info->total_green=0.0;
  node_info->total_blue=0.0;
  return(node_info);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   P r u n e C h i l d                                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function PruneChild deletes the given node and merges its statistics into
%  its parent.
%
%  The format of the PruneSubtree routine is:
%
%      PruneChild(color_cube,node_info)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o node_info: pointer to node in color cube tree that is to be pruned.
%
%
*/
static void PruneChild(CubeInfo *color_cube,const NodeInfo *node_info)
{
  NodeInfo
    *parent;

  register int
    id;

  /*
    Traverse any children.
  */
  if (node_info->census != 0)
    for (id=0; id < 8; id++)
      if (node_info->census & (1 << id))
        PruneChild(color_cube,node_info->child[id]);
  /*
    Merge color statistics into parent.
  */
  parent=node_info->parent;
  parent->census&=~(1 << node_info->id);
  parent->number_unique+=node_info->number_unique;
  parent->total_red+=node_info->total_red;
  parent->total_green+=node_info->total_green;
  parent->total_blue+=node_info->total_blue;
  color_cube->nodes--;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+  P r u n e L e v e l                                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure PruneLevel deletes all nodes at the bottom level of the color
%  tree merging their color statistics into their parent node.
%
%  The format of the PruneLevel routine is:
%
%      PruneLevel(color_cube,node_info)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o node_info: pointer to node in color cube tree that is to be pruned.
%
%
*/
static void PruneLevel(CubeInfo *color_cube,const NodeInfo *node_info)
{
  register int
    id;

  /*
    Traverse any children.
  */
  if (node_info->census != 0)
    for (id=0; id < 8; id++)
      if (node_info->census & (1 << id))
        PruneLevel(color_cube,node_info->child[id]);
  if (node_info->level == color_cube->depth)
    PruneChild(color_cube,node_info);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e d u c e                                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function Reduce traverses the color cube tree and prunes any node whose
%  quantization error falls below a particular threshold.
%
%  The format of the Reduce routine is:
%
%      Reduce(color_cube,node_info)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o node_info: pointer to node in color cube tree that is to be pruned.
%
%
*/
static void Reduce(CubeInfo *color_cube,const NodeInfo *node_info)
{
  register unsigned int
    id;

  /*
    Traverse any children.
  */
  if (node_info->census != 0)
    for (id=0; id < 8; id++)
      if (node_info->census & (1 << id))
        Reduce(color_cube,node_info->child[id]);
  if (node_info->quantization_error <= color_cube->pruning_threshold)
    PruneChild(color_cube,node_info);
  else
    {
      /*
        Find minimum pruning threshold.
      */
      if (node_info->number_unique > 0)
        color_cube->colors++;
      if (node_info->quantization_error < color_cube->next_pruning_threshold)
        color_cube->next_pruning_threshold=node_info->quantization_error;
    }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e d u c t i o n                                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function Reduction repeatedly prunes the tree until the number of nodes
%  with n2 > 0 is less than or equal to the maximum number of colors allowed
%  in the output image.  On any given iteration over the tree, it selects
%  those nodes whose E value is minimal for pruning and merges their
%  color statistics upward. It uses a pruning threshold, Ep, to govern
%  node selection as follows:
%
%    Ep = 0
%    while number of nodes with (n2 > 0) > required maximum number of colors
%      prune all nodes such that E <= Ep
%      Set Ep to minimum E in remaining nodes
%
%  This has the effect of minimizing any quantization error when merging
%  two nodes together.
%
%  When a node to be pruned has offspring, the pruning procedure invokes
%  itself recursively in order to prune the tree from the leaves upward.
%  n2,  Sr, Sg,  and  Sb in a node being pruned are always added to the
%  corresponding data in that node's parent.  This retains the pruned
%  node's color characteristics for later averaging.
%
%  For each node, n2 pixels exist for which that node represents the
%  smallest volume in RGB space containing those pixel's colors.  When n2
%  > 0 the node will uniquely define a color in the output image. At the
%  beginning of reduction,  n2 = 0  for all nodes except a the leaves of
%  the tree which represent colors present in the input image.
%
%  The other pixel count, n1, indicates the total number of colors
%  within the cubic volume which the node represents.  This includes n1 -
%  n2  pixels whose colors should be defined by nodes at a lower level in
%  the tree.
%
%  The format of the Reduction routine is:
%
%      Reduction(color_cube,number_colors)
%
%  A description of each parameter follows.
%
%    o cube: A pointer to the Cube structure.
%
%    o number_colors: This integer value indicates the maximum number of
%      colors in the quantized image or colormap.  The actual number of
%      colors allocated to the colormap may be less than this value, but
%      never more.
%
%
*/
static void Reduction(CubeInfo *color_cube,const unsigned int number_colors)
{
#define ReduceImageText  "  Reducing image colors...  "

  unsigned int
    span;

  span=(unsigned int) color_cube->colors;
  color_cube->next_pruning_threshold=0.0;
  while (color_cube->colors > number_colors)
  {
    color_cube->pruning_threshold=color_cube->next_pruning_threshold;
    color_cube->next_pruning_threshold=color_cube->root->quantization_error;
    color_cube->colors=0;
    Reduce(color_cube,color_cube->root);
 }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   M a p I m a g e                                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  vigraImpexMapImage replaces the colors of an image with the closest color from
%  a reference image.
%
%  The format of the vigraImpexMapImage routine is:
%
%      vigraImpexMapImage(image,map_image,dither)
%
%  A description of each parameter follows:
%
%    o image: Specifies a pointer to an VigraImpexImage structure.
%
%    o map_image: Specifies a pointer to a VigraImpexImage structure.  Reduce
%      image to a set of colors represented by this image.
%
%    o dither: Set this integer value to something other than zero to
%      dither the quantized image.
%
%
*/
void vigraImpexMapImage(VigraImpexImage *image,VigraImpexImage *map_image,const unsigned int dither)
{
  CubeInfo
    color_cube;

  Nodes
    *nodes;

  VigraImpexQuantizeInfo
    quantize_info;

  assert(image != (VigraImpexImage *) NULL);
  if (map_image == (VigraImpexImage *) NULL)
    return;
  /*
    Classify image colors from the reference image.
  */
  InitializeCube(&color_cube,8);
  Classification(&color_cube,map_image);
  quantize_info.number_colors=color_cube.colors;
  quantize_info.dither=dither;
  quantize_info.colorspace=image->matte ? VigraImpexTransparentColorspace : VigraImpexRGBColorspace;
  Assignment(&color_cube,&quantize_info,image);
  /*
    Release color cube tree storage.
  */
  do
  {
    nodes=color_cube.node_queue->next;
    free((char *) color_cube.node_queue);
    color_cube.node_queue=nodes;
  } while (color_cube.node_queue != (Nodes *) NULL);
  color_cube.squares-=MaxRGB;
  free((char *) color_cube.squares);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   O r d e r e d D i t h e r I m a g e                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure OrderedDitherImage uses the ordered dithering technique of
%  reducing color images to monochrome using positional information to retain
%  as much information as possible.
%
%  The format of the OrderedDitherImage routine is:
%
%      OrderedDitherImage(image)
%
%  A description of each parameter follows.
%
%    o image: Specifies a pointer to an VigraImpexImage structure;  returned from
%      ReadImage.
%
%
*/
static void OrderedDitherImage(VigraImpexImage *image)
{
  static Quantum
    DitherMatrix[8][8] =
    {
      { UpScale(  0), UpScale(192), UpScale( 48), UpScale(240),
        UpScale( 12), UpScale(204), UpScale( 60), UpScale(252) },
      { UpScale(128), UpScale( 64), UpScale(176), UpScale(112),
        UpScale(140), UpScale( 76), UpScale(188), UpScale(124) },
      { UpScale( 32), UpScale(224), UpScale( 16), UpScale(208),
        UpScale( 44), UpScale(236), UpScale( 28), UpScale(220) },
      { UpScale(160), UpScale( 96), UpScale(144), UpScale( 80),
        UpScale(172), UpScale(108), UpScale(156), UpScale( 92) },
      { UpScale(  8), UpScale(200), UpScale( 56), UpScale(248),
        UpScale(  4), UpScale(196), UpScale( 52), UpScale(244) },
      { UpScale(136), UpScale( 72), UpScale(184), UpScale(120),
        UpScale(132), UpScale( 68), UpScale(180), UpScale(116) },
      { UpScale( 40), UpScale(232), UpScale( 24), UpScale(216),
        UpScale( 36), UpScale(228), UpScale( 20), UpScale(212) },
      { UpScale(168), UpScale(104), UpScale(152), UpScale( 88),
        UpScale(164), UpScale(100), UpScale(148), UpScale( 84) }
    };

  register int
    x,
    y;

  register VigraImpexRunlengthPacket
    *p;

  /*
    Transform image to uncompressed normalized grayscale.
  */
  vigraImpexRGBTransformImage(image,VigraImpexGRAYColorspace);
  vigraImpexNormalizeImage(image);
  if (!vigraImpexUncondenseImage(image))
    return;
  /*
    Initialize colormap.
  */
  image->c_class=VigraImpexPseudoClass;
  if (image->colormap != (VigraImpexColorPacket *) NULL)
    free((char *) image->colormap);
  image->colors=2;
  image->colormap=(VigraImpexColorPacket *) malloc(image->colors*sizeof(VigraImpexColorPacket));
  if (image->colormap == (VigraImpexColorPacket *) NULL)
    vigraImpexMagickError(ResourceLimitError,"Unable to quantize image",
      "Memory allocation failed");
  image->colormap[0].red=0;
  image->colormap[0].green=0;
  image->colormap[0].blue=0;
  image->colormap[1].red=MaxRGB;
  image->colormap[1].green=MaxRGB;
  image->colormap[1].blue=MaxRGB;
  /*
    Dither image with the ordered dithering technique.
  */
  p=image->pixels;
  for (y=0; y < image->rows; y++)
  {
    for (x=0; x < image->columns; x++)
    {
      p->index=p->red > DitherMatrix[y & 0x07][x & 0x07] ? 1 : 0;
      p++;
    }
  }
  vigraImpexSyncImage(image);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  Q u a n t i z a t i o n E r r o r                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexQuantizationError measures the difference between the original
%  and quantized images.  This difference is the total quantization error.
%  The error is computed by summing over all pixels in an image the distance
%  squared in RGB space between each reference pixel value and its quantized
%  value.  These values are computed:
%
%    o mean_error_per_pixel:  This value is the mean error for any single
%      pixel in the image.
%
%    o normalized_mean_square_error:  This value is the normalized mean
%      quantization error for any single pixel in the image.  This distance
%      measure is normalized to a range between 0 and 1.  It is independent
%      of the range of red, green, and blue values in the image.
%
%    o normalized_maximum_square_error:  Thsi value is the normalized
%      maximum quantization error for any single pixel in the image.  This
%      distance measure is normalized to a range between 0 and 1.  It is
%      independent of the range of red, green, and blue values in your image.
%
%
%  The format of the vigraImpexQuantizationError routine is:
%
%      vigraImpexQuantizationError(image)
%
%  A description of each parameter follows.
%
%    o image: The address of a byte (8 bits) array of run-length
%      encoded pixel data of your reference image.  The sum of the
%      run-length counts in the reference image must be equal to or exceed
%      the number of pixels.
%
%
*/
void vigraImpexQuantizationError(VigraImpexImage *image)
{
  CubeInfo
    color_cube;

  double
    total_error;

  float
    distance_squared,
    maximum_error_per_pixel;

  int
    distance;

  register int
    i;

  register VigraImpexRunlengthPacket
    *p;

  /*
    Initialize measurement.
  */
  assert(image != (VigraImpexImage *) NULL);
  image->mean_error_per_pixel=0;
  image->normalized_mean_error=0.0;
  image->normalized_maximum_error=0.0;
  vigraImpexNumberColors(image,(FILE *) NULL);
  if (image->c_class == VigraImpexDirectClass)
    return;
  color_cube.squares=(unsigned int *)
    malloc((MaxRGB+MaxRGB+1)*sizeof(unsigned int));
  if (color_cube.squares == (unsigned int *) NULL)
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to measure error",
        "Memory allocation failed");
      return;
    }
  color_cube.squares+=MaxRGB;
  for (i=(-MaxRGB); i <= MaxRGB; i++)
    color_cube.squares[i]=i*i;
  /*
    For each pixel, collect error statistics.
  */
  maximum_error_per_pixel=0;
  total_error=0;
  p=image->pixels;
  for (i=0; i < image->packets; i++)
  {
    distance=(int) p->red-(int) image->colormap[p->index].red;
    distance_squared=color_cube.squares[distance];
    distance=(int) p->green-(int) image->colormap[p->index].green;
    distance_squared+=color_cube.squares[distance];
    distance=(int) p->blue-(int) image->colormap[p->index].blue;
    distance_squared+=color_cube.squares[distance];
    total_error+=(distance_squared*((double) p->length+1.0));
    if (distance_squared > maximum_error_per_pixel)
      maximum_error_per_pixel=distance_squared;
    p++;
  }
  /*
    Compute final error statistics.
  */
  image->mean_error_per_pixel=(unsigned int)
    total_error/(image->columns*image->rows);
  image->normalized_mean_error=
    (image->mean_error_per_pixel)/(3.0*(MaxRGB+1)*(MaxRGB+1));
  image->normalized_maximum_error=
    maximum_error_per_pixel/(3.0*(MaxRGB+1)*(MaxRGB+1));
  color_cube.squares-=MaxRGB;
  free((char *) color_cube.squares);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  Q u a n t i z e I m a g e                                                  %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexQuantizeImage analyzes the colors within a reference image and
%  chooses a fixed number of colors to represent the image.  The goal of the
%  algorithm is to minimize the difference between the input and output image
%  while minimizing the processing time.
%
%  The format of the vigraImpexQuantizeImage routine is:
%
%      vigraImpexQuantizeImage(quantize_info,image)
%
%  A description of each parameter follows:
%
%    o quantize_info: Specifies a pointer to an VigraImpexQuantizeInfo structure.
%
%    o image: Specifies a pointer to a VigraImpexImage structure.
%
*/
Export void vigraImpexQuantizeImage(VigraImpexQuantizeInfo *quantize_info,VigraImpexImage *image)
{
  CubeInfo
    color_cube;

  int
    depth;

  Nodes
    *nodes;

  unsigned int
    number_colors;

  assert(image != (VigraImpexImage *) NULL);
  if ((quantize_info->number_colors == 2) && quantize_info->dither &&
      (quantize_info->colorspace == VigraImpexGRAYColorspace))
    {
      OrderedDitherImage(image);
      return;
    }
  number_colors=quantize_info->number_colors;
  if (number_colors == 0)
    number_colors=MaxRGB+1;
  if (number_colors > MaxColormapSize)
    number_colors=MaxColormapSize;
  if (image->packets == (image->columns*image->rows))
    vigraImpexCondenseImage(image);
  depth=quantize_info->tree_depth;
  if (depth == 0)
    {
      unsigned int
        colors;

      /*
        Depth of color classification tree is: Log4(colormap size)+2.
      */
      colors=number_colors;
      for (depth=1; colors != 0; depth++)
        colors>>=2;
      if (quantize_info->dither)
        depth--;
      if (image->c_class == VigraImpexPseudoClass)
        depth+=2;
    }
  /*
    Reduce the number of colors in the continuous tone image.
  */
  InitializeCube(&color_cube,depth);
  if (quantize_info->colorspace != VigraImpexRGBColorspace)
    vigraImpexRGBTransformImage(image,quantize_info->colorspace);
  Classification(&color_cube,image);
  if ((color_cube.colors >> 1) < number_colors)
    quantize_info->dither=False;
  if (number_colors < color_cube.colors)
    Reduction(&color_cube,number_colors);
  Assignment(&color_cube,quantize_info,image);
  if (quantize_info->colorspace != VigraImpexRGBColorspace)
    vigraImpexTransformRGBImage(image,quantize_info->colorspace);
  /*
    Release color cube tree storage.
  */
  do
  {
    nodes=color_cube.node_queue->next;
    free((char *) color_cube.node_queue);
    color_cube.node_queue=nodes;
  } while (color_cube.node_queue != (Nodes *) NULL);
  color_cube.squares-=MaxRGB;
  free((char *) color_cube.squares);
}
