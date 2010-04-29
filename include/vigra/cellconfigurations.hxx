/************************************************************************/
/*                                                                      */
/*             Copyright 2009-2010 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_CELLCONFIGURATIONS_HXX
#define VIGRA_CELLCONFIGURATIONS_HXX

#include "cellimage.hxx"

namespace vigra {

namespace cellimage {

static CellType cellConfigurations[] = {
  /*     */
  /*  #  */
  /*     */
CellTypeVertex, /* 0 */

  /*     */
  /*  ## */
  /*     */
CellTypeVertex, /* 1 */

  /*   # */
  /*  #  */
  /*     */
CellTypeVertex, /* 2 */

  /*   # */
  /*  ## */
  /*     */
CellTypeVertex, /* 3 */

  /*  #  */
  /*  #  */
  /*     */
CellTypeVertex, /* 4 */

  /*  #  */
  /*  ## */
  /*     */
CellTypeVertexOrLine, /* 5 */

  /*  ## */
  /*  #  */
  /*     */
CellTypeVertex, /* 6 */

  /*  ## */
  /*  ## */
  /*     */
CellTypeVertex, /* 7 */

  /* #   */
  /*  #  */
  /*     */
CellTypeVertex, /* 8 */

  /* #   */
  /*  ## */
  /*     */
CellTypeLine,   /* 9 */

  /* # # */
  /*  #  */
  /*     */
CellTypeLine,   /* 10 */

  /* # # */
  /*  ## */
  /*     */
CellTypeLine,   /* 11 */

  /* ##  */
  /*  #  */
  /*     */
CellTypeVertex, /* 12 */

  /* ##  */
  /*  ## */
  /*     */
CellTypeVertexOrLine, /* 13 */

  /* ### */
  /*  #  */
  /*     */
CellTypeVertex,  /* 14 */

  /* ### */
  /*  ## */
  /*     */
CellTypeError,  /* 15 */

  /*     */
  /* ##  */
  /*     */
CellTypeVertex, /* 16 */

  /*     */
  /* ### */
  /*     */
CellTypeLine,   /* 17 */

  /*   # */
  /* ##  */
  /*     */
CellTypeLine,   /* 18 */

  /*   # */
  /* ### */
  /*     */
CellTypeLine,   /* 19 */

  /*  #  */
  /* ##  */
  /*     */
CellTypeVertexOrLine, /* 20 */

  /*  #  */
  /* ### */
  /*     */
CellTypeVertex, /* 21 */

  /*  ## */
  /* ##  */
  /*     */
CellTypeVertexOrLine, /* 22 */

  /*  ## */
  /* ### */
  /*     */
CellTypeVertex, /* 23 */

  /* #   */
  /* ##  */
  /*     */
CellTypeVertex, /* 24 */

  /* #   */
  /* ### */
  /*     */
CellTypeLine,   /* 25 */

  /* # # */
  /* ##  */
  /*     */
CellTypeLine,   /* 26 */

  /* # # */
  /* ### */
  /*     */
CellTypeLine,   /* 27 */

  /* ##  */
  /* ##  */
  /*     */
CellTypeVertex, /* 28 */

  /* ##  */
  /* ### */
  /*     */
CellTypeVertex, /* 29 */

  /* ### */
  /* ##  */
  /*     */
CellTypeError,  /* 30 */

  /* ### */
  /* ### */
  /*     */
CellTypeVertex, /* 31 */

  /*     */
  /*  #  */
  /* #   */
CellTypeVertex, /* 32 */

  /*     */
  /*  ## */
  /* #   */
CellTypeLine,   /* 33 */

  /*   # */
  /*  #  */
  /* #   */
CellTypeLine,   /* 34 */

  /*   # */
  /*  ## */
  /* #   */
CellTypeLine,   /* 35 */

  /*  #  */
  /*  #  */
  /* #   */
CellTypeLine,   /* 36 */

  /*  #  */
  /*  ## */
  /* #   */
CellTypeVertex, /* 37 */

  /*  ## */
  /*  #  */
  /* #   */
CellTypeLine,   /* 38 */

  /*  ## */
  /*  ## */
  /* #   */
CellTypeVertex, /* 39 */

  /* #   */
  /*  #  */
  /* #   */
CellTypeLine,   /* 40 */

  /* #   */
  /*  ## */
  /* #   */
CellTypeVertex, /* 41 */

  /* # # */
  /*  #  */
  /* #   */
CellTypeVertex, /* 42 */

  /* # # */
  /*  ## */
  /* #   */
CellTypeVertex, /* 43 */

  /* ##  */
  /*  #  */
  /* #   */
CellTypeLine,   /* 44 */

  /* ##  */
  /*  ## */
  /* #   */
CellTypeVertex, /* 45 */

  /* ### */
  /*  #  */
  /* #   */
CellTypeLine,   /* 46 */

  /* ### */
  /*  ## */
  /* #   */
CellTypeVertex, /* 47 */

  /*     */
  /* ##  */
  /* #   */
CellTypeVertex, /* 48 */

  /*     */
  /* ### */
  /* #   */
CellTypeLine,   /* 49 */

  /*   # */
  /* ##  */
  /* #   */
CellTypeLine,   /* 50 */

  /*   # */
  /* ### */
  /* #   */
CellTypeLine,   /* 51 */

  /*  #  */
  /* ##  */
  /* #   */
CellTypeVertexOrLine, /* 52 */

  /*  #  */
  /* ### */
  /* #   */
CellTypeVertex, /* 53 */

  /*  ## */
  /* ##  */
  /* #   */
CellTypeErrorOrLine,  /* 54 */

  /*  ## */
  /* ### */
  /* #   */
CellTypeError,  /* 55 */

  /* #   */
  /* ##  */
  /* #   */
CellTypeVertex,  /* 56 */

  /* #   */
  /* ### */
  /* #   */
CellTypeLine,   /* 57 */

  /* # # */
  /* ##  */
  /* #   */
CellTypeLine,   /* 58 */

  /* # # */
  /* ### */
  /* #   */
CellTypeLine,   /* 59 */

  /* ##  */
  /* ##  */
  /* #   */
CellTypeError,  /* 60 */

  /* ##  */
  /* ### */
  /* #   */
CellTypeVertex, /* 61 */

  /* ### */
  /* ##  */
  /* #   */
CellTypeError,  /* 62 */

  /* ### */
  /* ### */
  /* #   */
CellTypeError,  /* 63 */

  /*     */
  /*  #  */
  /*  #  */
CellTypeVertex, /* 64 */

  /*     */
  /*  ## */
  /*  #  */
CellTypeVertexOrLine, /* 65 */

  /*   # */
  /*  #  */
  /*  #  */
CellTypeLine,   /* 66 */

  /*   # */
  /*  ## */
  /*  #  */
CellTypeVertexOrLine, /* 67 */

  /*  #  */
  /*  #  */
  /*  #  */
CellTypeLine,   /* 68 */

  /*  #  */
  /*  ## */
  /*  #  */
CellTypeVertex, /* 69 */

  /*  ## */
  /*  #  */
  /*  #  */
CellTypeLine,   /* 70 */

  /*  ## */
  /*  ## */
  /*  #  */
CellTypeVertex, /* 71 */

  /* #   */
  /*  #  */
  /*  #  */
CellTypeLine,   /* 72 */

  /* #   */
  /*  ## */
  /*  #  */
CellTypeVertex, /* 73 */

  /* # # */
  /*  #  */
  /*  #  */
CellTypeVertex, /* 74 */

  /* # # */
  /*  ## */
  /*  #  */
CellTypeVertex, /* 75 */

  /* ##  */
  /*  #  */
  /*  #  */
CellTypeLine,   /* 76 */

  /* ##  */
  /*  ## */
  /*  #  */
CellTypeVertex, /* 77 */

  /* ### */
  /*  #  */
  /*  #  */
CellTypeLine,   /* 78 */

  /* ### */
  /*  ## */
  /*  #  */
CellTypeVertex, /* 79 */

  /*     */
  /* ##  */
  /*  #  */
CellTypeVertexOrLine, /* 80 */

  /*     */
  /* ### */
  /*  #  */
CellTypeVertex, /* 81 */

  /*   # */
  /* ##  */
  /*  #  */
CellTypeVertex, /* 82 */

  /*   # */
  /* ### */
  /*  #  */
CellTypeVertex, /* 83 */

  /*  #  */
  /* ##  */
  /*  #  */
CellTypeVertex, /* 84 */

  /*  #  */
  /* ### */
  /*  #  */
CellTypeVertex, /* 85 */

  /*  ## */
  /* ##  */
  /*  #  */
CellTypeVertex, /* 86 */

  /*  ## */
  /* ### */
  /*  #  */
CellTypeVertex, /* 87 */

  /* #   */
  /* ##  */
  /*  #  */
CellTypeVertexOrLine, /* 88 */

  /* #   */
  /* ### */
  /*  #  */
CellTypeVertex, /* 89 */

  /* # # */
  /* ##  */
  /*  #  */
CellTypeVertex, /* 90 */

  /* # # */
  /* ### */
  /*  #  */
CellTypeVertex,  /* 91 */

  /* ##  */
  /* ##  */
  /*  #  */
CellTypeVertex, /* 92 */

  /* ##  */
  /* ### */
  /*  #  */
CellTypeVertex, /* 93 */

  /* ### */
  /* ##  */
  /*  #  */
CellTypeVertex, /* 94 */

  /* ### */
  /* ### */
  /*  #  */
CellTypeVertex, /* 95 */

  /*     */
  /*  #  */
  /* ##  */
CellTypeVertex, /* 96 */

  /*     */
  /*  ## */
  /* ##  */
CellTypeVertexOrLine, /* 97 */

  /*   # */
  /*  #  */
  /* ##  */
CellTypeLine,   /* 98 */

  /*   # */
  /*  ## */
  /* ##  */
CellTypeErrorOrLine,  /* 99 */

  /*  #  */
  /*  #  */
  /* ##  */
CellTypeLine,   /* 100 */

  /*  #  */
  /*  ## */
  /* ##  */
CellTypeVertex, /* 101 */

  /*  ## */
  /*  #  */
  /* ##  */
CellTypeLine,   /* 102 */

  /*  ## */
  /*  ## */
  /* ##  */
CellTypeError,  /* 103 */

  /* #   */
  /*  #  */
  /* ##  */
CellTypeLine,   /* 104 */

  /* #   */
  /*  ## */
  /* ##  */
CellTypeVertex, /* 105 */

  /* # # */
  /*  #  */
  /* ##  */
CellTypeVertex, /* 106 */

  /* # # */
  /*  ## */
  /* ##  */
CellTypeVertex, /* 107 */

  /* ##  */
  /*  #  */
  /* ##  */
CellTypeLine,   /* 108 */

  /* ##  */
  /*  ## */
  /* ##  */
CellTypeVertex,  /* 109 */

  /* ### */
  /*  #  */
  /* ##  */
CellTypeLine,   /* 110 */

  /* ### */
  /*  ## */
  /* ##  */
CellTypeError,  /* 111 */

  /*     */
  /* ##  */
  /* ##  */
CellTypeVertex, /* 112 */

  /*     */
  /* ### */
  /* ##  */
CellTypeVertex, /* 113 */

  /*   # */
  /* ##  */
  /* ##  */
CellTypeVertex, /* 114 */

  /*   # */
  /* ### */
  /* ##  */
CellTypeError,  /* 115 */

  /*  #  */
  /* ##  */
  /* ##  */
CellTypeVertex, /* 116 */

  /*  #  */
  /* ### */
  /* ##  */
CellTypeVertex, /* 117 */

  /*  ## */
  /* ##  */
  /* ##  */
CellTypeError,  /* 118 */

  /*  ## */
  /* ### */
  /* ##  */
CellTypeVertex, /* 119 */

  /* #   */
  /* ##  */
  /* ##  */
CellTypeError,  /* 120 */

  /* #   */
  /* ### */
  /* ##  */
CellTypeVertex, /* 121 */

  /* # # */
  /* ##  */
  /* ##  */
CellTypeVertex, /* 122 */

  /* # # */
  /* ### */
  /* ##  */
CellTypeError,  /* 123 */

  /* ##  */
  /* ##  */
  /* ##  */
CellTypeVertex, /* 124 */

  /* ##  */
  /* ### */
  /* ##  */
CellTypeVertex, /* 125 */

  /* ### */
  /* ##  */
  /* ##  */
CellTypeError,  /* 126 */

  /* ### */
  /* ### */
  /* ##  */
CellTypeVertex, /* 127 */

  /*     */
  /*  #  */
  /*   # */
CellTypeVertex, /* 128 */

  /*     */
  /*  ## */
  /*   # */
CellTypeVertex, /* 129 */

  /*   # */
  /*  #  */
  /*   # */
CellTypeLine,   /* 130 */

  /*   # */
  /*  ## */
  /*   # */
CellTypeVertex,  /* 131 */

  /*  #  */
  /*  #  */
  /*   # */
CellTypeLine,   /* 132 */

  /*  #  */
  /*  ## */
  /*   # */
CellTypeVertexOrLine, /* 133 */

  /*  ## */
  /*  #  */
  /*   # */
CellTypeLine,   /* 134 */

  /*  ## */
  /*  ## */
  /*   # */
CellTypeError,  /* 135 */

  /* #   */
  /*  #  */
  /*   # */
CellTypeLine,   /* 136 */

  /* #   */
  /*  ## */
  /*   # */
CellTypeLine,   /* 137 */

  /* # # */
  /*  #  */
  /*   # */
CellTypeVertex, /* 138 */

  /* # # */
  /*  ## */
  /*   # */
CellTypeLine,   /* 139 */

  /* ##  */
  /*  #  */
  /*   # */
CellTypeLine,   /* 140 */

  /* ##  */
  /*  ## */
  /*   # */
CellTypeErrorOrLine,  /* 141 */

  /* ### */
  /*  #  */
  /*   # */
CellTypeLine,   /* 142 */

  /* ### */
  /*  ## */
  /*   # */
CellTypeError,  /* 143 */

  /*     */
  /* ##  */
  /*   # */
CellTypeLine,   /* 144 */

  /*     */
  /* ### */
  /*   # */
CellTypeLine,   /* 145 */

  /*   # */
  /* ##  */
  /*   # */
CellTypeVertex, /* 146 */

  /*   # */
  /* ### */
  /*   # */
CellTypeLine,   /* 147 */

  /*  #  */
  /* ##  */
  /*   # */
CellTypeVertex, /* 148 */

  /*  #  */
  /* ### */
  /*   # */
CellTypeVertex, /* 149 */

  /*  ## */
  /* ##  */
  /*   # */
CellTypeVertex, /* 150 */

  /*  ## */
  /* ### */
  /*   # */
CellTypeVertex, /* 151 */

  /* #   */
  /* ##  */
  /*   # */
CellTypeLine,   /* 152 */

  /* #   */
  /* ### */
  /*   # */
CellTypeLine,   /* 153 */

  /* # # */
  /* ##  */
  /*   # */
CellTypeVertex, /* 154 */

  /* # # */
  /* ### */
  /*   # */
CellTypeLine,   /* 155 */

  /* ##  */
  /* ##  */
  /*   # */
CellTypeVertex, /* 156 */

  /* ##  */
  /* ### */
  /*   # */
CellTypeError,  /* 157 */

  /* ### */
  /* ##  */
  /*   # */
CellTypeVertex, /* 158 */

  /* ### */
  /* ### */
  /*   # */
CellTypeError,  /* 159 */

  /*     */
  /*  #  */
  /* # # */
CellTypeLine,   /* 160 */

  /*     */
  /*  ## */
  /* # # */
CellTypeLine,   /* 161 */

  /*   # */
  /*  #  */
  /* # # */
CellTypeVertex, /* 162 */

  /*   # */
  /*  ## */
  /* # # */
CellTypeLine,   /* 163 */

  /*  #  */
  /*  #  */
  /* # # */
CellTypeVertex, /* 164 */

  /*  #  */
  /*  ## */
  /* # # */
CellTypeVertex, /* 165 */

  /*  ## */
  /*  #  */
  /* # # */
CellTypeVertex, /* 166 */

  /*  ## */
  /*  ## */
  /* # # */
CellTypeVertex, /* 167 */

  /* #   */
  /*  #  */
  /* # # */
CellTypeVertex, /* 168 */

  /* #   */
  /*  ## */
  /* # # */
CellTypeVertex, /* 169 */

  /* # # */
  /*  #  */
  /* # # */
CellTypeVertex, /* 170 */

  /* # # */
  /*  ## */
  /* # # */
CellTypeVertex, /* 171 */

  /* ##  */
  /*  #  */
  /* # # */
CellTypeVertex, /* 172 */

  /* ##  */
  /*  ## */
  /* # # */
CellTypeVertex, /* 173 */

  /* ### */
  /*  #  */
  /* # # */
CellTypeVertex, /* 174 */

  /* ### */
  /*  ## */
  /* # # */
CellTypeVertex, /* 175 */

  /*     */
  /* ##  */
  /* # # */
CellTypeLine,   /* 176 */

  /*     */
  /* ### */
  /* # # */
CellTypeLine,   /* 177 */

  /*   # */
  /* ##  */
  /* # # */
CellTypeVertex, /* 178 */

  /*   # */
  /* ### */
  /* # # */
CellTypeLine,   /* 179 */

  /*  #  */
  /* ##  */
  /* # # */
CellTypeVertex, /* 180 */

  /*  #  */
  /* ### */
  /* # # */
CellTypeVertex,  /* 181 */

  /*  ## */
  /* ##  */
  /* # # */
CellTypeVertex, /* 182 */

  /*  ## */
  /* ### */
  /* # # */
CellTypeError,  /* 183 */

  /* #   */
  /* ##  */
  /* # # */
CellTypeLine,   /* 184 */

  /* #   */
  /* ### */
  /* # # */
CellTypeLine,   /* 185 */

  /* # # */
  /* ##  */
  /* # # */
CellTypeVertex, /* 186 */

  /* # # */
  /* ### */
  /* # # */
CellTypeLine,   /* 187 */

  /* ##  */
  /* ##  */
  /* # # */
CellTypeVertex, /* 188 */

  /* ##  */
  /* ### */
  /* # # */
CellTypeError,  /* 189 */

  /* ### */
  /* ##  */
  /* # # */
CellTypeVertex, /* 190 */

  /* ### */
  /* ### */
  /* # # */
CellTypeError,  /* 191 */

  /*     */
  /*  #  */
  /*  ## */
CellTypeVertex, /* 192 */

  /*     */
  /*  ## */
  /*  ## */
CellTypeVertex, /* 193 */

  /*   # */
  /*  #  */
  /*  ## */
CellTypeLine,   /* 194 */

  /*   # */
  /*  ## */
  /*  ## */
CellTypeError,  /* 195 */

  /*  #  */
  /*  #  */
  /*  ## */
CellTypeLine,   /* 196 */

  /*  #  */
  /*  ## */
  /*  ## */
CellTypeVertex, /* 197 */

  /*  ## */
  /*  #  */
  /*  ## */
CellTypeLine,   /* 198 */

  /*  ## */
  /*  ## */
  /*  ## */
CellTypeVertex, /* 199 */

  /* #   */
  /*  #  */
  /*  ## */
CellTypeLine,   /* 200 */

  /* #   */
  /*  ## */
  /*  ## */
CellTypeVertex, /* 201 */

  /* # # */
  /*  #  */
  /*  ## */
CellTypeVertex, /* 202 */

  /* # # */
  /*  ## */
  /*  ## */
CellTypeVertex, /* 203 */

  /* ##  */
  /*  #  */
  /*  ## */
CellTypeLine,   /* 204 */

  /* ##  */
  /*  ## */
  /*  ## */
CellTypeError,  /* 205 */

  /* ### */
  /*  #  */
  /*  ## */
CellTypeLine,   /* 206 */

  /* ### */
  /*  ## */
  /*  ## */
CellTypeError,  /* 207 */

  /*     */
  /* ##  */
  /*  ## */
CellTypeVertexOrLine, /* 208 */

  /*     */
  /* ### */
  /*  ## */
CellTypeVertex, /* 209 */

  /*   # */
  /* ##  */
  /*  ## */
CellTypeVertex, /* 210 */

  /*   # */
  /* ### */
  /*  ## */
CellTypeVertex, /* 211 */

  /*  #  */
  /* ##  */
  /*  ## */
CellTypeVertex, /* 212 */

  /*  #  */
  /* ### */
  /*  ## */
CellTypeVertex, /* 213 */

  /*  ## */
  /* ##  */
  /*  ## */
CellTypeVertex,  /* 214 */

  /*  ## */
  /* ### */
  /*  ## */
CellTypeVertex, /* 215 */

  /* #   */
  /* ##  */
  /*  ## */
CellTypeErrorOrLine,  /* 216 */

  /* #   */
  /* ### */
  /*  ## */
CellTypeError,  /* 217 */

  /* # # */
  /* ##  */
  /*  ## */
CellTypeVertex, /* 218 */

  /* # # */
  /* ### */
  /*  ## */
CellTypeError,  /* 219 */

  /* ##  */
  /* ##  */
  /*  ## */
CellTypeError,  /* 220 */

  /* ##  */
  /* ### */
  /*  ## */
CellTypeVertex, /* 221 */

  /* ### */
  /* ##  */
  /*  ## */
CellTypeError,  /* 222 */

  /* ### */
  /* ### */
  /*  ## */
CellTypeVertex, /* 223 */

  /*     */
  /*  #  */
  /* ### */
CellTypeVertex,  /* 224 */

  /*     */
  /*  ## */
  /* ### */
CellTypeError,  /* 225 */

  /*   # */
  /*  #  */
  /* ### */
CellTypeLine,   /* 226 */

  /*   # */
  /*  ## */
  /* ### */
CellTypeError,  /* 227 */

  /*  #  */
  /*  #  */
  /* ### */
CellTypeLine,   /* 228 */

  /*  #  */
  /*  ## */
  /* ### */
CellTypeVertex, /* 229 */

  /*  ## */
  /*  #  */
  /* ### */
CellTypeLine,   /* 230 */

  /*  ## */
  /*  ## */
  /* ### */
CellTypeError,  /* 231 */

  /* #   */
  /*  #  */
  /* ### */
CellTypeLine,   /* 232 */

  /* #   */
  /*  ## */
  /* ### */
CellTypeVertex, /* 233 */

  /* # # */
  /*  #  */
  /* ### */
CellTypeVertex, /* 234 */

  /* # # */
  /*  ## */
  /* ### */
CellTypeVertex, /* 235 */

  /* ##  */
  /*  #  */
  /* ### */
CellTypeLine,   /* 236 */

  /* ##  */
  /*  ## */
  /* ### */
CellTypeError,  /* 237 */

  /* ### */
  /*  #  */
  /* ### */
CellTypeLine,   /* 238 */

  /* ### */
  /*  ## */
  /* ### */
CellTypeError,  /* 239 */

  /*     */
  /* ##  */
  /* ### */
CellTypeError,  /* 240 */

  /*     */
  /* ### */
  /* ### */
CellTypeVertex, /* 241 */

  /*   # */
  /* ##  */
  /* ### */
CellTypeVertex, /* 242 */

  /*   # */
  /* ### */
  /* ### */
CellTypeError,  /* 243 */

  /*  #  */
  /* ##  */
  /* ### */
CellTypeVertex, /* 244 */

  /*  #  */
  /* ### */
  /* ### */
CellTypeVertex, /* 245 */

  /*  ## */
  /* ##  */
  /* ### */
CellTypeError,  /* 246 */

  /*  ## */
  /* ### */
  /* ### */
CellTypeVertex, /* 247 */

  /* #   */
  /* ##  */
  /* ### */
CellTypeError,  /* 248 */

  /* #   */
  /* ### */
  /* ### */
CellTypeError,  /* 249 */

  /* # # */
  /* ##  */
  /* ### */
CellTypeVertex, /* 250 */

  /* # # */
  /* ### */
  /* ### */
CellTypeError,  /* 251 */

  /* ##  */
  /* ##  */
  /* ### */
CellTypeError,  /* 252 */

  /* ##  */
  /* ### */
  /* ### */
CellTypeVertex, /* 253 */

  /* ### */
  /* ##  */
  /* ### */
CellTypeError,  /* 254 */

  /* ### */
  /* ### */
  /* ### */
CellTypeVertex  /* 255 */
};

} // namespace cellimage

} // namespace vigra

#endif /* VIGRA_CELLCONFIGURATIONS_HXX */
