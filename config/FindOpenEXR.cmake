# Copyright (c) 2011, Lukas Jirkovsky
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

FIND_PATH(OPENEXR_INCLUDE_DIR OpenEXRConfig.h PATH_SUFFIXES OpenEXR)

FIND_LIBRARY(OPENEXR_HALF_LIBRARY NAMES Half)
FIND_LIBRARY(OPENEXR_IEX_LIBRARY NAMES Iex)
FIND_LIBRARY(OPENEXR_ILMTHREAD_LIBRARY NAMES IlmThread)
FIND_LIBRARY(OPENEXR_IMATH_LIBRARY NAMES Imath)
FIND_LIBRARY(OPENEXR_ILMIMF_LIBRARY NAMES IlmImf)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENEXR DEFAULT_MSG
    OPENEXR_HALF_LIBRARY OPENEXR_IEX_LIBRARY OPENEXR_IMATH_LIBRARY
    OPENEXR_ILMIMF_LIBRARY OPENEXR_INCLUDE_DIR
)

IF(OPENEXR_FOUND)
    SET(OPENEXR_LIBRARIES ${OPENEXR_HALF_LIBRARY} ${OPENEXR_IEX_LIBRARY}
        ${OPENEXR_ILMTHREAD_LIBRARY} ${OPENEXR_IMATH_LIBRARY} ${OPENEXR_ILMIMF_LIBRARY})
ENDIF(OPENEXR_FOUND)
