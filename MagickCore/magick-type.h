/*
  Copyright 1999-2011 ImageMagick Studio LLC, a non-profit organization
  dedicated to making software imaging solutions freely available.
  
  You may not use this file except in compliance with the License.
  obtain a copy of the License at
  
    http://www.imagemagick.org/script/license.php
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  MagickCore types.
*/
#ifndef _MAGICKCORE_MAGICK_TYPE_H
#define _MAGICKCORE_MAGICK_TYPE_H

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include "MagickCore/magick-config.h"

#if !defined(MAGICKCORE_QUANTUM_DEPTH)
#define MAGICKCORE_QUANTUM_DEPTH  16
#endif

#if defined(MAGICKCORE_WINDOWS_SUPPORT) && !defined(__MINGW32__)
#  define MagickLLConstant(c)  (MagickOffsetType) (c ## i64)
#  define MagickULLConstant(c)  (MagickSizeType) (c ## ui64)
#else
#  define MagickLLConstant(c)  (MagickOffsetType) (c ## LL)
#  define MagickULLConstant(c)  (MagickSizeType) (c ## ULL)
#endif

#if (MAGICKCORE_QUANTUM_DEPTH == 8)
#define MagickEpsilon  1.0e-6
#define MagickHuge     1.0e6
#define MaxColormapSize  256UL
#define MaxMap  255UL

typedef double MagickRealType;
#if defined(MAGICKCORE_HDRI_SUPPORT)
typedef float Quantum;
#define QuantumRange  255.0
#define QuantumFormat  "%g"
#else
typedef unsigned char Quantum;
#define QuantumRange  255
#define QuantumFormat  "%u"
#endif
#elif (MAGICKCORE_QUANTUM_DEPTH == 16)
#define MagickEpsilon  1.0e-10
#define MagickHuge     1.0e12
#define MaxColormapSize  65536UL
#define MaxMap  65535UL

typedef double MagickRealType;
#if defined(MAGICKCORE_HDRI_SUPPORT)
typedef float Quantum;
#define QuantumRange  65535.0
#define QuantumFormat  "%g"
#else
typedef unsigned short Quantum;
#define QuantumRange  65535
#define QuantumFormat  "%u"
#endif
#elif (MAGICKCORE_QUANTUM_DEPTH == 32)
#define MagickEpsilon  1.0e-10
#define MagickHuge     1.0e12
#define MaxColormapSize  65536UL
#define MaxMap  65535UL

typedef double MagickRealType;
#if defined(MAGICKCORE_HDRI_SUPPORT)
typedef float Quantum;
#define QuantumRange  4294967295.0
#define QuantumFormat  "%g"
#else
typedef unsigned int Quantum;
#define QuantumRange  4294967295
#define QuantumFormat  "%u"
#endif
#elif (MAGICKCORE_QUANTUM_DEPTH == 64) && defined(MAGICKCORE_HAVE_LONG_DOUBLE_WIDER)
#define MagickEpsilon  1.0e-10
#define MagickHuge     1.0e12
#define MaxColormapSize  65536UL
#define MaxMap  65535UL

typedef long double MagickRealType;
typedef double Quantum;
#define QuantumRange  18446744073709551615.0
#define QuantumFormat  "%g"
#else
#if !defined(_CH_)
# error "MAGICKCORE_QUANTUM_DEPTH must be one of 8, 16, 32, or 64"
#endif
#endif
#define MagickPI     3.14159265358979323846264338327950288419716939937510L
#define QuantumScale  ((double) 1.0/(double) QuantumRange)

/*
  Typedef declarations.
*/
typedef unsigned int MagickStatusType;
#if !defined(MAGICKCORE_WINDOWS_SUPPORT)
#if (MAGICKCORE_SIZEOF_UNSIGNED_LONG_LONG == 8)
typedef long long MagickOffsetType;
typedef unsigned long long MagickSizeType;
#define MagickOffsetFormat  "lld"
#define MagickSizeFormat  "llu"
#else
typedef ssize_t MagickOffsetType;
typedef size_t MagickSizeType;
#define MagickOffsetFormat  "ld"
#define MagickSizeFormat  "lu"
#endif
#else
typedef __int64 MagickOffsetType;
typedef unsigned __int64 MagickSizeType;
#define MagickOffsetFormat  "I64i"
#define MagickSizeFormat  "I64u"
#endif

#if defined(_MSC_VER) && (_MSC_VER == 1200)
typedef MagickOffsetType QuantumAny;
#else
typedef MagickSizeType QuantumAny;
#endif

#if defined(macintosh)
#define ExceptionInfo  MagickExceptionInfo
#endif

typedef enum
{
  UndefinedChannel = 0x0000,
  RedChannel = 0x0001,
  GrayChannel = 0x0001,
  CyanChannel = 0x0001,
  GreenChannel = 0x0002,
  MagentaChannel = 0x0002,
  BlueChannel = 0x0004,
  YellowChannel = 0x0004,
  AlphaChannel = 0x0008,
  OpacityChannel = 0x0008,
  BlackChannel = 0x0010,
  IndexChannel = 0x0020,
  CompositeChannels = 0x002F,
  AllChannels = 0x7ffffff,
  /*
    Special purpose channel types.
  */
  TrueAlphaChannel = 0x0040, /* extract actual alpha channel from opacity */
  RGBChannels = 0x0080,      /* set alpha from grayscale mask in RGB */
  GrayChannels = 0x0080,
  SyncChannels = 0x0100,     /* channels should be modified equally */
  DefaultChannels = ((AllChannels | SyncChannels) &~ AlphaChannel)
} ChannelType;

typedef enum
{
  UndefinedClass,
  DirectClass,
  PseudoClass
} ClassType;

typedef enum
{
  MagickFalse = 0,
  MagickTrue = 1
} MagickBooleanType;

typedef struct _BlobInfo BlobInfo;

typedef struct _ExceptionInfo ExceptionInfo;

typedef struct _Image Image;

typedef struct _ImageInfo ImageInfo;

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif
