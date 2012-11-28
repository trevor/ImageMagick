/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%    M   M    OOO    RRRR   PPPP   H   H   OOO   L       OOO    GGGG  Y   Y   %
%    MM MM   O   O   R   R  P   P  H   H  O   O  L      O   O  G       Y Y    %
%    M M M   O   O   RRRR   PPPP   HHHHH  O   O  L      O   O  G GGG    Y     %
%    M   M   O   O   R R    P      H   H  O   O  L      O   O  G   G    Y     %
%    M   M    OOO    R  R   P      H   H   OOO   LLLLL   OOO    GGG     Y     %
%                                                                             %
%                                                                             %
%                        MagickCore Morphology Methods                        %
%                                                                             %
%                              Software Design                                %
%                              Anthony Thyssen                                %
%                               January 2010                                  %
%                                                                             %
%                                                                             %
%  Copyright 1999-2011 ImageMagick Studio LLC, a non-profit organization      %
%  dedicated to making software imaging solutions freely available.           %
%                                                                             %
%  You may not use this file except in compliance with the License.  You may  %
%  obtain a copy of the License at                                            %
%                                                                             %
%    http://www.imagemagick.org/script/license.php                            %
%                                                                             %
%  Unless required by applicable law or agreed to in writing, software        %
%  distributed under the License is distributed on an "AS IS" BASIS,          %
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   %
%  See the License for the specific language governing permissions and        %
%  limitations under the License.                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Morpology is the the application of various kernels, of any size and even
% shape, to a image in various ways (typically binary, but not always).
%
% Convolution (weighted sum or average) is just one specific type of
% morphology. Just one that is very common for image bluring and sharpening
% effects.  Not only 2D Gaussian blurring, but also 2-pass 1D Blurring.
%
% This module provides not only a general morphology function, and the ability
% to apply more advanced or iterative morphologies, but also functions for the
% generation of many different types of kernel arrays from user supplied
% arguments. Prehaps even the generation of a kernel from a small image.
*/

/*
  Include declarations.
*/
#include "magick/studio.h"
#include "magick/artifact.h"
#include "magick/cache-view.h"
#include "magick/color-private.h"
#include "magick/enhance.h"
#include "magick/exception.h"
#include "magick/exception-private.h"
#include "magick/gem.h"
#include "magick/hashmap.h"
#include "magick/image.h"
#include "magick/image-private.h"
#include "magick/list.h"
#include "magick/magick.h"
#include "magick/memory_.h"
#include "magick/monitor-private.h"
#include "magick/morphology.h"
#include "magick/morphology-private.h"
#include "magick/option.h"
#include "magick/pixel-private.h"
#include "magick/prepress.h"
#include "magick/quantize.h"
#include "magick/registry.h"
#include "magick/semaphore.h"
#include "magick/splay-tree.h"
#include "magick/statistic.h"
#include "magick/string_.h"
#include "magick/string-private.h"
#include "magick/token.h"
#include "magick/utility.h"


/*
** The following test is for special floating point numbers of value NaN (not
** a number), that may be used within a Kernel Definition.  NaN's are defined
** as part of the IEEE standard for floating point number representation.
**
** These are used as a Kernel value to mean that this kernel position is not
** part of the kernel neighbourhood for convolution or morphology processing,
** and thus should be ignored.  This allows the use of 'shaped' kernels.
**
** The special properity that two NaN's are never equal, even if they are from
** the same variable allow you to test if a value is special NaN value.
**
** This macro  IsNaN() is thus is only true if the value given is NaN.
*/
#define IsNan(a)   ((a)!=(a))

/*
  Other global definitions used by module.
*/
static inline double MagickMin(const double x,const double y)
{
  return( x < y ? x : y);
}
static inline double MagickMax(const double x,const double y)
{
  return( x > y ? x : y);
}
#define Minimize(assign,value) assign=MagickMin(assign,value)
#define Maximize(assign,value) assign=MagickMax(assign,value)

/* Currently these are only internal to this module */
static void
  CalcKernelMetaData(KernelInfo *),
  ExpandMirrorKernelInfo(KernelInfo *),
  ExpandRotateKernelInfo(KernelInfo *, const double),
  RotateKernelInfo(KernelInfo *, double);


/* Quick function to find last kernel in a kernel list */
static inline KernelInfo *LastKernelInfo(KernelInfo *kernel)
{
  while (kernel->next != (KernelInfo *) NULL)
    kernel = kernel->next;
  return(kernel);
}


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     A c q u i r e K e r n e l I n f o                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  AcquireKernelInfo() takes the given string (generally supplied by the
%  user) and converts it into a Morphology/Convolution Kernel.  This allows
%  users to specify a kernel from a number of pre-defined kernels, or to fully
%  specify their own kernel for a specific Convolution or Morphology
%  Operation.
%
%  The kernel so generated can be any rectangular array of floating point
%  values (doubles) with the 'control point' or 'pixel being affected'
%  anywhere within that array of values.
%
%  Previously IM was restricted to a square of odd size using the exact
%  center as origin, this is no longer the case, and any rectangular kernel
%  with any value being declared the origin. This in turn allows the use of
%  highly asymmetrical kernels.
%
%  The floating point values in the kernel can also include a special value
%  known as 'nan' or 'not a number' to indicate that this value is not part
%  of the kernel array. This allows you to shaped the kernel within its
%  rectangular area. That is 'nan' values provide a 'mask' for the kernel
%  shape.  However at least one non-nan value must be provided for correct
%  working of a kernel.
%
%  The returned kernel should be freed using the DestroyKernelInfo() when you
%  are finished with it.  Do not free this memory yourself.
%
%  Input kernel defintion strings can consist of any of three types.
%
%    "name:args[[@><]"
%         Select from one of the built in kernels, using the name and
%         geometry arguments supplied.  See AcquireKernelBuiltIn()
%
%    "WxH[+X+Y][@><]:num, num, num ..."
%         a kernel of size W by H, with W*H floating point numbers following.
%         the 'center' can be optionally be defined at +X+Y (such that +0+0
%         is top left corner). If not defined the pixel in the center, for
%         odd sizes, or to the immediate top or left of center for even sizes
%         is automatically selected.
%
%    "num, num, num, num, ..."
%         list of floating point numbers defining an 'old style' odd sized
%         square kernel.  At least 9 values should be provided for a 3x3
%         square kernel, 25 for a 5x5 square kernel, 49 for 7x7, etc.
%         Values can be space or comma separated.  This is not recommended.
%
%  You can define a 'list of kernels' which can be used by some morphology
%  operators A list is defined as a semi-colon seperated list kernels.
%
%     " kernel ; kernel ; kernel ; "
%
%  Any extra ';' characters, at start, end or between kernel defintions are
%  simply ignored.
%
%  The special flags will expand a single kernel, into a list of rotated
%  kernels. A '@' flag will expand a 3x3 kernel into a list of 45-degree
%  cyclic rotations, while a '>' will generate a list of 90-degree rotations.
%  The '<' also exands using 90-degree rotates, but giving a 180-degree
%  reflected kernel before the +/- 90-degree rotations, which can be important
%  for Thinning operations.
%
%  Note that 'name' kernels will start with an alphabetic character while the
%  new kernel specification has a ':' character in its specification string.
%  If neither is the case, it is assumed an old style of a simple list of
%  numbers generating a odd-sized square kernel has been given.
%
%  The format of the AcquireKernal method is:
%
%      KernelInfo *AcquireKernelInfo(const char *kernel_string)
%
%  A description of each parameter follows:
%
%    o kernel_string: the Morphology/Convolution kernel wanted.
%
*/

/* This was separated so that it could be used as a separate
** array input handling function, such as for -color-matrix
*/
static KernelInfo *ParseKernelArray(const char *kernel_string)
{
  KernelInfo
    *kernel;

  char
    token[MaxTextExtent];

  const char
    *p,
    *end;

  register ssize_t
    i;

  double
    nan = sqrt((double)-1.0);  /* Special Value : Not A Number */

  MagickStatusType
    flags;

  GeometryInfo
    args;

  kernel=(KernelInfo *) AcquireMagickMemory(sizeof(*kernel));
  if (kernel == (KernelInfo *)NULL)
    return(kernel);
  (void) ResetMagickMemory(kernel,0,sizeof(*kernel));
  kernel->minimum = kernel->maximum = kernel->angle = 0.0;
  kernel->negative_range = kernel->positive_range = 0.0;
  kernel->type = UserDefinedKernel;
  kernel->next = (KernelInfo *) NULL;
  kernel->signature = MagickSignature;

  /* find end of this specific kernel definition string */
  end = strchr(kernel_string, ';');
  if ( end == (char *) NULL )
    end = strchr(kernel_string, '\0');

  /* clear flags - for Expanding kernal lists thorugh rotations */
   flags = NoValue;

  /* Has a ':' in argument - New user kernel specification */
  p = strchr(kernel_string, ':');
  if ( p != (char *) NULL && p < end)
    {
      /* ParseGeometry() needs the geometry separated! -- Arrgghh */
      memcpy(token, kernel_string, (size_t) (p-kernel_string));
      token[p-kernel_string] = '\0';
      SetGeometryInfo(&args);
      flags = ParseGeometry(token, &args);

      /* Size handling and checks of geometry settings */
      if ( (flags & WidthValue) == 0 ) /* if no width then */
        args.rho = args.sigma;         /* then  width = height */
      if ( args.rho < 1.0 )            /* if width too small */
         args.rho = 1.0;               /* then  width = 1 */
      if ( args.sigma < 1.0 )          /* if height too small */
        args.sigma = args.rho;         /* then  height = width */
      kernel->width = (size_t)args.rho;
      kernel->height = (size_t)args.sigma;

      /* Offset Handling and Checks */
      if ( args.xi  < 0.0 || args.psi < 0.0 )
        return(DestroyKernelInfo(kernel));
      kernel->x = ((flags & XValue)!=0) ? (ssize_t)args.xi
                                               : (ssize_t) (kernel->width-1)/2;
      kernel->y = ((flags & YValue)!=0) ? (ssize_t)args.psi
                                               : (ssize_t) (kernel->height-1)/2;
      if ( kernel->x >= (ssize_t) kernel->width ||
           kernel->y >= (ssize_t) kernel->height )
        return(DestroyKernelInfo(kernel));

      p++; /* advance beyond the ':' */
    }
  else
    { /* ELSE - Old old specification, forming odd-square kernel */
      /* count up number of values given */
      p=(const char *) kernel_string;
      while ((isspace((int) ((unsigned char) *p)) != 0) || (*p == '\''))
        p++;  /* ignore "'" chars for convolve filter usage - Cristy */
      for (i=0; p < end; i++)
      {
        GetMagickToken(p,&p,token);
        if (*token == ',')
          GetMagickToken(p,&p,token);
      }
      /* set the size of the kernel - old sized square */
      kernel->width = kernel->height= (size_t) sqrt((double) i+1.0);
      kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;
      p=(const char *) kernel_string;
      while ((isspace((int) ((unsigned char) *p)) != 0) || (*p == '\''))
        p++;  /* ignore "'" chars for convolve filter usage - Cristy */
    }

  /* Read in the kernel values from rest of input string argument */
  kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                        kernel->height*sizeof(double));
  if (kernel->values == (double *) NULL)
    return(DestroyKernelInfo(kernel));

  kernel->minimum = +MagickHuge;
  kernel->maximum = -MagickHuge;
  kernel->negative_range = kernel->positive_range = 0.0;

  for (i=0; (i < (ssize_t) (kernel->width*kernel->height)) && (p < end); i++)
  {
    GetMagickToken(p,&p,token);
    if (*token == ',')
      GetMagickToken(p,&p,token);
    if (    LocaleCompare("nan",token) == 0
        || LocaleCompare("-",token) == 0 ) {
      kernel->values[i] = nan; /* do not include this value in kernel */
    }
    else {
      kernel->values[i] = StringToDouble(token);
      ( kernel->values[i] < 0)
          ?  ( kernel->negative_range += kernel->values[i] )
          :  ( kernel->positive_range += kernel->values[i] );
      Minimize(kernel->minimum, kernel->values[i]);
      Maximize(kernel->maximum, kernel->values[i]);
    }
  }

  /* sanity check -- no more values in kernel definition */
  GetMagickToken(p,&p,token);
  if ( *token != '\0' && *token != ';' && *token != '\'' )
    return(DestroyKernelInfo(kernel));

#if 0
  /* this was the old method of handling a incomplete kernel */
  if ( i < (ssize_t) (kernel->width*kernel->height) ) {
    Minimize(kernel->minimum, kernel->values[i]);
    Maximize(kernel->maximum, kernel->values[i]);
    for ( ; i < (ssize_t) (kernel->width*kernel->height); i++)
      kernel->values[i]=0.0;
  }
#else
  /* Number of values for kernel was not enough - Report Error */
  if ( i < (ssize_t) (kernel->width*kernel->height) )
    return(DestroyKernelInfo(kernel));
#endif

  /* check that we recieved at least one real (non-nan) value! */
  if ( kernel->minimum == MagickHuge )
    return(DestroyKernelInfo(kernel));

  if ( (flags & AreaValue) != 0 )         /* '@' symbol in kernel size */
    ExpandRotateKernelInfo(kernel, 45.0); /* cyclic rotate 3x3 kernels */
  else if ( (flags & GreaterValue) != 0 ) /* '>' symbol in kernel args */
    ExpandRotateKernelInfo(kernel, 90.0); /* 90 degree rotate of kernel */
  else if ( (flags & LessValue) != 0 )    /* '<' symbol in kernel args */
    ExpandMirrorKernelInfo(kernel);       /* 90 degree mirror rotate */

  return(kernel);
}

static KernelInfo *ParseKernelName(const char *kernel_string)
{
  KernelInfo
    *kernel;

  char
    token[MaxTextExtent];

  ssize_t
    type;

  const char
    *p,
    *end;

  MagickStatusType
    flags;

  GeometryInfo
    args;

  /* Parse special 'named' kernel */
  GetMagickToken(kernel_string,&p,token);
  type=ParseMagickOption(MagickKernelOptions,MagickFalse,token);
  if ( type < 0 || type == UserDefinedKernel )
    return((KernelInfo *)NULL);  /* not a valid named kernel */

  while (((isspace((int) ((unsigned char) *p)) != 0) ||
          (*p == ',') || (*p == ':' )) && (*p != '\0') && (*p != ';'))
    p++;

  end = strchr(p, ';'); /* end of this kernel defintion */
  if ( end == (char *) NULL )
    end = strchr(p, '\0');

  /* ParseGeometry() needs the geometry separated! -- Arrgghh */
  memcpy(token, p, (size_t) (end-p));
  token[end-p] = '\0';
  SetGeometryInfo(&args);
  flags = ParseGeometry(token, &args);

#if 0
  /* For Debugging Geometry Input */
  fprintf(stderr, "Geometry = 0x%04X : %lg x %lg %+lg %+lg\n",
       flags, args.rho, args.sigma, args.xi, args.psi );
#endif

  /* special handling of missing values in input string */
  switch( type ) {
    case RectangleKernel:
      if ( (flags & WidthValue) == 0 ) /* if no width then */
        args.rho = args.sigma;         /* then  width = height */
      if ( args.rho < 1.0 )            /* if width too small */
          args.rho = 3;                 /* then  width = 3 */
      if ( args.sigma < 1.0 )          /* if height too small */
        args.sigma = args.rho;         /* then  height = width */
      if ( (flags & XValue) == 0 )     /* center offset if not defined */
        args.xi = (double)(((ssize_t)args.rho-1)/2);
      if ( (flags & YValue) == 0 )
        args.psi = (double)(((ssize_t)args.sigma-1)/2);
      break;
    case SquareKernel:
    case DiamondKernel:
    case DiskKernel:
    case PlusKernel:
    case CrossKernel:
      /* If no scale given (a 0 scale is valid! - set it to 1.0 */
      if ( (flags & HeightValue) == 0 )
        args.sigma = 1.0;
      break;
    case RingKernel:
      if ( (flags & XValue) == 0 )
        args.xi = 1.0;
      break;
    case ChebyshevKernel:
    case ManhattanKernel:
    case EuclideanKernel:
      if ( (flags & HeightValue) == 0 )           /* no distance scale */
        args.sigma = 100.0;                       /* default distance scaling */
      else if ( (flags & AspectValue ) != 0 )     /* '!' flag */
        args.sigma = QuantumRange/(args.sigma+1); /* maximum pixel distance */
      else if ( (flags & PercentValue ) != 0 )    /* '%' flag */
        args.sigma *= QuantumRange/100.0;         /* percentage of color range */
      break;
    default:
      break;
  }

  kernel = AcquireKernelBuiltIn((KernelInfoType)type, &args);

  /* global expand to rotated kernel list - only for single kernels */
  if ( kernel->next == (KernelInfo *) NULL ) {
    if ( (flags & AreaValue) != 0 )         /* '@' symbol in kernel args */
      ExpandRotateKernelInfo(kernel, 45.0);
    else if ( (flags & GreaterValue) != 0 ) /* '>' symbol in kernel args */
      ExpandRotateKernelInfo(kernel, 90.0);
    else if ( (flags & LessValue) != 0 )    /* '<' symbol in kernel args */
      ExpandMirrorKernelInfo(kernel);
  }

  return(kernel);
}

MagickExport KernelInfo *AcquireKernelInfo(const char *kernel_string)
{

  KernelInfo
    *kernel,
    *new_kernel;

  char
    token[MaxTextExtent];

  const char
    *p;

  size_t
    kernel_number;

  p = kernel_string;
  kernel = NULL;
  kernel_number = 0;

  while ( GetMagickToken(p,NULL,token),  *token != '\0' ) {

    /* ignore extra or multiple ';' kernel seperators */
    if ( *token != ';' ) {

      /* tokens starting with alpha is a Named kernel */
      if (isalpha((int) *token) != 0)
        new_kernel = ParseKernelName(p);
      else /* otherwise a user defined kernel array */
        new_kernel = ParseKernelArray(p);

      /* Error handling -- this is not proper error handling! */
      if ( new_kernel == (KernelInfo *) NULL ) {
        fprintf(stderr, "Failed to parse kernel number #%.20g\n",(double)
          kernel_number);
        if ( kernel != (KernelInfo *) NULL )
          kernel=DestroyKernelInfo(kernel);
        return((KernelInfo *) NULL);
      }

      /* initialise or append the kernel list */
      if ( kernel == (KernelInfo *) NULL )
        kernel = new_kernel;
      else
        LastKernelInfo(kernel)->next = new_kernel;
    }

    /* look for the next kernel in list */
    p = strchr(p, ';');
    if ( p == (char *) NULL )
      break;
    p++;

  }
  return(kernel);
}


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     A c q u i r e K e r n e l B u i l t I n                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  AcquireKernelBuiltIn() returned one of the 'named' built-in types of
%  kernels used for special purposes such as gaussian blurring, skeleton
%  pruning, and edge distance determination.
%
%  They take a KernelType, and a set of geometry style arguments, which were
%  typically decoded from a user supplied string, or from a more complex
%  Morphology Method that was requested.
%
%  The format of the AcquireKernalBuiltIn method is:
%
%      KernelInfo *AcquireKernelBuiltIn(const KernelInfoType type,
%           const GeometryInfo args)
%
%  A description of each parameter follows:
%
%    o type: the pre-defined type of kernel wanted
%
%    o args: arguments defining or modifying the kernel
%
%  Convolution Kernels
%
%    Unity
%       the No-Op kernel, also requivelent to  Gaussian of sigma zero.
%       Basically a 3x3 kernel of a 1 surrounded by zeros.
%
%    Gaussian:{radius},{sigma}
%       Generate a two-dimentional gaussian kernel, as used by -gaussian.
%       The sigma for the curve is required.  The resulting kernel is
%       normalized,
%
%       If 'sigma' is zero, you get a single pixel on a field of zeros.
%
%       NOTE: that the 'radius' is optional, but if provided can limit (clip)
%       the final size of the resulting kernel to a square 2*radius+1 in size.
%       The radius should be at least 2 times that of the sigma value, or
%       sever clipping and aliasing may result.  If not given or set to 0 the
%       radius will be determined so as to produce the best minimal error
%       result, which is usally much larger than is normally needed.
%
%    LoG:{radius},{sigma}
%        "Laplacian of a Gaussian" or "Mexician Hat" Kernel.
%        The supposed ideal edge detection, zero-summing kernel.
%
%        An alturnative to this kernel is to use a "DoG" with a sigma ratio of
%        approx 1.6 (according to wikipedia).
%
%    DoG:{radius},{sigma1},{sigma2}
%        "Difference of Gaussians" Kernel.
%        As "Gaussian" but with a gaussian produced by 'sigma2' subtracted
%        from the gaussian produced by 'sigma1'. Typically sigma2 > sigma1.
%        The result is a zero-summing kernel.
%
%    Blur:{radius},{sigma}[,{angle}]
%       Generates a 1 dimensional or linear gaussian blur, at the angle given
%       (current restricted to orthogonal angles).  If a 'radius' is given the
%       kernel is clipped to a width of 2*radius+1.  Kernel can be rotated
%       by a 90 degree angle.
%
%       If 'sigma' is zero, you get a single pixel on a field of zeros.
%
%       Note that two convolutions with two "Blur" kernels perpendicular to
%       each other, is equivelent to a far larger "Gaussian" kernel with the
%       same sigma value, However it is much faster to apply. This is how the
%       "-blur" operator actually works.
%
%    Comet:{width},{sigma},{angle}
%       Blur in one direction only, much like how a bright object leaves
%       a comet like trail.  The Kernel is actually half a gaussian curve,
%       Adding two such blurs in opposite directions produces a Blur Kernel.
%       Angle can be rotated in multiples of 90 degrees.
%
%       Note that the first argument is the width of the kernel and not the
%       radius of the kernel.
%
%    # Still to be implemented...
%    #
%    # Filter2D
%    # Filter1D
%    #    Set kernel values using a resize filter, and given scale (sigma)
%    #    Cylindrical or Linear.   Is this posible with an image?
%    #
%
%  Named Constant Convolution Kernels
%
%  All these are unscaled, zero-summing kernels by default. As such for
%  non-HDRI version of ImageMagick some form of normalization, user scaling,
%  and biasing the results is recommended, to prevent the resulting image
%  being 'clipped'.
%
%  The 3x3 kernels (most of these) can be circularly rotated in multiples of
%  45 degrees to generate the 8 angled varients of each of the kernels.
%
%    Laplacian:{type}
%      Discrete Lapacian Kernels, (without normalization)
%        Type 0 :  3x3 with center:8 surounded by -1  (8 neighbourhood)
%        Type 1 :  3x3 with center:4 edge:-1 corner:0 (4 neighbourhood)
%        Type 2 :  3x3 with center:4 edge:1 corner:-2
%        Type 3 :  3x3 with center:4 edge:-2 corner:1
%        Type 5 :  5x5 laplacian
%        Type 7 :  7x7 laplacian
%        Type 15 : 5x5 LoG (sigma approx 1.4)
%        Type 19 : 9x9 LoG (sigma approx 1.4)
%
%    Sobel:{angle}
%      Sobel 'Edge' convolution kernel (3x3)
%          | -1, 0, 1 |
%          | -2, 0,-2 |
%          | -1, 0, 1 |
%
%    Sobel:{type},{angle}
%      Type 0:  default un-nomalized version shown above.
%
%      Type 1:  As default but pre-normalized
%          | 1, 0, -1 |
%          | 2, 0, -2 |  / 4
%          | 1, 0, -1 |
%
%      Type 2:  Diagonal version with same normalization as 1
%          | 1, 0, -1 |
%          | 2, 0, -2 |  / 4
%          | 1, 0, -1 |
%
%    Roberts:{angle}
%      Roberts convolution kernel (3x3)
%          |  0, 0, 0 |
%          | -1, 1, 0 |
%          |  0, 0, 0 |
%
%    Prewitt:{angle}
%      Prewitt Edge convolution kernel (3x3)
%          | -1, 0, 1 |
%          | -1, 0, 1 |
%          | -1, 0, 1 |
%
%    Compass:{angle}
%      Prewitt's "Compass" convolution kernel (3x3)
%          | -1, 1, 1 |
%          | -1,-2, 1 |
%          | -1, 1, 1 |
%
%    Kirsch:{angle}
%      Kirsch's "Compass" convolution kernel (3x3)
%          | -3,-3, 5 |
%          | -3, 0, 5 |
%          | -3,-3, 5 |
%
%    FreiChen:{angle}
%      Frei-Chen Edge Detector is based on a kernel that is similar to
%      the Sobel Kernel, but is designed to be isotropic. That is it takes
%      into account the distance of the diagonal in the kernel.
%
%          |   1,     0,   -1     |
%          | sqrt(2), 0, -sqrt(2) |
%          |   1,     0,   -1     |
%
%    FreiChen:{type},{angle}
%
%      Frei-Chen Pre-weighted kernels...
%
%        Type 0:  default un-nomalized version shown above.
%
%        Type 1: Orthogonal Kernel (same as type 11 below)
%          |   1,     0,   -1     |
%          | sqrt(2), 0, -sqrt(2) | / 2*sqrt(2)
%          |   1,     0,   -1     |
%
%        Type 2: Diagonal form of Kernel...
%          |   1,     sqrt(2),    0     |
%          | sqrt(2),   0,     -sqrt(2) | / 2*sqrt(2)
%          |   0,    -sqrt(2)    -1     |
%
%      However this kernel is als at the heart of the FreiChen Edge Detection
%      Process which uses a set of 9 specially weighted kernel.  These 9
%      kernels not be normalized, but directly applied to the image. The
%      results is then added together, to produce the intensity of an edge in
%      a specific direction.  The square root of the pixel value can then be
%      taken as the cosine of the edge, and at least 2 such runs at 90 degrees
%      from each other, both the direction and the strength of the edge can be
%      determined.
%
%        Type 10: All 9 of the following pre-weighted kernels...
%
%        Type 11: |   1,     0,   -1     |
%                 | sqrt(2), 0, -sqrt(2) | / 2*sqrt(2)
%                 |   1,     0,   -1     |
%
%        Type 12: | 1, sqrt(2), 1 |
%                 | 0,   0,     0 | / 2*sqrt(2)
%                 | 1, sqrt(2), 1 |
%
%        Type 13: | sqrt(2), -1,    0     |
%                 |  -1,      0,    1     | / 2*sqrt(2)
%                 |   0,      1, -sqrt(2) |
%
%        Type 14: |    0,     1, -sqrt(2) |
%                 |   -1,     0,     1    | / 2*sqrt(2)
%                 | sqrt(2), -1,     0    |
%
%        Type 15: | 0, -1, 0 |
%                 | 1,  0, 1 | / 2
%                 | 0, -1, 0 |
%
%        Type 16: |  1, 0, -1 |
%                 |  0, 0,  0 | / 2
%                 | -1, 0,  1 |
%
%        Type 17: |  1, -2,  1 |
%                 | -2,  4, -2 | / 6
%                 | -1, -2,  1 |
%
%        Type 18: | -2, 1, -2 |
%                 |  1, 4,  1 | / 6
%                 | -2, 1, -2 |
%
%        Type 19: | 1, 1, 1 |
%                 | 1, 1, 1 | / 3
%                 | 1, 1, 1 |
%
%      The first 4 are for edge detection, the next 4 are for line detection
%      and the last is to add a average component to the results.
%
%      Using a special type of '-1' will return all 9 pre-weighted kernels
%      as a multi-kernel list, so that you can use them directly (without
%      normalization) with the special "-set option:morphology:compose Plus"
%      setting to apply the full FreiChen Edge Detection Technique.
%
%      If 'type' is large it will be taken to be an actual rotation angle for
%      the default FreiChen (type 0) kernel.  As such  FreiChen:45  will look
%      like a  Sobel:45  but with 'sqrt(2)' instead of '2' values.
%
%      WARNING: The above was layed out as per
%          http://www.math.tau.ac.il/~turkel/notes/edge_detectors.pdf
%      But rotated 90 degrees so direction is from left rather than the top.
%      I have yet to find any secondary confirmation of the above. The only
%      other source found was actual source code at
%          http://ltswww.epfl.ch/~courstiv/exos_labos/sol3.pdf
%      Neigher paper defineds the kernels in a way that looks locical or
%      correct when taken as a whole.
%
%  Boolean Kernels
%
%    Diamond:[{radius}[,{scale}]]
%       Generate a diamond shaped kernel with given radius to the points.
%       Kernel size will again be radius*2+1 square and defaults to radius 1,
%       generating a 3x3 kernel that is slightly larger than a square.
%
%    Square:[{radius}[,{scale}]]
%       Generate a square shaped kernel of size radius*2+1, and defaulting
%       to a 3x3 (radius 1).
%
%       Note that using a larger radius for the "Square" or the "Diamond" is
%       also equivelent to iterating the basic morphological method that many
%       times. However iterating with the smaller radius is actually faster
%       than using a larger kernel radius.
%
%    Rectangle:{geometry}
%       Simply generate a rectangle of 1's with the size given. You can also
%       specify the location of the 'control point', otherwise the closest
%       pixel to the center of the rectangle is selected.
%
%       Properly centered and odd sized rectangles work the best.
%
%    Disk:[{radius}[,{scale}]]
%       Generate a binary disk of the radius given, radius may be a float.
%       Kernel size will be ceil(radius)*2+1 square.
%       NOTE: Here are some disk shapes of specific interest
%          "Disk:1"    => "diamond" or "cross:1"
%          "Disk:1.5"  => "square"
%          "Disk:2"    => "diamond:2"
%          "Disk:2.5"  => a general disk shape of radius 2
%          "Disk:2.9"  => "square:2"
%          "Disk:3.5"  => default - octagonal/disk shape of radius 3
%          "Disk:4.2"  => roughly octagonal shape of radius 4
%          "Disk:4.3"  => a general disk shape of radius 4
%       After this all the kernel shape becomes more and more circular.
%
%       Because a "disk" is more circular when using a larger radius, using a
%       larger radius is preferred over iterating the morphological operation.
%
%  Symbol Dilation Kernels
%
%    These kernel is not a good general morphological kernel, but is used
%    more for highlighting and marking any single pixels in an image using,
%    a "Dilate" method as appropriate.
%
%    For the same reasons iterating these kernels does not produce the
%    same result as using a larger radius for the symbol.
%
%    Plus:[{radius}[,{scale}]]
%    Cross:[{radius}[,{scale}]]
%       Generate a kernel in the shape of a 'plus' or a 'cross' with
%       a each arm the length of the given radius (default 2).
%
%       NOTE: "plus:1" is equivelent to a "Diamond" kernel.
%
%    Ring:{radius1},{radius2}[,{scale}]
%       A ring of the values given that falls between the two radii.
%       Defaults to a ring of approximataly 3 radius in a 7x7 kernel.
%       This is the 'edge' pixels of the default "Disk" kernel,
%       More specifically, "Ring" -> "Ring:2.5,3.5,1.0"
%
%  Hit and Miss Kernels
%
%    Peak:radius1,radius2
%       Find any peak larger than the pixels the fall between the two radii.
%       The default ring of pixels is as per "Ring".
%    Edges
%       Find flat orthogonal edges of a binary shape
%    Corners
%       Find 90 degree corners of a binary shape
%    LineEnds:type
%       Find end points of lines (for pruning a skeletion)
%       Two types of lines ends (default to both) can be searched for
%         Type 0: All line ends
%         Type 1: single kernel for 4-conneected line ends
%         Type 2: single kernel for simple line ends
%    LineJunctions
%       Find three line junctions (within a skeletion)
%         Type 0: all line junctions
%         Type 1: Y Junction kernel
%         Type 2: Diagonal T Junction kernel
%         Type 3: Orthogonal T Junction kernel
%         Type 4: Diagonal X Junction kernel
%         Type 5: Orthogonal + Junction kernel
%    Ridges:type
%       Find single pixel ridges or thin lines
%         Type 1: Fine single pixel thick lines and ridges
%         Type 2: Find two pixel thick lines and ridges
%    ConvexHull
%       Octagonal thicken kernel, to generate convex hulls of 45 degrees
%    Skeleton:type
%       Traditional skeleton generating kernels.
%         Type 1: Tradional Skeleton kernel (4 connected skeleton)
%         Type 2: HIPR2 Skeleton kernel (8 connected skeleton)
%         Type 3: Experimental Variation to try to present left-right symmetry
%         Type 4: Experimental Variation to preserve left-right symmetry
%
%  Distance Measuring Kernels
%
%    Different types of distance measuring methods, which are used with the
%    a 'Distance' morphology method for generating a gradient based on
%    distance from an edge of a binary shape, though there is a technique
%    for handling a anti-aliased shape.
%
%    See the 'Distance' Morphological Method, for information of how it is
%    applied.
%
%    Chebyshev:[{radius}][x{scale}[%!]]
%       Chebyshev Distance (also known as Tchebychev Distance) is a value of
%       one to any neighbour, orthogonal or diagonal. One why of thinking of
%       it is the number of squares a 'King' or 'Queen' in chess needs to
%       traverse reach any other position on a chess board.  It results in a
%       'square' like distance function, but one where diagonals are closer
%       than expected.
%
%    Manhattan:[{radius}][x{scale}[%!]]
%       Manhattan Distance (also known as Rectilinear Distance, or the Taxi
%       Cab metric), is the distance needed when you can only travel in
%       orthogonal (horizontal or vertical) only.  It is the distance a 'Rook'
%       in chess would travel. It results in a diamond like distances, where
%       diagonals are further than expected.
%
%    Euclidean:[{radius}][x{scale}[%!]]
%       Euclidean Distance is the 'direct' or 'as the crow flys distance.
%       However by default the kernel size only has a radius of 1, which
%       limits the distance to 'Knight' like moves, with only orthogonal and
%       diagonal measurements being correct.  As such for the default kernel
%       you will get octagonal like distance function, which is reasonally
%       accurate.
%
%       However if you use a larger radius such as "Euclidean:4" you will
%       get a much smoother distance gradient from the edge of the shape.
%       Of course a larger kernel is slower to use, and generally not needed.
%
%       To allow the use of fractional distances that you get with diagonals
%       the actual distance is scaled by a fixed value which the user can
%       provide.  This is not actually nessary for either ""Chebyshev" or
%       "Manhattan" distance kernels, but is done for all three distance
%       kernels.  If no scale is provided it is set to a value of 100,
%       allowing for a maximum distance measurement of 655 pixels using a Q16
%       version of IM, from any edge.  However for small images this can
%       result in quite a dark gradient.
%
*/

MagickExport KernelInfo *AcquireKernelBuiltIn(const KernelInfoType type,
   const GeometryInfo *args)
{
  KernelInfo
    *kernel;

  register ssize_t
    i;

  register ssize_t
    u,
    v;

  double
    nan = sqrt((double)-1.0);  /* Special Value : Not A Number */

  /* Generate a new empty kernel if needed */
  kernel=(KernelInfo *) NULL;
  switch(type) {
    case UndefinedKernel:    /* These should not call this function */
    case UserDefinedKernel:
      break;
    case UnityKernel:      /* Named Descrete Convolution Kernels */
    case LaplacianKernel:
    case SobelKernel:
    case RobertsKernel:
    case PrewittKernel:
    case CompassKernel:
    case KirschKernel:
    case FreiChenKernel:
    case EdgesKernel:       /* Hit and Miss kernels */
    case CornersKernel:
    case ThinDiagonalsKernel:
    case LineEndsKernel:
    case LineJunctionsKernel:
    case RidgesKernel:
    case ConvexHullKernel:
    case SkeletonKernel:
      break;               /* A pre-generated kernel is not needed */
#if 0
    /* set to 1 to do a compile-time check that we haven't missed anything */
    case GaussianKernel:
    case DoGKernel:
    case LoGKernel:
    case BlurKernel:
    case CometKernel:
    case DiamondKernel:
    case SquareKernel:
    case RectangleKernel:
    case DiskKernel:
    case PlusKernel:
    case CrossKernel:
    case RingKernel:
    case PeaksKernel:
    case ChebyshevKernel:
    case ManhattanKernel:
    case EuclideanKernel:
#else
    default:
#endif
      /* Generate the base Kernel Structure */
      kernel=(KernelInfo *) AcquireMagickMemory(sizeof(*kernel));
      if (kernel == (KernelInfo *) NULL)
        return(kernel);
      (void) ResetMagickMemory(kernel,0,sizeof(*kernel));
      kernel->minimum = kernel->maximum = kernel->angle = 0.0;
      kernel->negative_range = kernel->positive_range = 0.0;
      kernel->type = type;
      kernel->next = (KernelInfo *) NULL;
      kernel->signature = MagickSignature;
      break;
  }

  switch(type) {
    /* Convolution Kernels */
    case GaussianKernel:
    case DoGKernel:
    case LoGKernel:
      { double
          sigma = fabs(args->sigma),
          sigma2 = fabs(args->xi),
          A, B, R;

        if ( args->rho >= 1.0 )
          kernel->width = (size_t)args->rho*2+1;
        else if ( (type != DoGKernel) || (sigma >= sigma2) )
          kernel->width = GetOptimalKernelWidth2D(args->rho,sigma);
        else
          kernel->width = GetOptimalKernelWidth2D(args->rho,sigma2);
        kernel->height = kernel->width;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;
        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* WARNING: The following generates a 'sampled gaussian' kernel.
         * What we really want is a 'discrete gaussian' kernel.
         *
         * How to do this is currently not known, but appears to be
         * basied on the Error Function 'erf()' (intergral of a gaussian)
         */

        if ( type == GaussianKernel || type == DoGKernel )
          { /* Calculate a Gaussian,  OR positive half of a DoG */
            if ( sigma > MagickEpsilon )
              { A = 1.0/(2.0*sigma*sigma);  /* simplify loop expressions */
                B = (double) (1.0/(Magick2PI*sigma*sigma));
                for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
                  for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
                      kernel->values[i] = exp(-((double)(u*u+v*v))*A)*B;
              }
            else /* limiting case - a unity (normalized Dirac) kernel */
              { (void) ResetMagickMemory(kernel->values,0, (size_t)
                            kernel->width*kernel->height*sizeof(double));
                kernel->values[kernel->x+kernel->y*kernel->width] = 1.0;
              }
          }

        if ( type == DoGKernel )
          { /* Subtract a Negative Gaussian for "Difference of Gaussian" */
            if ( sigma2 > MagickEpsilon )
              { sigma = sigma2;                /* simplify loop expressions */
                A = 1.0/(2.0*sigma*sigma);
                B = (double) (1.0/(Magick2PI*sigma*sigma));
                for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
                  for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
                    kernel->values[i] -= exp(-((double)(u*u+v*v))*A)*B;
              }
            else /* limiting case - a unity (normalized Dirac) kernel */
              kernel->values[kernel->x+kernel->y*kernel->width] -= 1.0;
          }

        if ( type == LoGKernel )
          { /* Calculate a Laplacian of a Gaussian - Or Mexician Hat */
            if ( sigma > MagickEpsilon )
              { A = 1.0/(2.0*sigma*sigma);  /* simplify loop expressions */
                B = (double) (1.0/(MagickPI*sigma*sigma*sigma*sigma));
                for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
                  for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
                    { R = ((double)(u*u+v*v))*A;
                      kernel->values[i] = (1-R)*exp(-R)*B;
                    }
              }
            else /* special case - generate a unity kernel */
              { (void) ResetMagickMemory(kernel->values,0, (size_t)
                            kernel->width*kernel->height*sizeof(double));
                kernel->values[kernel->x+kernel->y*kernel->width] = 1.0;
              }
          }

        /* Note the above kernels may have been 'clipped' by a user defined
        ** radius, producing a smaller (darker) kernel.  Also for very small
        ** sigma's (> 0.1) the central value becomes larger than one, and thus
        ** producing a very bright kernel.
        **
        ** Normalization will still be needed.
        */

        /* Normalize the 2D Gaussian Kernel
        **
        ** NB: a CorrelateNormalize performs a normal Normalize if
        ** there are no negative values.
        */
        CalcKernelMetaData(kernel);  /* the other kernel meta-data */
        ScaleKernelInfo(kernel, 1.0, CorrelateNormalizeValue);

        break;
      }
    case BlurKernel:
      { double
          sigma = fabs(args->sigma),
          alpha, beta;

        if ( args->rho >= 1.0 )
          kernel->width = (size_t)args->rho*2+1;
        else
          kernel->width = GetOptimalKernelWidth1D(args->rho,sigma);
        kernel->height = 1;
        kernel->x = (ssize_t) (kernel->width-1)/2;
        kernel->y = 0;
        kernel->negative_range = kernel->positive_range = 0.0;
        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

#if 1
#define KernelRank 3
        /* Formula derived from GetBlurKernel() in "effect.c" (plus bug fix).
        ** It generates a gaussian 3 times the width, and compresses it into
        ** the expected range.  This produces a closer normalization of the
        ** resulting kernel, especially for very low sigma values.
        ** As such while wierd it is prefered.
        **
        ** I am told this method originally came from Photoshop.
        **
        ** A properly normalized curve is generated (apart from edge clipping)
        ** even though we later normalize the result (for edge clipping)
        ** to allow the correct generation of a "Difference of Blurs".
        */

        /* initialize */
        v = (ssize_t) (kernel->width*KernelRank-1)/2; /* start/end points to fit range */
        (void) ResetMagickMemory(kernel->values,0, (size_t)
                     kernel->width*kernel->height*sizeof(double));
        /* Calculate a Positive 1D Gaussian */
        if ( sigma > MagickEpsilon )
          { sigma *= KernelRank;               /* simplify loop expressions */
            alpha = 1.0/(2.0*sigma*sigma);
            beta= (double) (1.0/(MagickSQ2PI*sigma ));
            for ( u=-v; u <= v; u++) {
              kernel->values[(u+v)/KernelRank] +=
                              exp(-((double)(u*u))*alpha)*beta;
            }
          }
        else /* special case - generate a unity kernel */
          kernel->values[kernel->x+kernel->y*kernel->width] = 1.0;
#else
        /* Direct calculation without curve averaging */

        /* Calculate a Positive Gaussian */
        if ( sigma > MagickEpsilon )
          { alpha = 1.0/(2.0*sigma*sigma);    /* simplify loop expressions */
            beta = 1.0/(MagickSQ2PI*sigma);
            for ( i=0, u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
              kernel->values[i] = exp(-((double)(u*u))*alpha)*beta;
          }
        else /* special case - generate a unity kernel */
          { (void) ResetMagickMemory(kernel->values,0, (size_t)
                         kernel->width*kernel->height*sizeof(double));
            kernel->values[kernel->x+kernel->y*kernel->width] = 1.0;
          }
#endif
        /* Note the above kernel may have been 'clipped' by a user defined
        ** radius, producing a smaller (darker) kernel.  Also for very small
        ** sigma's (> 0.1) the central value becomes larger than one, and thus
        ** producing a very bright kernel.
        **
        ** Normalization will still be needed.
        */

        /* Normalize the 1D Gaussian Kernel
        **
        ** NB: a CorrelateNormalize performs a normal Normalize if
        ** there are no negative values.
        */
        CalcKernelMetaData(kernel);  /* the other kernel meta-data */
        ScaleKernelInfo(kernel, 1.0, CorrelateNormalizeValue);

        /* rotate the 1D kernel by given angle */
        RotateKernelInfo(kernel, args->xi );
        break;
      }
    case CometKernel:
      { double
          sigma = fabs(args->sigma),
          A;

        if ( args->rho < 1.0 )
          kernel->width = (GetOptimalKernelWidth1D(args->rho,sigma)-1)/2+1;
        else
          kernel->width = (size_t)args->rho;
        kernel->x = kernel->y = 0;
        kernel->height = 1;
        kernel->negative_range = kernel->positive_range = 0.0;
        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* A comet blur is half a 1D gaussian curve, so that the object is
        ** blurred in one direction only.  This may not be quite the right
        ** curve to use so may change in the future. The function must be
        ** normalised after generation, which also resolves any clipping.
        **
        ** As we are normalizing and not subtracting gaussians,
        ** there is no need for a divisor in the gaussian formula
        **
        ** It is less comples
        */
        if ( sigma > MagickEpsilon )
          {
#if 1
#define KernelRank 3
            v = (ssize_t) kernel->width*KernelRank; /* start/end points */
            (void) ResetMagickMemory(kernel->values,0, (size_t)
                          kernel->width*sizeof(double));
            sigma *= KernelRank;            /* simplify the loop expression */
            A = 1.0/(2.0*sigma*sigma);
            /* B = 1.0/(MagickSQ2PI*sigma); */
            for ( u=0; u < v; u++) {
              kernel->values[u/KernelRank] +=
                  exp(-((double)(u*u))*A);
              /*  exp(-((double)(i*i))/2.0*sigma*sigma)/(MagickSQ2PI*sigma); */
            }
            for (i=0; i < (ssize_t) kernel->width; i++)
              kernel->positive_range += kernel->values[i];
#else
            A = 1.0/(2.0*sigma*sigma);     /* simplify the loop expression */
            /* B = 1.0/(MagickSQ2PI*sigma); */
            for ( i=0; i < (ssize_t) kernel->width; i++)
              kernel->positive_range +=
                kernel->values[i] =
                  exp(-((double)(i*i))*A);
                /* exp(-((double)(i*i))/2.0*sigma*sigma)/(MagickSQ2PI*sigma); */
#endif
          }
        else /* special case - generate a unity kernel */
          { (void) ResetMagickMemory(kernel->values,0, (size_t)
                         kernel->width*kernel->height*sizeof(double));
            kernel->values[kernel->x+kernel->y*kernel->width] = 1.0;
            kernel->positive_range = 1.0;
          }

        kernel->minimum = 0.0;
        kernel->maximum = kernel->values[0];
        kernel->negative_range = 0.0;

        ScaleKernelInfo(kernel, 1.0, NormalizeValue); /* Normalize */
        RotateKernelInfo(kernel, args->xi); /* Rotate by angle */
        break;
      }

    /* Convolution Kernels - Well Known Constants */
    case LaplacianKernel:
      { switch ( (int) args->rho ) {
          case 0:
          default: /* laplacian square filter -- default */
            kernel=ParseKernelArray("3: -1,-1,-1  -1,8,-1  -1,-1,-1");
            break;
          case 1:  /* laplacian diamond filter */
            kernel=ParseKernelArray("3: 0,-1,0  -1,4,-1  0,-1,0");
            break;
          case 2:
            kernel=ParseKernelArray("3: -2,1,-2  1,4,1  -2,1,-2");
            break;
          case 3:
            kernel=ParseKernelArray("3: 1,-2,1  -2,4,-2  1,-2,1");
            break;
          case 5:   /* a 5x5 laplacian */
            kernel=ParseKernelArray(
              "5: -4,-1,0,-1,-4  -1,2,3,2,-1  0,3,4,3,0  -1,2,3,2,-1  -4,-1,0,-1,-4");
            break;
          case 7:   /* a 7x7 laplacian */
            kernel=ParseKernelArray(
              "7:-10,-5,-2,-1,-2,-5,-10 -5,0,3,4,3,0,-5 -2,3,6,7,6,3,-2 -1,4,7,8,7,4,-1 -2,3,6,7,6,3,-2 -5,0,3,4,3,0,-5 -10,-5,-2,-1,-2,-5,-10" );
            break;
          case 15:  /* a 5x5 LoG (sigma approx 1.4) */
            kernel=ParseKernelArray(
              "5: 0,0,-1,0,0  0,-1,-2,-1,0  -1,-2,16,-2,-1  0,-1,-2,-1,0  0,0,-1,0,0");
            break;
          case 19:  /* a 9x9 LoG (sigma approx 1.4) */
            /* http://www.cscjournals.org/csc/manuscript/Journals/IJIP/volume3/Issue1/IJIP-15.pdf */
            kernel=ParseKernelArray(
              "9: 0,-1,-1,-2,-2,-2,-1,-1,0  -1,-2,-4,-5,-5,-5,-4,-2,-1  -1,-4,-5,-3,-0,-3,-5,-4,-1  -2,-5,-3,12,24,12,-3,-5,-2  -2,-5,-0,24,40,24,-0,-5,-2  -2,-5,-3,12,24,12,-3,-5,-2  -1,-4,-5,-3,-0,-3,-5,-4,-1  -1,-2,-4,-5,-5,-5,-4,-2,-1  0,-1,-1,-2,-2,-2,-1,-1,0");
            break;
        }
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        break;
      }
    case SobelKernel:
#if 0
      { /* Sobel with optional 'sub-types' */
        switch ( (int) args->rho ) {
          default:
          case 0:
            kernel=ParseKernelArray("3: 1,0,-1  2,0,-2  1,0,-1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            break;
          case 1:
            kernel=ParseKernelArray("3: 1,0,-1  2,0,-2  1,0,-1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 0.25, NoValue);
            break;
          case 2:
            kernel=ParseKernelArray("3: 1,2,0  2,0,-2  0,-2,-1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 0.25, NoValue);
            break;
        }
        if ( fabs(args->sigma) > MagickEpsilon )
          /* Rotate by correctly supplied 'angle' */
          RotateKernelInfo(kernel, args->sigma);
        else if ( args->rho > 30.0 || args->rho < -30.0 )
          /* Rotate by out of bounds 'type' */
          RotateKernelInfo(kernel, args->rho);
        break;
      }
#else
      { /* Simple Sobel Kernel */
        kernel=ParseKernelArray("3: 1,0,-1  2,0,-2  1,0,-1");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        RotateKernelInfo(kernel, args->rho);
        break;
      }
#endif
    case RobertsKernel:
      {
        kernel=ParseKernelArray("3: 0,0,0  1,-1,0  0,0,0");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        RotateKernelInfo(kernel, args->rho);
        break;
      }
    case PrewittKernel:
      {
        kernel=ParseKernelArray("3: 1,0,-1  1,0,-1  1,0,-1");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        RotateKernelInfo(kernel, args->rho);
        break;
      }
    case CompassKernel:
      {
        kernel=ParseKernelArray("3: 1,1,-1  1,-2,-1  1,1,-1");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        RotateKernelInfo(kernel, args->rho);
        break;
      }
    case KirschKernel:
      {
        kernel=ParseKernelArray("3: 5,-3,-3  5,0,-3  5,-3,-3");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        RotateKernelInfo(kernel, args->rho);
        break;
      }
    case FreiChenKernel:
      /* Direction is set to be left to right positive */
      /* http://www.math.tau.ac.il/~turkel/notes/edge_detectors.pdf -- RIGHT? */
      /* http://ltswww.epfl.ch/~courstiv/exos_labos/sol3.pdf -- WRONG? */
      { switch ( (int) args->rho ) {
          default:
          case 0:
            kernel=ParseKernelArray("3: 1,0,-1  2,0,-2  1,0,-1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            kernel->values[3] = +MagickSQ2;
            kernel->values[5] = -MagickSQ2;
            CalcKernelMetaData(kernel);     /* recalculate meta-data */
            break;
          case 2:
            kernel=ParseKernelArray("3: 1,2,0  2,0,-2  0,-2,-1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            kernel->values[1] = kernel->values[3] = +MagickSQ2;
            kernel->values[5] = kernel->values[7] = -MagickSQ2;
            CalcKernelMetaData(kernel);     /* recalculate meta-data */
            ScaleKernelInfo(kernel, 1.0/2.0*MagickSQ2, NoValue);
            break;
          case 10:
            kernel=AcquireKernelInfo("FreiChen:11;FreiChen:12;FreiChen:13;FreiChen:14;FreiChen:15;FreiChen:16;FreiChen:17;FreiChen:18;FreiChen:19");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            break;
          case 1:
          case 11:
            kernel=ParseKernelArray("3: 1,0,-1  2,0,-2  1,0,-1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            kernel->values[3] = +MagickSQ2;
            kernel->values[5] = -MagickSQ2;
            CalcKernelMetaData(kernel);     /* recalculate meta-data */
            ScaleKernelInfo(kernel, (double) (1.0/2.0*MagickSQ2), NoValue);
            break;
          case 12:
            kernel=ParseKernelArray("3: 1,2,1  0,0,0  1,2,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            kernel->values[1] = +MagickSQ2;
            kernel->values[7] = +MagickSQ2;
            CalcKernelMetaData(kernel);
            ScaleKernelInfo(kernel, (double) (1.0/2.0*MagickSQ2), NoValue);
            break;
          case 13:
            kernel=ParseKernelArray("3: 2,-1,0  -1,0,1  0,1,-2");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            kernel->values[0] = +MagickSQ2;
            kernel->values[8] = -MagickSQ2;
            CalcKernelMetaData(kernel);
            ScaleKernelInfo(kernel, (double) (1.0/2.0*MagickSQ2), NoValue);
            break;
          case 14:
            kernel=ParseKernelArray("3: 0,1,-2  -1,0,1  2,-1,0");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            kernel->values[2] = -MagickSQ2;
            kernel->values[6] = +MagickSQ2;
            CalcKernelMetaData(kernel);
            ScaleKernelInfo(kernel, (double) (1.0/2.0*MagickSQ2), NoValue);
            break;
          case 15:
            kernel=ParseKernelArray("3: 0,-1,0  1,0,1  0,-1,0");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 1.0/2.0, NoValue);
            break;
          case 16:
            kernel=ParseKernelArray("3: 1,0,-1  0,0,0  -1,0,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 1.0/2.0, NoValue);
            break;
          case 17:
            kernel=ParseKernelArray("3: 1,-2,1  -2,4,-2  -1,-2,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 1.0/6.0, NoValue);
            break;
          case 18:
            kernel=ParseKernelArray("3: -2,1,-2  1,4,1  -2,1,-2");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 1.0/6.0, NoValue);
            break;
          case 19:
            kernel=ParseKernelArray("3: 1,1,1  1,1,1  1,1,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ScaleKernelInfo(kernel, 1.0/3.0, NoValue);
            break;
        }
        if ( fabs(args->sigma) > MagickEpsilon )
          /* Rotate by correctly supplied 'angle' */
          RotateKernelInfo(kernel, args->sigma);
        else if ( args->rho > 30.0 || args->rho < -30.0 )
          /* Rotate by out of bounds 'type' */
          RotateKernelInfo(kernel, args->rho);
        break;
      }

    /* Boolean Kernels */
    case DiamondKernel:
      {
        if (args->rho < 1.0)
          kernel->width = kernel->height = 3;  /* default radius = 1 */
        else
          kernel->width = kernel->height = ((size_t)args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* set all kernel values within diamond area to scale given */
        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            if ( (labs((long) u)+labs((long) v)) <= (long) kernel->x)
              kernel->positive_range += kernel->values[i] = args->sigma;
            else
              kernel->values[i] = nan;
        kernel->minimum = kernel->maximum = args->sigma;   /* a flat shape */
        break;
      }
    case SquareKernel:
    case RectangleKernel:
      { double
          scale;
        if ( type == SquareKernel )
          {
            if (args->rho < 1.0)
              kernel->width = kernel->height = 3;  /* default radius = 1 */
            else
              kernel->width = kernel->height = (size_t) (2*args->rho+1);
            kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;
            scale = args->sigma;
          }
        else {
            /* NOTE: user defaults set in "AcquireKernelInfo()" */
            if ( args->rho < 1.0 || args->sigma < 1.0 )
              return(DestroyKernelInfo(kernel));    /* invalid args given */
            kernel->width = (size_t)args->rho;
            kernel->height = (size_t)args->sigma;
            if ( args->xi  < 0.0 || args->xi  > (double)kernel->width ||
                 args->psi < 0.0 || args->psi > (double)kernel->height )
              return(DestroyKernelInfo(kernel));    /* invalid args given */
            kernel->x = (ssize_t) args->xi;
            kernel->y = (ssize_t) args->psi;
            scale = 1.0;
          }
        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* set all kernel values to scale given */
        u=(ssize_t) (kernel->width*kernel->height);
        for ( i=0; i < u; i++)
            kernel->values[i] = scale;
        kernel->minimum = kernel->maximum = scale;   /* a flat shape */
        kernel->positive_range = scale*u;
        break;
      }
    case DiskKernel:
      {
        ssize_t
         limit = (ssize_t)(args->rho*args->rho);

        if (args->rho < 0.4)           /* default radius approx 3.5 */
          kernel->width = kernel->height = 7L, limit = 10L;
        else
           kernel->width = kernel->height = (size_t)fabs(args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* set all kernel values within disk area to scale given */
        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            if ((u*u+v*v) <= limit)
              kernel->positive_range += kernel->values[i] = args->sigma;
            else
              kernel->values[i] = nan;
        kernel->minimum = kernel->maximum = args->sigma;   /* a flat shape */
        break;
      }
    case PlusKernel:
      {
        if (args->rho < 1.0)
          kernel->width = kernel->height = 5;  /* default radius 2 */
        else
           kernel->width = kernel->height = ((size_t)args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* set all kernel values along axises to given scale */
        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            kernel->values[i] = (u == 0 || v == 0) ? args->sigma : nan;
        kernel->minimum = kernel->maximum = args->sigma;   /* a flat shape */
        kernel->positive_range = args->sigma*(kernel->width*2.0 - 1.0);
        break;
      }
    case CrossKernel:
      {
        if (args->rho < 1.0)
          kernel->width = kernel->height = 5;  /* default radius 2 */
        else
           kernel->width = kernel->height = ((size_t)args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* set all kernel values along axises to given scale */
        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            kernel->values[i] = (u == v || u == -v) ? args->sigma : nan;
        kernel->minimum = kernel->maximum = args->sigma;   /* a flat shape */
        kernel->positive_range = args->sigma*(kernel->width*2.0 - 1.0);
        break;
      }
    /* HitAndMiss Kernels */
    case RingKernel:
    case PeaksKernel:
      {
        ssize_t
          limit1,
          limit2,
          scale;

        if (args->rho < args->sigma)
          {
            kernel->width = ((size_t)args->sigma)*2+1;
            limit1 = (ssize_t)(args->rho*args->rho);
            limit2 = (ssize_t)(args->sigma*args->sigma);
          }
        else
          {
            kernel->width = ((size_t)args->rho)*2+1;
            limit1 = (ssize_t)(args->sigma*args->sigma);
            limit2 = (ssize_t)(args->rho*args->rho);
          }
        if ( limit2 <= 0 )
          kernel->width = 7L, limit1 = 7L, limit2 = 11L;

        kernel->height = kernel->width;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;
        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        /* set a ring of points of 'scale' ( 0.0 for PeaksKernel ) */
        scale = (ssize_t) (( type == PeaksKernel) ? 0.0 : args->xi);
        for ( i=0, v= -kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            { ssize_t radius=u*u+v*v;
              if (limit1 < radius && radius <= limit2)
                kernel->positive_range += kernel->values[i] = (double) scale;
              else
                kernel->values[i] = nan;
            }
        kernel->minimum = kernel->maximum = (double) scale;
        if ( type == PeaksKernel ) {
          /* set the central point in the middle */
          kernel->values[kernel->x+kernel->y*kernel->width] = 1.0;
          kernel->positive_range = 1.0;
          kernel->maximum = 1.0;
        }
        break;
      }
    case EdgesKernel:
      {
        kernel=ParseKernelArray("3: 0,0,0  -,1,-  1,1,1");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        ExpandMirrorKernelInfo(kernel); /* mirror expansion of other kernels */
        break;
      }
    case CornersKernel:
      {
        kernel=ParseKernelArray("3: 0,0,-  0,1,1  -,1,-");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        ExpandRotateKernelInfo(kernel, 90.0); /* Expand 90 degree rotations */
        break;
      }
    case ThinDiagonalsKernel:
      {
        switch ( (int) args->rho ) {
          case 0:
          default:
            { KernelInfo
                *new_kernel;
              kernel=ParseKernelArray("3: 0,0,0  0,1,1  1,1,-");
              if (kernel == (KernelInfo *) NULL)
                return(kernel);
              kernel->type = type;
              new_kernel=ParseKernelArray("3: 0,0,1  0,1,1  0,1,-");
              if (new_kernel == (KernelInfo *) NULL)
                return(DestroyKernelInfo(kernel));
              new_kernel->type = type;
              LastKernelInfo(kernel)->next = new_kernel;
              ExpandMirrorKernelInfo(kernel);
              break;
            }
          case 1:
            kernel=ParseKernelArray("3: 0,0,0  0,1,1  1,1,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
          case 2:
            kernel=ParseKernelArray("3: 0,0,1  0,1,1  0,1,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
        }
        break;
      }
    case LineEndsKernel:
      { /* Kernels for finding the end of thin lines */
        switch ( (int) args->rho ) {
          case 0:
          default:
            /* set of kernels to find all end of lines */
            kernel=AcquireKernelInfo("LineEnds:1>;LineEnds:2>");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            break;
          case 1:
            /* kernel for 4-connected line ends - no rotation */
            kernel=ParseKernelArray("3: 0,0,-  0,1,1  0,0,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
         case 2:
            /* kernel to add for 8-connected lines - no rotation */
            kernel=ParseKernelArray("3: 0,0,0  0,1,0  0,0,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
         case 3:
            /* kernel to add for orthogonal line ends - does not find corners */
            kernel=ParseKernelArray("3: 0,0,0  0,1,1  0,0,0");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
         case 4:
            /* traditional line end - fails on last T end */
            kernel=ParseKernelArray("3: 0,0,0  0,1,-  0,0,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
        }
        break;
      }
    case LineJunctionsKernel:
      { /* kernels for finding the junctions of multiple lines */
        switch ( (int) args->rho ) {
          case 0:
          default:
            /* set of kernels to find all line junctions */
            kernel=AcquireKernelInfo("LineJunctions:1@;LineJunctions:2>");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            break;
          case 1:
            /* Y Junction */
            kernel=ParseKernelArray("3: 1,-,1  -,1,-  -,1,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
          case 2:
            /* Diagonal T Junctions */
            kernel=ParseKernelArray("3: 1,-,-  -,1,-  1,-,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
          case 3:
            /* Orthogonal T Junctions */
            kernel=ParseKernelArray("3: -,-,-  1,1,1  -,1,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
          case 4:
            /* Diagonal X Junctions */
            kernel=ParseKernelArray("3: 1,-,1  -,1,-  1,-,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
          case 5:
            /* Orthogonal X Junctions - minimal diamond kernel */
            kernel=ParseKernelArray("3: -,1,-  1,1,1  -,1,-");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            RotateKernelInfo(kernel, args->sigma);
            break;
        }
        break;
      }
    case RidgesKernel:
      { /* Ridges - Ridge finding kernels */
        KernelInfo
          *new_kernel;
        switch ( (int) args->rho ) {
          case 1:
          default:
            kernel=ParseKernelArray("3x1:0,1,0");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ExpandRotateKernelInfo(kernel, 90.0); /* 2 rotated kernels (symmetrical) */
            break;
          case 2:
            kernel=ParseKernelArray("4x1:0,1,1,0");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ExpandRotateKernelInfo(kernel, 90.0); /* 4 rotated kernels */

            /* Kernels to find a stepped 'thick' line, 4 rotates + mirrors */
            /* Unfortunatally we can not yet rotate a non-square kernel */
            /* But then we can't flip a non-symetrical kernel either */
            new_kernel=ParseKernelArray("4x3+1+1:0,1,1,- -,1,1,- -,1,1,0");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("4x3+2+1:0,1,1,- -,1,1,- -,1,1,0");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("4x3+1+1:-,1,1,0 -,1,1,- 0,1,1,-");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("4x3+2+1:-,1,1,0 -,1,1,- 0,1,1,-");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("3x4+1+1:0,-,- 1,1,1 1,1,1 -,-,0");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("3x4+1+2:0,-,- 1,1,1 1,1,1 -,-,0");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("3x4+1+1:-,-,0 1,1,1 1,1,1 0,-,-");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            new_kernel=ParseKernelArray("3x4+1+2:-,-,0 1,1,1 1,1,1 0,-,-");
            if (new_kernel == (KernelInfo *) NULL)
              return(DestroyKernelInfo(kernel));
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            break;
        }
        break;
      }
    case ConvexHullKernel:
      {
        KernelInfo
          *new_kernel;
        /* first set of 8 kernels */
        kernel=ParseKernelArray("3: 1,1,-  1,0,-  1,-,0");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = type;
        ExpandRotateKernelInfo(kernel, 90.0);
        /* append the mirror versions too - no flip function yet */
        new_kernel=ParseKernelArray("3: 1,1,1  1,0,-  -,-,0");
        if (new_kernel == (KernelInfo *) NULL)
          return(DestroyKernelInfo(kernel));
        new_kernel->type = type;
        ExpandRotateKernelInfo(new_kernel, 90.0);
        LastKernelInfo(kernel)->next = new_kernel;
        break;
      }
    case SkeletonKernel:
      {
        KernelInfo
          *new_kernel;
        switch ( (int) args->rho ) {
          case 1:
          default:
            /* Traditional Skeleton...
            ** A cyclically rotated single kernel
            */
            kernel=ParseKernelArray("3: 0,0,0  -,1,-  1,1,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            ExpandRotateKernelInfo(kernel, 45.0); /* 8 rotations */
            break;
          case 2:
            /* HIPR Variation of the cyclic skeleton
            ** Corners of the traditional method made more forgiving,
            ** but the retain the same cyclic order.
            */
            kernel=ParseKernelArray("3: 0,0,0  -,1,-  1,1,1");
            if (kernel == (KernelInfo *) NULL)
              return(kernel);
            kernel->type = type;
            new_kernel=ParseKernelArray("3: -,0,0  1,1,0  -,1,-");
            if (new_kernel == (KernelInfo *) NULL)
              return(new_kernel);
            new_kernel->type = type;
            LastKernelInfo(kernel)->next = new_kernel;
            ExpandRotateKernelInfo(kernel, 90.0); /* 4 rotations of the 2 kernels */
            break;
        }
        break;
      }
    /* Distance Measuring Kernels */
    case ChebyshevKernel:
      {
        if (args->rho < 1.0)
          kernel->width = kernel->height = 3;  /* default radius = 1 */
        else
          kernel->width = kernel->height = ((size_t)args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            kernel->positive_range += ( kernel->values[i] =
                 args->sigma*((labs((long) u)>labs((long) v)) ? labs((long) u) : labs((long) v)) );
        kernel->maximum = kernel->values[0];
        break;
      }
    case ManhattanKernel:
      {
        if (args->rho < 1.0)
          kernel->width = kernel->height = 3;  /* default radius = 1 */
        else
           kernel->width = kernel->height = ((size_t)args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            kernel->positive_range += ( kernel->values[i] =
                 args->sigma*(labs((long) u)+labs((long) v)) );
        kernel->maximum = kernel->values[0];
        break;
      }
    case EuclideanKernel:
      {
        if (args->rho < 1.0)
          kernel->width = kernel->height = 3;  /* default radius = 1 */
        else
           kernel->width = kernel->height = ((size_t)args->rho)*2+1;
        kernel->x = kernel->y = (ssize_t) (kernel->width-1)/2;

        kernel->values=(double *) AcquireQuantumMemory(kernel->width,
                              kernel->height*sizeof(double));
        if (kernel->values == (double *) NULL)
          return(DestroyKernelInfo(kernel));

        for ( i=0, v=-kernel->y; v <= (ssize_t)kernel->y; v++)
          for ( u=-kernel->x; u <= (ssize_t)kernel->x; u++, i++)
            kernel->positive_range += ( kernel->values[i] =
                 args->sigma*sqrt((double)(u*u+v*v)) );
        kernel->maximum = kernel->values[0];
        break;
      }
    case UnityKernel:
    default:
      {
        /* Unity or No-Op Kernel - Basically just a single pixel on its own */
        kernel=ParseKernelArray("1:1");
        if (kernel == (KernelInfo *) NULL)
          return(kernel);
        kernel->type = ( type == UnityKernel ) ? UnityKernel : UndefinedKernel;
        break;
      }
      break;
  }

  return(kernel);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     C l o n e K e r n e l I n f o                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CloneKernelInfo() creates a new clone of the given Kernel List so that its
%  can be modified without effecting the original.  The cloned kernel should
%  be destroyed using DestoryKernelInfo() when no longer needed.
%
%  The format of the CloneKernelInfo method is:
%
%      KernelInfo *CloneKernelInfo(const KernelInfo *kernel)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel to be cloned
%
*/
MagickExport KernelInfo *CloneKernelInfo(const KernelInfo *kernel)
{
  register ssize_t
    i;

  KernelInfo
    *new_kernel;

  assert(kernel != (KernelInfo *) NULL);
  new_kernel=(KernelInfo *) AcquireMagickMemory(sizeof(*kernel));
  if (new_kernel == (KernelInfo *) NULL)
    return(new_kernel);
  *new_kernel=(*kernel); /* copy values in structure */

  /* replace the values with a copy of the values */
  new_kernel->values=(double *) AcquireQuantumMemory(kernel->width,
    kernel->height*sizeof(double));
  if (new_kernel->values == (double *) NULL)
    return(DestroyKernelInfo(new_kernel));
  for (i=0; i < (ssize_t) (kernel->width*kernel->height); i++)
    new_kernel->values[i]=kernel->values[i];

  /* Also clone the next kernel in the kernel list */
  if ( kernel->next != (KernelInfo *) NULL ) {
    new_kernel->next = CloneKernelInfo(kernel->next);
    if ( new_kernel->next == (KernelInfo *) NULL )
      return(DestroyKernelInfo(new_kernel));
  }

  return(new_kernel);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     D e s t r o y K e r n e l I n f o                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DestroyKernelInfo() frees the memory used by a Convolution/Morphology
%  kernel.
%
%  The format of the DestroyKernelInfo method is:
%
%      KernelInfo *DestroyKernelInfo(KernelInfo *kernel)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel to be destroyed
%
*/
MagickExport KernelInfo *DestroyKernelInfo(KernelInfo *kernel)
{
  assert(kernel != (KernelInfo *) NULL);

  if ( kernel->next != (KernelInfo *) NULL )
    kernel->next = DestroyKernelInfo(kernel->next);

  kernel->values = (double *)RelinquishMagickMemory(kernel->values);
  kernel = (KernelInfo *) RelinquishMagickMemory(kernel);
  return(kernel);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+     E x p a n d M i r r o r K e r n e l I n f o                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ExpandMirrorKernelInfo() takes a single kernel, and expands it into a
%  sequence of 90-degree rotated kernels but providing a reflected 180
%  rotatation, before the -/+ 90-degree rotations.
%
%  This special rotation order produces a better, more symetrical thinning of
%  objects.
%
%  The format of the ExpandMirrorKernelInfo method is:
%
%      void ExpandMirrorKernelInfo(KernelInfo *kernel)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
% This function is only internel to this module, as it is not finalized,
% especially with regard to non-orthogonal angles, and rotation of larger
% 2D kernels.
*/

#if 0
static void FlopKernelInfo(KernelInfo *kernel)
    { /* Do a Flop by reversing each row. */
      size_t
        y;
      register ssize_t
        x,r;
      register double
        *k,t;

      for ( y=0, k=kernel->values; y < kernel->height; y++, k+=kernel->width)
        for ( x=0, r=kernel->width-1; x<kernel->width/2; x++, r--)
          t=k[x],  k[x]=k[r],  k[r]=t;

      kernel->x = kernel->width - kernel->x - 1;
      angle = fmod(angle+180.0, 360.0);
    }
#endif

static void ExpandMirrorKernelInfo(KernelInfo *kernel)
{
  KernelInfo
    *clone,
    *last;

  last = kernel;

  clone = CloneKernelInfo(last);
  RotateKernelInfo(clone, 180);   /* flip */
  LastKernelInfo(last)->next = clone;
  last = clone;

  clone = CloneKernelInfo(last);
  RotateKernelInfo(clone, 90);   /* transpose */
  LastKernelInfo(last)->next = clone;
  last = clone;

  clone = CloneKernelInfo(last);
  RotateKernelInfo(clone, 180);  /* flop */
  LastKernelInfo(last)->next = clone;

  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+     E x p a n d R o t a t e K e r n e l I n f o                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ExpandRotateKernelInfo() takes a kernel list, and expands it by rotating
%  incrementally by the angle given, until the first kernel repeats.
%
%  WARNING: 45 degree rotations only works for 3x3 kernels.
%  While 90 degree roatations only works for linear and square kernels
%
%  The format of the ExpandRotateKernelInfo method is:
%
%      void ExpandRotateKernelInfo(KernelInfo *kernel, double angle)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
%    o angle: angle to rotate in degrees
%
% This function is only internel to this module, as it is not finalized,
% especially with regard to non-orthogonal angles, and rotation of larger
% 2D kernels.
*/

/* Internal Routine - Return true if two kernels are the same */
static MagickBooleanType SameKernelInfo(const KernelInfo *kernel1,
     const KernelInfo *kernel2)
{
  register size_t
    i;

  /* check size and origin location */
  if (    kernel1->width != kernel2->width
       || kernel1->height != kernel2->height
       || kernel1->x != kernel2->x
       || kernel1->y != kernel2->y )
    return MagickFalse;

  /* check actual kernel values */
  for (i=0; i < (kernel1->width*kernel1->height); i++) {
    /* Test for Nan equivelence */
    if ( IsNan(kernel1->values[i]) && !IsNan(kernel2->values[i]) )
      return MagickFalse;
    if ( IsNan(kernel2->values[i]) && !IsNan(kernel1->values[i]) )
      return MagickFalse;
    /* Test actual values are equivelent */
    if ( fabs(kernel1->values[i] - kernel2->values[i]) > MagickEpsilon )
      return MagickFalse;
  }

  return MagickTrue;
}

static void ExpandRotateKernelInfo(KernelInfo *kernel, const double angle)
{
  KernelInfo
    *clone,
    *last;

  last = kernel;
  while(1) {
    clone = CloneKernelInfo(last);
    RotateKernelInfo(clone, angle);
    if ( SameKernelInfo(kernel, clone) == MagickTrue )
      break;
    LastKernelInfo(last)->next = clone;
    last = clone;
  }
  clone = DestroyKernelInfo(clone); /* kernel has repeated - junk the clone */
  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+     C a l c M e t a K e r n a l I n f o                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CalcKernelMetaData() recalculate the KernelInfo meta-data of this kernel only,
%  using the kernel values.  This should only ne used if it is not posible to
%  calculate that meta-data in some easier way.
%
%  It is important that the meta-data is correct before ScaleKernelInfo() is
%  used to perform kernel normalization.
%
%  The format of the CalcKernelMetaData method is:
%
%      void CalcKernelMetaData(KernelInfo *kernel, const double scale )
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel to modify
%
%  WARNING: Minimum and Maximum values are assumed to include zero, even if
%  zero is not part of the kernel (as in Gaussian Derived kernels). This
%  however is not true for flat-shaped morphological kernels.
%
%  WARNING: Only the specific kernel pointed to is modified, not a list of
%  multiple kernels.
%
% This is an internal function and not expected to be useful outside this
% module.  This could change however.
*/
static void CalcKernelMetaData(KernelInfo *kernel)
{
  register size_t
    i;

  kernel->minimum = kernel->maximum = 0.0;
  kernel->negative_range = kernel->positive_range = 0.0;
  for (i=0; i < (kernel->width*kernel->height); i++)
    {
      if ( fabs(kernel->values[i]) < MagickEpsilon )
        kernel->values[i] = 0.0;
      ( kernel->values[i] < 0)
          ?  ( kernel->negative_range += kernel->values[i] )
          :  ( kernel->positive_range += kernel->values[i] );
      Minimize(kernel->minimum, kernel->values[i]);
      Maximize(kernel->maximum, kernel->values[i]);
    }

  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     M o r p h o l o g y A p p l y                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MorphologyApply() applies a morphological method, multiple times using
%  a list of multiple kernels.
%
%  It is basically equivelent to as MorphologyImageChannel() (see below) but
%  without any user controls.  This allows internel programs to use this
%  function, to actually perform a specific task without posible interference
%  by any API user supplied settings.
%
%  It is MorphologyImageChannel() task to extract any such user controls, and
%  pass them to this function for processing.
%
%  More specifically kernels are not normalized/scaled/blended by the
%  'convolve:scale' Image Artifact (setting), nor is the convolve bias
%  (-bias setting or image->bias) loooked at, but must be supplied from the
%  function arguments.
%
%  The format of the MorphologyApply method is:
%
%      Image *MorphologyApply(const Image *image,MorphologyMethod method,
%        const ssize_t iterations,const KernelInfo *kernel,
%        const CompositeMethod compose, const double bias,
%        ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o image: the source image
%
%    o method: the morphology method to be applied.
%
%    o iterations: apply the operation this many times (or no change).
%                  A value of -1 means loop until no change found.
%                  How this is applied may depend on the morphology method.
%                  Typically this is a value of 1.
%
%    o channel: the channel type.
%
%    o kernel: An array of double representing the morphology kernel.
%
%    o compose: How to handle or merge multi-kernel results.
%          If 'UndefinedCompositeOp' use default for the Morphology method.
%          If 'NoCompositeOp' force image to be re-iterated by each kernel.
%          Otherwise merge the results using the compose method given.
%
%    o bias: Convolution Output Bias.
%
%    o exception: return any errors or warnings in this structure.
%
*/


/* Apply a Morphology Primative to an image using the given kernel.
** Two pre-created images must be provided, no image is created.
** It returns the number of pixels that changed betwene the images
** for convergence determination.
*/
static size_t MorphologyPrimitive(const Image *image, Image
     *result_image, const MorphologyMethod method, const ChannelType channel,
     const KernelInfo *kernel,const double bias,ExceptionInfo *exception)
{
#define MorphologyTag  "Morphology/Image"

  CacheView
    *p_view,
    *q_view;

  ssize_t
    y, offx, offy,
    changed;

  MagickBooleanType
    status;

  MagickOffsetType
    progress;

  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  assert(result_image != (Image *) NULL);
  assert(result_image->signature == MagickSignature);
  assert(kernel != (KernelInfo *) NULL);
  assert(kernel->signature == MagickSignature);
  assert(exception != (ExceptionInfo *) NULL);
  assert(exception->signature == MagickSignature);

  status=MagickTrue;
  changed=0;
  progress=0;

  p_view=AcquireCacheView(image);
  q_view=AcquireCacheView(result_image);

  /* Some methods (including convolve) needs use a reflected kernel.
   * Adjust 'origin' offsets to loop though kernel as a reflection.
   */
  offx = kernel->x;
  offy = kernel->y;
  switch(method) {
    case ConvolveMorphology:
    case DilateMorphology:
    case DilateIntensityMorphology:
    case DistanceMorphology:
      /* kernel needs to used with reflection about origin */
      offx = (ssize_t) kernel->width-offx-1;
      offy = (ssize_t) kernel->height-offy-1;
      break;
    case ErodeMorphology:
    case ErodeIntensityMorphology:
    case HitAndMissMorphology:
    case ThinningMorphology:
    case ThickenMorphology:
      /* kernel is used as is, without reflection */
      break;
    default:
      assert("Not a Primitive Morphology Method" != (char *) NULL);
      break;
  }


  if ( method == ConvolveMorphology && kernel->width == 1 )
  { /* Special handling (for speed) of vertical (blur) kernels.
    ** This performs its handling in columns rather than in rows.
    ** This is only done fo convolve as it is the only method that
    ** generates very large 1-D vertical kernels (such as a 'BlurKernel')
    **
    ** Timing tests (on single CPU laptop)
    ** Using a vertical 1-d Blue with normal row-by-row (below)
    **   time convert logo: -morphology Convolve Blur:0x10+90 null:
    **      0.807u
    ** Using this column method
    **   time convert logo: -morphology Convolve Blur:0x10+90 null:
    **      0.620u
    **
    ** Anthony Thyssen, 14 June 2010
    */
    register ssize_t
      x;

#if defined(MAGICKCORE_OPENMP_SUPPORT)
#pragma omp parallel for schedule(dynamic,4) shared(progress,status)
#endif
    for (x=0; x < (ssize_t) image->columns; x++)
    {
      register const PixelPacket
        *restrict p;

      register const IndexPacket
        *restrict p_indexes;

      register PixelPacket
        *restrict q;

      register IndexPacket
        *restrict q_indexes;

      register ssize_t
        y;

      ssize_t
        r;

      if (status == MagickFalse)
        continue;
      p=GetCacheViewVirtualPixels(p_view, x,  -offy,1,
          image->rows+kernel->height, exception);
      q=GetCacheViewAuthenticPixels(q_view,x,0,1,result_image->rows,exception);
      if ((p == (const PixelPacket *) NULL) || (q == (PixelPacket *) NULL))
        {
          status=MagickFalse;
          continue;
        }
      p_indexes=GetCacheViewVirtualIndexQueue(p_view);
      q_indexes=GetCacheViewAuthenticIndexQueue(q_view);
      r = offy;  /* offset to the origin pixel in 'p' */

      for (y=0; y < (ssize_t) image->rows; y++)
      {
        register ssize_t
          v;

        register const double
          *restrict k;

        register const PixelPacket
          *restrict k_pixels;

        register const IndexPacket
          *restrict k_indexes;

        MagickPixelPacket
          result;

        /* Copy input image to the output image for unused channels
        * This removes need for 'cloning' a new image every iteration
        */
        *q = p[r];
        if (image->colorspace == CMYKColorspace)
          q_indexes[y] = p_indexes[r];

        /* Set the bias of the weighted average output */
        result.red     =
        result.green   =
        result.blue    =
        result.opacity =
        result.index   = bias;


        /* Weighted Average of pixels using reflected kernel
        **
        ** NOTE for correct working of this operation for asymetrical
        ** kernels, the kernel needs to be applied in its reflected form.
        ** That is its values needs to be reversed.
        */
        k = &kernel->values[ kernel->height-1 ];
        k_pixels = p;
        k_indexes = p_indexes;
        if ( ((channel & SyncChannels) == 0 ) ||
                             (image->matte == MagickFalse) )
          { /* No 'Sync' involved.
            ** Convolution is simple greyscale channel operation
            */
            for (v=0; v < (ssize_t) kernel->height; v++) {
              if ( IsNan(*k) ) continue;
              result.red     += (*k)*k_pixels->red;
              result.green   += (*k)*k_pixels->green;
              result.blue    += (*k)*k_pixels->blue;
              result.opacity += (*k)*k_pixels->opacity;
              if ( image->colorspace == CMYKColorspace)
                result.index += (*k)*(*k_indexes);
              k--;
              k_pixels++;
              k_indexes++;
            }
            if ((channel & RedChannel) != 0)
              q->red = ClampToQuantum(result.red);
            if ((channel & GreenChannel) != 0)
              q->green = ClampToQuantum(result.green);
            if ((channel & BlueChannel) != 0)
              q->blue = ClampToQuantum(result.blue);
            if ((channel & OpacityChannel) != 0
                && image->matte == MagickTrue )
              q->opacity = ClampToQuantum(result.opacity);
            if ((channel & IndexChannel) != 0
                && image->colorspace == CMYKColorspace)
              q_indexes[x] = ClampToQuantum(result.index);
          }
        else
          { /* Channel 'Sync' Flag, and Alpha Channel enabled.
            ** Weight the color channels with Alpha Channel so that
            ** transparent pixels are not part of the results.
            */
            MagickRealType
              alpha,  /* alpha weighting of colors : kernel*alpha  */
              gamma;  /* divisor, sum of color weighting values */

            gamma=0.0;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              if ( IsNan(*k) ) continue;
              alpha=(*k)*(QuantumScale*(QuantumRange-k_pixels->opacity));
              gamma += alpha;
              result.red     += alpha*k_pixels->red;
              result.green   += alpha*k_pixels->green;
              result.blue    += alpha*k_pixels->blue;
              result.opacity += (*k)*k_pixels->opacity;
              if ( image->colorspace == CMYKColorspace)
                result.index += alpha*(*k_indexes);
              k--;
              k_pixels++;
              k_indexes++;
            }
            /* Sync'ed channels, all channels are modified */
            gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
            q->red = ClampToQuantum(gamma*result.red);
            q->green = ClampToQuantum(gamma*result.green);
            q->blue = ClampToQuantum(gamma*result.blue);
            q->opacity = ClampToQuantum(result.opacity);
            if (image->colorspace == CMYKColorspace)
              q_indexes[x] = ClampToQuantum(gamma*result.index);
          }

        /* Count up changed pixels */
        if (   ( p[r].red != q->red )
            || ( p[r].green != q->green )
            || ( p[r].blue != q->blue )
            || ( p[r].opacity != q->opacity )
            || ( image->colorspace == CMYKColorspace &&
                    p_indexes[r] != q_indexes[x] ) )
          changed++;  /* The pixel was changed in some way! */
        p++;
        q++;
      } /* y */
      if ( SyncCacheViewAuthenticPixels(q_view,exception) == MagickFalse)
        status=MagickFalse;
      if (image->progress_monitor != (MagickProgressMonitor) NULL)
        {
          MagickBooleanType
            proceed;

#if defined(MAGICKCORE_OPENMP_SUPPORT)
  #pragma omp critical (MagickCore_MorphologyImage)
#endif
          proceed=SetImageProgress(image,MorphologyTag,progress++,image->rows);
          if (proceed == MagickFalse)
            status=MagickFalse;
        }
    } /* x */
    result_image->type=image->type;
    q_view=DestroyCacheView(q_view);
    p_view=DestroyCacheView(p_view);
    return(status ? (size_t) changed : 0);
  }

  /*
  ** Normal handling of horizontal or rectangular kernels (row by row)
  */
#if defined(MAGICKCORE_OPENMP_SUPPORT)
  #pragma omp parallel for schedule(dynamic,4) shared(progress,status)
#endif
  for (y=0; y < (ssize_t) image->rows; y++)
  {
    register const PixelPacket
      *restrict p;

    register const IndexPacket
      *restrict p_indexes;

    register PixelPacket
      *restrict q;

    register IndexPacket
      *restrict q_indexes;

    register ssize_t
      x;

    size_t
      r;

    if (status == MagickFalse)
      continue;
    p=GetCacheViewVirtualPixels(p_view, -offx,  y-offy,
         image->columns+kernel->width,  kernel->height,  exception);
    q=GetCacheViewAuthenticPixels(q_view,0,y,result_image->columns,1,
         exception);
    if ((p == (const PixelPacket *) NULL) || (q == (PixelPacket *) NULL))
      {
        status=MagickFalse;
        continue;
      }
    p_indexes=GetCacheViewVirtualIndexQueue(p_view);
    q_indexes=GetCacheViewAuthenticIndexQueue(q_view);
    r = (image->columns+kernel->width)*offy+offx; /* offset to origin in 'p' */

    for (x=0; x < (ssize_t) image->columns; x++)
    {
       ssize_t
        v;

      register ssize_t
        u;

      register const double
        *restrict k;

      register const PixelPacket
        *restrict k_pixels;

      register const IndexPacket
        *restrict k_indexes;

      MagickPixelPacket
        result,
        min,
        max;

      /* Copy input image to the output image for unused channels
       * This removes need for 'cloning' a new image every iteration
       */
      *q = p[r];
      if (image->colorspace == CMYKColorspace)
        q_indexes[x] = p_indexes[r];

      /* Defaults */
      min.red     =
      min.green   =
      min.blue    =
      min.opacity =
      min.index   = (MagickRealType) QuantumRange;
      max.red     =
      max.green   =
      max.blue    =
      max.opacity =
      max.index   = (MagickRealType) 0;
      /* default result is the original pixel value */
      result.red     = (MagickRealType) p[r].red;
      result.green   = (MagickRealType) p[r].green;
      result.blue    = (MagickRealType) p[r].blue;
      result.opacity = QuantumRange - (MagickRealType) p[r].opacity;
      result.index   = 0.0;
      if ( image->colorspace == CMYKColorspace)
         result.index   = (MagickRealType) p_indexes[r];

      switch (method) {
        case ConvolveMorphology:
          /* Set the bias of the weighted average output */
          result.red     =
          result.green   =
          result.blue    =
          result.opacity =
          result.index   = bias;
          break;
        case DilateIntensityMorphology:
        case ErodeIntensityMorphology:
          /* use a boolean flag indicating when first match found */
          result.red = 0.0;  /* result is not used otherwise */
          break;
        default:
          break;
      }

      switch ( method ) {
        case ConvolveMorphology:
            /* Weighted Average of pixels using reflected kernel
            **
            ** NOTE for correct working of this operation for asymetrical
            ** kernels, the kernel needs to be applied in its reflected form.
            ** That is its values needs to be reversed.
            **
            ** Correlation is actually the same as this but without reflecting
            ** the kernel, and thus 'lower-level' that Convolution.  However
            ** as Convolution is the more common method used, and it does not
            ** really cost us much in terms of processing to use a reflected
            ** kernel, so it is Convolution that is implemented.
            **
            ** Correlation will have its kernel reflected before calling
            ** this function to do a Convolve.
            **
            ** For more details of Correlation vs Convolution see
            **   http://www.cs.umd.edu/~djacobs/CMSC426/Convolution.pdf
            */
            k = &kernel->values[ kernel->width*kernel->height-1 ];
            k_pixels = p;
            k_indexes = p_indexes;
            if ( ((channel & SyncChannels) == 0 ) ||
                                 (image->matte == MagickFalse) )
              { /* No 'Sync' involved.
                ** Convolution is simple greyscale channel operation
                */
                for (v=0; v < (ssize_t) kernel->height; v++) {
                  for (u=0; u < (ssize_t) kernel->width; u++, k--) {
                    if ( IsNan(*k) ) continue;
                    result.red     += (*k)*k_pixels[u].red;
                    result.green   += (*k)*k_pixels[u].green;
                    result.blue    += (*k)*k_pixels[u].blue;
                    result.opacity += (*k)*k_pixels[u].opacity;
                    if ( image->colorspace == CMYKColorspace)
                      result.index   += (*k)*k_indexes[u];
                  }
                  k_pixels += image->columns+kernel->width;
                  k_indexes += image->columns+kernel->width;
                }
                if ((channel & RedChannel) != 0)
                  q->red = ClampToQuantum(result.red);
                if ((channel & GreenChannel) != 0)
                  q->green = ClampToQuantum(result.green);
                if ((channel & BlueChannel) != 0)
                  q->blue = ClampToQuantum(result.blue);
                if ((channel & OpacityChannel) != 0
                    && image->matte == MagickTrue )
                  q->opacity = ClampToQuantum(result.opacity);
                if ((channel & IndexChannel) != 0
                    && image->colorspace == CMYKColorspace)
                  q_indexes[x] = ClampToQuantum(result.index);
              }
            else
              { /* Channel 'Sync' Flag, and Alpha Channel enabled.
                ** Weight the color channels with Alpha Channel so that
                ** transparent pixels are not part of the results.
                */
                MagickRealType
                  alpha,  /* alpha weighting of colors : kernel*alpha  */
                  gamma;  /* divisor, sum of color weighting values */

                gamma=0.0;
                for (v=0; v < (ssize_t) kernel->height; v++) {
                  for (u=0; u < (ssize_t) kernel->width; u++, k--) {
                    if ( IsNan(*k) ) continue;
                    alpha=(*k)*(QuantumScale*(QuantumRange-
                                          k_pixels[u].opacity));
                    gamma += alpha;
                    result.red     += alpha*k_pixels[u].red;
                    result.green   += alpha*k_pixels[u].green;
                    result.blue    += alpha*k_pixels[u].blue;
                    result.opacity += (*k)*k_pixels[u].opacity;
                    if ( image->colorspace == CMYKColorspace)
                      result.index   += alpha*k_indexes[u];
                  }
                  k_pixels += image->columns+kernel->width;
                  k_indexes += image->columns+kernel->width;
                }
                /* Sync'ed channels, all channels are modified */
                gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
                q->red = ClampToQuantum(gamma*result.red);
                q->green = ClampToQuantum(gamma*result.green);
                q->blue = ClampToQuantum(gamma*result.blue);
                q->opacity = ClampToQuantum(result.opacity);
                if (image->colorspace == CMYKColorspace)
                  q_indexes[x] = ClampToQuantum(gamma*result.index);
              }
            break;

        case ErodeMorphology:
            /* Minimum Value within kernel neighbourhood
            **
            ** NOTE that the kernel is not reflected for this operation!
            **
            ** NOTE: in normal Greyscale Morphology, the kernel value should
            ** be added to the real value, this is currently not done, due to
            ** the nature of the boolean kernels being used.
            */
            k = kernel->values;
            k_pixels = p;
            k_indexes = p_indexes;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              for (u=0; u < (ssize_t) kernel->width; u++, k++) {
                if ( IsNan(*k) || (*k) < 0.5 ) continue;
                Minimize(min.red,     (double) k_pixels[u].red);
                Minimize(min.green,   (double) k_pixels[u].green);
                Minimize(min.blue,    (double) k_pixels[u].blue);
                Minimize(min.opacity,
                            QuantumRange-(double) k_pixels[u].opacity);
                if ( image->colorspace == CMYKColorspace)
                  Minimize(min.index,   (double) k_indexes[u]);
              }
              k_pixels += image->columns+kernel->width;
              k_indexes += image->columns+kernel->width;
            }
            break;

        case DilateMorphology:
            /* Maximum Value within kernel neighbourhood
            **
            ** NOTE for correct working of this operation for asymetrical
            ** kernels, the kernel needs to be applied in its reflected form.
            ** That is its values needs to be reversed.
            **
            ** NOTE: in normal Greyscale Morphology, the kernel value should
            ** be added to the real value, this is currently not done, due to
            ** the nature of the boolean kernels being used.
            **
            */
            k = &kernel->values[ kernel->width*kernel->height-1 ];
            k_pixels = p;
            k_indexes = p_indexes;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              for (u=0; u < (ssize_t) kernel->width; u++, k--) {
                if ( IsNan(*k) || (*k) < 0.5 ) continue;
                Maximize(max.red,     (double) k_pixels[u].red);
                Maximize(max.green,   (double) k_pixels[u].green);
                Maximize(max.blue,    (double) k_pixels[u].blue);
                Maximize(max.opacity,
                            QuantumRange-(double) k_pixels[u].opacity);
                if ( image->colorspace == CMYKColorspace)
                  Maximize(max.index,   (double) k_indexes[u]);
              }
              k_pixels += image->columns+kernel->width;
              k_indexes += image->columns+kernel->width;
            }
            break;

        case HitAndMissMorphology:
        case ThinningMorphology:
        case ThickenMorphology:
            /* Minimum of Foreground Pixel minus Maxumum of Background Pixels
            **
            ** NOTE that the kernel is not reflected for this operation,
            ** and consists of both foreground and background pixel
            ** neighbourhoods, 0.0 for background, and 1.0 for foreground
            ** with either Nan or 0.5 values for don't care.
            **
            ** Note that this will never produce a meaningless negative
            ** result.  Such results can cause Thinning/Thicken to not work
            ** correctly when used against a greyscale image.
            */
            k = kernel->values;
            k_pixels = p;
            k_indexes = p_indexes;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              for (u=0; u < (ssize_t) kernel->width; u++, k++) {
                if ( IsNan(*k) ) continue;
                if ( (*k) > 0.7 )
                { /* minimim of foreground pixels */
                  Minimize(min.red,     (double) k_pixels[u].red);
                  Minimize(min.green,   (double) k_pixels[u].green);
                  Minimize(min.blue,    (double) k_pixels[u].blue);
                  Minimize(min.opacity,
                              QuantumRange-(double) k_pixels[u].opacity);
                  if ( image->colorspace == CMYKColorspace)
                    Minimize(min.index,   (double) k_indexes[u]);
                }
                else if ( (*k) < 0.3 )
                { /* maximum of background pixels */
                  Maximize(max.red,     (double) k_pixels[u].red);
                  Maximize(max.green,   (double) k_pixels[u].green);
                  Maximize(max.blue,    (double) k_pixels[u].blue);
                  Maximize(max.opacity,
                              QuantumRange-(double) k_pixels[u].opacity);
                  if ( image->colorspace == CMYKColorspace)
                    Maximize(max.index,   (double) k_indexes[u]);
                }
              }
              k_pixels += image->columns+kernel->width;
              k_indexes += image->columns+kernel->width;
            }
            /* Pattern Match if difference is positive */
            min.red     -= max.red;     Maximize( min.red,     0.0 );
            min.green   -= max.green;   Maximize( min.green,   0.0 );
            min.blue    -= max.blue;    Maximize( min.blue,    0.0 );
            min.opacity -= max.opacity; Maximize( min.opacity, 0.0 );
            min.index   -= max.index;   Maximize( min.index,   0.0 );
            break;

        case ErodeIntensityMorphology:
            /* Select Pixel with Minimum Intensity within kernel neighbourhood
            **
            ** WARNING: the intensity test fails for CMYK and does not
            ** take into account the moderating effect of the alpha channel
            ** on the intensity.
            **
            ** NOTE that the kernel is not reflected for this operation!
            */
            k = kernel->values;
            k_pixels = p;
            k_indexes = p_indexes;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              for (u=0; u < (ssize_t) kernel->width; u++, k++) {
                if ( IsNan(*k) || (*k) < 0.5 ) continue;
                if ( result.red == 0.0 ||
                     PixelIntensity(&(k_pixels[u])) < PixelIntensity(q) ) {
                  /* copy the whole pixel - no channel selection */
                  *q = k_pixels[u];
                  if ( result.red > 0.0 ) changed++;
                  result.red = 1.0;
                }
              }
              k_pixels += image->columns+kernel->width;
              k_indexes += image->columns+kernel->width;
            }
            break;

        case DilateIntensityMorphology:
            /* Select Pixel with Maximum Intensity within kernel neighbourhood
            **
            ** WARNING: the intensity test fails for CMYK and does not
            ** take into account the moderating effect of the alpha channel
            ** on the intensity (yet).
            **
            ** NOTE for correct working of this operation for asymetrical
            ** kernels, the kernel needs to be applied in its reflected form.
            ** That is its values needs to be reversed.
            */
            k = &kernel->values[ kernel->width*kernel->height-1 ];
            k_pixels = p;
            k_indexes = p_indexes;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              for (u=0; u < (ssize_t) kernel->width; u++, k--) {
                if ( IsNan(*k) || (*k) < 0.5 ) continue; /* boolean kernel */
                if ( result.red == 0.0 ||
                     PixelIntensity(&(k_pixels[u])) > PixelIntensity(q) ) {
                  /* copy the whole pixel - no channel selection */
                  *q = k_pixels[u];
                  if ( result.red > 0.0 ) changed++;
                  result.red = 1.0;
                }
              }
              k_pixels += image->columns+kernel->width;
              k_indexes += image->columns+kernel->width;
            }
            break;


        case DistanceMorphology:
            /* Add kernel Value and select the minimum value found.
            ** The result is a iterative distance from edge of image shape.
            **
            ** All Distance Kernels are symetrical, but that may not always
            ** be the case. For example how about a distance from left edges?
            ** To work correctly with asymetrical kernels the reflected kernel
            ** needs to be applied.
            **
            ** Actually this is really a GreyErode with a negative kernel!
            **
            */
            k = &kernel->values[ kernel->width*kernel->height-1 ];
            k_pixels = p;
            k_indexes = p_indexes;
            for (v=0; v < (ssize_t) kernel->height; v++) {
              for (u=0; u < (ssize_t) kernel->width; u++, k--) {
                if ( IsNan(*k) ) continue;
                Minimize(result.red,     (*k)+k_pixels[u].red);
                Minimize(result.green,   (*k)+k_pixels[u].green);
                Minimize(result.blue,    (*k)+k_pixels[u].blue);
                Minimize(result.opacity, (*k)+QuantumRange-k_pixels[u].opacity);
                if ( image->colorspace == CMYKColorspace)
                  Minimize(result.index,   (*k)+k_indexes[u]);
              }
              k_pixels += image->columns+kernel->width;
              k_indexes += image->columns+kernel->width;
            }
            break;

        case UndefinedMorphology:
        default:
            break; /* Do nothing */
      }
      /* Final mathematics of results (combine with original image?)
      **
      ** NOTE: Difference Morphology operators Edge* and *Hat could also
      ** be done here but works better with iteration as a image difference
      ** in the controling function (below).  Thicken and Thinning however
      ** should be done here so thay can be iterated correctly.
      */
      switch ( method ) {
        case HitAndMissMorphology:
        case ErodeMorphology:
          result = min;    /* minimum of neighbourhood */
          break;
        case DilateMorphology:
          result = max;    /* maximum of neighbourhood */
          break;
        case ThinningMorphology:
          /* subtract pattern match from original */
          result.red     -= min.red;
          result.green   -= min.green;
          result.blue    -= min.blue;
          result.opacity -= min.opacity;
          result.index   -= min.index;
          break;
        case ThickenMorphology:
          /* Add the pattern matchs to the original */
          result.red     += min.red;
          result.green   += min.green;
          result.blue    += min.blue;
          result.opacity += min.opacity;
          result.index   += min.index;
          break;
        default:
          /* result directly calculated or assigned */
          break;
      }
      /* Assign the resulting pixel values - Clamping Result */
      switch ( method ) {
        case UndefinedMorphology:
        case ConvolveMorphology:
        case DilateIntensityMorphology:
        case ErodeIntensityMorphology:
          break;  /* full pixel was directly assigned - not a channel method */
        default:
          if ((channel & RedChannel) != 0)
            q->red = ClampToQuantum(result.red);
          if ((channel & GreenChannel) != 0)
            q->green = ClampToQuantum(result.green);
          if ((channel & BlueChannel) != 0)
            q->blue = ClampToQuantum(result.blue);
          if ((channel & OpacityChannel) != 0
              && image->matte == MagickTrue )
            q->opacity = ClampToQuantum(QuantumRange-result.opacity);
          if ((channel & IndexChannel) != 0
              && image->colorspace == CMYKColorspace)
            q_indexes[x] = ClampToQuantum(result.index);
          break;
      }
      /* Count up changed pixels */
      if (   ( p[r].red != q->red )
          || ( p[r].green != q->green )
          || ( p[r].blue != q->blue )
          || ( p[r].opacity != q->opacity )
          || ( image->colorspace == CMYKColorspace &&
                  p_indexes[r] != q_indexes[x] ) )
        changed++;  /* The pixel was changed in some way! */
      p++;
      q++;
    } /* x */
    if ( SyncCacheViewAuthenticPixels(q_view,exception) == MagickFalse)
      status=MagickFalse;
    if (image->progress_monitor != (MagickProgressMonitor) NULL)
      {
        MagickBooleanType
          proceed;

#if defined(MAGICKCORE_OPENMP_SUPPORT)
  #pragma omp critical (MagickCore_MorphologyImage)
#endif
        proceed=SetImageProgress(image,MorphologyTag,progress++,image->rows);
        if (proceed == MagickFalse)
          status=MagickFalse;
      }
  } /* y */
  result_image->type=image->type;
  q_view=DestroyCacheView(q_view);
  p_view=DestroyCacheView(p_view);
  return(status ? (size_t) changed : 0);
}


MagickExport Image *MorphologyApply(const Image *image, const ChannelType
     channel,const MorphologyMethod method, const ssize_t iterations,
     const KernelInfo *kernel, const CompositeOperator compose,
     const double bias, ExceptionInfo *exception)
{
  Image
    *curr_image,    /* Image we are working with or iterating */
    *work_image,    /* secondary image for primative iteration */
    *save_image,    /* saved image - for 'edge' method only */
    *rslt_image;    /* resultant image - after multi-kernel handling */

  KernelInfo
    *reflected_kernel, /* A reflected copy of the kernel (if needed) */
    *norm_kernel,      /* the current normal un-reflected kernel */
    *rflt_kernel,      /* the current reflected kernel (if needed) */
    *this_kernel;      /* the kernel being applied */

  MorphologyMethod
    primative;      /* the current morphology primative being applied */

  CompositeOperator
    rslt_compose;   /* multi-kernel compose method for results to use */

  MagickBooleanType
    verbose;        /* verbose output of results */

  size_t
    method_loop,    /* Loop 1: number of compound method iterations */
    method_limit,   /*         maximum number of compound method iterations */
    kernel_number,  /* Loop 2: the kernel number being applied */
    stage_loop,     /* Loop 3: primative loop for compound morphology */
    stage_limit,    /*         how many primatives in this compound */
    kernel_loop,    /* Loop 4: iterate the kernel (basic morphology) */
    kernel_limit,   /*         number of times to iterate kernel */
    count,          /* total count of primative steps applied */
    changed,        /* number pixels changed by last primative operation */
    kernel_changed, /* total count of changed using iterated kernel */
    method_changed; /* total count of changed over method iteration */

  char
    v_info[80];

  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  assert(kernel != (KernelInfo *) NULL);
  assert(kernel->signature == MagickSignature);
  assert(exception != (ExceptionInfo *) NULL);
  assert(exception->signature == MagickSignature);

  count = 0;      /* number of low-level morphology primatives performed */
  if ( iterations == 0 )
    return((Image *)NULL);   /* null operation - nothing to do! */

  kernel_limit = (size_t) iterations;
  if ( iterations < 0 )  /* negative interations = infinite (well alomst) */
     kernel_limit = image->columns > image->rows ? image->columns : image->rows;


  verbose = IsMagickTrue(GetImageArtifact(image,"verbose"));

  /* initialise for cleanup */
  curr_image = (Image *) image;
  work_image = save_image = rslt_image = (Image *) NULL;
  reflected_kernel = (KernelInfo *) NULL;

  /* Initialize specific methods
   * + which loop should use the given iteratations
   * + how many primatives make up the compound morphology
   * + multi-kernel compose method to use (by default)
   */
  method_limit = 1;       /* just do method once, unless otherwise set */
  stage_limit = 1;        /* assume method is not a compount */
  rslt_compose = compose; /* and we are composing multi-kernels as given */
  switch( method ) {
    case SmoothMorphology:  /* 4 primative compound morphology */
      stage_limit = 4;
      break;
    case OpenMorphology:    /* 2 primative compound morphology */
    case OpenIntensityMorphology:
    case TopHatMorphology:
    case CloseMorphology:
    case CloseIntensityMorphology:
    case BottomHatMorphology:
    case EdgeMorphology:
      stage_limit = 2;
      break;
    case HitAndMissMorphology:
      rslt_compose = LightenCompositeOp;  /* Union of multi-kernel results */
      /* FALL THUR */
    case ThinningMorphology:
    case ThickenMorphology:
      method_limit = kernel_limit;  /* iterate the whole method */
      kernel_limit = 1;             /* do not do kernel iteration  */
      break;
    default:
      break;
  }

  /* Handle user (caller) specified multi-kernel composition method */
  if ( compose != UndefinedCompositeOp )
    rslt_compose = compose;  /* override default composition for method */
  if ( rslt_compose == UndefinedCompositeOp )
    rslt_compose = NoCompositeOp; /* still not defined! Then re-iterate */

  /* Some methods require a reflected kernel to use with primatives.
   * Create the reflected kernel for those methods. */
  switch ( method ) {
    case CorrelateMorphology:
    case CloseMorphology:
    case CloseIntensityMorphology:
    case BottomHatMorphology:
    case SmoothMorphology:
      reflected_kernel = CloneKernelInfo(kernel);
      if (reflected_kernel == (KernelInfo *) NULL)
        goto error_cleanup;
      RotateKernelInfo(reflected_kernel,180);
      break;
    default:
      break;
  }

  /* Loop 1:  iterate the compound method */
  method_loop = 0;
  method_changed = 1;
  while ( method_loop < method_limit && method_changed > 0 ) {
    method_loop++;
    method_changed = 0;

    /* Loop 2:  iterate over each kernel in a multi-kernel list */
    norm_kernel = (KernelInfo *) kernel;
    this_kernel = (KernelInfo *) kernel;
    rflt_kernel = reflected_kernel;

    kernel_number = 0;
    while ( norm_kernel != NULL ) {

      /* Loop 3: Compound Morphology Staging - Select Primative to apply */
      stage_loop = 0;          /* the compound morphology stage number */
      while ( stage_loop < stage_limit ) {
        stage_loop++;   /* The stage of the compound morphology */

        /* Select primative morphology for this stage of compound method */
        this_kernel = norm_kernel; /* default use unreflected kernel */
        primative = method;        /* Assume method is a primative */
        switch( method ) {
          case ErodeMorphology:      /* just erode */
          case EdgeInMorphology:     /* erode and image difference */
            primative = ErodeMorphology;
            break;
          case DilateMorphology:     /* just dilate */
          case EdgeOutMorphology:    /* dilate and image difference */
            primative = DilateMorphology;
            break;
          case OpenMorphology:       /* erode then dialate */
          case TopHatMorphology:     /* open and image difference */
            primative = ErodeMorphology;
            if ( stage_loop == 2 )
              primative = DilateMorphology;
            break;
          case OpenIntensityMorphology:
            primative = ErodeIntensityMorphology;
            if ( stage_loop == 2 )
              primative = DilateIntensityMorphology;
            break;
          case CloseMorphology:      /* dilate, then erode */
          case BottomHatMorphology:  /* close and image difference */
            this_kernel = rflt_kernel; /* use the reflected kernel */
            primative = DilateMorphology;
            if ( stage_loop == 2 )
              primative = ErodeMorphology;
            break;
          case CloseIntensityMorphology:
            this_kernel = rflt_kernel; /* use the reflected kernel */
            primative = DilateIntensityMorphology;
            if ( stage_loop == 2 )
              primative = ErodeIntensityMorphology;
            break;
          case SmoothMorphology:         /* open, close */
            switch ( stage_loop ) {
              case 1: /* start an open method, which starts with Erode */
                primative = ErodeMorphology;
                break;
              case 2:  /* now Dilate the Erode */
                primative = DilateMorphology;
                break;
              case 3:  /* Reflect kernel a close */
                this_kernel = rflt_kernel; /* use the reflected kernel */
                primative = DilateMorphology;
                break;
              case 4:  /* Finish the Close */
                this_kernel = rflt_kernel; /* use the reflected kernel */
                primative = ErodeMorphology;
                break;
            }
            break;
          case EdgeMorphology:        /* dilate and erode difference */
            primative = DilateMorphology;
            if ( stage_loop == 2 ) {
              save_image = curr_image;      /* save the image difference */
              curr_image = (Image *) image;
              primative = ErodeMorphology;
            }
            break;
          case CorrelateMorphology:
            /* A Correlation is a Convolution with a reflected kernel.
            ** However a Convolution is a weighted sum using a reflected
            ** kernel.  It may seem stange to convert a Correlation into a
            ** Convolution as the Correlation is the simplier method, but
            ** Convolution is much more commonly used, and it makes sense to
            ** implement it directly so as to avoid the need to duplicate the
            ** kernel when it is not required (which is typically the
            ** default).
            */
            this_kernel = rflt_kernel; /* use the reflected kernel */
            primative = ConvolveMorphology;
            break;
          default:
            break;
        }
        assert( this_kernel != (KernelInfo *) NULL );

        /* Extra information for debugging compound operations */
        if ( verbose == MagickTrue ) {
          if ( stage_limit > 1 )
            (void) FormatMagickString(v_info,MaxTextExtent,"%s:%.20g.%.20g -> ",
             MagickOptionToMnemonic(MagickMorphologyOptions,method),(double)
             method_loop,(double) stage_loop);
          else if ( primative != method )
            (void) FormatMagickString(v_info, MaxTextExtent, "%s:%.20g -> ",
              MagickOptionToMnemonic(MagickMorphologyOptions, method),(double)
              method_loop);
          else
            v_info[0] = '\0';
        }

        /* Loop 4: Iterate the kernel with primative */
        kernel_loop = 0;
        kernel_changed = 0;
        changed = 1;
        while ( kernel_loop < kernel_limit && changed > 0 ) {
          kernel_loop++;     /* the iteration of this kernel */

          /* Create a destination image, if not yet defined */
          if ( work_image == (Image *) NULL )
            {
              work_image=CloneImage(image,0,0,MagickTrue,exception);
              if (work_image == (Image *) NULL)
                goto error_cleanup;
              if (SetImageStorageClass(work_image,DirectClass) == MagickFalse)
                {
                  InheritException(exception,&work_image->exception);
                  goto error_cleanup;
                }
            }

          /* APPLY THE MORPHOLOGICAL PRIMITIVE (curr -> work) */
          count++;
          changed = MorphologyPrimitive(curr_image, work_image, primative,
                        channel, this_kernel, bias, exception);
          kernel_changed += changed;
          method_changed += changed;

          if ( verbose == MagickTrue ) {
            if ( kernel_loop > 1 )
              fprintf(stderr, "\n"); /* add end-of-line from previous */
            (void) fprintf(stderr, "%s%s%s:%.20g.%.20g #%.20g => Changed %.20g",
              v_info,MagickOptionToMnemonic(MagickMorphologyOptions,
              primative),(this_kernel == rflt_kernel ) ? "*" : "",
              (double) (method_loop+kernel_loop-1),(double) kernel_number,
              (double) count,(double) changed);
          }
          /* prepare next loop */
          { Image *tmp = work_image;   /* swap images for iteration */
            work_image = curr_image;
            curr_image = tmp;
          }
          if ( work_image == image )
            work_image = (Image *) NULL; /* replace input 'image' */

        } /* End Loop 4: Iterate the kernel with primative */

        if ( verbose == MagickTrue && kernel_changed != changed )
          fprintf(stderr, "   Total %.20g",(double) kernel_changed);
        if ( verbose == MagickTrue && stage_loop < stage_limit )
          fprintf(stderr, "\n"); /* add end-of-line before looping */

#if 0
    fprintf(stderr, "--E-- image=0x%lx\n", (unsigned long)image);
    fprintf(stderr, "      curr =0x%lx\n", (unsigned long)curr_image);
    fprintf(stderr, "      work =0x%lx\n", (unsigned long)work_image);
    fprintf(stderr, "      save =0x%lx\n", (unsigned long)save_image);
    fprintf(stderr, "      union=0x%lx\n", (unsigned long)rslt_image);
#endif

      } /* End Loop 3: Primative (staging) Loop for Coumpound Methods */

      /*  Final Post-processing for some Compound Methods
      **
      ** The removal of any 'Sync' channel flag in the Image Compositon
      ** below ensures the methematical compose method is applied in a
      ** purely mathematical way, and only to the selected channels.
      ** Turn off SVG composition 'alpha blending'.
      */
      switch( method ) {
        case EdgeOutMorphology:
        case EdgeInMorphology:
        case TopHatMorphology:
        case BottomHatMorphology:
          if ( verbose == MagickTrue )
            fprintf(stderr, "\n%s: Difference with original image",
                 MagickOptionToMnemonic(MagickMorphologyOptions, method) );
          (void) CompositeImageChannel(curr_image,
                  (ChannelType) (channel & ~SyncChannels),
                  DifferenceCompositeOp, image, 0, 0);
          break;
        case EdgeMorphology:
          if ( verbose == MagickTrue )
            fprintf(stderr, "\n%s: Difference of Dilate and Erode",
                 MagickOptionToMnemonic(MagickMorphologyOptions, method) );
          (void) CompositeImageChannel(curr_image,
                  (ChannelType) (channel & ~SyncChannels),
                  DifferenceCompositeOp, save_image, 0, 0);
          save_image = DestroyImage(save_image); /* finished with save image */
          break;
        default:
          break;
      }

      /* multi-kernel handling:  re-iterate, or compose results */
      if ( kernel->next == (KernelInfo *) NULL )
        rslt_image = curr_image;   /* just return the resulting image */
      else if ( rslt_compose == NoCompositeOp )
        { if ( verbose == MagickTrue ) {
            if ( this_kernel->next != (KernelInfo *) NULL )
              fprintf(stderr, " (re-iterate)");
            else
              fprintf(stderr, " (done)");
          }
          rslt_image = curr_image; /* return result, and re-iterate */
        }
      else if ( rslt_image == (Image *) NULL)
        { if ( verbose == MagickTrue )
            fprintf(stderr, " (save for compose)");
          rslt_image = curr_image;
          curr_image = (Image *) image;  /* continue with original image */
        }
      else
        { /* add the new 'current' result to the composition
          **
          ** The removal of any 'Sync' channel flag in the Image Compositon
          ** below ensures the methematical compose method is applied in a
          ** purely mathematical way, and only to the selected channels.
          ** Turn off SVG composition 'alpha blending'.
          **
          ** The compose image order is specifically so that the new image can
          ** be subtarcted 'Minus' from the collected result, to allow you to
          ** convert a HitAndMiss methd into a Thinning method.
          */
          if ( verbose == MagickTrue )
            fprintf(stderr, " (compose \"%s\")",
                 MagickOptionToMnemonic(MagickComposeOptions, rslt_compose) );
          (void) CompositeImageChannel(curr_image,
               (ChannelType) (channel & ~SyncChannels), rslt_compose,
               rslt_image, 0, 0);
          rslt_image = DestroyImage(rslt_image);
          rslt_image = curr_image;
          curr_image = (Image *) image;  /* continue with original image */
        }
      if ( verbose == MagickTrue )
        fprintf(stderr, "\n");

      /* loop to the next kernel in a multi-kernel list */
      norm_kernel = norm_kernel->next;
      if ( rflt_kernel != (KernelInfo *) NULL )
        rflt_kernel = rflt_kernel->next;
      kernel_number++;
    } /* End Loop 2: Loop over each kernel */

  } /* End Loop 1: compound method interation */

  goto exit_cleanup;

  /* Yes goto's are bad, but it makes cleanup lot more efficient */
error_cleanup:
  if ( curr_image != (Image *) NULL &&
       curr_image != rslt_image &&
       curr_image != image )
    curr_image = DestroyImage(curr_image);
  if ( rslt_image != (Image *) NULL )
    rslt_image = DestroyImage(rslt_image);
exit_cleanup:
  if ( curr_image != (Image *) NULL &&
       curr_image != rslt_image &&
       curr_image != image )
    curr_image = DestroyImage(curr_image);
  if ( work_image != (Image *) NULL )
    work_image = DestroyImage(work_image);
  if ( save_image != (Image *) NULL )
    save_image = DestroyImage(save_image);
  if ( reflected_kernel != (KernelInfo *) NULL )
    reflected_kernel = DestroyKernelInfo(reflected_kernel);
  return(rslt_image);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     M o r p h o l o g y I m a g e C h a n n e l                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MorphologyImageChannel() applies a user supplied kernel to the image
%  according to the given mophology method.
%
%  This function applies any and all user defined settings before calling
%  the above internal function MorphologyApply().
%
%  User defined settings include...
%    * Output Bias for Convolution and correlation   ("-bias")
%    * Kernel Scale/normalize settings     ("-set 'option:convolve:scale'")
%      This can also includes the addition of a scaled unity kernel.
%    * Show Kernel being applied           ("-set option:showkernel 1")
%
%  The format of the MorphologyImage method is:
%
%      Image *MorphologyImage(const Image *image,MorphologyMethod method,
%        const ssize_t iterations,KernelInfo *kernel,ExceptionInfo *exception)
%
%      Image *MorphologyImageChannel(const Image *image, const ChannelType
%        channel,MorphologyMethod method,const ssize_t iterations,
%        KernelInfo *kernel,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o image: the image.
%
%    o method: the morphology method to be applied.
%
%    o iterations: apply the operation this many times (or no change).
%                  A value of -1 means loop until no change found.
%                  How this is applied may depend on the morphology method.
%                  Typically this is a value of 1.
%
%    o channel: the channel type.
%
%    o kernel: An array of double representing the morphology kernel.
%              Warning: kernel may be normalized for the Convolve method.
%
%    o exception: return any errors or warnings in this structure.
%
*/

MagickExport Image *MorphologyImageChannel(const Image *image,
  const ChannelType channel,const MorphologyMethod method,
  const ssize_t iterations,const KernelInfo *kernel,ExceptionInfo *exception)
{
  KernelInfo
    *curr_kernel;

  CompositeOperator
    compose;

  Image
    *morphology_image;


  /* Apply Convolve/Correlate Normalization and Scaling Factors.
   * This is done BEFORE the ShowKernelInfo() function is called so that
   * users can see the results of the 'option:convolve:scale' option.
   */
  curr_kernel = (KernelInfo *) kernel;
  if ( method == ConvolveMorphology ||  method == CorrelateMorphology )
    {
      const char
        *artifact;
      artifact = GetImageArtifact(image,"convolve:scale");
      if ( artifact != (const char *)NULL ) {
        if ( curr_kernel == kernel )
          curr_kernel = CloneKernelInfo(kernel);
        if (curr_kernel == (KernelInfo *) NULL) {
          curr_kernel=DestroyKernelInfo(curr_kernel);
          return((Image *) NULL);
        }
        ScaleGeometryKernelInfo(curr_kernel, artifact);
      }
    }

  /* display the (normalized) kernel via stderr */
  if ( IsMagickTrue(GetImageArtifact(image,"showkernel"))
    || IsMagickTrue(GetImageArtifact(image,"convolve:showkernel"))
    || IsMagickTrue(GetImageArtifact(image,"morphology:showkernel")) )
    ShowKernelInfo(curr_kernel);

  /* Override the default handling of multi-kernel morphology results
   * If 'Undefined' use the default method
   * If 'None' (default for 'Convolve') re-iterate previous result
   * Otherwise merge resulting images using compose method given.
   * Default for 'HitAndMiss' is 'Lighten'.
   */
  { const char
      *artifact;
    artifact = GetImageArtifact(image,"morphology:compose");
    compose = UndefinedCompositeOp;  /* use default for method */
    if ( artifact != (const char *) NULL)
      compose = (CompositeOperator) ParseMagickOption(
                             MagickComposeOptions,MagickFalse,artifact);
  }
  /* Apply the Morphology */
  morphology_image = MorphologyApply(image, channel, method, iterations,
                         curr_kernel, compose, image->bias, exception);

  /* Cleanup and Exit */
  if ( curr_kernel != kernel )
    curr_kernel=DestroyKernelInfo(curr_kernel);
  return(morphology_image);
}

MagickExport Image *MorphologyImage(const Image *image, const MorphologyMethod
  method, const ssize_t iterations,const KernelInfo *kernel, ExceptionInfo
  *exception)
{
  Image
    *morphology_image;

  morphology_image=MorphologyImageChannel(image,DefaultChannels,method,
    iterations,kernel,exception);
  return(morphology_image);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+     R o t a t e K e r n e l I n f o                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  RotateKernelInfo() rotates the kernel by the angle given.
%
%  Currently it is restricted to 90 degree angles, of either 1D kernels
%  or square kernels. And 'circular' rotations of 45 degrees for 3x3 kernels.
%  It will ignore usless rotations for specific 'named' built-in kernels.
%
%  The format of the RotateKernelInfo method is:
%
%      void RotateKernelInfo(KernelInfo *kernel, double angle)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
%    o angle: angle to rotate in degrees
%
% This function is currently internal to this module only, but can be exported
% to other modules if needed.
*/
static void RotateKernelInfo(KernelInfo *kernel, double angle)
{
  /* angle the lower kernels first */
  if ( kernel->next != (KernelInfo *) NULL)
    RotateKernelInfo(kernel->next, angle);

  /* WARNING: Currently assumes the kernel (rightly) is horizontally symetrical
  **
  ** TODO: expand beyond simple 90 degree rotates, flips and flops
  */

  /* Modulus the angle */
  angle = fmod(angle, 360.0);
  if ( angle < 0 )
    angle += 360.0;

  if ( 337.5 < angle || angle <= 22.5 )
    return;   /* Near zero angle - no change! - At least not at this time */

  /* Handle special cases */
  switch (kernel->type) {
    /* These built-in kernels are cylindrical kernels, rotating is useless */
    case GaussianKernel:
    case DoGKernel:
    case LoGKernel:
    case DiskKernel:
    case PeaksKernel:
    case LaplacianKernel:
    case ChebyshevKernel:
    case ManhattanKernel:
    case EuclideanKernel:
      return;

    /* These may be rotatable at non-90 angles in the future */
    /* but simply rotating them in multiples of 90 degrees is useless */
    case SquareKernel:
    case DiamondKernel:
    case PlusKernel:
    case CrossKernel:
      return;

    /* These only allows a +/-90 degree rotation (by transpose) */
    /* A 180 degree rotation is useless */
    case BlurKernel:
    case RectangleKernel:
      if ( 135.0 < angle && angle <= 225.0 )
        return;
      if ( 225.0 < angle && angle <= 315.0 )
        angle -= 180;
      break;

    default:
      break;
  }
  /* Attempt rotations by 45 degrees */
  if ( 22.5 < fmod(angle,90.0) && fmod(angle,90.0) <= 67.5 )
    {
      if ( kernel->width == 3 && kernel->height == 3 )
        { /* Rotate a 3x3 square by 45 degree angle */
          MagickRealType t  = kernel->values[0];
          kernel->values[0] = kernel->values[3];
          kernel->values[3] = kernel->values[6];
          kernel->values[6] = kernel->values[7];
          kernel->values[7] = kernel->values[8];
          kernel->values[8] = kernel->values[5];
          kernel->values[5] = kernel->values[2];
          kernel->values[2] = kernel->values[1];
          kernel->values[1] = t;
          /* rotate non-centered origin */
          if ( kernel->x != 1 || kernel->y != 1 ) {
            ssize_t x,y;
            x = (ssize_t) kernel->x-1;
            y = (ssize_t) kernel->y-1;
                 if ( x == y  ) x = 0;
            else if ( x == 0  ) x = -y;
            else if ( x == -y ) y = 0;
            else if ( y == 0  ) y = x;
            kernel->x = (ssize_t) x+1;
            kernel->y = (ssize_t) y+1;
          }
          angle = fmod(angle+315.0, 360.0);  /* angle reduced 45 degrees */
          kernel->angle = fmod(kernel->angle+45.0, 360.0);
        }
      else
        perror("Unable to rotate non-3x3 kernel by 45 degrees");
    }
  if ( 45.0 < fmod(angle, 180.0)  && fmod(angle,180.0) <= 135.0 )
    {
      if ( kernel->width == 1 || kernel->height == 1 )
        { /* Do a transpose of a 1 dimentional kernel,
          ** which results in a fast 90 degree rotation of some type.
          */
          ssize_t
            t;
          t = (ssize_t) kernel->width;
          kernel->width = kernel->height;
          kernel->height = (size_t) t;
          t = kernel->x;
          kernel->x = kernel->y;
          kernel->y = t;
          if ( kernel->width == 1 ) {
            angle = fmod(angle+270.0, 360.0);     /* angle reduced 90 degrees */
            kernel->angle = fmod(kernel->angle+90.0, 360.0);
          } else {
            angle = fmod(angle+90.0, 360.0);   /* angle increased 90 degrees */
            kernel->angle = fmod(kernel->angle+270.0, 360.0);
          }
        }
      else if ( kernel->width == kernel->height )
        { /* Rotate a square array of values by 90 degrees */
          { register size_t
              i,j,x,y;
            register MagickRealType
              *k,t;
            k=kernel->values;
            for( i=0, x=kernel->width-1;  i<=x;   i++, x--)
              for( j=0, y=kernel->height-1;  j<y;   j++, y--)
                { t                    = k[i+j*kernel->width];
                  k[i+j*kernel->width] = k[j+x*kernel->width];
                  k[j+x*kernel->width] = k[x+y*kernel->width];
                  k[x+y*kernel->width] = k[y+i*kernel->width];
                  k[y+i*kernel->width] = t;
                }
          }
          /* rotate the origin - relative to center of array */
          { register ssize_t x,y;
            x = (ssize_t) (kernel->x*2-kernel->width+1);
            y = (ssize_t) (kernel->y*2-kernel->height+1);
            kernel->x = (ssize_t) ( -y +(ssize_t) kernel->width-1)/2;
            kernel->y = (ssize_t) ( +x +(ssize_t) kernel->height-1)/2;
          }
          angle = fmod(angle+270.0, 360.0);     /* angle reduced 90 degrees */
          kernel->angle = fmod(kernel->angle+90.0, 360.0);
        }
      else
        perror("Unable to rotate a non-square, non-linear kernel 90 degrees");
    }
  if ( 135.0 < angle && angle <= 225.0 )
    {
      /* For a 180 degree rotation - also know as a reflection
       * This is actually a very very common operation!
       * Basically all that is needed is a reversal of the kernel data!
       * And a reflection of the origon
       */
      size_t
        i,j;
      register double
        *k,t;

      k=kernel->values;
      for ( i=0, j=kernel->width*kernel->height-1;  i<j;  i++, j--)
        t=k[i],  k[i]=k[j],  k[j]=t;

      kernel->x = (ssize_t) kernel->width  - kernel->x - 1;
      kernel->y = (ssize_t) kernel->height - kernel->y - 1;
      angle = fmod(angle-180.0, 360.0);   /* angle+180 degrees */
      kernel->angle = fmod(kernel->angle+180.0, 360.0);
    }
  /* At this point angle should at least between -45 (315) and +45 degrees
   * In the future some form of non-orthogonal angled rotates could be
   * performed here, posibily with a linear kernel restriction.
   */

  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     S c a l e G e o m e t r y K e r n e l I n f o                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ScaleGeometryKernelInfo() takes a geometry argument string, typically
%  provided as a  "-set option:convolve:scale {geometry}" user setting,
%  and modifies the kernel according to the parsed arguments of that setting.
%
%  The first argument (and any normalization flags) are passed to
%  ScaleKernelInfo() to scale/normalize the kernel.  The second argument
%  is then passed to UnityAddKernelInfo() to add a scled unity kernel
%  into the scaled/normalized kernel.
%
%  The format of the ScaleKernelInfo method is:
%
%      void ScaleKernelInfo(KernelInfo *kernel, const double scaling_factor,
%               const MagickStatusType normalize_flags )
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel to modify
%
%    o geometry:
%             The geometry string to parse, typically from the user provided
%             "-set option:convolve:scale {geometry}" setting.
%
*/
MagickExport void ScaleGeometryKernelInfo (KernelInfo *kernel,
     const char *geometry)
{
  GeometryFlags
    flags;
  GeometryInfo
    args;

  SetGeometryInfo(&args);
  flags = (GeometryFlags) ParseGeometry(geometry, &args);

#if 0
  /* For Debugging Geometry Input */
  fprintf(stderr, "Geometry = 0x%04X : %lg x %lg %+lg %+lg\n",
       flags, args.rho, args.sigma, args.xi, args.psi );
#endif

  if ( (flags & PercentValue) != 0 )      /* Handle Percentage flag*/
    args.rho *= 0.01,  args.sigma *= 0.01;

  if ( (flags & RhoValue) == 0 )          /* Set Defaults for missing args */
    args.rho = 1.0;
  if ( (flags & SigmaValue) == 0 )
    args.sigma = 0.0;

  /* Scale/Normalize the input kernel */
  ScaleKernelInfo(kernel, args.rho, flags);

  /* Add Unity Kernel, for blending with original */
  if ( (flags & SigmaValue) != 0 )
    UnityAddKernelInfo(kernel, args.sigma);

  return;
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     S c a l e K e r n e l I n f o                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ScaleKernelInfo() scales the given kernel list by the given amount, with or
%  without normalization of the sum of the kernel values (as per given flags).
%
%  By default (no flags given) the values within the kernel is scaled
%  directly using given scaling factor without change.
%
%  If either of the two 'normalize_flags' are given the kernel will first be
%  normalized and then further scaled by the scaling factor value given.
%
%  Kernel normalization ('normalize_flags' given) is designed to ensure that
%  any use of the kernel scaling factor with 'Convolve' or 'Correlate'
%  morphology methods will fall into -1.0 to +1.0 range.  Note that for
%  non-HDRI versions of IM this may cause images to have any negative results
%  clipped, unless some 'bias' is used.
%
%  More specifically.  Kernels which only contain positive values (such as a
%  'Gaussian' kernel) will be scaled so that those values sum to +1.0,
%  ensuring a 0.0 to +1.0 output range for non-HDRI images.
%
%  For Kernels that contain some negative values, (such as 'Sharpen' kernels)
%  the kernel will be scaled by the absolute of the sum of kernel values, so
%  that it will generally fall within the +/- 1.0 range.
%
%  For kernels whose values sum to zero, (such as 'Laplician' kernels) kernel
%  will be scaled by just the sum of the postive values, so that its output
%  range will again fall into the  +/- 1.0 range.
%
%  For special kernels designed for locating shapes using 'Correlate', (often
%  only containing +1 and -1 values, representing foreground/brackground
%  matching) a special normalization method is provided to scale the positive
%  values seperatally to those of the negative values, so the kernel will be
%  forced to become a zero-sum kernel better suited to such searches.
%
%  WARNING: Correct normalization of the kernel assumes that the '*_range'
%  attributes within the kernel structure have been correctly set during the
%  kernels creation.
%
%  NOTE: The values used for 'normalize_flags' have been selected specifically
%  to match the use of geometry options, so that '!' means NormalizeValue, '^'
%  means CorrelateNormalizeValue.  All other GeometryFlags values are ignored.
%
%  The format of the ScaleKernelInfo method is:
%
%      void ScaleKernelInfo(KernelInfo *kernel, const double scaling_factor,
%               const MagickStatusType normalize_flags )
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
%    o scaling_factor:
%             multiply all values (after normalization) by this factor if not
%             zero.  If the kernel is normalized regardless of any flags.
%
%    o normalize_flags:
%             GeometryFlags defining normalization method to use.
%             specifically: NormalizeValue, CorrelateNormalizeValue,
%                           and/or PercentValue
%
*/
MagickExport void ScaleKernelInfo(KernelInfo *kernel,
  const double scaling_factor,const GeometryFlags normalize_flags)
{
  register ssize_t
    i;

  register double
    pos_scale,
    neg_scale;

  /* do the other kernels in a multi-kernel list first */
  if ( kernel->next != (KernelInfo *) NULL)
    ScaleKernelInfo(kernel->next, scaling_factor, normalize_flags);

  /* Normalization of Kernel */
  pos_scale = 1.0;
  if ( (normalize_flags&NormalizeValue) != 0 ) {
    if ( fabs(kernel->positive_range + kernel->negative_range) > MagickEpsilon )
      /* non-zero-summing kernel (generally positive) */
      pos_scale = fabs(kernel->positive_range + kernel->negative_range);
    else
      /* zero-summing kernel */
      pos_scale = kernel->positive_range;
  }
  /* Force kernel into a normalized zero-summing kernel */
  if ( (normalize_flags&CorrelateNormalizeValue) != 0 ) {
    pos_scale = ( fabs(kernel->positive_range) > MagickEpsilon )
                 ? kernel->positive_range : 1.0;
    neg_scale = ( fabs(kernel->negative_range) > MagickEpsilon )
                 ? -kernel->negative_range : 1.0;
  }
  else
    neg_scale = pos_scale;

  /* finialize scaling_factor for positive and negative components */
  pos_scale = scaling_factor/pos_scale;
  neg_scale = scaling_factor/neg_scale;

  for (i=0; i < (ssize_t) (kernel->width*kernel->height); i++)
    if ( ! IsNan(kernel->values[i]) )
      kernel->values[i] *= (kernel->values[i] >= 0) ? pos_scale : neg_scale;

  /* convolution output range */
  kernel->positive_range *= pos_scale;
  kernel->negative_range *= neg_scale;
  /* maximum and minimum values in kernel */
  kernel->maximum *= (kernel->maximum >= 0.0) ? pos_scale : neg_scale;
  kernel->minimum *= (kernel->minimum >= 0.0) ? pos_scale : neg_scale;

  /* swap kernel settings if user's scaling factor is negative */
  if ( scaling_factor < MagickEpsilon ) {
    double t;
    t = kernel->positive_range;
    kernel->positive_range = kernel->negative_range;
    kernel->negative_range = t;
    t = kernel->maximum;
    kernel->maximum = kernel->minimum;
    kernel->minimum = 1;
  }

  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     S h o w K e r n e l I n f o                                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ShowKernelInfo() outputs the details of the given kernel defination to
%  standard error, generally due to a users 'showkernel' option request.
%
%  The format of the ShowKernel method is:
%
%      void ShowKernelInfo(KernelInfo *kernel)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
*/
MagickExport void ShowKernelInfo(KernelInfo *kernel)
{
  KernelInfo
    *k;

  size_t
    c, i, u, v;

  for (c=0, k=kernel;  k != (KernelInfo *) NULL;  c++, k=k->next ) {

    fprintf(stderr, "Kernel");
    if ( kernel->next != (KernelInfo *) NULL )
      fprintf(stderr, " #%lu", (unsigned long) c );
    fprintf(stderr, " \"%s",
          MagickOptionToMnemonic(MagickKernelOptions, k->type) );
    if ( fabs(k->angle) > MagickEpsilon )
      fprintf(stderr, "@%lg", k->angle);
    fprintf(stderr, "\" of size %lux%lu%+ld%+ld",(unsigned long) k->width,
      (unsigned long) k->height,(long) k->x,(long) k->y);
    fprintf(stderr,
          " with values from %.*lg to %.*lg\n",
          GetMagickPrecision(), k->minimum,
          GetMagickPrecision(), k->maximum);
    fprintf(stderr, "Forming a output range from %.*lg to %.*lg",
          GetMagickPrecision(), k->negative_range,
          GetMagickPrecision(), k->positive_range);
    if ( fabs(k->positive_range+k->negative_range) < MagickEpsilon )
      fprintf(stderr, " (Zero-Summing)\n");
    else if ( fabs(k->positive_range+k->negative_range-1.0) < MagickEpsilon )
      fprintf(stderr, " (Normalized)\n");
    else
      fprintf(stderr, " (Sum %.*lg)\n",
          GetMagickPrecision(), k->positive_range+k->negative_range);
    for (i=v=0; v < k->height; v++) {
      fprintf(stderr, "%2lu:", (unsigned long) v );
      for (u=0; u < k->width; u++, i++)
        if ( IsNan(k->values[i]) )
          fprintf(stderr," %*s", GetMagickPrecision()+3, "nan");
        else
          fprintf(stderr," %*.*lg", GetMagickPrecision()+3,
              GetMagickPrecision(), k->values[i]);
      fprintf(stderr,"\n");
    }
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     U n i t y A d d K e r n a l I n f o                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  UnityAddKernelInfo() Adds a given amount of the 'Unity' Convolution Kernel
%  to the given pre-scaled and normalized Kernel.  This in effect adds that
%  amount of the original image into the resulting convolution kernel.  This
%  value is usually provided by the user as a percentage value in the
%  'convolve:scale' setting.
%
%  The resulting effect is to convert the defined kernels into blended
%  soft-blurs, unsharp kernels or into sharpening kernels.
%
%  The format of the UnityAdditionKernelInfo method is:
%
%      void UnityAdditionKernelInfo(KernelInfo *kernel, const double scale )
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
%    o scale:
%             scaling factor for the unity kernel to be added to
%             the given kernel.
%
*/
MagickExport void UnityAddKernelInfo(KernelInfo *kernel,
  const double scale)
{
  /* do the other kernels in a multi-kernel list first */
  if ( kernel->next != (KernelInfo *) NULL)
    UnityAddKernelInfo(kernel->next, scale);

  /* Add the scaled unity kernel to the existing kernel */
  kernel->values[kernel->x+kernel->y*kernel->width] += scale;
  CalcKernelMetaData(kernel);  /* recalculate the meta-data */

  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%     Z e r o K e r n e l N a n s                                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ZeroKernelNans() replaces any special 'nan' value that may be present in
%  the kernel with a zero value.  This is typically done when the kernel will
%  be used in special hardware (GPU) convolution processors, to simply
%  matters.
%
%  The format of the ZeroKernelNans method is:
%
%      void ZeroKernelNans (KernelInfo *kernel)
%
%  A description of each parameter follows:
%
%    o kernel: the Morphology/Convolution kernel
%
*/
MagickExport void ZeroKernelNans(KernelInfo *kernel)
{
  register size_t
    i;

  /* do the other kernels in a multi-kernel list first */
  if ( kernel->next != (KernelInfo *) NULL)
    ZeroKernelNans(kernel->next);

  for (i=0; i < (kernel->width*kernel->height); i++)
    if ( IsNan(kernel->values[i]) )
      kernel->values[i] = 0.0;

  return;
}