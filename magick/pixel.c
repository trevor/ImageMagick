/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%                      PPPP   IIIII  X   X  EEEEE  L                          %
%                      P   P    I     X X   E      L                          %
%                      PPPP     I      X    EEE    L                          %
%                      P        I     X X   E      L                          %
%                      P      IIIII  X   X  EEEEE  LLLLL                      %
%                                                                             %
%                  MagickCore Methods to Import/Export Pixels                 %
%                                                                             %
%                             Software Design                                 %
%                               John Cristy                                   %
%                               October 1998                                  %
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
%
*/

/*
  Include declarations.
*/
#include "magick/studio.h"
#include "magick/property.h"
#include "magick/blob.h"
#include "magick/blob-private.h"
#include "magick/color-private.h"
#include "magick/draw.h"
#include "magick/exception.h"
#include "magick/exception-private.h"
#include "magick/cache.h"
#include "magick/constitute.h"
#include "magick/delegate.h"
#include "magick/geometry.h"
#include "magick/list.h"
#include "magick/magick.h"
#include "magick/memory_.h"
#include "magick/monitor.h"
#include "magick/option.h"
#include "magick/pixel.h"
#include "magick/pixel-private.h"
#include "magick/quantum.h"
#include "magick/resource_.h"
#include "magick/semaphore.h"
#include "magick/statistic.h"
#include "magick/stream.h"
#include "magick/string_.h"
#include "magick/transform.h"
#include "magick/utility.h"

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   E x p o r t I m a g e P i x e l s                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ExportImagePixels() extracts pixel data from an image and returns it to you.
%  The method returns MagickTrue on success otherwise MagickFalse if an error is
%  encountered.  The data is returned as char, short int, int, ssize_t, float,
%  or double in the order specified by map.
%
%  Suppose you want to extract the first scanline of a 640x480 image as
%  character data in red-green-blue order:
%
%      ExportImagePixels(image,0,0,640,1,"RGB",CharPixel,pixels,exception);
%
%  The format of the ExportImagePixels method is:
%
%      MagickBooleanType ExportImagePixels(const Image *image,
%        const ssize_t x_offset,const ssize_t y_offset,const size_t columns,
%        const size_t rows,const char *map,const StorageType type,
%        void *pixels,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o image: the image.
%
%    o x_offset,y_offset,columns,rows:  These values define the perimeter
%      of a region of pixels you want to extract.
%
%    o map:  This string reflects the expected ordering of the pixel array.
%      It can be any combination or order of R = red, G = green, B = blue,
%      A = alpha (0 is transparent), O = opacity (0 is opaque), C = cyan,
%      Y = yellow, M = magenta, K = black, I = intensity (for grayscale),
%      P = pad.
%
%    o type: Define the data type of the pixels.  Float and double types are
%      normalized to [0..1] otherwise [0..QuantumRange].  Choose from these
%      types: CharPixel, DoublePixel, FloatPixel, IntegerPixel, LongPixel,
%      QuantumPixel, or ShortPixel.
%
%    o pixels: This array of values contain the pixel components as defined by
%      map and type.  You must preallocate this array where the expected
%      length varies depending on the values of width, height, map, and type.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport MagickBooleanType ExportImagePixels(const Image *image,
  const ssize_t x_offset,const ssize_t y_offset,const size_t columns,
  const size_t rows,const char *map,const StorageType type,void *pixels,
  ExceptionInfo *exception)
{
  QuantumType
    *quantum_map;

  register ssize_t
    i,
    x;

  register const IndexPacket
    *indexes;

  register const PixelPacket
    *p;

  size_t
    length;

  ssize_t
    y;

  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  if (image->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",image->filename);
  length=strlen(map);
  quantum_map=(QuantumType *) AcquireQuantumMemory(length,sizeof(*quantum_map));
  if (quantum_map == (QuantumType *) NULL)
    {
      (void) ThrowMagickException(exception,GetMagickModule(),
        ResourceLimitError,"MemoryAllocationFailed","`%s'",image->filename);
      return(MagickFalse);
    }
  for (i=0; i < (ssize_t) length; i++)
  {
    switch (map[i])
    {
      case 'A':
      case 'a':
      {
        quantum_map[i]=AlphaQuantum;
        break;
      }
      case 'B':
      case 'b':
      {
        quantum_map[i]=BlueQuantum;
        break;
      }
      case 'C':
      case 'c':
      {
        quantum_map[i]=CyanQuantum;
        if (image->colorspace == CMYKColorspace)
          break;
        quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
        (void) ThrowMagickException(exception,GetMagickModule(),ImageError,
          "ColorSeparatedImageRequired","`%s'",map);
        return(MagickFalse);
      }
      case 'g':
      case 'G':
      {
        quantum_map[i]=GreenQuantum;
        break;
      }
      case 'I':
      case 'i':
      {
        quantum_map[i]=IndexQuantum;
        break;
      }
      case 'K':
      case 'k':
      {
        quantum_map[i]=BlackQuantum;
        if (image->colorspace == CMYKColorspace)
          break;
        quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
        (void) ThrowMagickException(exception,GetMagickModule(),ImageError,
          "ColorSeparatedImageRequired","`%s'",map);
        return(MagickFalse);
      }
      case 'M':
      case 'm':
      {
        quantum_map[i]=MagentaQuantum;
        if (image->colorspace == CMYKColorspace)
          break;
        quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
        (void) ThrowMagickException(exception,GetMagickModule(),ImageError,
          "ColorSeparatedImageRequired","`%s'",map);
        return(MagickFalse);
      }
      case 'o':
      case 'O':
      {
        quantum_map[i]=OpacityQuantum;
        break;
      }
      case 'P':
      case 'p':
      {
        quantum_map[i]=UndefinedQuantum;
        break;
      }
      case 'R':
      case 'r':
      {
        quantum_map[i]=RedQuantum;
        break;
      }
      case 'Y':
      case 'y':
      {
        quantum_map[i]=YellowQuantum;
        if (image->colorspace == CMYKColorspace)
          break;
        quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
        (void) ThrowMagickException(exception,GetMagickModule(),ImageError,
          "ColorSeparatedImageRequired","`%s'",map);
        return(MagickFalse);
      }
      default:
      {
        quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
        (void) ThrowMagickException(exception,GetMagickModule(),OptionError,
          "UnrecognizedPixelMap","`%s'",map);
        return(MagickFalse);
      }
    }
  }
  switch (type)
  {
    case CharPixel:
    {
      register unsigned char
        *q;

      q=(unsigned char *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(GetBluePixelComponent(p));
              *q++=ScaleQuantumToChar(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToChar(GetRedPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(GetBluePixelComponent(p));
              *q++=ScaleQuantumToChar(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToChar(GetRedPixelComponent(p));
              *q++=ScaleQuantumToChar((Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(GetBluePixelComponent(p));
              *q++=ScaleQuantumToChar(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToChar(GetRedPixelComponent(p));
              *q++=ScaleQuantumToChar((Quantum) 0);
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(PixelIntensityToQuantum(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(GetRedPixelComponent(p));
              *q++=ScaleQuantumToChar(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToChar(GetBluePixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(GetRedPixelComponent(p));
              *q++=ScaleQuantumToChar(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToChar(GetBluePixelComponent(p));
              *q++=ScaleQuantumToChar((Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToChar(GetRedPixelComponent(p));
              *q++=ScaleQuantumToChar(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToChar(GetBluePixelComponent(p));
              *q++=ScaleQuantumToChar((Quantum) 0);
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=ScaleQuantumToChar(GetRedPixelComponent(p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=ScaleQuantumToChar(GetGreenPixelComponent(p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=ScaleQuantumToChar(GetBluePixelComponent(p));
                break;
              }
              case AlphaQuantum:
              {
                *q=ScaleQuantumToChar((Quantum) (GetAlphaPixelComponent(p)));
                break;
              }
              case OpacityQuantum:
              {
                *q=ScaleQuantumToChar(GetOpacityPixelComponent(p));
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=ScaleQuantumToChar(indexes[x]);
                break;
              }
              case IndexQuantum:
              {
                *q=ScaleQuantumToChar(PixelIntensityToQuantum(p));
                break;
              }
              default:
                break;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    case DoublePixel:
    {
      register double
        *q;

      q=(double *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*GetBluePixelComponent(p));
              *q++=(double) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(double) (QuantumScale*GetRedPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*GetBluePixelComponent(p));
              *q++=(double) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(double) (QuantumScale*GetRedPixelComponent(p));
              *q++=(double) (QuantumScale*((Quantum) (QuantumRange-
                GetOpacityPixelComponent(p))));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*GetBluePixelComponent(p));
              *q++=(double) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(double) (QuantumScale*GetRedPixelComponent(p));
              *q++=0.0;
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*PixelIntensityToQuantum(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*GetRedPixelComponent(p));
              *q++=(double) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(double) (QuantumScale*GetBluePixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*GetRedPixelComponent(p));
              *q++=(double) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(double) (QuantumScale*GetBluePixelComponent(p));
              *q++=(double) (QuantumScale*((Quantum) (QuantumRange-
                GetOpacityPixelComponent(p))));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(double) (QuantumScale*GetRedPixelComponent(p));
              *q++=(double) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(double) (QuantumScale*GetBluePixelComponent(p));
              *q++=0.0;
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=(double) (QuantumScale*GetRedPixelComponent(p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=(double) (QuantumScale*GetGreenPixelComponent(p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=(double) (QuantumScale*GetBluePixelComponent(p));
                break;
              }
              case AlphaQuantum:
              {
                *q=(double) (QuantumScale*((Quantum) (QuantumRange-
                  GetOpacityPixelComponent(p))));
                break;
              }
              case OpacityQuantum:
              {
                *q=(double) (QuantumScale*GetOpacityPixelComponent(p));
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=(double) (QuantumScale*indexes[x]);
                break;
              }
              case IndexQuantum:
              {
                *q=(double) (QuantumScale*PixelIntensityToQuantum(p));
                break;
              }
              default:
                *q=0;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    case FloatPixel:
    {
      register float
        *q;

      q=(float *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*GetBluePixelComponent(p));
              *q++=(float) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(float) (QuantumScale*GetRedPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*GetBluePixelComponent(p));
              *q++=(float) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(float) (QuantumScale*GetRedPixelComponent(p));
              *q++=(float) (QuantumScale*(Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*GetBluePixelComponent(p));
              *q++=(float) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(float) (QuantumScale*GetRedPixelComponent(p));
              *q++=0.0;
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*PixelIntensityToQuantum(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*GetRedPixelComponent(p));
              *q++=(float) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(float) (QuantumScale*GetBluePixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*GetRedPixelComponent(p));
              *q++=(float) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(float) (QuantumScale*GetBluePixelComponent(p));
              *q++=(float) (QuantumScale*((Quantum) (GetAlphaPixelComponent(p))));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(float) (QuantumScale*GetRedPixelComponent(p));
              *q++=(float) (QuantumScale*GetGreenPixelComponent(p));
              *q++=(float) (QuantumScale*GetBluePixelComponent(p));
              *q++=0.0;
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=(float) (QuantumScale*GetRedPixelComponent(p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=(float) (QuantumScale*GetGreenPixelComponent(p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=(float) (QuantumScale*GetBluePixelComponent(p));
                break;
              }
              case AlphaQuantum:
              {
                *q=(float) (QuantumScale*((Quantum) (GetAlphaPixelComponent(p))));
                break;
              }
              case OpacityQuantum:
              {
                *q=(float) (QuantumScale*GetOpacityPixelComponent(p));
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=(float) (QuantumScale*indexes[x]);
                break;
              }
              case IndexQuantum:
              {
                *q=(float) (QuantumScale*PixelIntensityToQuantum(p));
                break;
              }
              default:
                *q=0;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    case IntegerPixel:
    {
      register unsigned int
        *q;

      q=(unsigned int *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong((Quantum) (QuantumRange-
                GetOpacityPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=0U;
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int)
                ScaleQuantumToLong(PixelIntensityToQuantum(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong((Quantum)
                (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=0U;
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=(unsigned int) ScaleQuantumToLong(GetRedPixelComponent(p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=(unsigned int) ScaleQuantumToLong(GetGreenPixelComponent(p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=(unsigned int) ScaleQuantumToLong(GetBluePixelComponent(p));
                break;
              }
              case AlphaQuantum:
              {
                *q=(unsigned int) ScaleQuantumToLong((Quantum) (QuantumRange-
                  GetOpacityPixelComponent(p)));
                break;
              }
              case OpacityQuantum:
              {
                *q=(unsigned int) ScaleQuantumToLong(GetOpacityPixelComponent(p));
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=(unsigned int) ScaleQuantumToLong(indexes[x]);
                break;
              }
              case IndexQuantum:
              {
                *q=(unsigned int)
                  ScaleQuantumToLong(PixelIntensityToQuantum(p));
                break;
              }
              default:
                *q=0;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    case LongPixel:
    {
      register size_t
        *q;

      q=(size_t *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToLong(GetRedPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=ScaleQuantumToLong((Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=0;
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(PixelIntensityToQuantum(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToLong(GetBluePixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=ScaleQuantumToLong((Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToLong(GetRedPixelComponent(p));
              *q++=ScaleQuantumToLong(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToLong(GetBluePixelComponent(p));
              *q++=0;
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=ScaleQuantumToLong(GetRedPixelComponent(p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=ScaleQuantumToLong(GetGreenPixelComponent(p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=ScaleQuantumToLong(GetBluePixelComponent(p));
                break;
              }
              case AlphaQuantum:
              {
                *q=ScaleQuantumToLong((Quantum) (GetAlphaPixelComponent(p)));
                break;
              }
              case OpacityQuantum:
              {
                *q=ScaleQuantumToLong(GetOpacityPixelComponent(p));
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=ScaleQuantumToLong(indexes[x]);
                break;
              }
              case IndexQuantum:
              {
                *q=ScaleQuantumToLong(PixelIntensityToQuantum(p));
                break;
              }
              default:
                break;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    case QuantumPixel:
    {
      register Quantum
        *q;

      q=(Quantum *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=GetBluePixelComponent(p);
              *q++=GetGreenPixelComponent(p);
              *q++=GetRedPixelComponent(p);
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=GetBluePixelComponent(p);
              *q++=GetGreenPixelComponent(p);
              *q++=GetRedPixelComponent(p);
              *q++=(Quantum) (GetAlphaPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=GetBluePixelComponent(p);
              *q++=GetGreenPixelComponent(p);
              *q++=GetRedPixelComponent(p);
              *q++=(Quantum) 0;
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=PixelIntensityToQuantum(p);
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=GetRedPixelComponent(p);
              *q++=GetGreenPixelComponent(p);
              *q++=GetBluePixelComponent(p);
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=GetRedPixelComponent(p);
              *q++=GetGreenPixelComponent(p);
              *q++=GetBluePixelComponent(p);
              *q++=(Quantum) (GetAlphaPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=GetRedPixelComponent(p);
              *q++=GetGreenPixelComponent(p);
              *q++=GetBluePixelComponent(p);
              *q++=(Quantum) 0;
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=(Quantum) 0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=GetRedPixelComponent(p);
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=GetGreenPixelComponent(p);
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=GetBluePixelComponent(p);
                break;
              }
              case AlphaQuantum:
              {
                *q=(Quantum) (GetAlphaPixelComponent(p));
                break;
              }
              case OpacityQuantum:
              {
                *q=GetOpacityPixelComponent(p);
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=indexes[x];
                break;
              }
              case IndexQuantum:
              {
                *q=(PixelIntensityToQuantum(p));
                break;
              }
              default:
                *q=(Quantum) 0;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    case ShortPixel:
    {
      register unsigned short
        *q;

      q=(unsigned short *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(GetBluePixelComponent(p));
              *q++=ScaleQuantumToShort(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToShort(GetRedPixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(GetBluePixelComponent(p));
              *q++=ScaleQuantumToShort(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToShort(GetRedPixelComponent(p));
              *q++=ScaleQuantumToShort((Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(GetBluePixelComponent(p));
              *q++=ScaleQuantumToShort(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToShort(GetRedPixelComponent(p));
              *q++=0;
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(PixelIntensityToQuantum(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(GetRedPixelComponent(p));
              *q++=ScaleQuantumToShort(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToShort(GetBluePixelComponent(p));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(GetRedPixelComponent(p));
              *q++=ScaleQuantumToShort(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToShort(GetBluePixelComponent(p));
              *q++=ScaleQuantumToShort((Quantum) (GetAlphaPixelComponent(p)));
              p++;
            }
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (p == (const PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              *q++=ScaleQuantumToShort(GetRedPixelComponent(p));
              *q++=ScaleQuantumToShort(GetGreenPixelComponent(p));
              *q++=ScaleQuantumToShort(GetBluePixelComponent(p));
              *q++=0;
              p++;
            }
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        p=GetVirtualPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          break;
        indexes=GetVirtualIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            *q=0;
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                *q=ScaleQuantumToShort(GetRedPixelComponent(p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                *q=ScaleQuantumToShort(GetGreenPixelComponent(p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                *q=ScaleQuantumToShort(GetBluePixelComponent(p));
                break;
              }
              case AlphaQuantum:
              {
                *q=ScaleQuantumToShort((Quantum) (GetAlphaPixelComponent(p)));
                break;
              }
              case OpacityQuantum:
              {
                *q=ScaleQuantumToShort(GetOpacityPixelComponent(p));
                break;
              }
              case BlackQuantum:
              {
                if (image->colorspace == CMYKColorspace)
                  *q=ScaleQuantumToShort(indexes[x]);
                break;
              }
              case IndexQuantum:
              {
                *q=ScaleQuantumToShort(PixelIntensityToQuantum(p));
                break;
              }
              default:
                break;
            }
            q++;
          }
          p++;
        }
      }
      break;
    }
    default:
    {
      quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
      (void) ThrowMagickException(exception,GetMagickModule(),OptionError,
        "UnrecognizedPixelMap","`%s'",map);
      break;
    }
  }
  quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t M a g i c k P i x e l P a c k e t                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetMagickPixelPacket() initializes the MagickPixelPacket structure.
%
%  The format of the GetMagickPixelPacket method is:
%
%      GetMagickPixelPacket(const Image *image,MagickPixelPacket *pixel)
%
%  A description of each parameter follows:
%
%    o image: the image.
%
%    o pixel: Specifies a pointer to a PixelPacket structure.
%
*/
MagickExport void GetMagickPixelPacket(const Image *image,
  MagickPixelPacket *pixel)
{
  pixel->storage_class=DirectClass;
  pixel->colorspace=RGBColorspace;
  pixel->matte=MagickFalse;
  pixel->fuzz=0.0;
  pixel->depth=MAGICKCORE_QUANTUM_DEPTH;
  pixel->red=0.0;
  pixel->green=0.0;
  pixel->blue=0.0;
  pixel->opacity=(MagickRealType) OpaqueOpacity;
  pixel->index=0.0;
  if (image == (const Image *) NULL)
    return;
  pixel->storage_class=image->storage_class;
  pixel->colorspace=image->colorspace;
  pixel->matte=image->matte;
  pixel->depth=image->depth;
  pixel->fuzz=image->fuzz;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   I m p o r t I m a g e P i x e l s                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ImportImagePixels() accepts pixel data and stores in the image at the
%  location you specify.  The method returns MagickTrue on success otherwise
%  MagickFalse if an error is encountered.  The pixel data can be either char,
%  short int, int, ssize_t, float, or double in the order specified by map.
%
%  Suppose your want to upload the first scanline of a 640x480 image from
%  character data in red-green-blue order:
%
%      ImportImagePixels(image,0,0,640,1,"RGB",CharPixel,pixels);
%
%  The format of the ImportImagePixels method is:
%
%      MagickBooleanType ImportImagePixels(Image *image,const ssize_t x_offset,
%        const ssize_t y_offset,const size_t columns,
%        const size_t rows,const char *map,const StorageType type,
%        const void *pixels)
%
%  A description of each parameter follows:
%
%    o image: the image.
%
%    o x_offset,y_offset,columns,rows:  These values define the perimeter
%      of a region of pixels you want to define.
%
%    o map:  This string reflects the expected ordering of the pixel array.
%      It can be any combination or order of R = red, G = green, B = blue,
%      A = alpha (0 is transparent), O = opacity (0 is opaque), C = cyan,
%      Y = yellow, M = magenta, K = black, I = intensity (for grayscale),
%      P = pad.
%
%    o type: Define the data type of the pixels.  Float and double types are
%      normalized to [0..1] otherwise [0..QuantumRange].  Choose from these
%      types: CharPixel, ShortPixel, IntegerPixel, LongPixel, FloatPixel, or
%      DoublePixel.
%
%    o pixels: This array of values contain the pixel components as defined by
%      map and type.  You must preallocate this array where the expected
%      length varies depending on the values of width, height, map, and type.
%
*/
MagickExport MagickBooleanType ImportImagePixels(Image *image,
  const ssize_t x_offset,const ssize_t y_offset,const size_t columns,
  const size_t rows,const char *map,const StorageType type,
  const void *pixels)
{
  ExceptionInfo
    *exception;

  PixelPacket
    *q;

  QuantumType
    *quantum_map;

  register IndexPacket
    *indexes;

  register ssize_t
    i,
    x;

  size_t
    length;

  ssize_t
    y;

  /*
    Allocate image structure.
  */
  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  if (image->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",image->filename);
  length=strlen(map);
  quantum_map=(QuantumType *) AcquireQuantumMemory(length,sizeof(*quantum_map));
  if (quantum_map == (QuantumType *) NULL)
    ThrowBinaryException(ResourceLimitError,"MemoryAllocationFailed",
      image->filename);
  for (i=0; i < (ssize_t) length; i++)
  {
    switch (map[i])
    {
      case 'a':
      case 'A':
      {
        quantum_map[i]=AlphaQuantum;
        image->matte=MagickTrue;
        break;
      }
      case 'B':
      case 'b':
      {
        quantum_map[i]=BlueQuantum;
        break;
      }
      case 'C':
      case 'c':
      {
        quantum_map[i]=CyanQuantum;
        (void) SetImageColorspace(image,CMYKColorspace);
        break;
      }
      case 'g':
      case 'G':
      {
        quantum_map[i]=GreenQuantum;
        break;
      }
      case 'K':
      case 'k':
      {
        quantum_map[i]=BlackQuantum;
        (void) SetImageColorspace(image,CMYKColorspace);
        break;
      }
      case 'I':
      case 'i':
      {
        quantum_map[i]=IndexQuantum;
        break;
      }
      case 'm':
      case 'M':
      {
        quantum_map[i]=MagentaQuantum;
        (void) SetImageColorspace(image,CMYKColorspace);
        break;
      }
      case 'O':
      case 'o':
      {
        quantum_map[i]=OpacityQuantum;
        image->matte=MagickTrue;
        break;
      }
      case 'P':
      case 'p':
      {
        quantum_map[i]=UndefinedQuantum;
        break;
      }
      case 'R':
      case 'r':
      {
        quantum_map[i]=RedQuantum;
        break;
      }
      case 'Y':
      case 'y':
      {
        quantum_map[i]=YellowQuantum;
        (void) SetImageColorspace(image,CMYKColorspace);
        break;
      }
      default:
      {
        quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
        (void) ThrowMagickException(&image->exception,GetMagickModule(),
          OptionError,"UnrecognizedPixelMap","`%s'",map);
        return(MagickFalse);
      }
    }
  }
  if (SetImageStorageClass(image,DirectClass) == MagickFalse)
    return(MagickFalse);
  /*
    Transfer the pixels from the pixel datarray to the image.
  */
  exception=(&image->exception);
  switch (type)
  {
    case CharPixel:
    {
      register const unsigned char
        *p;

      p=(const unsigned char *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->red=ScaleCharToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->red=ScaleCharToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleCharToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRO") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->red=ScaleCharToQuantum(*p++);
              q->opacity=ScaleCharToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->red=ScaleCharToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleCharToQuantum(*p++);
              q->green=q->red;
              q->blue=q->red;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->blue=ScaleCharToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->blue=ScaleCharToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleCharToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBO") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->blue=ScaleCharToQuantum(*p++);
              q->opacity=ScaleCharToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleCharToQuantum(*p++);
              q->green=ScaleCharToQuantum(*p++);
              q->blue=ScaleCharToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=ScaleCharToQuantum(*p);
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=ScaleCharToQuantum(*p);
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=ScaleCharToQuantum(*p);
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-ScaleCharToQuantum(*p);
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=ScaleCharToQuantum(*p);
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=ScaleCharToQuantum(*p);
                break;
              }
              case IndexQuantum:
              {
                q->red=ScaleCharToQuantum(*p);
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    case DoublePixel:
    {
      register const double
        *p;

      p=(const double *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->opacity=(Quantum) QuantumRange-ClampToQuantum((MagickRealType)
                QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              q->green=q->red;
              q->blue=q->red;
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->opacity=(Quantum) QuantumRange-ClampToQuantum((MagickRealType)
                QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-ClampToQuantum(
                  (MagickRealType) QuantumRange*(*p));
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case IndexQuantum:
              {
                q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    case FloatPixel:
    {
      register const float
        *p;

      p=(const float *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->opacity=(Quantum) QuantumRange-ClampToQuantum((MagickRealType)
                QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              q->green=q->red;
              q->blue=q->red;
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->opacity=(Quantum) QuantumRange-ClampToQuantum((MagickRealType)
                QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-ClampToQuantum(
                  (MagickRealType) QuantumRange*(*p));
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                break;
              }
              case IndexQuantum:
              {
                q->red=ClampToQuantum((MagickRealType) QuantumRange*(*p));
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    case IntegerPixel:
    {
      register const unsigned int
        *p;

      p=(const unsigned int *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->red=ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->red=ScaleLongToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->red=ScaleLongToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=q->red;
              q->blue=q->red;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->blue=ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->blue=ScaleLongToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->blue=ScaleLongToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=ScaleLongToQuantum(*p);
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=ScaleLongToQuantum(*p);
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=ScaleLongToQuantum(*p);
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-ScaleLongToQuantum(*p);
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=ScaleLongToQuantum(*p);
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=ScaleLongToQuantum(*p);
                break;
              }
              case IndexQuantum:
              {
                q->red=ScaleLongToQuantum(*p);
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    case LongPixel:
    {
      register const unsigned int
        *p;

      p=(const unsigned int *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->red=ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->red=ScaleLongToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->red=ScaleLongToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=q->red;
              q->blue=q->red;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->blue=ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->blue=ScaleLongToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleLongToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleLongToQuantum(*p++);
              q->green=ScaleLongToQuantum(*p++);
              q->blue=ScaleLongToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=ScaleLongToQuantum(*p);
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=ScaleLongToQuantum(*p);
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=ScaleLongToQuantum(*p);
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-ScaleLongToQuantum(*p);
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=ScaleLongToQuantum(*p);
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=ScaleLongToQuantum(*p);
                break;
              }
              case IndexQuantum:
              {
                q->red=ScaleLongToQuantum(*p);
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    case QuantumPixel:
    {
      register const Quantum
        *p;

      p=(const Quantum *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=(*p++);
              q->green=(*p++);
              q->red=(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=(*p++);
              q->green=(*p++);
              q->red=(*p++);
              q->opacity=(Quantum) QuantumRange-(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=(*p++);
              q->green=(*p++);
              q->red=(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=(*p++);
              q->green=q->red;
              q->blue=q->red;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=(*p++);
              q->green=(*p++);
              q->blue=(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=(*p++);
              q->green=(*p++);
              q->blue=(*p++);
              q->opacity=(Quantum) QuantumRange-(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=(*p++);
              q->green=(*p++);
              q->blue=(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=(*p);
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=(*p);
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=(*p);
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-(*p);
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=(*p);
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=(*p);
                break;
              }
              case IndexQuantum:
              {
                q->red=(*p);
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    case ShortPixel:
    {
      register const unsigned short
        *p;

      p=(const unsigned short *) pixels;
      if (LocaleCompare(map,"BGR") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleShortToQuantum(*p++);
              q->green=ScaleShortToQuantum(*p++);
              q->red=ScaleShortToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleShortToQuantum(*p++);
              q->green=ScaleShortToQuantum(*p++);
              q->red=ScaleShortToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleShortToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"BGRP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->blue=ScaleShortToQuantum(*p++);
              q->green=ScaleShortToQuantum(*p++);
              q->red=ScaleShortToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"I") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleShortToQuantum(*p++);
              q->green=q->red;
              q->blue=q->red;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGB") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleShortToQuantum(*p++);
              q->green=ScaleShortToQuantum(*p++);
              q->blue=ScaleShortToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBA") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleShortToQuantum(*p++);
              q->green=ScaleShortToQuantum(*p++);
              q->blue=ScaleShortToQuantum(*p++);
              q->opacity=(Quantum) QuantumRange-ScaleShortToQuantum(*p++);
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      if (LocaleCompare(map,"RGBP") == 0)
        {
          for (y=0; y < (ssize_t) rows; y++)
          {
            q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
            if (q == (PixelPacket *) NULL)
              break;
            for (x=0; x < (ssize_t) columns; x++)
            {
              q->red=ScaleShortToQuantum(*p++);
              q->green=ScaleShortToQuantum(*p++);
              q->blue=ScaleShortToQuantum(*p++);
              p++;
              q++;
            }
            if (SyncAuthenticPixels(image,exception) == MagickFalse)
              break;
          }
          break;
        }
      for (y=0; y < (ssize_t) rows; y++)
      {
        q=GetAuthenticPixels(image,x_offset,y_offset+y,columns,1,exception);
        if (q == (PixelPacket *) NULL)
          break;
        indexes=GetAuthenticIndexQueue(image);
        for (x=0; x < (ssize_t) columns; x++)
        {
          for (i=0; i < (ssize_t) length; i++)
          {
            switch (quantum_map[i])
            {
              case RedQuantum:
              case CyanQuantum:
              {
                q->red=ScaleShortToQuantum(*p);
                break;
              }
              case GreenQuantum:
              case MagentaQuantum:
              {
                q->green=ScaleShortToQuantum(*p);
                break;
              }
              case BlueQuantum:
              case YellowQuantum:
              {
                q->blue=ScaleShortToQuantum(*p);
                break;
              }
              case AlphaQuantum:
              {
                q->opacity=(Quantum) QuantumRange-ScaleShortToQuantum(*p);
                break;
              }
              case OpacityQuantum:
              {
                q->opacity=ScaleShortToQuantum(*p);
                break;
              }
              case BlackQuantum:
              {
                indexes[x]=ScaleShortToQuantum(*p);
                break;
              }
              case IndexQuantum:
              {
                q->red=ScaleShortToQuantum(*p);
                q->green=q->red;
                q->blue=q->red;
                break;
              }
              default:
                break;
            }
            p++;
          }
          q++;
        }
        if (SyncAuthenticPixels(image,exception) == MagickFalse)
          break;
      }
      break;
    }
    default:
    {
      quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
      (void) ThrowMagickException(&image->exception,GetMagickModule(),
        OptionError,"UnrecognizedPixelMap","`%s'",map);
      break;
    }
  }
  quantum_map=(QuantumType *) RelinquishMagickMemory(quantum_map);
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   I n t e r p o l a t e P i x e l P a c k e t                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  InterpolatePixelPacket() applies bi-linear or tri-linear interpolation
%  between a floating point coordinate and the pixels surrounding that
%  coordinate.  No pixel area resampling, or scaling of the result is
%  performed.
%
%  The format of the InterpolatePixelPacket method is:
%
%      MagickBooleanType InterpolatePixelPacket(const Image *image,
%        CacheView *image_view,const InterpolatePixelMethod method,
%        const double x,const double y,MagickPixelPacket *pixel,
%        ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o image: the image.
%
%    o image_view: the image view.
%
%    o method: the pixel color interpolation method.
%
%    o x,y: A double representing the current (x,y) position of the pixel.
%
%    o pixel: return the interpolated pixel here.
%
%    o exception: return any errors or warnings in this structure.
%
*/

static inline double MagickMax(const double x,const double y)
{
  if (x > y)
    return(x);
  return(y);
}

static void BicubicInterpolate(const MagickPixelPacket *pixels,const double dx,
  MagickPixelPacket *pixel)
{
  MagickRealType
    dx2,
    p,
    q,
    r,
    s;

  dx2=dx*dx;
  p=(pixels[3].red-pixels[2].red)-(pixels[0].red-pixels[1].red);
  q=(pixels[0].red-pixels[1].red)-p;
  r=pixels[2].red-pixels[0].red;
  s=pixels[1].red;
  pixel->red=(dx*dx2*p)+(dx2*q)+(dx*r)+s;
  p=(pixels[3].green-pixels[2].green)-(pixels[0].green-pixels[1].green);
  q=(pixels[0].green-pixels[1].green)-p;
  r=pixels[2].green-pixels[0].green;
  s=pixels[1].green;
  pixel->green=(dx*dx2*p)+(dx2*q)+(dx*r)+s;
  p=(pixels[3].blue-pixels[2].blue)-(pixels[0].blue-pixels[1].blue);
  q=(pixels[0].blue-pixels[1].blue)-p;
  r=pixels[2].blue-pixels[0].blue;
  s=pixels[1].blue;
  pixel->blue=(dx*dx2*p)+(dx2*q)+(dx*r)+s;
  p=(pixels[3].opacity-pixels[2].opacity)-(pixels[0].opacity-pixels[1].opacity);
  q=(pixels[0].opacity-pixels[1].opacity)-p;
  r=pixels[2].opacity-pixels[0].opacity;
  s=pixels[1].opacity;
  pixel->opacity=(dx*dx2*p)+(dx2*q)+(dx*r)+s;
  if (pixel->colorspace == CMYKColorspace)
    {
      p=(pixels[3].index-pixels[2].index)-(pixels[0].index-pixels[1].index);
      q=(pixels[0].index-pixels[1].index)-p;
      r=pixels[2].index-pixels[0].index;
      s=pixels[1].index;
      pixel->index=(dx*dx2*p)+(dx2*q)+(dx*r)+s;
    }
}

static inline MagickRealType CubicWeightingFunction(const MagickRealType x)
{
  MagickRealType
    alpha,
    gamma;

  alpha=MagickMax(x+2.0,0.0);
  gamma=1.0*alpha*alpha*alpha;
  alpha=MagickMax(x+1.0,0.0);
  gamma-=4.0*alpha*alpha*alpha;
  alpha=MagickMax(x+0.0,0.0);
  gamma+=6.0*alpha*alpha*alpha;
  alpha=MagickMax(x-1.0,0.0);
  gamma-=4.0*alpha*alpha*alpha;
  return(gamma/6.0);
}

static inline double MeshInterpolate(const PointInfo *delta,const double p,
  const double x,const double y)
{
  return(delta->x*x+delta->y*y+(1.0-delta->x-delta->y)*p);
}

static inline ssize_t NearestNeighbor(MagickRealType x)
{
  if (x >= 0.0)
    return((ssize_t) (x+0.5));
  return((ssize_t) (x-0.5));
}

MagickExport MagickBooleanType InterpolatePixelPacket(const Image *image,
  CacheView *image_view,const InterpolatePixelMethod method,const double x,
  const double y,MagickPixelPacket *pixel,ExceptionInfo *exception)
{
  MagickBooleanType
    status;

  register const IndexPacket
    *indexes;

  register const PixelPacket
    *p;

  register ssize_t
    i;

  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  assert(image_view != (CacheView *) NULL);
  status=MagickTrue;
  switch (method)
  {
    case AverageInterpolatePixel:
    {
      MagickPixelPacket
        pixels[16];

      MagickRealType
        alpha[16],
        gamma;

      p=GetCacheViewVirtualPixels(image_view,(ssize_t) floor(x)-1,
        (ssize_t) floor(y)-1,4,4,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      for (i=0; i < 16L; i++)
      {
        GetMagickPixelPacket(image,pixels+i);
        SetMagickPixelPacket(image,p,indexes+i,pixels+i);
        alpha[i]=1.0;
        if (pixels[i].matte != MagickFalse)
          {
            alpha[i]=QuantumScale*((MagickRealType) GetAlphaPixelComponent(p));
            pixels[i].red*=alpha[i];
            pixels[i].green*=alpha[i];
            pixels[i].blue*=alpha[i];
            if (pixels[i].colorspace == CMYKColorspace)
              pixels[i].index*=alpha[i];
          }
        gamma=alpha[i];
        gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
        pixel->red+=gamma*0.0625*pixels[i].red;
        pixel->green+=gamma*0.0625*pixels[i].green;
        pixel->blue+=gamma*0.0625*pixels[i].blue;
        pixel->opacity+=0.0625*pixels[i].opacity;
        if (pixel->colorspace == CMYKColorspace)
          pixel->index+=gamma*0.0625*pixels[i].index;
        p++;
      }
      break;
    }
    case BicubicInterpolatePixel:
    {
      MagickPixelPacket
        pixels[16],
        u[4];

      MagickRealType
        alpha[16];

      PointInfo
        delta;

      p=GetCacheViewVirtualPixels(image_view,(ssize_t) floor(x)-1,(ssize_t)
        floor(y)-1,4,4,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      for (i=0; i < 16L; i++)
      {
        GetMagickPixelPacket(image,pixels+i);
        SetMagickPixelPacket(image,p,indexes+i,pixels+i);
        alpha[i]=1.0;
        if (pixels[i].matte != MagickFalse)
          {
            alpha[i]=QuantumScale*((MagickRealType) GetAlphaPixelComponent(p));
            pixels[i].red*=alpha[i];
            pixels[i].green*=alpha[i];
            pixels[i].blue*=alpha[i];
            if (pixels[i].colorspace == CMYKColorspace)
              pixels[i].index*=alpha[i];
          }
        p++;
      }
      delta.x=x-floor(x);
      for (i=0; i < 4L; i++)
        BicubicInterpolate(pixels+4*i,delta.x,u+i);
      delta.y=y-floor(y);
      BicubicInterpolate(u,delta.y,pixel);
      break;
    }
    case BilinearInterpolatePixel:
    default:
    {
      MagickPixelPacket
        pixels[4];

      MagickRealType
        alpha[4],
        gamma;

      PointInfo
        delta,
        epsilon;

      p=GetCacheViewVirtualPixels(image_view,(ssize_t) floor(x),(ssize_t)
        floor(y),2,2,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      for (i=0; i < 4L; i++)
      {
        GetMagickPixelPacket(image,pixels+i);
        pixels[i].red=(MagickRealType) p[i].red;
        pixels[i].green=(MagickRealType) p[i].green;
        pixels[i].blue=(MagickRealType) p[i].blue;
        pixels[i].opacity=(MagickRealType) p[i].opacity;
        alpha[i]=1.0;
        if (pixels[i].matte != MagickFalse)
          {
            alpha[i]=QuantumScale*((MagickRealType) QuantumRange-p[i].opacity);
            pixels[i].red*=alpha[i];
            pixels[i].green*=alpha[i];
            pixels[i].blue*=alpha[i];
          }
        if (indexes != (IndexPacket *) NULL)
          {
            pixels[i].index=(MagickRealType) indexes[i];
            if (pixels[i].colorspace == CMYKColorspace)
              pixels[i].index*=alpha[i];
          }
      }
      delta.x=x-floor(x);
      delta.y=y-floor(y);
      epsilon.x=1.0-delta.x;
      epsilon.y=1.0-delta.y;
      gamma=((epsilon.y*(epsilon.x*alpha[0]+delta.x*alpha[1])+delta.y*
        (epsilon.x*alpha[2]+delta.x*alpha[3])));
      gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
      pixel->red=gamma*(epsilon.y*(epsilon.x*pixels[0].red+delta.x*
        pixels[1].red)+delta.y*(epsilon.x*pixels[2].red+delta.x*pixels[3].red));
      pixel->green=gamma*(epsilon.y*(epsilon.x*pixels[0].green+delta.x*
        pixels[1].green)+delta.y*(epsilon.x*pixels[2].green+delta.x*
        pixels[3].green));
      pixel->blue=gamma*(epsilon.y*(epsilon.x*pixels[0].blue+delta.x*
        pixels[1].blue)+delta.y*(epsilon.x*pixels[2].blue+delta.x*
        pixels[3].blue));
      pixel->opacity=(epsilon.y*(epsilon.x*pixels[0].opacity+delta.x*
        pixels[1].opacity)+delta.y*(epsilon.x*pixels[2].opacity+delta.x*
        pixels[3].opacity));
      if (pixel->colorspace == CMYKColorspace)
        pixel->index=gamma*(epsilon.y*(epsilon.x*pixels[0].index+delta.x*
          pixels[1].index)+delta.y*(epsilon.x*pixels[2].index+delta.x*
          pixels[3].index));
      break;
    }
    case FilterInterpolatePixel:
    {
      CacheView
        *filter_view;

      Image
        *excerpt_image,
        *filter_image;

      MagickPixelPacket
        pixels[1];

      RectangleInfo
        geometry;

      geometry.width=4L;
      geometry.height=4L;
      geometry.x=(ssize_t) floor(x)-1L;
      geometry.y=(ssize_t) floor(y)-1L;
      excerpt_image=ExcerptImage(image,&geometry,exception);
      if (excerpt_image == (Image *) NULL)
        {
          status=MagickFalse;
          break;
        }
      filter_image=ResizeImage(excerpt_image,1,1,image->filter,image->blur,
        exception);
      excerpt_image=DestroyImage(excerpt_image);
      if (filter_image == (Image *) NULL)
        break;
      filter_view=AcquireCacheView(filter_image);
      p=GetCacheViewVirtualPixels(filter_view,0,0,1,1,exception);
      if (p != (const PixelPacket *) NULL)
        {
          indexes=GetVirtualIndexQueue(filter_image);
          GetMagickPixelPacket(image,pixels);
          SetMagickPixelPacket(image,p,indexes,pixel);
        }
      filter_view=DestroyCacheView(filter_view);
      filter_image=DestroyImage(filter_image);
      break;
    }
    case IntegerInterpolatePixel:
    {
      MagickPixelPacket
        pixels[1];

      p=GetCacheViewVirtualPixels(image_view,(ssize_t) floor(x),(ssize_t)
        floor(y),1,1,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      GetMagickPixelPacket(image,pixels);
      SetMagickPixelPacket(image,p,indexes,pixel);
      break;
    }
    case MeshInterpolatePixel:
    {
      MagickPixelPacket
        pixels[4];

      MagickRealType
        alpha[4],
        gamma;

      PointInfo
        delta,
        luminance;

      p=GetCacheViewVirtualPixels(image_view,(ssize_t) floor(x),(ssize_t)
        floor(y),2,2,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      for (i=0; i < 4L; i++)
      {
        GetMagickPixelPacket(image,pixels+i);
        SetMagickPixelPacket(image,p,indexes+i,pixels+i);
        alpha[i]=1.0;
        if (pixels[i].matte != MagickFalse)
          {
            alpha[i]=QuantumScale*((MagickRealType) GetAlphaPixelComponent(p));
            pixels[i].red*=alpha[i];
            pixels[i].green*=alpha[i];
            pixels[i].blue*=alpha[i];
            if (pixels[i].colorspace == CMYKColorspace)
              pixels[i].index*=alpha[i];
          }
        p++;
      }
      delta.x=x-floor(x);
      delta.y=y-floor(y);
      luminance.x=MagickPixelLuminance(pixels+0)-MagickPixelLuminance(pixels+3);
      luminance.y=MagickPixelLuminance(pixels+1)-MagickPixelLuminance(pixels+2);
      if (fabs(luminance.x) < fabs(luminance.y))
        {
          /*
            Diagonal 0-3 NW-SE.
          */
          if (delta.x <= delta.y)
            {
              /*
                Bottom-left triangle  (pixel:2, diagonal: 0-3).
              */
              delta.y=1.0-delta.y;
              gamma=MeshInterpolate(&delta,alpha[2],alpha[3],alpha[0]);
              gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
              pixel->red=gamma*MeshInterpolate(&delta,pixels[2].red,
                pixels[3].red,pixels[0].red);
              pixel->green=gamma*MeshInterpolate(&delta,pixels[2].green,
                pixels[3].green,pixels[0].green);
              pixel->blue=gamma*MeshInterpolate(&delta,pixels[2].blue,
                pixels[3].blue,pixels[0].blue);
              pixel->opacity=gamma*MeshInterpolate(&delta,pixels[2].opacity,
                pixels[3].opacity,pixels[0].opacity);
              if (pixel->colorspace == CMYKColorspace)
                pixel->index=gamma*MeshInterpolate(&delta,pixels[2].index,
                  pixels[3].index,pixels[0].index);
            }
          else
            {
              /*
                Top-right triangle (pixel:1, diagonal: 0-3).
              */
              delta.x=1.0-delta.x;
              gamma=MeshInterpolate(&delta,alpha[1],alpha[0],alpha[3]);
              gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
              pixel->red=gamma*MeshInterpolate(&delta,pixels[1].red,
                pixels[0].red,pixels[3].red);
              pixel->green=gamma*MeshInterpolate(&delta,pixels[1].green,
                pixels[0].green,pixels[3].green);
              pixel->blue=gamma*MeshInterpolate(&delta,pixels[1].blue,
                pixels[0].blue,pixels[3].blue);
              pixel->opacity=gamma*MeshInterpolate(&delta,pixels[1].opacity,
                pixels[0].opacity,pixels[3].opacity);
              if (pixel->colorspace == CMYKColorspace)
                pixel->index=gamma*MeshInterpolate(&delta,pixels[1].index,
                  pixels[0].index,pixels[3].index);
            }
        }
      else
        {
          /*
            Diagonal 1-2 NE-SW.
          */
          if (delta.x <= (1.0-delta.y))
            {
              /*
                Top-left triangle (pixel 0, diagonal: 1-2).
              */
              gamma=MeshInterpolate(&delta,alpha[0],alpha[1],alpha[2]);
              gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
              pixel->red=gamma*MeshInterpolate(&delta,pixels[0].red,
                pixels[1].red,pixels[2].red);
              pixel->green=gamma*MeshInterpolate(&delta,pixels[0].green,
                pixels[1].green,pixels[2].green);
              pixel->blue=gamma*MeshInterpolate(&delta,pixels[0].blue,
                pixels[1].blue,pixels[2].blue);
              pixel->opacity=gamma*MeshInterpolate(&delta,pixels[0].opacity,
                pixels[1].opacity,pixels[2].opacity);
              if (pixel->colorspace == CMYKColorspace)
                pixel->index=gamma*MeshInterpolate(&delta,pixels[0].index,
                  pixels[1].index,pixels[2].index);
            }
          else
            {
              /*
                Bottom-right triangle (pixel: 3, diagonal: 1-2).
              */
              delta.x=1.0-delta.x;
              delta.y=1.0-delta.y;
              gamma=MeshInterpolate(&delta,alpha[3],alpha[2],alpha[1]);
              gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
              pixel->red=gamma*MeshInterpolate(&delta,pixels[3].red,
                pixels[2].red,pixels[1].red);
              pixel->green=gamma*MeshInterpolate(&delta,pixels[3].green,
                pixels[2].green,pixels[1].green);
              pixel->blue=gamma*MeshInterpolate(&delta,pixels[3].blue,
                pixels[2].blue,pixels[1].blue);
              pixel->opacity=gamma*MeshInterpolate(&delta,pixels[3].opacity,
                pixels[2].opacity,pixels[1].opacity);
              if (pixel->colorspace == CMYKColorspace)
                pixel->index=gamma*MeshInterpolate(&delta,pixels[3].index,
                  pixels[2].index,pixels[1].index);
            }
        }
      break;
    }
    case NearestNeighborInterpolatePixel:
    {
      MagickPixelPacket
        pixels[1];

      p=GetCacheViewVirtualPixels(image_view,NearestNeighbor(x),
        NearestNeighbor(y),1,1,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      GetMagickPixelPacket(image,pixels);
      SetMagickPixelPacket(image,p,indexes,pixel);
      break;
    }
    case SplineInterpolatePixel:
    {
      MagickPixelPacket
        pixels[16];

      MagickRealType
        alpha[16],
        dx,
        dy,
        gamma;

      PointInfo
        delta;

      ssize_t
        j,
        n;

      p=GetCacheViewVirtualPixels(image_view,(ssize_t) floor(x)-1,(ssize_t)
        floor(y)-1,4,4,exception);
      if (p == (const PixelPacket *) NULL)
        {
          status=MagickFalse;
          break;
        }
      indexes=GetCacheViewVirtualIndexQueue(image_view);
      n=0;
      delta.x=x-floor(x);
      delta.y=y-floor(y);
      for (i=(-1); i < 3L; i++)
      {
        dy=CubicWeightingFunction((MagickRealType) i-delta.y);
        for (j=(-1); j < 3L; j++)
        {
          GetMagickPixelPacket(image,pixels+n);
          SetMagickPixelPacket(image,p,indexes+n,pixels+n);
          alpha[n]=1.0;
          if (pixels[i].matte != MagickFalse)
            {
              alpha[n]=QuantumScale*((MagickRealType)
                GetAlphaPixelComponent(p));
              pixels[n].red*=alpha[n];
              pixels[n].green*=alpha[n];
              pixels[n].blue*=alpha[n];
              if (pixels[i].colorspace == CMYKColorspace)
                pixels[n].index*=alpha[n];
            }
          dx=CubicWeightingFunction(delta.x-(MagickRealType) j);
          gamma=alpha[n];
          gamma=1.0/(fabs((double) gamma) <= MagickEpsilon ? 1.0 : gamma);
          pixel->red+=gamma*dx*dy*pixels[n].red;
          pixel->green+=gamma*dx*dy*pixels[n].green;
          pixel->blue+=gamma*dx*dy*pixels[n].blue;
          if (pixel->matte != MagickFalse)
            pixel->opacity+=dx*dy*pixels[n].opacity;
          if (pixel->colorspace == CMYKColorspace)
            pixel->index+=gamma*dx*dy*pixels[n].index;
          n++;
          p++;
        }
      }
      break;
    }
  }
  return(status);
}
