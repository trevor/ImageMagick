/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%                        FFFFF  IIIII  TTTTT  SSSSS                           %
%                        F        I      T    SS                              %
%                        FFF      I      T     SSS                            %
%                        F        I      T       SS                           %
%                        F      IIIII    T    SSSSS                           %
%                                                                             %
%                                                                             %
%            Read/Write Flexible Image Transport System Images.               %
%                                                                             %
%                              Software Design                                %
%                                John Cristy                                  %
%                                 July 1992                                   %
%                                                                             %
%                                                                             %
%  Copyright 1999-2009 ImageMagick Studio LLC, a non-profit organization      %
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
#include "magick/cache.h"
#include "magick/color-private.h"
#include "magick/colorspace.h"
#include "magick/constitute.h"
#include "magick/exception.h"
#include "magick/exception-private.h"
#include "magick/image.h"
#include "magick/image-private.h"
#include "magick/list.h"
#include "magick/magick.h"
#include "magick/memory_.h"
#include "magick/monitor.h"
#include "magick/monitor-private.h"
#include "magick/pixel-private.h"
#include "magick/pixel-private.h"
#include "magick/static.h"
#include "magick/statistic.h"
#include "magick/string_.h"
#include "magick/module.h"
#include "magick/module.h"

/*
  Forward declarations.
*/
#define FITSBlocksize  2880UL

/*
  Forward declarations.
*/
static MagickBooleanType
  WriteFITSImage(const ImageInfo *,Image *);

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   I s F I T S                                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  IsFITS() returns MagickTrue if the image format type, identified by the
%  magick string, is FITS.
%
%  The format of the IsFITS method is:
%
%      MagickBooleanType IsFITS(const unsigned char *magick,const size_t length)
%
%  A description of each parameter follows:
%
%    o magick: compare image format pattern against these bytes.
%
%    o length: Specifies the length of the magick string.
%
*/
static MagickBooleanType IsFITS(const unsigned char *magick,const size_t length)
{
  if (length < 6)
    return(MagickFalse);
  if (LocaleNCompare((char *) magick,"IT0",3) == 0)
    return(MagickTrue);
  if (LocaleNCompare((char *) magick,"SIMPLE",6) == 0)
    return(MagickTrue);
  return(MagickFalse);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   R e a d F I T S I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ReadFITSImage() reads a FITS image file and returns it.  It allocates the
%  memory necessary for the new Image structure and returns a pointer to the
%  new image.
%
%  The format of the ReadFITSImage method is:
%
%      Image *ReadFITSImage(const ImageInfo *image_info,
%        ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o image_info: the image info.
%
%    o exception: return any errors or warnings in this structure.
%
*/

static inline double GetFITSPixel(Image *image,int bits_per_pixel)
{
  switch (image->depth >> 3)
  {
    case 1:
      return((double) ReadBlobByte(image));
    case 2:
      return((double) ((short) ReadBlobShort(image)));
    case 4:
    {
      if (bits_per_pixel > 0)
        return((double) ((long) ReadBlobLong(image)));
      return((double) ReadBlobFloat(image));
    }
    case 8:
    {
      if (bits_per_pixel > 0)
        return((double) ((MagickOffsetType) ReadBlobLongLong(image)));
    }
    default:
      break;
  }
  return(ReadBlobDouble(image));
}

static void GetFITSPixelExtrema(Image *image,const int bits_per_pixel,
  double *minima,double *maxima)
{
  double
    pixel;

  MagickOffsetType
    offset;

  MagickSizeType
    number_pixels;

  register MagickOffsetType
    i;

  offset=TellBlob(image);
  number_pixels=(MagickSizeType) image->columns*image->rows;
  *minima=GetFITSPixel(image,bits_per_pixel);
  *maxima=(*minima);
  for (i=1; i < (MagickOffsetType) number_pixels; i++)
  {
    pixel=GetFITSPixel(image,bits_per_pixel);
    if (pixel < *minima)
      *minima=pixel;
    if (pixel > *maxima)
      *maxima=pixel;
  }
  (void) SeekBlob(image,offset,SEEK_SET);
}

static inline double GetFITSPixelRange(const unsigned long depth)
{
  return((double) ((MagickOffsetType) GetQuantumRange(depth)));
}

static void SetFITSUnsignedPixels(const size_t length,
  const unsigned long bits_per_pixel,unsigned char *pixels)
{
  register long
    i;

  for (i=0; i < (long) length; i++)
  {
    *pixels^=0x80;
    pixels+=bits_per_pixel >> 3;
  }
}

static Image *ReadFITSImage(const ImageInfo *image_info,
  ExceptionInfo *exception)
{
  typedef struct _FITSInfo
  {
    EndianType
      endian;

    int
      simple,
      bits_per_pixel,
      columns,
      rows,
      number_axes,
      number_planes;

    double
      min_data,
      max_data,
      zero,
      scale;
  } FITSInfo;

  char
    keyword[FITSBlocksize],
    value[FITSBlocksize];

  double
    pixel,
    scale;

  FITSInfo
    fits_info;

  Image
    *image;

  int
    c;

  long
    scene,
    y;

  MagickBooleanType
    status,
    value_expected;

  MagickSizeType
    number_pixels;

  register long
    x;

  register PixelPacket
    *q;

  /*
    Open image file.
  */
  assert(image_info != (const ImageInfo *) NULL);
  assert(image_info->signature == MagickSignature);
  if (image_info->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      image_info->filename);
  assert(exception != (ExceptionInfo *) NULL);
  assert(exception->signature == MagickSignature);
  image=AcquireImage(image_info);
  status=OpenBlob(image_info,image,ReadBinaryBlobMode,exception);
  if (status == MagickFalse)
    {
      image=DestroyImageList(image);
      return((Image *) NULL);
    }
  /*
    Initialize image header.
  */
  (void) ResetMagickMemory(&fits_info,0,sizeof(fits_info));
  fits_info.bits_per_pixel=8;
  fits_info.columns=1;
  fits_info.rows=1;
  fits_info.rows=1;
  fits_info.number_planes=1;
  fits_info.min_data=0.0;
  fits_info.max_data=0.0;
  fits_info.zero=0.0;
  fits_info.scale=1.0;
  fits_info.endian=MSBEndian;
  /*
    Decode image header.
  */
  c=ReadBlobByte(image);
  if (c == EOF)
    {
      image=DestroyImage(image);
      return((Image *) NULL);
    }
  for ( ; c != EOF; )
  {
    if (isalnum((int) ((unsigned char) c)) == 0)
      c=ReadBlobByte(image);
    else
      {
        register char
          *p;

        /*
          Determine a keyword and its value.
        */
        p=keyword;
        do
        {
          if ((size_t) (p-keyword) < (FITSBlocksize-1))
            *p++=(char) c;
          c=ReadBlobByte(image);
        } while ((isalnum((int) ((unsigned char) c)) != MagickFalse) ||
                 ((int) c == '_'));
        *p='\0';
        if (LocaleCompare(keyword,"END") == 0)
          break;
        value_expected=MagickFalse;
        while ((isspace((int) ((unsigned char) c)) != 0) || ((int) c == '='))
        {
          if ((int) c == '=')
            value_expected=MagickTrue;
          c=ReadBlobByte(image);
        }
        if (value_expected == MagickFalse)
          continue;
        p=value;
        while (isalnum(c) || (c == '-') || (c == '+') || (c == '.'))
        {
          if ((size_t) (p-value) < (FITSBlocksize-1))
            *p++=c;
          c=ReadBlobByte(image);
        }
        *p='\0';
        /*
          Assign a value to the specified keyword.
        */
        if (LocaleCompare(keyword,"SIMPLE") == 0)
          fits_info.simple=(int) ((*value == 'T') || (*value == 't'));
        if (LocaleCompare(keyword,"BITPIX") == 0)
          fits_info.bits_per_pixel=atoi(value);
        if (LocaleCompare(keyword,"NAXIS") == 0)
          fits_info.number_axes=atoi(value);
        if (LocaleCompare(keyword,"NAXIS1") == 0)
          fits_info.columns=atoi(value);
        if (LocaleCompare(keyword,"NAXIS2") == 0)
          fits_info.rows=atoi(value);
        if (LocaleCompare(keyword,"NAXIS3") == 0)
          fits_info.number_planes=atoi(value);
        if (LocaleCompare(keyword,"DATAMAX") == 0)
          fits_info.max_data=atof(value);
        if (LocaleCompare(keyword,"DATAMIN") == 0)
          fits_info.min_data=atof(value);
        if (LocaleCompare(keyword,"BZERO") == 0)
          fits_info.zero=atof(value);
        if (LocaleCompare(keyword,"BSCALE") == 0)
          fits_info.scale=atof(value);
        if (LocaleCompare(keyword,"XENDIAN") == 0)
          {
            if (LocaleCompare(keyword,"BIG") == 0)
              fits_info.endian=MSBEndian;
            else
              fits_info.endian=LSBEndian;
          }
        (void) SetImageProperty(image,keyword,value);
      }
    while (((TellBlob(image) % 80) != 0) && (c != EOF))
      c=ReadBlobByte(image);
    c=ReadBlobByte(image);
  }
  while (((TellBlob(image) % FITSBlocksize) != 0) && (c != EOF))
    c=ReadBlobByte(image);
  /*
    Verify that required image information is defined.
  */
  number_pixels=(MagickSizeType) fits_info.columns*fits_info.rows;
  if ((fits_info.simple == 0) || (fits_info.number_axes < 1) ||
      (fits_info.number_axes > 4) || (number_pixels == 0))
    ThrowReaderException(CorruptImageError,"ImageTypeNotSupported");
  for (scene=0; scene < (long) fits_info.number_planes; scene++)
  {
    image->columns=(unsigned long) fits_info.columns;
    image->rows=(unsigned long) fits_info.rows;
    image->depth=(unsigned long) (fits_info.bits_per_pixel < 0 ? -1 : 1)*
      fits_info.bits_per_pixel;
    image->endian=fits_info.endian;
    image->scene=(unsigned long) scene;
    if ((image_info->ping != MagickFalse) && (image_info->number_scenes != 0))
      if (image->scene >= (image_info->scene+image_info->number_scenes-1))
        break;
    /*
      Initialize image structure.
    */
    if ((fits_info.min_data != 0.0) || (fits_info.max_data != 0.0))
      {
        if ((fits_info.bits_per_pixel != 0) && (fits_info.max_data == 0.0))
          fits_info.max_data=GetFITSPixelRange((unsigned long)
            fits_info.bits_per_pixel);
      }
    else
      GetFITSPixelExtrema(image,fits_info.bits_per_pixel,&fits_info.min_data,
        &fits_info.max_data);
    /*
      Convert FITS pixels to pixel packets.
    */
    scale=(double) QuantumRange/(fits_info.scale*(fits_info.max_data-
      fits_info.min_data)+fits_info.zero);
    for (y=(long) image->rows-1; y >= 0; y--)
    {
      q=QueueAuthenticPixels(image,0,y,image->columns,1,exception);
      if (q == (PixelPacket *) NULL)
        break;
      for (x=0; x < (long) image->columns; x++)
      {
        pixel=GetFITSPixel(image,fits_info.bits_per_pixel);
        q->red=(Quantum) RoundToQuantum(scale*(fits_info.scale*(pixel-
          fits_info.min_data)+fits_info.zero));
        q->green=q->red;
        q->blue=q->red;
        q++;
      }
      if (SyncAuthenticPixels(image,exception) == MagickFalse)
        break;
      if (image->previous == (Image *) NULL)
        {
          status=SetImageProgress(image,LoadImageTag,y,image->rows);
          if (status == MagickFalse)
            break;
        }
    }
    if (EOFBlob(image) != MagickFalse)
      {
        ThrowFileException(exception,CorruptImageError,"UnexpectedEndOfFile",
          image->filename);
        break;
      }
    /*
      Proceed to next image.
    */
    if (image_info->number_scenes != 0)
      if (image->scene >= (image_info->scene+image_info->number_scenes-1))
        break;
    if (scene < (long) (fits_info.number_planes-1))
      {
        /*
          Allocate next image structure.
        */
        AcquireNextImage(image_info,image);
        if (GetNextImageInList(image) == (Image *) NULL)
          {
            image=DestroyImageList(image);
            return((Image *) NULL);
          }
        image=SyncNextImageInList(image);
        status=SetImageProgress(image,LoadImagesTag,TellBlob(image),
          GetBlobSize(image));
        if (status == MagickFalse)
          break;
      }
  }
  (void) CloseBlob(image);
  return(GetFirstImageInList(image));
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   R e g i s t e r F I T S I m a g e                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  RegisterFITSImage() adds attributes for the FITS image format to
%  the list of supported formats.  The attributes include the image format
%  tag, a method to read and/or write the format, whether the format
%  supports the saving of more than one frame to the same file or blob,
%  whether the format supports native in-memory I/O, and a brief
%  description of the format.
%
%  The format of the RegisterFITSImage method is:
%
%      unsigned long RegisterFITSImage(void)
%
*/
ModuleExport unsigned long RegisterFITSImage(void)
{
  MagickInfo
    *entry;

  entry=SetMagickInfo("FITS");
  entry->decoder=(DecodeImageHandler *) ReadFITSImage;
  entry->encoder=(EncodeImageHandler *) WriteFITSImage;
  entry->magick=(IsImageFormatHandler *) IsFITS;
  entry->adjoin=MagickFalse;
  entry->seekable_stream=MagickTrue;
  entry->description=ConstantString("Flexible Image Transport System");
  entry->module=ConstantString("FITS");
  (void) RegisterMagickInfo(entry);
  entry=SetMagickInfo("FTS");
  entry->decoder=(DecodeImageHandler *) ReadFITSImage;
  entry->encoder=(EncodeImageHandler *) WriteFITSImage;
  entry->magick=(IsImageFormatHandler *) IsFITS;
  entry->adjoin=MagickFalse;
  entry->seekable_stream=MagickTrue;
  entry->description=ConstantString("Flexible Image Transport System");
  entry->module=ConstantString("FTS");
  (void) RegisterMagickInfo(entry);
  return(MagickImageCoderSignature);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   U n r e g i s t e r F I T S I m a g e                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  UnregisterFITSImage() removes format registrations made by the
%  FITS module from the list of supported formats.
%
%  The format of the UnregisterFITSImage method is:
%
%      UnregisterFITSImage(void)
%
*/
ModuleExport void UnregisterFITSImage(void)
{
  (void) UnregisterMagickInfo("FITS");
  (void) UnregisterMagickInfo("FTS");
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   W r i t e F I T S I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  WriteFITSImage() writes a Flexible Image Transport System image to a
%  file as gray scale intensities [0..255].
%
%  The format of the WriteFITSImage method is:
%
%      MagickBooleanType WriteFITSImage(const ImageInfo *image_info,
%        Image *image)
%
%  A description of each parameter follows.
%
%    o image_info: the image info.
%
%    o image:  The image.
%
*/
static MagickBooleanType WriteFITSImage(const ImageInfo *image_info,
  Image *image)
{
  char
    header[FITSBlocksize],
    *fits_info;

  long
    y;

  MagickBooleanType
    status;

  QuantumInfo
    *quantum_info;

  register const PixelPacket
    *p;

  size_t
    length;

  ssize_t
    count,
    offset;

  unsigned char
    *pixels;

  /*
    Open output image file.
  */
  assert(image_info != (const ImageInfo *) NULL);
  assert(image_info->signature == MagickSignature);
  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  if (image->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",image->filename);
  status=OpenBlob(image_info,image,WriteBinaryBlobMode,&image->exception);
  if (status == MagickFalse)
    return(status);
  if (image->colorspace != RGBColorspace)
    (void) TransformImageColorspace(image,RGBColorspace);
  /*
    Allocate image memory.
  */
  fits_info=(char *) AcquireQuantumMemory(FITSBlocksize,sizeof(*fits_info));
  if (fits_info == (char *) NULL)
    ThrowWriterException(ResourceLimitError,"MemoryAllocationFailed");
  (void) ResetMagickMemory(fits_info,' ',FITSBlocksize*sizeof(*fits_info));
  /*
    Initialize image header.
  */
  image->depth=GetImageQuantumDepth(image,MagickFalse);
  quantum_info=AcquireQuantumInfo((const ImageInfo *) NULL,image);
  if (quantum_info == (QuantumInfo *) NULL)
    ThrowWriterException(ResourceLimitError,"MemoryAllocationFailed");
  offset=0;
  (void) FormatMagickString(header,FITSBlocksize,
    "SIMPLE  =                    T");
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"BITPIX  =           %10ld",
    (quantum_info->format == FloatingPointQuantumFormat ? -1 : 1)*image->depth);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"NAXIS   =           %10lu",
    2UL);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"NAXIS1  =           %10lu",
    image->columns);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"NAXIS2  =           %10lu",
    image->rows);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"BSCALE  =         %E",1.0);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"BZERO   =         %E",
    image->depth > 8 ? GetFITSPixelRange(image->depth) : 0.0);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"DATAMAX =         %E",
    1.0*((MagickOffsetType) GetQuantumRange(image->depth)));
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) FormatMagickString(header,FITSBlocksize,"DATAMIN =         %E",0.0);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  if (image->endian == LSBEndian)
    {
      (void) FormatMagickString(header,FITSBlocksize,"XENDIAN = 'SMALL'");
      (void) strncpy(fits_info+offset,header,strlen(header));
      offset+=80;
    }
  (void) FormatMagickString(header,FITSBlocksize,"HISTORY %.72s",
    GetMagickVersion((unsigned long *) NULL));
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) strncpy(header,"END",FITSBlocksize);
  (void) strncpy(fits_info+offset,header,strlen(header));
  offset+=80;
  (void) WriteBlob(image,FITSBlocksize,(unsigned char *) fits_info);
  /*
    Convert image to fits scale PseudoColor class.
  */
  pixels=GetQuantumPixels(quantum_info);
  length=GetQuantumExtent(image,quantum_info,GrayQuantum);
  for (y=(long) image->rows-1; y >= 0; y--)
  {
    p=GetVirtualPixels(image,0,y,image->columns,1,&image->exception);
    if (p == (const PixelPacket *) NULL)
      break;
    length=ExportQuantumPixels(image,(const ViewInfo *) NULL,quantum_info,
      GrayQuantum,pixels,&image->exception);
    if (image->depth == 16)
      SetFITSUnsignedPixels(image->columns,image->depth,pixels);
    if (((image->depth == 32) || (image->depth == 64)) &&
        (quantum_info->format != FloatingPointQuantumFormat))
      SetFITSUnsignedPixels(image->columns,image->depth,pixels);
    count=WriteBlob(image,length,pixels);
    if (count != (ssize_t) length)
      break;
    status=SetImageProgress(image,SaveImageTag,y,image->rows);
    if (status == MagickFalse)
      break;
  }
  quantum_info=DestroyQuantumInfo(quantum_info);
  length=(size_t) (FITSBlocksize-TellBlob(image) % FITSBlocksize);
  if (length != 0)
    {
      (void) ResetMagickMemory(fits_info,0,length*sizeof(*fits_info));
      (void) WriteBlob(image,length,(unsigned char *) fits_info);
    }
  fits_info=DestroyString(fits_info);
  (void) CloseBlob(image);
  return(MagickTrue);
}
