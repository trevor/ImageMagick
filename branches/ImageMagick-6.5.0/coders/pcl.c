/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%                            PPPP    CCCC  L                                  %
%                            P   P  C      L                                  %
%                            PPPP   C      L                                  %
%                            P      C      L                                  %
%                            P       CCCC  LLLLL                              %
%                                                                             %
%                                                                             %
%                      Read/Write HP PCL Printer Format                       %
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
#include "magick/color.h"
#include "magick/color-private.h"
#include "magick/colorspace.h"
#include "magick/constitute.h"
#include "magick/delegate.h"
#include "magick/draw.h"
#include "magick/exception.h"
#include "magick/exception-private.h"
#include "magick/geometry.h"
#include "magick/image.h"
#include "magick/image-private.h"
#include "magick/list.h"
#include "magick/magick.h"
#include "magick/memory_.h"
#include "magick/monitor.h"
#include "magick/monitor-private.h"
#include "magick/profile.h"
#include "magick/resource_.h"
#include "magick/quantum-private.h"
#include "magick/static.h"
#include "magick/string_.h"
#include "magick/module.h"
#include "magick/token.h"
#include "magick/transform.h"
#include "magick/utility.h"

/*
  Forward declarations.
*/
static MagickBooleanType
  WritePCLImage(const ImageInfo *,Image *);

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   I s P C L                                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  IsPCL() returns MagickTrue if the image format type, identified by the
%  magick string, is PCL.
%
%  The format of the IsPCL method is:
%
%      MagickBooleanType IsPCL(const unsigned char *magick,const size_t length)
%
%  A description of each parameter follows:
%
%    o magick: This string is generally the first few bytes of an image file
%      or blob.
%
%    o length: Specifies the length of the magick string.
%
%
*/
static MagickBooleanType IsPCL(const unsigned char *magick,const size_t length)
{
  if (length < 4)
    return(MagickFalse);
  if (memcmp(magick,"\033E\033",3) == 0)
    return(MagickTrue);
  if (memcmp(magick,"\033E\033&",4) == 0)
    return(MagickFalse);
  return(MagickFalse);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   R e a d P C L I m a g e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ReadPCLImage() reads a Printer Control Language image file and returns it.
%  It allocates the memory necessary for the new Image structure and returns a
%  pointer to the new image.
%
%  The format of the ReadPCLImage method is:
%
%      Image *ReadPCLImage(const ImageInfo *image_info,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o image_info: the image info.
%
%    o exception: return any errors or warnings in this structure.
%
*/
static Image *ReadPCLImage(const ImageInfo *image_info,ExceptionInfo *exception)
{
#define CropBox  "CropBox"
#define DeviceCMYK  "DeviceCMYK"
#define MediaBox  "MediaBox"
#define RenderPCLText  "  Rendering PCL...  "

  char
    command[MaxTextExtent],
    density[MaxTextExtent],
    filename[MaxTextExtent],
    geometry[MaxTextExtent],
    options[MaxTextExtent],
    input_filename[MaxTextExtent];

  const DelegateInfo
    *delegate_info;

  Image
    *image,
    *next_image;

  ImageInfo
    *read_info;

  MagickBooleanType
    cmyk,
    status;

  PointInfo
    delta;

  RectangleInfo
    bounding_box,
    page;

  register char
    *p;

  register long
    c;

  SegmentInfo
    bounds;

  ssize_t
    count;

  unsigned long
    height,
    width;

  assert(image_info != (const ImageInfo *) NULL);
  assert(image_info->signature == MagickSignature);
  if (image_info->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      image_info->filename);
  assert(exception != (ExceptionInfo *) NULL);
  assert(exception->signature == MagickSignature);
  /*
    Open image file.
  */
  image=AcquireImage(image_info);
  status=OpenBlob(image_info,image,ReadBinaryBlobMode,exception);
  if (status == MagickFalse)
    {
      image=DestroyImageList(image);
      return((Image *) NULL);
    }
  status=AcquireUniqueSymbolicLink(image_info->filename,input_filename);
  if (status == MagickFalse)
    {
      ThrowFileException(exception,FileOpenError,"UnableToCreateTemporaryFile",
        image_info->filename);
      image=DestroyImageList(image);
      return((Image *) NULL);
    }
  /*
    Set the page density.
  */
  delta.x=DefaultResolution;
  delta.y=DefaultResolution;
  if ((image->x_resolution == 0.0) || (image->y_resolution == 0.0))
    {
      GeometryInfo
        geometry_info;

      MagickStatusType
        flags;

      flags=ParseGeometry(PSDensityGeometry,&geometry_info);
      image->x_resolution=geometry_info.rho;
      image->y_resolution=geometry_info.sigma;
      if ((flags & SigmaValue) == 0)
        image->y_resolution=image->x_resolution;
    }
  /*
    Determine page geometry from the PCL media box.
  */
  cmyk=image->colorspace == CMYKColorspace ? MagickTrue : MagickFalse;
  count=0;
  (void) ResetMagickMemory(&bounding_box,0,sizeof(bounding_box));
  (void) ResetMagickMemory(&bounds,0,sizeof(bounds));
  (void) ResetMagickMemory(&page,0,sizeof(page));
  (void) ResetMagickMemory(command,0,sizeof(command));
  p=command;
  for (c=ReadBlobByte(image); c != EOF; c=ReadBlobByte(image))
  {
    if (image_info->page != (char *) NULL)
      continue;
    /*
      Note PCL elements.
    */
    *p++=(char) c;
    if ((c != (int) '/') && (c != '\n') &&
        ((size_t) (p-command) < (MaxTextExtent-1)))
      continue;
    *p='\0';
    p=command;
    /*
      Is this a CMYK document?
    */
    if (LocaleNCompare(DeviceCMYK,command,strlen(DeviceCMYK)) == 0)
      cmyk=MagickTrue;
    if (LocaleNCompare(CropBox,command,strlen(CropBox)) == 0)
      {
        /*
          Note region defined by crop box.
        */
        count=(ssize_t) sscanf(command,"CropBox [%lf %lf %lf %lf",
          &bounds.x1,&bounds.y1,&bounds.x2,&bounds.y2);
        if (count != 4)
          count=(ssize_t) sscanf(command,"CropBox[%lf %lf %lf %lf",
            &bounds.x1,&bounds.y1,&bounds.x2,&bounds.y2);
      }
    if (LocaleNCompare(MediaBox,command,strlen(MediaBox)) == 0)
      {
        /*
          Note region defined by media box.
        */
        count=(ssize_t) sscanf(command,"MediaBox [%lf %lf %lf %lf",
          &bounds.x1,&bounds.y1,&bounds.x2,&bounds.y2);
        if (count != 4)
          count=(ssize_t) sscanf(command,"MediaBox[%lf %lf %lf %lf",
            &bounds.x1,&bounds.y1,&bounds.x2,&bounds.y2);
      }
    if (count != 4)
      continue;
    /*
      Set PCL render geometry.
    */
    width=(unsigned long) (bounds.x2-bounds.x1+0.5);
    height=(unsigned long) (bounds.y2-bounds.y1+0.5);
    if (width > page.width)
      page.width=width;
    if (height > page.height)
      page.height=height;
  }
  (void) CloseBlob(image);
  /*
    Render PCL with the GhostPCL delegate.
  */
  if ((page.width == 0) || (page.height == 0))
    (void) ParseAbsoluteGeometry(PSPageGeometry,&page);
  if (image_info->page != (char *) NULL)
    (void) ParseAbsoluteGeometry(image_info->page,&page);
  (void) FormatMagickString(geometry,MaxTextExtent,"%lux%lu",
    page.width,page.height);
  if (image_info->monochrome != MagickFalse)
    delegate_info=GetDelegateInfo("pcl:mono",(char *) NULL,exception);
  else
     if (cmyk != MagickFalse)
       delegate_info=GetDelegateInfo("pcl:cmyk",(char *) NULL,exception);
     else
       delegate_info=GetDelegateInfo("pcl:color",(char *) NULL,exception);
  if (delegate_info == (const DelegateInfo *) NULL)
    return((Image *) NULL);
  *options='\0';
  if ((page.width == 0) || (page.height == 0))
    (void) ParseAbsoluteGeometry(PSPageGeometry,&page);
  if (image_info->page != (char *) NULL)
    (void) ParseAbsoluteGeometry(image_info->page,&page);
  (void) FormatMagickString(density,MaxTextExtent,"%gx%g",
    image->x_resolution,image->y_resolution);
  page.width=(unsigned long) (page.width*image->x_resolution/delta.x+0.5);
  page.height=(unsigned long) (page.height*image->y_resolution/delta.y+0.5);
  (void) FormatMagickString(options,MaxTextExtent,"-g%lux%lu ",
    page.width,page.height);
  image=DestroyImage(image);
  read_info=CloneImageInfo(image_info);
  *read_info->magick='\0';
  if (read_info->number_scenes != 0)
    {
      if (read_info->number_scenes != 1)
        (void) FormatMagickString(options,MaxTextExtent,"-dLastPage=%lu",
          read_info->scene+read_info->number_scenes);
      else
        (void) FormatMagickString(options,MaxTextExtent,
          "-dFirstPage=%lu -dLastPage=%lu",read_info->scene+1,read_info->scene+
          read_info->number_scenes);
      read_info->number_scenes=0;
      if (read_info->scenes != (char *) NULL)
        *read_info->scenes='\0';
    }
  if (read_info->authenticate != (char *) NULL)
    (void) FormatMagickString(options+strlen(options),MaxTextExtent,
      " -sPCLPassword=%s",read_info->authenticate);
  (void) CopyMagickString(filename,read_info->filename,MaxTextExtent);
  (void) AcquireUniqueFilename(read_info->filename);
  (void) FormatMagickString(command,MaxTextExtent,
    GetDelegateCommands(delegate_info),
    read_info->antialias != MagickFalse ? 4 : 1,
    read_info->antialias != MagickFalse ? 4 : 1,density,options,
    read_info->filename,input_filename);
  status=SystemCommand(read_info->verbose,command) != 0 ? MagickTrue :
    MagickFalse;
  image=ReadImage(read_info,exception);
  (void) RelinquishUniqueFileResource(read_info->filename);
  (void) RelinquishUniqueFileResource(input_filename);
  read_info=DestroyImageInfo(read_info);
  if (image == (Image *) NULL)
    ThrowReaderException(DelegateError,"PCLDelegateFailed");
  if (LocaleCompare(image->magick,"BMP") == 0)
    {
      Image
        *cmyk_image;

      cmyk_image=ConsolidateCMYKImages(image,&image->exception);
      if (cmyk_image != (Image *) NULL)
        {
          image=DestroyImageList(image);
          image=cmyk_image;
        }
    }
  do
  {
    (void) CopyMagickString(image->filename,filename,MaxTextExtent);
    image->page=page;
    next_image=SyncNextImageInList(image);
    if (next_image != (Image *) NULL)
      image=next_image;
  } while (next_image != (Image *) NULL);
  return(GetFirstImageInList(image));
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   R e g i s t e r P C L I m a g e                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  RegisterPCLImage() adds attributes for the PCL image format to
%  the list of supported formats.  The attributes include the image format
%  tag, a method to read and/or write the format, whether the format
%  supports the saving of more than one frame to the same file or blob,
%  whether the format supports native in-memory I/O, and a brief
%  description of the format.
%
%  The format of the RegisterPCLImage method is:
%
%      unsigned long RegisterPCLImage(void)
%
*/
ModuleExport unsigned long RegisterPCLImage(void)
{
  MagickInfo
    *entry;

  entry=SetMagickInfo("PCL");
  entry->decoder=(DecodeImageHandler *) ReadPCLImage;
  entry->encoder=(EncodeImageHandler *) WritePCLImage;
  entry->magick=(IsImageFormatHandler *) IsPCL;
  entry->adjoin=MagickFalse;
  entry->blob_support=MagickFalse;
  entry->seekable_stream=MagickTrue;
  entry->thread_support=EncoderThreadSupport;
  entry->description=ConstantString("Printer Control Language");
  entry->module=ConstantString("PCL");
  (void) RegisterMagickInfo(entry);
  return(MagickImageCoderSignature);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   U n r e g i s t e r P C L I m a g e                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  UnregisterPCLImage() removes format registrations made by the
%  PCL module from the list of supported formats.
%
%  The format of the UnregisterPCLImage method is:
%
%      UnregisterPCLImage(void)
%
*/
ModuleExport void UnregisterPCLImage(void)
{
  (void) UnregisterMagickInfo("PCL");
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   W r i t e P C L I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  WritePCLImage() writes an image in the Page Control Language encoded
%  image format.
%
%  The format of the WritePCLImage method is:
%
%      MagickBooleanType WritePCLImage(const ImageInfo *image_info,Image *image)
%
%  A description of each parameter follows.
%
%    o image_info: the image info.
%
%    o image:  The image.
%
*/
static MagickBooleanType WritePCLImage(const ImageInfo *image_info,Image *image)
{
  char
    buffer[MaxTextExtent];

  long
    y;

  MagickBooleanType
    status;

  register const PixelPacket
    *p;

  register long
    x;

  register unsigned char
    *q;

  unsigned long
    density;

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
    (void) SetImageColorspace(image,RGBColorspace);
  /*
    Initialize the printer.
  */
  (void) WriteBlobString(image,"\033E");  /* printer reset */
  (void) WriteBlobString(image,"\033&l0E");  /* top margin 0 */
  density=75;
  if (image_info->density != (char *) NULL)
    {
      GeometryInfo
        geometry;

      (void) ParseGeometry(image_info->density,&geometry);
      density=(unsigned long) geometry.rho;
    }
  if (IsMonochromeImage(image,&image->exception) != MagickFalse)
    {
      register unsigned char
        bit,
        byte;

      /*
        Write PCL monochrome image.
      */
      (void) FormatMagickString(buffer,MaxTextExtent,"\033*t%ldR",density);
      (void) WriteBlobString(image,buffer);
      (void) WriteBlobString(image,"\033*r1A");  /* start graphics */
      (void) WriteBlobString(image,"\033*b0M");  /* no compression */
      for (y=0; y < (long) image->rows; y++)
      {
        p=GetVirtualPixels(image,0,y,image->columns,1,&image->exception);
        if (p == (const PixelPacket *) NULL)
          break;
        bit=0;
        byte=0;
        (void) FormatMagickString(buffer,MaxTextExtent,"\033*b%luW",
          (image->columns+7)/8);
        (void) WriteBlobString(image,buffer);
        for (x=0; x < (long) image->columns; x++)
        {
          byte<<=1;
          if (PixelIntensity(p) < ((MagickRealType) QuantumRange/2.0))
            byte|=0x01;
          bit++;
          if (bit == 8)
            {
              (void) WriteBlobByte(image,byte);
              bit=0;
              byte=0;
            }
          p++;
        }
        if (bit != 0)
          (void) WriteBlobByte(image,(unsigned char) (byte << (8-bit)));
        status=SetImageProgress(image,SaveImageTag,y,image->rows);
        if (status == MagickFalse)
          break;
      }
      (void) WriteBlobString(image,"\033*rB");  /* end graphics */
    }
  else
    {
      static char
        color_mode[6] = { 0, 3, 0, 8, 8, 8 };

      unsigned char
        *pixels;

      /*
        Allocate pixel buffers.
      */
      pixels=(unsigned char *) AcquireQuantumMemory(image->columns,
        3UL*sizeof(*pixels));
      if (pixels == (unsigned char *) NULL)
        ThrowWriterException(ResourceLimitError,"MemoryAllocationFailed");
      /*
        Write PCL color image.
      */
      (void) WriteBlobString(image,"\033*r3F");  /* set presentation mode */
      (void) FormatMagickString(buffer,MaxTextExtent,"\033*t%luR",density);  /* set resolution */
      (void) WriteBlobString(image,buffer);
      (void) FormatMagickString(buffer,MaxTextExtent,"\033*r%luT",image->rows);
      (void) WriteBlobString(image,buffer);
      (void) FormatMagickString(buffer,MaxTextExtent,"\033*r%luS",
        image->columns);
      (void) WriteBlobString(image,buffer);
      (void) WriteBlobString(image,"\033*v6W");  /* set color mode */
      (void) WriteBlob(image,6,(unsigned char *) color_mode);
      (void) WriteBlobString(image,"\033*r1A");  /* start raster graphics */
      (void) WriteBlobString(image,"\033*b0Y");  /* set y offset */
      (void) WriteBlobString(image,"\033*b0M");  /* no compression */
      for (y=0; y < (long) image->rows; y++)
      {
        p=GetVirtualPixels(image,0,y,image->columns,1,&image->exception);
        if (p == (const PixelPacket *) NULL)
          break;
        q=pixels;
        for (x=0; x < (long) image->columns; x++)
        {
          *q++=ScaleQuantumToChar(p->red);
          *q++=ScaleQuantumToChar(p->green);
          *q++=ScaleQuantumToChar(p->blue);
          p++;
        }
        (void) FormatMagickString(buffer,MaxTextExtent,"\033*b%luW",
          3*image->columns);  /* send row */
        (void) WriteBlobString(image,buffer);
        (void) WriteBlob(image,3*image->columns,pixels);
        status=SetImageProgress(image,SaveImageTag,y,image->rows);
        if (status == MagickFalse)
          break;
      }
      (void) WriteBlobString(image,"\033*r0C");  /* end graphics */
      pixels=(unsigned char *) RelinquishMagickMemory(pixels);
    }
  (void) WriteBlobString(image,"\033E");
  (void) CloseBlob(image);
  return(MagickTrue);
}
