/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%                      CCCC   AAA    CCCC  H   H  EEEEE                       %
%                     C      A   A  C      H   H  E                           %
%                     C      AAAAA  C      HHHHH  EEE                         %
%                     C      A   A  C      H   H  E                           %
%                      CCCC  A   A   CCCC  H   H  EEEEE                       %
%                                                                             %
%                        V   V  IIIII  EEEEE  W   W                           %
%                        V   V    I    E      W   W                           %
%                        V   V    I    EEE    W W W                           %
%                         V V     I    E      WW WW                           %
%                          V    IIIII  EEEEE  W   W                           %
%                                                                             %
%                                                                             %
%                        MagickCore Cache View Methods                        %
%                                                                             %
%                              Software Design                                %
%                                John Cristy                                  %
%                               February 2000                                 %
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
%
*/

/*
  Include declarations.
*/
#include "magick/studio.h"
#include "magick/cache.h"
#include "magick/cache-private.h"
#include "magick/cache-view.h"
#include "magick/memory_.h"
#include "magick/exception.h"
#include "magick/exception-private.h"
#include "magick/string_.h"

/*
  Typedef declarations.
*/
struct _ViewInfo
{
  Image
    *image;

  VirtualPixelMethod
    virtual_pixel_method;

  unsigned long
    number_threads;

  NexusInfo
    **nexus_info;

  MagickBooleanType
    debug;

  unsigned long
    signature;
};

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   A c q u i r e C a c h e V i e w                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  AcquireCacheView() acquires a view into the pixel cache, using the
%  VirtualPixelMethod that is defined within the given image itself.
%
%  The format of the AcquireCacheView method is:
%
%      ViewInfo *AcquireCacheView(const Image *image)
%
%  A description of each parameter follows:
%
%    o image: the image.
%
*/
MagickExport ViewInfo *AcquireCacheView(const Image *image)
{
  ViewInfo
    *cache_view;

  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  if (image->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",image->filename);
  cache_view=(ViewInfo *) AcquireMagickMemory(sizeof(*cache_view));
  if (cache_view == (ViewInfo *) NULL)
    ThrowFatalException(ResourceLimitFatalError,"MemoryAllocationFailed");
  (void) ResetMagickMemory(cache_view,0,sizeof(*cache_view));
  cache_view->image=ReferenceImage((Image *) image);
  cache_view->number_threads=GetPixelCacheMaximumThreads();
  cache_view->nexus_info=AcquirePixelCacheNexus(cache_view->number_threads);
  cache_view->virtual_pixel_method=GetImageVirtualPixelMethod(image);
  cache_view->debug=IsEventLogging();
  cache_view->signature=MagickSignature;
  if (cache_view->nexus_info == (NexusInfo **) NULL)
    ThrowFatalException(CacheFatalError,"UnableToAcquireCacheView");
  return(cache_view);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C l o n e C a c h e V i e w                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CloneCacheView()  makes an exact copy of the specified cache view.
%
%  The format of the CloneCacheView method is:
%
%      ViewInfo *CloneCacheView(const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport ViewInfo *CloneCacheView(const ViewInfo *cache_view)
{
  ViewInfo
    *clone_view;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  clone_view=(ViewInfo *) AcquireMagickMemory(sizeof(*clone_view));
  if (clone_view == (ViewInfo *) NULL)
    ThrowFatalException(ResourceLimitFatalError,"MemoryAllocationFailed");
  (void) ResetMagickMemory(clone_view,0,sizeof(*clone_view));
  clone_view->image=ReferenceImage(cache_view->image);
  clone_view->number_threads=cache_view->number_threads;
  clone_view->nexus_info=AcquirePixelCacheNexus(cache_view->number_threads);
  clone_view->virtual_pixel_method=cache_view->virtual_pixel_method;
  clone_view->debug=cache_view->debug;
  clone_view->signature=MagickSignature;
  return(clone_view);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   D e s t r o y C a c h e V i e w                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DestroyCacheView() destroys the specified view returned by a previous call
%  to AcquireCacheView().
%
%  The format of the DestroyCacheView method is:
%
%      ViewInfo *DestroyCacheView(ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport ViewInfo *DestroyCacheView(ViewInfo *cache_view)
{
  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  if (cache_view->nexus_info != (NexusInfo **) NULL)
    cache_view->nexus_info=DestroyPixelCacheNexus(cache_view->nexus_info,
      cache_view->number_threads);
  cache_view->image=DestroyImage(cache_view->image);
  cache_view->signature=(~MagickSignature);
  cache_view=(ViewInfo *) RelinquishMagickMemory(cache_view);
  return(cache_view);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w C o l o r s p a c e                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewColorspace() returns the image colorspace associated with the
%  specified view.
%
%  The format of the GetCacheViewColorspace method is:
%
%      ColorspaceType GetCacheViewColorspace(const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport ColorspaceType GetCacheViewColorspace(const ViewInfo *cache_view)
{
  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  return(cache_view->image->colorspace);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w E x c e p t i o n                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewException() returns the image exception associated with the
%  specified view.
%
%  The format of the GetCacheViewException method is:
%
%      ExceptionInfo GetCacheViewException(const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport ExceptionInfo *GetCacheViewException(const ViewInfo *cache_view)
{
  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  return(&cache_view->image->exception);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   G e t C a c h e V i e w E x t e n t                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewExtent() returns the extent of the pixels associated with the
%  last call to QueueCacheViewAuthenticPixels() or
%  GetCacheViewAuthenticPixels().
%
%  The format of the GetCacheViewExtent() method is:
%
%      MagickSizeType GetCacheViewExtent(const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport MagickSizeType GetCacheViewExtent(const ViewInfo *cache_view)
{
  MagickSizeType
    extent;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  assert(cache_view->image->cache != (Cache) NULL);
  extent=GetPixelCacheNexusExtent(cache_view->image->cache,
    cache_view->nexus_info[GetPixelCacheThreadId()]);
  return(extent);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w S t o r a g e C l a s s                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewStorageClass() returns the image storage class  associated with
%  the specified view.
%
%  The format of the GetCacheViewStorageClass method is:
%
%      ClassType GetCacheViewStorageClass(const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport ClassType GetCacheViewStorageClass(const ViewInfo *cache_view)
{
  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  return(cache_view->image->storage_class);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w A u t h e n t i c P i x e l s                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewAuthenticPixels() gets pixels from the in-memory or disk pixel
%  cache as defined by the geometry parameters.   A pointer to the pixels is
%  returned if the pixels are transferred, otherwise a NULL is returned.
%
%  The format of the GetCacheViewAuthenticPixels method is:
%
%      PixelPacket *GetCacheViewAuthenticPixels(ViewInfo *cache_view,
%        const long x,const long y,const unsigned long columns,
%        const unsigned long rows,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o x,y,columns,rows:  These values define the perimeter of a region of
%      pixels.
%
*/
MagickExport PixelPacket *GetCacheViewAuthenticPixels(ViewInfo *cache_view,
  const long x,const long y,const unsigned long columns,
  const unsigned long rows,ExceptionInfo *exception)
{
  Cache
    cache;

  PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  cache=GetImagePixelCache(cache_view->image,MagickTrue,exception);
  if (cache == (Cache) NULL)
    return((PixelPacket *) NULL);
  pixels=GetAuthenticPixelCacheNexus(cache_view->image,x,y,columns,rows,
    cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  return(pixels);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t O n e C a c h e V i e w A u t h e n t i c P i x e l                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetOneCacheViewAuthenticPixel() returns a single pixel at the specified (x,y)
%  location.  The image background color is returned if an error occurs.
%
%  The format of the GetOneCacheViewAuthenticPixel method is:
%
%      MagickBooleaNType GetOneCacheViewAuthenticPixel(
%        const ViewInfo *cache_view,const long x,const long y,
%        Pixelpacket *pixel,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o x,y:  These values define the offset of the pixel.
%
%    o pixel: return a pixel at the specified (x,y) location.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport MagickBooleanType GetOneCacheViewAuthenticPixel(
  const ViewInfo *cache_view,const long x,const long y,PixelPacket *pixel,
  ExceptionInfo *exception)
{
  Cache
    cache;

  PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  cache=GetImagePixelCache(cache_view->image,MagickTrue,exception);
  if (cache == (Cache) NULL)
    return(MagickFalse);
  *pixel=cache_view->image->background_color;
  pixels=GetAuthenticPixelCacheNexus(cache_view->image,x,y,1,1,
    cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  if (pixels == (const PixelPacket *) NULL)
    return(MagickFalse);
  *pixel=(*pixels);
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w A u t h e n t i c I n d e x Q u e u e             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewAuthenticIndexQueue() returns the indexes associated with the
%  last call to SetCacheViewIndexes() or GetCacheViewAuthenticIndexQueue().  The
%  indexes are authentic and can be updated.
%
%  The format of the GetCacheViewAuthenticIndexQueue() method is:
%
%      IndexPacket *GetCacheViewAuthenticIndexQueue(ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport IndexPacket *GetCacheViewAuthenticIndexQueue(ViewInfo *cache_view)
{
  IndexPacket
    *indexes;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  assert(cache_view->image->cache != (Cache) NULL);
  indexes=GetPixelCacheNexusIndexes(cache_view->image->cache,
    cache_view->nexus_info[GetPixelCacheThreadId()]);
  return(indexes);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w A u t h e n t i c P i x e l Q u e u e             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewAuthenticPixelQueue() returns the pixels associated with the
%  last call to QueueCacheViewAuthenticPixels() or
%  GetCacheViewAuthenticPixels().  The pixels are authentic and therefore can be
%  updated.
%
%  The format of the GetCacheViewAuthenticPixelQueue() method is:
%
%      PixelPacket *GetCacheViewAuthenticPixelQueue(ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport PixelPacket *GetCacheViewAuthenticPixelQueue(ViewInfo *cache_view)
{
  PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  assert(cache_view->image->cache != (Cache) NULL);
  pixels=GetPixelCacheNexusPixels(cache_view->image->cache,
    cache_view->nexus_info[GetPixelCacheThreadId()]);
  return(pixels);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w V i r t u a l I n d e x Q u e u e                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewVirtualIndexQueue() returns the indexes associated with the
%  last call to GetCacheViewVirtualIndexQueue().  The indexes are virtual and
%  therefore cannot be updated.
%
%  The format of the GetCacheViewVirtualIndexQueue() method is:
%
%      const IndexPacket *GetCacheViewVirtualIndexQueue(
%        const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport const IndexPacket *GetCacheViewVirtualIndexQueue(
  const ViewInfo *cache_view)
{
  const IndexPacket
    *indexes;

  assert(cache_view != (const ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  assert(cache_view->image->cache != (Cache) NULL);
  indexes=GetVirtualIndexesFromNexus(cache_view->image->cache,
    cache_view->nexus_info[GetPixelCacheThreadId()]);
  return(indexes);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w V i r t u a l P i x e l Q u e u e                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewVirtualPixelQueue() returns the the pixels associated with
%  the last call to GetCacheViewVirtualPixels().  The pixels are virtual
%  and therefore cannot be updated.
%
%  The format of the GetCacheViewVirtualPixelQueue() method is:
%
%      const PixelPacket *GetCacheViewVirtualPixelQueue(
%        const ViewInfo *cache_view)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
*/
MagickExport const PixelPacket *GetCacheViewVirtualPixelQueue(
  const ViewInfo *cache_view)
{
  const PixelPacket
    *pixels;

  assert(cache_view != (const ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  assert(cache_view->image->cache != (Cache) NULL);
  pixels=GetVirtualPixelsNexus(cache_view->image->cache,
    cache_view->nexus_info[GetPixelCacheThreadId()]);
  return(pixels);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t C a c h e V i e w V i r t u a l P i x e l s                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetCacheViewVirtualPixels() gets virtual pixels from the in-memory or
%  disk pixel cache as defined by the geometry parameters.   A pointer to the
%  pixels is returned if the pixels are transferred, otherwise a NULL is
%  returned.
%
%  The format of the GetCacheViewVirtualPixels method is:
%
%      const PixelPacket *GetCacheViewVirtualPixels(
%        const ViewInfo *cache_view,const long x,const long y,
%        const unsigned long columns,const unsigned long rows,
%        ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o x,y,columns,rows:  These values define the perimeter of a region of
%      pixels.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport const PixelPacket *GetCacheViewVirtualPixels(
  const ViewInfo *cache_view,const long x,const long y,
  const unsigned long columns,const unsigned long rows,ExceptionInfo *exception)
{
  const PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  pixels=GetVirtualPixelsFromNexus(cache_view->image,
    cache_view->virtual_pixel_method,x,y,columns,rows,
    cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  return(pixels);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t O n e C a c h e V i e w V i r t u a l P i x e l                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetOneCacheViewVirtualPixel() returns a single pixel at the specified (x,y)
%  location.  The image background color is returned if an error occurs.  If
%  you plan to modify the pixel, use GetOneCacheViewAuthenticPixel() instead.
%
%  The format of the GetOneCacheViewVirtualPixel method is:
%
%      MagickBooleanType GetOneCacheViewVirtualPixel(
%        const ViewInfo *cache_view,const long x,const long y,
%        PixelPacket *pixel,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o x,y:  These values define the offset of the pixel.
%
%    o pixel: return a pixel at the specified (x,y) location.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport MagickBooleanType GetOneCacheViewVirtualPixel(
  const ViewInfo *cache_view,const long x,const long y,PixelPacket *pixel,
  ExceptionInfo *exception)
{
  const PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  *pixel=cache_view->image->background_color;
  pixels=GetVirtualPixelsFromNexus(cache_view->image,
    cache_view->virtual_pixel_method,x,y,1,1,
    cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  if (pixels == (const PixelPacket *) NULL)
    return(MagickFalse);
  *pixel=(*pixels);
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t O n e C a c h e V i e w V i r t u a l P i x e l                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GetOneCacheViewVirtualMethodPixel() returns a single virtual pixel at
%  the specified (x,y) location.  The image background color is returned if an
%  error occurs.  If you plan to modify the pixel, use
%  GetOneCacheViewAuthenticPixel() instead.
%
%  The format of the GetOneCacheViewVirtualPixel method is:
%
%      MagickBooleanType GetOneCacheViewVirtualMethodPixel(
%        const ViewInfo *cache_view,
%        const VirtualPixelMethod virtual_pixel_method,const long x,
%        const long y,PixelPacket *pixel,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o virtual_pixel_method: the virtual pixel method.
%
%    o x,y:  These values define the offset of the pixel.
%
%    o pixel: return a pixel at the specified (x,y) location.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport MagickBooleanType GetOneCacheViewVirtualMethodPixel(
  const ViewInfo *cache_view,const VirtualPixelMethod virtual_pixel_method,
  const long x,const long y,PixelPacket *pixel,ExceptionInfo *exception)
{
  const PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  *pixel=cache_view->image->background_color;
  pixels=GetVirtualPixelsFromNexus(cache_view->image,virtual_pixel_method,x,y,1,
    1,cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  if (pixels == (const PixelPacket *) NULL)
    return(MagickFalse);
  *pixel=(*pixels);
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   Q u e u e C a c h e V i e w A u t h e n t i c P i x e l s                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  QueueCacheViewAuthenticPixels() queues authentic pixels from the in-memory or
%  disk pixel cache as defined by the geometry parameters.   A pointer to the
%  pixels is returned if the pixels are transferred, otherwise a NULL is
%  returned.
%
%  The format of the QueueCacheViewAuthenticPixels method is:
%
%      PixelPacket *QueueCacheViewAuthenticPixels(ViewInfo *cache_view,
%        const long x,const long y,const unsigned long columns,
%        const unsigned long rows,ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o x,y,columns,rows:  These values define the perimeter of a region of
%      pixels.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport PixelPacket *QueueCacheViewAuthenticPixels(ViewInfo *cache_view,
  const long x,const long y,const unsigned long columns,
  const unsigned long rows,ExceptionInfo *exception)
{
  Cache
    cache;

  PixelPacket
    *pixels;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  cache=GetImagePixelCache(cache_view->image,MagickFalse,exception);
  if (cache == (Cache) NULL)
    return((PixelPacket *) NULL);
  pixels=QueueAuthenticNexus(cache_view->image,x,y,columns,rows,
    cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  return(pixels);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S e t C a c h e V i e w S t o r a g e C l a s s                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SetCacheViewStorageClass() sets the image storage class associated with
%  the specified view.
%
%  The format of the SetCacheViewStorageClass method is:
%
%      MagickBooleanType SetCacheViewStorageClass(ViewInfo *cache_view,
%        const ClassType storage_class)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o storage_class: the image storage class: PseudoClass or DirectClass.
%
*/
MagickExport MagickBooleanType SetCacheViewStorageClass(ViewInfo *cache_view,
  const ClassType storage_class)
{
  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  return(SetImageStorageClass(cache_view->image,storage_class));
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S e t C a c h e V i e w V i r t u a l P i x e l M e t h o d               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SetCacheViewVirtualPixelMethod() sets the virtual pixel method associated
%  with the specified cache view.
%
%  The format of the SetCacheViewVirtualPixelMethod method is:
%
%      MagickBooleanType SetCacheViewVirtualPixelMethod(ViewInfo *cache_view,
%        const VirtualPixelMethod virtual_pixel_method)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o virtual_pixel_method: the virtual pixel method.
%
*/
MagickExport MagickBooleanType SetCacheViewVirtualPixelMethod(
  ViewInfo *cache_view,const VirtualPixelMethod virtual_pixel_method)
{
  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  cache_view->virtual_pixel_method=virtual_pixel_method;
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S y n c C a c h e V i e w A u t h e n t i c P i x e l s                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SyncCacheViewAuthenticPixels() saves the cache view pixels to the in-memory
%  or disk cache.  It returns MagickTrue if the pixel region is flushed,
%  otherwise MagickFalse.
%
%  The format of the SyncCacheViewAuthenticPixels method is:
%
%      MagickBooleanType SyncCacheViewAuthenticPixels(ViewInfo *cache_view,
%        ExceptionInfo *exception)
%
%  A description of each parameter follows:
%
%    o cache_view: the cache view.
%
%    o exception: return any errors or warnings in this structure.
%
*/
MagickExport MagickBooleanType SyncCacheViewAuthenticPixels(
  ViewInfo *cache_view,ExceptionInfo *exception)
{
  MagickBooleanType
    status;

  assert(cache_view != (ViewInfo *) NULL);
  assert(cache_view->signature == MagickSignature);
  if (cache_view->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",
      cache_view->image->filename);
  status=SyncAuthenticPixelCacheNexus(cache_view->image,
    cache_view->nexus_info[GetPixelCacheThreadId()],exception);
  return(status);
}
