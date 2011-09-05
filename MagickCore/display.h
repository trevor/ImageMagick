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

  MagickCore methods to interactively display and edit an image.
*/
#ifndef _MAGICKCORE_DISPLAY_H
#define _MAGICKCORE_DISPLAY_H

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

extern MagickExport MagickBooleanType
  DisplayImages(const ImageInfo *,Image *,ExceptionInfo *),
  RemoteDisplayCommand(const ImageInfo *,const char *,const char *,
    ExceptionInfo *);

#if defined(MAGICKCORE_X11_DELEGATE)
#include "MagickCore/xwindow.h"

extern MagickExport Image
  *XDisplayImage(Display *,XResourceInfo *,char **,int,Image **,size_t *,
    ExceptionInfo *);

extern MagickExport MagickBooleanType XDisplayBackgroundImage(Display *,
  XResourceInfo *,Image *,ExceptionInfo *);
#endif

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif
