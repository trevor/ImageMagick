/*
  Copyright 1999-2009 ImageMagick Studio LLC, a non-profit organization
  dedicated to making software imaging solutions freely available.
  
  You may not use this file except in compliance with the License.
  obtain a copy of the License at
  
    http://www.imagemagick.org/script/license.php
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  MagickCore option methods.
*/
#ifndef _MAGICKCORE_OPTION_H
#define _MAGICKCORE_OPTION_H

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

typedef enum
{
  MagickUndefinedOptions = -1,
  MagickAlignOptions = 0,
  MagickAlphaOptions,
  MagickBooleanOptions,
  MagickChannelOptions,
  MagickClassOptions,
  MagickClipPathOptions,
  MagickCoderOptions,
  MagickColorOptions,
  MagickColorspaceOptions,
  MagickCommandOptions,
  MagickComposeOptions,
  MagickCompressOptions,
  MagickConfigureOptions,
  MagickDataTypeOptions,
  MagickDebugOptions,
  MagickDecorateOptions,
  MagickDelegateOptions,
  MagickDisposeOptions,
  MagickDistortOptions,
  MagickDitherOptions,
  MagickEndianOptions,
  MagickEvaluateOptions,
  MagickFillRuleOptions,
  MagickFilterOptions,
  MagickFontOptions,
  MagickFontsOptions,
  MagickFormatOptions,
  MagickFunctionOptions,
  MagickGravityOptions,
  MagickImageListOptions,
  MagickIntentOptions,
  MagickInterlaceOptions,
  MagickInterpolateOptions,
  MagickLayerOptions,
  MagickLineCapOptions,
  MagickLineJoinOptions,
  MagickListOptions,
  MagickLocaleOptions,
  MagickLogEventOptions,
  MagickLogOptions,
  MagickMagicOptions,
  MagickMethodOptions,
  MagickMetricOptions,
  MagickMimeOptions,
  MagickModeOptions,
  MagickModuleOptions,
  MagickNoiseOptions,
  MagickOrientationOptions,
  MagickPreviewOptions,
  MagickPrimitiveOptions,
  MagickQuantumFormatOptions,
  MagickResolutionOptions,
  MagickResourceOptions,
  MagickSparseColorOptions,
  MagickStorageOptions,
  MagickStretchOptions,
  MagickStyleOptions,
  MagickThresholdOptions,
  MagickTypeOptions,
  MagickVirtualPixelOptions
} MagickOption;

typedef struct _OptionInfo
{
  const char
    *mnemonic;

  long
    type;

  MagickBooleanType
    stealth;
} OptionInfo;

extern MagickExport char
  **GetMagickOptions(const MagickOption),
  *GetNextImageOption(const ImageInfo *),
  *RemoveImageOption(ImageInfo *,const char *);

extern MagickExport const char
  *GetImageOption(const ImageInfo *,const char *),
  *MagickOptionToMnemonic(const MagickOption,const long);

extern MagickExport long
  ParseChannelOption(const char *),
  ParseMagickOption(const MagickOption,const MagickBooleanType,const char *);

extern MagickExport MagickBooleanType
  CloneImageOptions(ImageInfo *,const ImageInfo *),
  DefineImageOption(ImageInfo *,const char *),
  DeleteImageOption(ImageInfo *,const char *),
  IsMagickOption(const char *),
  ListMagickOptions(FILE *,const MagickOption,ExceptionInfo *),
  SetImageOption(ImageInfo *,const char *,const char *);

extern MagickExport void
  DestroyImageOptions(ImageInfo *),
  ResetImageOptionIterator(const ImageInfo *);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif
