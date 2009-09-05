/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%        CCCC   OOO   M   M  PPPP    OOO   SSSSS  IIIII  TTTTT  EEEEE         %
%       C      O   O  MM MM  P   P  O   O  SS       I      T    E             %
%       C      O   O  M M M  PPPP   O   O   SSS     I      T    EEE           %
%       C      O   O  M   M  P      O   O     SS    I      T    E             %
%        CCCC   OOO   M   M  P       OOO   SSSSS  IIIII    T    EEEEE         %
%                                                                             %
%                                                                             %
%                     MagickCore Image Composite Methods                      %
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
%
*/

/*
  Include declarations.
*/
#include "magick/studio.h"
#include "magick/artifact.h"
#include "magick/cache-view.h"
#include "magick/client.h"
#include "magick/color.h"
#include "magick/color-private.h"
#include "magick/colorspace.h"
#include "magick/colorspace-private.h"
#include "magick/composite.h"
#include "magick/composite-private.h"
#include "magick/constitute.h"
#include "magick/draw.h"
#include "magick/fx.h"
#include "magick/gem.h"
#include "magick/geometry.h"
#include "magick/image.h"
#include "magick/image-private.h"
#include "magick/list.h"
#include "magick/log.h"
#include "magick/monitor.h"
#include "magick/monitor-private.h"
#include "magick/memory_.h"
#include "magick/option.h"
#include "magick/pixel-private.h"
#include "magick/property.h"
#include "magick/quantum.h"
#include "magick/resample.h"
#include "magick/resource_.h"
#include "magick/string_.h"
#include "magick/utility.h"
#include "magick/version.h"

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C o m p o s i t e I m a g e C h a n n e l                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CompositeImageChannel() returns the second image composited onto the first
%  at the specified offset, using the specified composite method.
%
%  The format of the CompositeImageChannel method is:
%
%      MagickBooleanType CompositeImage(Image *image,
%        const CompositeOperator compose,Image *composite_image,
%        const long x_offset,const long y_offset)
%      MagickBooleanType CompositeImageChannel(Image *image,
%        const ChannelType channel,const CompositeOperator compose,
%        Image *composite_image,const long x_offset,const long y_offset)
%
%  A description of each parameter follows:
%
%    o image: the destination image, modified by he composition
%
%    o channel: the channel.
%
%    o compose: This operator affects how the composite is applied to
%      the image.  The operators and how they are utilized are listed here
%      http://www.w3.org/TR/SVG12/#compositing.
%
%    o composite_image: the composite (source) image.
%
%    o x_offset: the column offset of the composited image.
%
%    o y_offset: the row offset of the composited image.
%
%  Extra Controls from Image meta-data in 'composite_image' (artifacts)
%
%    o "compose:args"
%        A string containing extra numerical arguments for specific compose
%        methods, generally expressed as a 'geometry' or a comma separated list
%        of numbers.
%
%        Compose methods needing such arguments include "BlendCompositeOp" and
%        "DisplaceCompositeOp".
%
%    o "compose:outside-overlay"
%        Modify how the composition is to effect areas not directly covered
%        by the 'composite_image' at the offset given.  Normally this is
%        dependant on the 'compose' method, especially Duff-Porter methods.
%
%        If set to "false" then disable all normal handling of pixels not
%        covered by the composite_image.  Typically used for repeated tiling
%        of the composite_image by the calling API.
%
%        FUTURE: If set to "virtual" then the whole  destination canvas is
%        composed with the 'composite_image' surrounded by its virtual pixels.
%
%        Previous to IM v6.5.3-3  this was called "modify-outside-overlay"
%
*/

static inline double MagickMin(const double x,const double y)
{
  if (x < y)
    return(x);
  return(y);
}
static inline double MagickMax(const double x,const double y)
{
  if (x > y)
    return(x);
  return(y);
}

static inline MagickRealType Add(const MagickRealType p,const MagickRealType q)
{
  MagickRealType
    pixel;

  pixel=p+q;
  if (pixel > QuantumRange)
    pixel-=(QuantumRange+1.0);
  return(pixel);
}

static inline void CompositeAdd(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  composite->red=Add(p->red,q->red);
  composite->green=Add(p->green,q->green);
  composite->blue=Add(p->blue,q->blue);
  composite->opacity=Add(alpha,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=Add(p->index,q->index);
}

static inline MagickRealType Atop(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=((1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)+
    (1.0-QuantumScale*beta)*q*QuantumScale*alpha);
  return(pixel);
}

static inline void CompositeAtop(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=(1.0-QuantumScale*beta);
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Atop(p->red,alpha,q->red,beta);
  composite->green=gamma*Atop(p->green,alpha,q->green,beta);
  composite->blue=gamma*Atop(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Atop(p->index,alpha,q->index,beta);
}

static inline void CompositeBumpmap(const MagickPixelPacket *p,
  const MagickRealType magick_unused(alpha),const MagickPixelPacket *q,
  const MagickRealType magick_unused(beta),MagickPixelPacket *composite)
{
  MagickRealType
    intensity;

  intensity=MagickPixelIntensity(p);
  composite->red=QuantumScale*intensity*q->red;
  composite->green=QuantumScale*intensity*q->green;
  composite->blue=QuantumScale*intensity*q->blue;
  composite->opacity=(MagickRealType) QuantumScale*intensity*p->opacity;
  if (q->colorspace == CMYKColorspace)
    composite->index=QuantumScale*intensity*q->index;
}

static inline void CompositeClear(const MagickPixelPacket *q,
  MagickPixelPacket *composite)
{
  composite->red=0.0;
  composite->green=0.0;
  composite->blue=0.0;
  composite->opacity=(MagickRealType) TransparentOpacity;
  if (q->colorspace == CMYKColorspace)
    composite->index=0.0;
}

static MagickRealType ColorBurn(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    delta,
    pixel;

  delta=QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha);
  if (delta <= ((1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta)))
    {
      pixel=QuantumScale*(1.0-QuantumScale*alpha)*p*
        (1.0-(1.0-QuantumScale*beta))+QuantumScale*(1.0-QuantumScale*beta)*q*
        (1.0-(1.0-QuantumScale*alpha));
      return(pixel);
    }
  pixel=QuantumScale*((1.0-QuantumScale*alpha)*(QuantumScale*(1.0-QuantumScale*
    alpha)*p*(1.0-QuantumScale*beta)+QuantumScale*(1.0-QuantumScale*beta)*q*
    (1.0-QuantumScale*alpha)-(1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta))/
    QuantumScale*(1.0-QuantumScale*alpha)*p+QuantumScale*(1.0-
    QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha)));
  return(pixel);
}

static inline void CompositeColorBurn(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*ColorBurn(p->red,alpha,q->red,beta);
  composite->green=gamma*ColorBurn(p->green,alpha,q->green,beta);
  composite->blue=gamma*ColorBurn(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*ColorBurn(p->index,alpha,q->index,beta);
}

static MagickRealType ColorDodge(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    delta,
    pixel;

  delta=QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha);
  if (delta >= ((1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta)))
    {
      pixel=(1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta)+QuantumScale*
        (1.0-QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+QuantumScale*
        (1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha));
      return((MagickRealType) QuantumRange*pixel);
    }
  pixel=QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha)/
    (1.0-QuantumScale*(1.0-QuantumScale*alpha)*p/(1.0-QuantumScale*alpha))+
    QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha));
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeColorDodge(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*ColorDodge(p->red,alpha,q->red,beta);
  composite->green=gamma*ColorDodge(p->green,alpha,q->green,beta);
  composite->blue=gamma*ColorDodge(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*ColorDodge(p->index,alpha,q->index,beta);
}

static inline MagickRealType Darken(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  if (((1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)) <
      ((1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha)))
    {
      pixel=(1.0-QuantumScale*alpha)*p+(1.0-QuantumScale*beta)*q*QuantumScale*
        alpha;
      return(pixel);
    }
  pixel=(1.0-QuantumScale*beta)*q+(1.0-QuantumScale*alpha)*p*QuantumScale*beta;
  return(pixel);
}

static inline void CompositeDarken(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Darken(p->red,alpha,q->red,beta);
  composite->green=gamma*Darken(p->green,alpha,q->green,beta);
  composite->blue=gamma*Darken(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Darken(p->index,alpha,q->index,beta);
}

static inline MagickRealType Difference(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=QuantumScale*(1.0-QuantumScale*alpha)*p+QuantumScale*
    (1.0-QuantumScale*beta)*q-2.0*MagickMin(QuantumScale*(1.0-QuantumScale*
    alpha)*p*(1.0-QuantumScale*beta),QuantumScale*(1.0-QuantumScale*beta)*q*
    (1.0-QuantumScale*alpha));
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeDifference(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Difference(p->red,alpha,q->red,beta);
  composite->green=gamma*Difference(p->green,alpha,q->green,beta);
  composite->blue=gamma*Difference(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Difference(p->index,alpha,q->index,beta);
}

static inline MagickRealType Divide(const MagickRealType p,
  const MagickRealType q)
{
  MagickRealType
    pixel;

  if (q == 0.0)
    return(0.0);
  pixel=p/q;
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeDivide(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  composite->red=Divide(p->red,q->red);
  composite->green=Divide(p->green,q->green);
  composite->blue=Divide(p->blue,q->blue);
  composite->opacity=Divide(alpha,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=Divide(p->index,q->index);
}

static inline MagickRealType Exclusion(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=(QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha)-2.0*
    QuantumScale*(1.0-QuantumScale*alpha)*p*QuantumScale*
    (1.0-QuantumScale*beta)*q)+QuantumScale*(1.0-QuantumScale*alpha)*p*
    (1.0-(1.0-QuantumScale*beta))+QuantumScale*(1.0-QuantumScale*beta)*q*
    (1.0 -(1.0-QuantumScale*alpha));
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeExclusion(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Exclusion(p->red,alpha,q->red,beta);
  composite->green=gamma*Exclusion(p->green,alpha,q->green,beta);
  composite->blue=gamma*Exclusion(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Exclusion(p->index,alpha,q->index,beta);
}

static MagickRealType HardLight(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  if ((2.0*QuantumScale*(1.0-QuantumScale*alpha)*p) < (1.0-QuantumScale*alpha))
    {
      pixel=2.0*QuantumScale*(1.0-QuantumScale*alpha)*p*QuantumScale*
        (1.0-QuantumScale*beta)*q+QuantumScale*(1.0-QuantumScale*alpha)*p*
        (1.0-(1.0-QuantumScale*beta))+QuantumScale*(1.0-QuantumScale*beta)*q*
        (1.0-(1.0-QuantumScale*alpha));
      return((MagickRealType) QuantumRange*pixel);
    }
  pixel=(1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta)-2.0*
    ((1.0-QuantumScale*beta)-QuantumScale*(1.0-QuantumScale*beta)*q)*
    ((1.0-QuantumScale*alpha)-QuantumScale*(1.0-QuantumScale*alpha)*p)+
    QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha));
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeHardLight(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*HardLight(p->red,alpha,q->red,beta);
  composite->green=gamma*HardLight(p->green,alpha,q->green,beta);
  composite->blue=gamma*HardLight(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*HardLight(p->index,alpha,q->index,beta);
}

static void CompositeHSB(const MagickRealType red,const MagickRealType green,
  const MagickRealType blue,double *hue,double *saturation,double *brightness)
{
  MagickRealType
    delta,
    max,
    min;

  /*
    Convert RGB to HSB colorspace.
  */
  assert(hue != (double *) NULL);
  assert(saturation != (double *) NULL);
  assert(brightness != (double *) NULL);
  max=(red > green ? red : green);
  if (blue > max)
    max=blue;
  min=(red < green ? red : green);
  if (blue < min)
    min=blue;
  *hue=0.0;
  *saturation=0.0;
  *brightness=(double) (QuantumScale*max);
  if (max == 0.0)
    return;
  *saturation=(double) (1.0-min/max);
  delta=max-min;
  if (delta == 0.0)
    return;
  if (red == max)
    *hue=(double) ((green-blue)/delta);
  else
    if (green == max)
      *hue=(double) (2.0+(blue-red)/delta);
    else
      if (blue == max)
        *hue=(double) (4.0+(red-green)/delta);
  *hue/=6.0;
  if (*hue < 0.0)
    *hue+=1.0;
}

static inline MagickRealType In(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType magick_unused(q),
  const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=(1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta);
  return(pixel);
}

static inline void CompositeIn(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=(1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta);
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*In(p->red,alpha,q->red,beta);
  composite->green=gamma*In(p->green,alpha,q->green,beta);
  composite->blue=gamma*In(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*In(p->index,alpha,q->index,beta);
}

static inline MagickRealType Lighten(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  if (((1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)) >
      ((1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha)))
    {
      pixel=(1.0-QuantumScale*alpha)*p+(1.0-QuantumScale*beta)*q*QuantumScale*
        alpha;
      return(pixel);
    }
  pixel=(1.0-QuantumScale*beta)*q+(1.0-QuantumScale*alpha)*p*QuantumScale*beta;
  return(pixel);
}

static inline void CompositeLighten(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Lighten(p->red,alpha,q->red,beta);
  composite->green=gamma*Lighten(p->green,alpha,q->green,beta);
  composite->blue=gamma*Lighten(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Lighten(p->index,alpha,q->index,beta);
}

static inline MagickRealType LinearLight(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  return((1.0-QuantumScale*beta)*q+2.0*(1.0-QuantumScale*alpha)*p-1.0);
}

static inline void CompositeLinearLight(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*beta)+2.0*(1.0-QuantumScale*alpha)-1.0);
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*LinearLight(p->red,alpha,q->red,beta);
  composite->green=gamma*LinearLight(p->green,alpha,q->green,beta);
  composite->blue=gamma*LinearLight(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*LinearLight(p->index,alpha,q->index,beta);
}

static inline MagickRealType Minus(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  return((1.0-QuantumScale*alpha)*p-(1.0-QuantumScale*beta)*q);
}

static inline void CompositeMinus(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)-(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Minus(p->red,alpha,q->red,beta);
  composite->green=gamma*Minus(p->green,alpha,q->green,beta);
  composite->blue=gamma*Minus(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Minus(p->index,alpha,q->index,beta);
}

static inline MagickRealType Multiply(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-QuantumScale*beta)*q+
    (1.0-QuantumScale*alpha)*p*QuantumScale*beta+(1.0-QuantumScale*beta)*q*
    QuantumScale*alpha;
  return(pixel);
}

static inline void CompositeMultiply(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Multiply(p->red,alpha,q->red,beta);
  composite->green=gamma*Multiply(p->green,alpha,q->green,beta);
  composite->blue=gamma*Multiply(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Multiply(p->index,alpha,q->index,beta);
}

static inline MagickRealType Out(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType magick_unused(q),
  const MagickRealType beta)
{
  return((1.0-QuantumScale*alpha)*p*QuantumScale*beta);
}

static inline void CompositeOut(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=(1.0-QuantumScale*alpha)*QuantumScale*beta;
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Out(p->red,alpha,q->red,beta);
  composite->green=gamma*Out(p->green,alpha,q->green,beta);
  composite->blue=gamma*Out(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Out(p->index,alpha,q->index,beta);
}

static inline void CompositeOver(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=1.0-QuantumScale*QuantumScale*alpha*beta;
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*MagickOver_(p->red,alpha,q->red,beta);
  composite->green=gamma*MagickOver_(p->green,alpha,q->green,beta);
  composite->blue=gamma*MagickOver_(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*MagickOver_(p->index,alpha,q->index,beta);
}

static MagickRealType Overlay(const MagickRealType p,const MagickRealType alpha,
  const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  if ((2.0*QuantumScale*(1.0-QuantumScale*beta)*q) < (1.0-QuantumScale*beta))
    {
      pixel=2.0*QuantumScale*(1.0-QuantumScale*alpha)*p*QuantumScale*
        (1.0-QuantumScale*beta)*q+QuantumScale*(1.0-QuantumScale*alpha)*p*
        (1.0-(1.0-QuantumScale*beta))+QuantumScale*(1.0-QuantumScale*beta)*q*
        (1.0-(1.0-QuantumScale*alpha));
      return((MagickRealType) QuantumRange*pixel);
    }
  pixel=(1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta)-2.0*
    ((1.0-QuantumScale*beta)-QuantumScale*(1.0-QuantumScale*beta)*q)*
    ((1.0-QuantumScale*alpha)-QuantumScale*(1.0-QuantumScale*alpha)*p)+
    QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha));
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeOverlay(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Overlay(p->red,alpha,q->red,beta);
  composite->green=gamma*Overlay(p->green,alpha,q->green,beta);
  composite->blue=gamma*Overlay(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Overlay(p->index,alpha,q->index,beta);
}

static inline MagickRealType Plus(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  return((1.0-QuantumScale*alpha)*p+(1.0-QuantumScale*beta)*q);
}

static inline void CompositePlus(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Plus(p->red,alpha,q->red,beta);
  composite->green=gamma*Plus(p->green,alpha,q->green,beta);
  composite->blue=gamma*Plus(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Plus(p->index,alpha,q->index,beta);
}

static inline MagickRealType Screen(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=QuantumScale*(1.0-QuantumScale*alpha)*p+QuantumScale*
    (1.0-QuantumScale*beta)*q-QuantumScale*(1.0-QuantumScale*alpha)*p*
    QuantumScale*(1.0-QuantumScale*beta)*q;
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeScreen(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Screen(p->red,alpha,q->red,beta);
  composite->green=gamma*Screen(p->green,alpha,q->green,beta);
  composite->blue=gamma*Screen(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Screen(p->index,alpha,q->index,beta);
}

static MagickRealType SoftLight(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  if ((2.0*QuantumScale*(1.0-QuantumScale*alpha)*p) < (1.0-QuantumScale*alpha))
    {
      pixel=QuantumScale*(1.0-QuantumScale*beta)*q*((1.0-QuantumScale*alpha)-
        (1.0-QuantumScale*(1.0-QuantumScale*beta)*q/(1.0-QuantumScale*beta))*
        (2.0*QuantumScale*(1.0-QuantumScale*alpha)*p-(1.0-QuantumScale*alpha)))+
        QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+
        QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha));
      return((MagickRealType) QuantumRange*pixel);
    }
  if ((8.0*QuantumScale*(1.0-QuantumScale*beta)*q) <= (1.0-QuantumScale*beta))
    {
      pixel=QuantumScale*(1.0-QuantumScale*beta)*q*((1.0-QuantumScale*alpha)-
        (1.0-QuantumScale*(1.0-QuantumScale*beta)*q/(1.0-QuantumScale*beta))*
        (2.0*QuantumScale*(1.0-QuantumScale*alpha)*p-(1.0-QuantumScale*alpha))*
        (3.0-8.0*QuantumScale*(1.0-QuantumScale*beta)*q/
        (1.0-QuantumScale*beta)))+QuantumScale*(1.0-QuantumScale*alpha)*p*
        (1.0-(1.0-QuantumScale*beta))+QuantumScale*(1.0-QuantumScale*beta)*q*
        (1.0-(1.0-QuantumScale*alpha));
      return((MagickRealType) QuantumRange*pixel);
    }
  pixel=(QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-QuantumScale*alpha)+
    (pow(QuantumScale*(1.0-QuantumScale*beta)*q/(1.0-QuantumScale*beta),0.5)*
    (1.0-QuantumScale*beta)-QuantumScale*(1.0-QuantumScale*beta)*q)*
    (2.0*QuantumScale*(1.0-QuantumScale*alpha)*p-(1.0-QuantumScale*alpha)))+
    QuantumScale*(1.0-QuantumScale*alpha)*p*(1.0-(1.0-QuantumScale*beta))+
    QuantumScale*(1.0-QuantumScale*beta)*q*(1.0-(1.0-QuantumScale*alpha));
  return((MagickRealType) QuantumRange*pixel);
}

static inline void CompositeSoftLight(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=RoundToUnity((1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    (1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta));
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*SoftLight(p->red,alpha,q->red,beta);
  composite->green=gamma*SoftLight(p->green,alpha,q->green,beta);
  composite->blue=gamma*SoftLight(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*SoftLight(p->index,alpha,q->index,beta);
}

static inline MagickRealType Subtract(const MagickRealType p,
  const MagickRealType magick_unused(alpha),const MagickRealType q,
  const MagickRealType magick_unused(beta))
{
  MagickRealType
    pixel;

  pixel=p-q;
  if (pixel < 0.0)
    pixel+=(QuantumRange+1.0);
  return(pixel);
}

static inline void CompositeSubtract(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  composite->red=Subtract(p->red,alpha,q->red,beta);
  composite->green=Subtract(p->green,alpha,q->green,beta);
  composite->blue=Subtract(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=Subtract(p->index,alpha,q->index,beta);
}

static inline MagickRealType Threshold(const MagickRealType p,
  const MagickRealType magick_unused(alpha),const MagickRealType q,
  const MagickRealType magick_unused(beta),const MagickRealType threshold,
  const MagickRealType amount)
{
  MagickRealType
    pixel;

  pixel=p-q;
  if ((MagickRealType) fabs((double) (2.0*pixel)) < threshold)
    {
      pixel=q;
      return(pixel);
    }
  pixel=q+(pixel*amount);
  return(pixel);
}

static inline void CompositeThreshold(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,const MagickRealType threshold,
  const MagickRealType amount,MagickPixelPacket *composite)
{
  composite->red=Threshold(p->red,alpha,q->red,beta,threshold,amount);
  composite->green=Threshold(p->green,alpha,q->green,beta,threshold,amount);
  composite->blue=Threshold(p->blue,alpha,q->blue,beta,threshold,amount);
  composite->opacity=(MagickRealType) QuantumRange-
    Threshold(p->opacity,alpha,q->opacity,beta,threshold,amount);
  if (q->colorspace == CMYKColorspace)
    composite->index=Threshold(p->index,alpha,q->index,beta,threshold,amount);
}

static inline MagickRealType Xor(const MagickRealType p,
  const MagickRealType alpha,const MagickRealType q,const MagickRealType beta)
{
  MagickRealType
    pixel;

  pixel=(1.0-QuantumScale*alpha)*p*QuantumScale*beta+(1.0-QuantumScale*beta)*q*
    QuantumScale*alpha;
  return(pixel);
}

static inline void CompositeXor(const MagickPixelPacket *p,
  const MagickRealType alpha,const MagickPixelPacket *q,
  const MagickRealType beta,MagickPixelPacket *composite)
{
  MagickRealType
    gamma;

  gamma=(1.0-QuantumScale*alpha)+(1.0-QuantumScale*beta)-
    2*(1.0-QuantumScale*alpha)*(1.0-QuantumScale*beta);
  composite->opacity=(MagickRealType) QuantumRange*(1.0-gamma);
  gamma=1.0/(fabs(gamma) <= MagickEpsilon ? 1.0 : gamma);
  composite->red=gamma*Xor(p->red,alpha,q->red,beta);
  composite->green=gamma*Xor(p->green,alpha,q->green,beta);
  composite->blue=gamma*Xor(p->blue,alpha,q->blue,beta);
  if (q->colorspace == CMYKColorspace)
    composite->index=gamma*Xor(p->index,alpha,q->index,beta);
}

static void HSBComposite(const double hue,const double saturation,
  const double brightness,MagickRealType *red,MagickRealType *green,
  MagickRealType *blue)
{
  MagickRealType
    f,
    h,
    p,
    q,
    t;

  /*
    Convert HSB to RGB colorspace.
  */
  assert(red != (MagickRealType *) NULL);
  assert(green != (MagickRealType *) NULL);
  assert(blue != (MagickRealType *) NULL);
  if (saturation == 0.0)
    {
      *red=(MagickRealType) QuantumRange*brightness;
      *green=(*red);
      *blue=(*red);
      return;
    }
  h=6.0*(hue-floor(hue));
  f=h-floor((double) h);
  p=brightness*(1.0-saturation);
  q=brightness*(1.0-saturation*f);
  t=brightness*(1.0-saturation*(1.0-f));
  switch ((int) h)
  {
    case 0:
    default:
    {
      *red=(MagickRealType) QuantumRange*brightness;
      *green=(MagickRealType) QuantumRange*t;
      *blue=(MagickRealType) QuantumRange*p;
      break;
    }
    case 1:
    {
      *red=(MagickRealType) QuantumRange*q;
      *green=(MagickRealType) QuantumRange*brightness;
      *blue=(MagickRealType) QuantumRange*p;
      break;
    }
    case 2:
    {
      *red=(MagickRealType) QuantumRange*p;
      *green=(MagickRealType) QuantumRange*brightness;
      *blue=(MagickRealType) QuantumRange*t;
      break;
    }
    case 3:
    {
      *red=(MagickRealType) QuantumRange*p;
      *green=(MagickRealType) QuantumRange*q;
      *blue=(MagickRealType) QuantumRange*brightness;
      break;
    }
    case 4:
    {
      *red=(MagickRealType) QuantumRange*t;
      *green=(MagickRealType) QuantumRange*p;
      *blue=(MagickRealType) QuantumRange*brightness;
      break;
    }
    case 5:
    {
      *red=(MagickRealType) QuantumRange*brightness;
      *green=(MagickRealType) QuantumRange*p;
      *blue=(MagickRealType) QuantumRange*q;
      break;
    }
  }
}

MagickExport MagickBooleanType CompositeImage(Image *image,
  const CompositeOperator compose,const Image *composite_image,
  const long x_offset,const long y_offset)
{
  MagickBooleanType
    status;

  status=CompositeImageChannel(image,DefaultChannels,compose,composite_image,
    x_offset,y_offset);
  return(status);
}

MagickExport MagickBooleanType CompositeImageChannel(Image *image,
  const ChannelType magick_unused(channel),const CompositeOperator compose,
  const Image *composite_image,const long x_offset,const long y_offset)
{
#define CompositeImageTag  "Composite/Image"

  const char
    *value;

  double
    sans;

  ExceptionInfo
    *exception;

  GeometryInfo
    geometry_info;

  Image
    *dest_image;

  long
    progress,
    y;

  MagickBooleanType
    modify_outside_overlay,
    status;

  MagickPixelPacket
    zero;

  MagickRealType
    amount,
    destination_dissolve,
    midpoint,
    percent_brightness,
    percent_saturation,
    source_dissolve,
    threshold;

  MagickStatusType
    flags;

  ViewInfo
    *composite_view,
    *image_view;

  /*
    Prepare composite image.
  */
  assert(image != (Image *) NULL);
  assert(image->signature == MagickSignature);
  if (image->debug != MagickFalse)
    (void) LogMagickEvent(TraceEvent,GetMagickModule(),"%s",image->filename);
  assert(composite_image != (Image *) NULL);
  assert(composite_image->signature == MagickSignature);
  if (SetImageStorageClass(image,DirectClass) == MagickFalse)
    return(MagickFalse);
  GetMagickPixelPacket(image,&zero);
  dest_image=(Image *) NULL;
  amount=0.5;
  destination_dissolve=1.0;
  modify_outside_overlay=MagickFalse;
  percent_brightness=100.0;
  percent_saturation=100.0;
  source_dissolve=1.0;
  threshold=0.05f;
  switch (compose)
  {
    case ClearCompositeOp:
    case SrcCompositeOp:
    case SrcInCompositeOp:
    case InCompositeOp:
    case SrcOutCompositeOp:
    case OutCompositeOp:
    case DstInCompositeOp:
    case DstAtopCompositeOp:
    {
      /*
        Modify destination outside the overlaid region.
      */
      modify_outside_overlay=MagickTrue;
      break;
    }
    case CopyOpacityCompositeOp:
    case ChangeMaskCompositeOp:
    {
      /*
        Modify destination outside the overlaid region and require an alpha
        channel to exist, to add transparency.
      */
      if (image->matte == MagickFalse)
        (void) SetImageAlphaChannel(image,OpaqueAlphaChannel);
      modify_outside_overlay=MagickTrue;
      break;
    }
    case BlurCompositeOp:
    {
      /* Blur Image based on overlay gradient map
         X = red_channel    Y = green_channel
         compose:args =  x_scale[,y_scale[,angle]]
       */
      MagickPixelPacket
        pixel;

      MagickRealType
        blur_xx,blur_xy,blur_yx,blur_yy;

      register IndexPacket
        *__restrict dest_indexes;

      register PixelPacket
        *__restrict r;

      ResampleFilter
        *resample_filter;

      ViewInfo
        *composite_view,
        *dest_view;

      /*
        Allocate the destination image.
      */
      dest_image=CloneImage(composite_image,composite_image->columns,
        composite_image->rows,MagickTrue,&image->exception);
      if (dest_image == (Image *) NULL)
        return(MagickFalse);
      /*
        Determine the horizontal and vertical maximim blur
      */
      SetGeometryInfo(&geometry_info);
      flags=NoValue;
      value=GetImageArtifact(composite_image,"compose:args");
      if (value != (char *) NULL)
        flags=ParseGeometry(value,&geometry_info);
      if ((flags & WidthValue) == 0 )
        {
          dest_image=DestroyImage(dest_image);
          return(MagickFalse);
        }
      blur_xx = geometry_info.rho;
      blur_yy = geometry_info.sigma;
      blur_xy = blur_yx = 0.0;
      if ((flags & HeightValue) == 0)
        blur_yy=blur_xx;
      if ((flags & XValue) != 0)
        {
/* The rotate part is not working properly...
 * I think I may have found a mistake in the ResamplePixelColor()
 */
          MagickRealType
             x = blur_xx,
             y = blur_yy,
             angle = DegreesToRadians(geometry_info.xi);
          blur_xx = x*cos(angle);
          blur_xy = x*sin(angle);
          blur_yx = -y*sin(angle);
          blur_yy = y*cos(angle);
        }

      /*
        Blur Image by Resampling
      */
      pixel=zero;
      exception=(&image->exception);
      resample_filter=AcquireResampleFilter(image,&image->exception);
      SetResampleFilter(resample_filter,GaussianFilter,1.0);
      dest_view=AcquireCacheView(dest_image);
      composite_view=AcquireCacheView(composite_image);
      for (y=0; y < (long) composite_image->rows; y++)
      {
        MagickBooleanType
          sync;

        register const PixelPacket
          *__restrict p;

        register long
          x;

        if (((y+y_offset) < 0) || ((y+y_offset) >= (long) image->rows))
          continue;
        p=GetCacheViewVirtualPixels(composite_view,0,y,composite_image->columns,
          1,exception);
        r=QueueCacheViewAuthenticPixels(dest_view,0,y,
          dest_image->columns,1,&image->exception);
        if ((p == (const PixelPacket *) NULL) || (r == (PixelPacket *) NULL))
          break;
        dest_indexes=GetCacheViewAuthenticIndexQueue(dest_view);
        for (x=0; x < (long) composite_image->columns; x++)
        {
          if (((x_offset+x) < 0) || ((x_offset+x) >= (long) image->columns))
            {
              p++;
              continue;
            }
          /* Get the blurred pixel */
          ScaleResampleFilter(resample_filter,
               blur_xx*QuantumScale*p->red,
               blur_xy*QuantumScale*p->red,
               blur_yx*QuantumScale*p->green,
               blur_yy*QuantumScale*p->green );
          (void) ResamplePixelColor(resample_filter,
               (double) x_offset+x,  (double) y_offset+y,   &pixel);
          /* Save the blurred lookup */
          SetPixelPacket(dest_image,&pixel,r,dest_indexes+x);
          p++;
          r++;
        }
        sync=SyncCacheViewAuthenticPixels(dest_view,exception);
        if (sync == MagickFalse)
          break;
      }
      resample_filter=DestroyResampleFilter(resample_filter);
      composite_view=DestroyCacheView(composite_view);
      dest_view=DestroyCacheView(dest_view);
      composite_image=dest_image;
      break;
    }
    case DisplaceCompositeOp:
    case DistortCompositeOp:
    {
      /* Displace/Distort based on overlay gradient map
         X = red_channel    Y = green_channel
         compose:args =  x_scale[,y_scale[,x_center,y_center]]
       */
      MagickPixelPacket
        pixel;

      MagickRealType
        horizontal_scale,
        vertical_scale,
        x_center,
        y_center,
        x_lookup,
        y_lookup;

      register IndexPacket
        *__restrict dest_indexes;

      register PixelPacket
        *__restrict r;

      ResampleFilter
        *resample_filter;

      ViewInfo
        *composite_view,
        *dest_view;
#if 0
NB: image_view and the related 'q' PixelPacket has been replaced by
ResamplePixelColor().  This also allows access to virtual pixels.
        *image_view;
#endif

      /*
        Allocate the destination image.
      */
      dest_image=CloneImage(composite_image,composite_image->columns,
        composite_image->rows,MagickTrue,&image->exception);
      if (dest_image == (Image *) NULL)
        return(MagickFalse);
      /*
        Determine the horizontal and vertical displacement scale.
      */
      SetGeometryInfo(&geometry_info);
      flags=NoValue;
      value=GetImageArtifact(composite_image,"compose:args");
      if (value != (char *) NULL)
        flags=ParseGeometry(value,&geometry_info);
      if ((flags & (WidthValue|HeightValue)) == 0 )
        if ((flags & AspectValue) == 0)
          {
            horizontal_scale=(MagickRealType) composite_image->columns/2;
            vertical_scale=(MagickRealType) composite_image->rows/2;
          }
        else
          {
            horizontal_scale=(MagickRealType) image->columns/2;
            vertical_scale=(MagickRealType) image->rows/2;
          }
      else
        {
          horizontal_scale = geometry_info.rho;
          vertical_scale = geometry_info.sigma;
          if ((flags & PercentValue) != 0)
            if ((flags & AspectValue) == 0)
              {
                horizontal_scale *= composite_image->columns / 200.0;
                vertical_scale *= composite_image->rows / 200.0;
              }
            else
              {
                horizontal_scale *= image->columns / 200.0;
                vertical_scale *= image->rows / 200.0;
              }
          if ((flags & HeightValue) == 0)
            vertical_scale=horizontal_scale;
        }
      /*
        Determine fixed center point for absolute distortion map
         Absolute distort ==
           Displace lookup relative to a fixed absolute point
           Select that point according to +X+Y user inputs.
           default = center of overlay image
           flag '!' = locations/percentage relative to background image
      */
      x_center = x_offset;
      y_center = y_offset;
      if ( compose == DistortCompositeOp )
        {
          if ((flags & XValue) == 0)
            if ((flags & AspectValue) == 0)
              x_center = x_offset + composite_image->columns/2;
            else
              x_center = image->columns/2;
          else
            if ((flags & AspectValue) == 0)
              x_center = x_offset + geometry_info.xi;
            else
              x_center = geometry_info.xi;

          if ((flags & YValue) == 0)
            if ((flags & AspectValue) == 0)
              y_center = y_offset + composite_image->rows/2;
            else
              y_center = image->rows/2;
          else
            if ((flags & AspectValue) == 0)
              y_center = y_offset + geometry_info.psi;
            else
              y_center = geometry_info.psi;
        }

      /*
        Shift the pixel lookup point as defined by the provided,
        displacement/distortion map.  -- Like a lens...
      */
      pixel=zero;
      exception=(&image->exception);
      resample_filter=AcquireResampleFilter(image,&image->exception);
#if 0
      image_view=AcquireCacheView(image);
#endif
      dest_view=AcquireCacheView(dest_image);
      composite_view=AcquireCacheView(composite_image);
      for (y=0; y < (long) composite_image->rows; y++)
      {
        MagickBooleanType
          sync;

        register const PixelPacket
          *__restrict p;

        register long
          x;

#if 0
        register PixelPacket
          *__restrict q;
#endif

        if (((y+y_offset) < 0) || ((y+y_offset) >= (long) image->rows))
          continue;
        p=GetCacheViewVirtualPixels(composite_view,0,y,composite_image->columns,
          1,exception);
        r=QueueCacheViewAuthenticPixels(dest_view,0,y,
          dest_image->columns,1,&image->exception);
        if ((p == (const PixelPacket *) NULL) || (r == (PixelPacket *) NULL))
          break;
        dest_indexes=GetCacheViewAuthenticIndexQueue(dest_view);
#if 0
        q=GetCacheViewAuthenticPixels(image_view,0,y+y_offset,image->columns,1,
          exception);
        if ((q == (const PixelPacket *) NULL))
          break;
        q+=x_offset;
#endif
        for (x=0; x < (long) composite_image->columns; x++)
        {
          if (((x_offset+x) < 0) || ((x_offset+x) >= (long) image->columns))
            {
              p++;
#if 0
              q++;
#endif
              continue;
            }
          /* Displace the lookup */
          x_lookup=(horizontal_scale*(p->red -
            (((MagickRealType) QuantumRange+1.0)/2.0)))/
            (((MagickRealType) QuantumRange+1.0)/2.0) +
            x_center + (( compose == DisplaceCompositeOp ) ? x : 0 );
          y_lookup=(vertical_scale*(p->green -
            (((MagickRealType) QuantumRange+1.0)/2.0)))/
            (((MagickRealType) QuantumRange+1.0)/2.0) +
            y_center + (( compose == DisplaceCompositeOp ) ? y : 0 );
          /* Get the pixel */
          (void) ResamplePixelColor(resample_filter,
               (double) x_lookup,  (double) y_lookup,   &pixel);
          /* Mask with 'invalid pixel mask' in alpha channel */
          pixel.opacity = (MagickRealType) QuantumRange*(1.0 -
            (1.0-QuantumScale*pixel.opacity)*(1.0-QuantumScale*p->opacity));
          /* save the resulting lookup */
          SetPixelPacket(dest_image,&pixel,r,dest_indexes+x);
          p++;
#if 0
          q++;
#endif
          r++;
        }
        sync=SyncCacheViewAuthenticPixels(dest_view,exception);
        if (sync == MagickFalse)
          break;
      }
      resample_filter=DestroyResampleFilter(resample_filter);
      composite_view=DestroyCacheView(composite_view);
      dest_view=DestroyCacheView(dest_view);
#if 0
      image_view=DestroyCacheView(image_view);
#endif
      composite_image=dest_image;
      break;
    }
    case DissolveCompositeOp:
    {
      /*
        Geometry arguments to dissolve factors.
      */
      value=GetImageArtifact(composite_image,"compose:args");
      if (value != (char *) NULL)
        {
          flags=ParseGeometry(value,&geometry_info);
          source_dissolve=geometry_info.rho/100.0;
          destination_dissolve=1.0;
          if ((source_dissolve-MagickEpsilon) < 0.0)
            source_dissolve=0.0;
          if ((source_dissolve+MagickEpsilon) > 1.0)
            {
              destination_dissolve=2.0-source_dissolve;
              source_dissolve=1.0;
            }
          if ((flags & SigmaValue) != 0)
            destination_dissolve=geometry_info.sigma/100.0;
          if ((destination_dissolve-MagickEpsilon) < 0.0)
            destination_dissolve=0.0;
          modify_outside_overlay=MagickTrue;
          if ((destination_dissolve+MagickEpsilon) > 1.0 )
            {
              destination_dissolve=1.0;
              modify_outside_overlay=MagickFalse;
            }
        }
      break;
    }
    case BlendCompositeOp:
    {
      value=GetImageArtifact(composite_image,"compose:args");
      if (value != (char *) NULL)
        {
          flags=ParseGeometry(value,&geometry_info);
          source_dissolve=geometry_info.rho/100.0;
          destination_dissolve=1.0-source_dissolve;
          if ((flags & SigmaValue) != 0)
            destination_dissolve=geometry_info.sigma/100.0;
          modify_outside_overlay=MagickTrue;
          if ((destination_dissolve+MagickEpsilon) > 1.0)
            modify_outside_overlay=MagickFalse;
        }
      break;
    }
    case ModulateCompositeOp:
    {
      /*
        Determine the brightness and saturation scale.
      */
      value=GetImageArtifact(composite_image,"compose:args");
      if (value != (char *) NULL)
        {
          flags=ParseGeometry(value,&geometry_info);
          percent_brightness=geometry_info.rho;
          if ((flags & SigmaValue) != 0)
            percent_saturation=geometry_info.sigma;
        }
      break;
    }
    case ThresholdCompositeOp:
    {
      /*
        Determine the amount and threshold.
      */
      value=GetImageArtifact(composite_image,"compose:args");
      if (value != (char *) NULL)
        {
          flags=ParseGeometry(value,&geometry_info);
          amount=geometry_info.rho;
          threshold=geometry_info.sigma;
          if ((flags & SigmaValue) == 0)
            threshold=0.05f;
        }
      threshold*=QuantumRange;
      break;
    }
    default:
      break;
  }
  value=GetImageArtifact(composite_image,"compose:outside-overlay");
  if (value != (const char *) NULL)
    modify_outside_overlay=MagickFalse;
  /*
    Composite image.
  */
  status=MagickTrue;
  progress=0;
  midpoint=((MagickRealType) QuantumRange+1.0)/2;
  exception=(&image->exception);
  image_view=AcquireCacheView(image);
  composite_view=AcquireCacheView(composite_image);
#if defined(MAGICKCORE_OPENMP_SUPPORT)
  #pragma omp parallel for schedule(dynamic,4) shared(progress,status)
#endif
  for (y=0; y < (long) image->rows; y++)
  {
    const PixelPacket
      *pixels;

    double
      brightness,
      hue,
      saturation;

    MagickPixelPacket
      composite,
      destination,
      source;

    register const IndexPacket
      *__restrict composite_indexes;

    register const PixelPacket
      *__restrict p;

    register IndexPacket
      *__restrict indexes;

    register long
      x;

    register PixelPacket
      *__restrict q;

    if (status == MagickFalse)
      continue;
    if (modify_outside_overlay == MagickFalse)
      {
        if (y < y_offset)
          continue;
        if ((y-y_offset) >= (long) composite_image->rows)
          continue;
      }
    /*
      If pixels is NULL, y is outside overlay region.
    */
    pixels=(PixelPacket *) NULL;
    p=(PixelPacket *) NULL;
    if ((y >= y_offset) && ((y-y_offset) < (long) composite_image->rows))
      {
        p=GetCacheViewVirtualPixels(composite_view,0,y-y_offset,
          composite_image->columns,1,exception);
        if (p == (const PixelPacket *) NULL)
          {
            status=MagickFalse;
            continue;
          }
        pixels=p;
        if (x_offset < 0)
          p-=x_offset;
      }
    q=GetCacheViewAuthenticPixels(image_view,0,y,image->columns,1,
      exception);
    if (q == (PixelPacket *) NULL)
      {
        status=MagickFalse;
        continue;
      }
    indexes=GetCacheViewAuthenticIndexQueue(image_view);
    composite_indexes=GetCacheViewVirtualIndexQueue(composite_view);
    GetMagickPixelPacket(composite_image,&source);
    GetMagickPixelPacket(image,&destination);
    hue=0.0;
    saturation=0.0;
    brightness=0.0;
    for (x=0; x < (long) image->columns; x++)
    {
      if (modify_outside_overlay == MagickFalse)
        {
          if (x < x_offset)
            {
              q++;
              continue;
            }
          if ((x-x_offset) >= (long) composite_image->columns)
            break;
        }
      destination.red=(MagickRealType) q->red;
      destination.green=(MagickRealType) q->green;
      destination.blue=(MagickRealType) q->blue;
      if (image->matte != MagickFalse)
        destination.opacity=(MagickRealType) q->opacity;
      if (image->colorspace == CMYKColorspace)
        {
          destination.red=(MagickRealType) QuantumRange-destination.red;
          destination.green=(MagickRealType) QuantumRange-destination.green;
          destination.blue=(MagickRealType) QuantumRange-destination.blue;
          destination.index=(MagickRealType) (QuantumRange-indexes[x]);
        }
      /*
        Handle destination modifications outside overlaid region.
      */
      composite=destination;
      if ((pixels == (PixelPacket *) NULL) || (x < x_offset) ||
          ((x-x_offset) >= (long) composite_image->columns))
        {
          switch (compose)
          {
            case DissolveCompositeOp:
            case BlendCompositeOp:
            {
              composite.opacity=(MagickRealType) (QuantumRange-
                destination_dissolve*(QuantumRange-composite.opacity));
              break;
            }
            case ClearCompositeOp:
            case SrcCompositeOp:
            {
              CompositeClear(&destination,&composite);
              break;
            }
            case InCompositeOp:
            case SrcInCompositeOp:
            case SrcOutCompositeOp:
            case DstInCompositeOp:
            case DstAtopCompositeOp:
            case CopyOpacityCompositeOp:
            case ChangeMaskCompositeOp:
            {
              composite.opacity=(MagickRealType) TransparentOpacity;
              break;
            }
            default:
              break;
          }
          if (image->colorspace == CMYKColorspace)
            {
              composite.red=(MagickRealType) QuantumRange-composite.red;
              composite.green=(MagickRealType) QuantumRange-composite.green;
              composite.blue=(MagickRealType) QuantumRange-composite.blue;
              composite.index=(MagickRealType) QuantumRange-composite.index;
            }
          q->red=RoundToQuantum(composite.red);
          q->green=RoundToQuantum(composite.green);
          q->blue=RoundToQuantum(composite.blue);
          if (image->matte != MagickFalse)
            q->opacity=RoundToQuantum(composite.opacity);
          if (image->colorspace == CMYKColorspace)
            indexes[x]=RoundToQuantum(composite.index);
          q++;
          continue;
        }
      /*
        Handle normal overlay of source onto destination.
      */
      source.red=(MagickRealType) p->red;
      source.green=(MagickRealType) p->green;
      source.blue=(MagickRealType) p->blue;
      if (composite_image->matte != MagickFalse)
        source.opacity=(MagickRealType) p->opacity;
      if (composite_image->colorspace == CMYKColorspace)
        {
          source.red=(MagickRealType) QuantumRange-source.red;
          source.green=(MagickRealType) QuantumRange-source.green;
          source.blue=(MagickRealType) QuantumRange-source.blue;
          source.index=(MagickRealType) QuantumRange-(MagickRealType)
            composite_indexes[x-x_offset];
        }
      switch (compose)
      {
        case AddCompositeOp:
        {
          CompositeAdd(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case ClearCompositeOp:
        {
          CompositeClear(&destination,&composite);
          break;
        }
        case SrcCompositeOp:
        case CopyCompositeOp:
        case ReplaceCompositeOp:
        {
          composite=source;
          break;
        }
        case ChangeMaskCompositeOp:
        {
          if ((composite.opacity > ((MagickRealType) QuantumRange/2.0)) ||
              IsMagickColorSimilar(&source,&destination))
            composite.opacity=(MagickRealType) TransparentOpacity;
          else
            composite.opacity=(MagickRealType) OpaqueOpacity;
          break;
        }
        case DivideCompositeOp:
        {
          CompositeDivide(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case DstCompositeOp:
          break;
        case OverCompositeOp:
        case SrcOverCompositeOp:
        {
          CompositeOver(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case DstOverCompositeOp:
        {
          CompositeOver(&destination,destination.opacity,&source,source.opacity,
            &composite);
          break;
        }
        case SrcInCompositeOp:
        case InCompositeOp:
        {
          CompositeIn(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case DstInCompositeOp:
        {
          CompositeIn(&destination,destination.opacity,&source,source.opacity,
            &composite);
          break;
        }
        case OutCompositeOp:
        case SrcOutCompositeOp:
        {
          CompositeOut(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case DstOutCompositeOp:
        {
          CompositeOut(&destination,destination.opacity,&source,source.opacity,
            &composite);
          break;
        }
        case AtopCompositeOp:
        case SrcAtopCompositeOp:
        {
          CompositeAtop(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case DstAtopCompositeOp:
        {
          CompositeAtop(&destination,destination.opacity,&source,source.opacity,
            &composite);
          break;
        }
        case XorCompositeOp:
        {
          CompositeXor(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case PlusCompositeOp:
        {
          CompositePlus(&source,source.opacity,&destination,destination.opacity,
            &composite);
          break;
        }
        case MultiplyCompositeOp:
        {
          CompositeMultiply(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case ScreenCompositeOp:
        {
          CompositeScreen(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case OverlayCompositeOp:
        {
          CompositeOverlay(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case DarkenCompositeOp:
        {
          CompositeDarken(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case LightenCompositeOp:
        {
          CompositeLighten(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case LinearLightCompositeOp:
        {
          CompositeLinearLight(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case ColorDodgeCompositeOp:
        {
          CompositeColorDodge(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case ColorBurnCompositeOp:
        {
          CompositeColorBurn(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case HardLightCompositeOp:
        {
          CompositeHardLight(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case SoftLightCompositeOp:
        {
          CompositeSoftLight(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case DifferenceCompositeOp:
        {
          CompositeDifference(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case ExclusionCompositeOp:
        {
          CompositeExclusion(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case MinusCompositeOp:
        {
          CompositeMinus(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case BumpmapCompositeOp:
        {
          if (source.opacity == TransparentOpacity)
            break;
          CompositeBumpmap(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case DissolveCompositeOp:
        {
          CompositeOver(&source,(MagickRealType) (QuantumRange-source_dissolve*
            (QuantumRange-source.opacity)),&destination,(MagickRealType)
            (QuantumRange-destination_dissolve*(QuantumRange-
            destination.opacity)),&composite);
          break;
        }
        case BlendCompositeOp:
        {
          CompositePlus(&source,(MagickRealType) (QuantumRange-source_dissolve*
            (QuantumRange-source.opacity)),&destination,(MagickRealType)
            (QuantumRange-destination_dissolve*(QuantumRange-
            destination.opacity)),&composite);
          break;
        }
        case BlurCompositeOp:
        case DisplaceCompositeOp:
        case DistortCompositeOp:
        {
          composite=source;
          break;
        }
        case ThresholdCompositeOp:
        {
          CompositeThreshold(&source,source.opacity,&destination,
            destination.opacity,threshold,amount,&composite);
          break;
        }
        case ModulateCompositeOp:
        {
          long
            offset;

          if (source.opacity == TransparentOpacity)
            break;
          offset=(long) (MagickPixelIntensityToQuantum(&source)-midpoint);
          if (offset == 0)
            break;
          CompositeHSB(destination.red,destination.green,destination.blue,&hue,
            &saturation,&brightness);
          brightness+=(0.01*percent_brightness*offset)/midpoint;
          saturation*=0.01*percent_saturation;
          HSBComposite(hue,saturation,brightness,&composite.red,
            &composite.green,&composite.blue);
          break;
        }
        case HueCompositeOp:
        {
          if (source.opacity == TransparentOpacity)
            break;
          if (destination.opacity == TransparentOpacity)
            {
              composite=source;
              break;
            }
          CompositeHSB(destination.red,destination.green,destination.blue,&hue,
            &saturation,&brightness);
          CompositeHSB(source.red,source.green,source.blue,&hue,&sans,&sans);
          HSBComposite(hue,saturation,brightness,&composite.red,&composite.green,
            &composite.blue);
          if (source.opacity < destination.opacity)
            composite.opacity=source.opacity;
          break;
        }
        case SaturateCompositeOp:
        {
          if (source.opacity == TransparentOpacity)
            break;
          if (destination.opacity == TransparentOpacity)
            {
              composite=source;
              break;
            }
          CompositeHSB(destination.red,destination.green,destination.blue,&hue,
            &saturation,&brightness);
          CompositeHSB(source.red,source.green,source.blue,&sans,&saturation,
            &sans);
          HSBComposite(hue,saturation,brightness,&composite.red,&composite.green,
            &composite.blue);
          if (source.opacity < destination.opacity)
            composite.opacity=source.opacity;
          break;
        }
        case SubtractCompositeOp:
        {
          CompositeSubtract(&source,source.opacity,&destination,
            destination.opacity,&composite);
          break;
        }
        case LuminizeCompositeOp:
        {
          if (source.opacity == TransparentOpacity)
            break;
          if (destination.opacity == TransparentOpacity)
            {
              composite=source;
              break;
            }
          CompositeHSB(destination.red,destination.green,destination.blue,&hue,
            &saturation,&brightness);
          CompositeHSB(source.red,source.green,source.blue,&sans,&sans,
            &brightness);
          HSBComposite(hue,saturation,brightness,&composite.red,&composite.green,
            &composite.blue);
          if (source.opacity < destination.opacity)
            composite.opacity=source.opacity;
          break;
        }
        case ColorizeCompositeOp:
        {
          if (source.opacity == TransparentOpacity)
            break;
          if (destination.opacity == TransparentOpacity)
            {
              composite=source;
              break;
            }
          CompositeHSB(destination.red,destination.green,destination.blue,&sans,
            &sans,&brightness);
          CompositeHSB(source.red,source.green,source.blue,&hue,&saturation,
            &sans);
          HSBComposite(hue,saturation,brightness,&composite.red,
            &composite.green,&composite.blue);
          if (source.opacity < destination.opacity)
            composite.opacity=source.opacity;
          break;
        }
        case CopyRedCompositeOp:
        case CopyCyanCompositeOp:
        {
          composite.red=source.red;
          break;
        }
        case CopyGreenCompositeOp:
        case CopyMagentaCompositeOp:
        {
          composite.green=source.green;
          break;
        }
        case CopyBlueCompositeOp:
        case CopyYellowCompositeOp:
        {
          composite.blue=source.blue;
          break;
        }
        case CopyOpacityCompositeOp:
        {
          if (source.matte == MagickFalse)
            {
              composite.opacity=(MagickRealType) (QuantumRange-
                MagickPixelIntensityToQuantum(&source));
              break;
            }
          composite.opacity=source.opacity;
          break;
        }
        case CopyBlackCompositeOp:
        {
          if (source.colorspace != CMYKColorspace)
            ConvertRGBToCMYK(&source);
          composite.index=source.index;
          break;
        }
        default:
          break;
      }
      if (image->colorspace == CMYKColorspace)
        {
          composite.red=(MagickRealType) QuantumRange-composite.red;
          composite.green=(MagickRealType) QuantumRange-composite.green;
          composite.blue=(MagickRealType) QuantumRange-composite.blue;
          composite.index=(MagickRealType) QuantumRange-composite.index;
        }
      q->red=RoundToQuantum(composite.red);
      q->green=RoundToQuantum(composite.green);
      q->blue=RoundToQuantum(composite.blue);
      q->opacity=RoundToQuantum(composite.opacity);
      if (image->colorspace == CMYKColorspace)
        indexes[x]=RoundToQuantum(composite.index);
      p++;
      if (p >= (pixels+composite_image->columns))
        p=pixels;
      q++;
    }
    if (SyncCacheViewAuthenticPixels(image_view,exception) == MagickFalse)
      status=MagickFalse;
    if (image->progress_monitor != (MagickProgressMonitor) NULL)
      {
        MagickBooleanType
          proceed;

#if defined(MAGICKCORE_OPENMP_SUPPORT)
  #pragma omp critical (MagickCore_CompositeImageChannel)
#endif
        proceed=SetImageProgress(image,CompositeImageTag,progress++,
          image->rows);
        if (proceed == MagickFalse)
          status=MagickFalse;
      }
  }
  composite_view=DestroyCacheView(composite_view);
  image_view=DestroyCacheView(image_view);
  if (dest_image != (Image * ) NULL)
    dest_image=DestroyImage(dest_image);
  return(status);
}
