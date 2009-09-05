/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%        SSSSS  EEEEE  M   M   AAA   PPPP   H   H   OOO   RRRR   EEEEE        %
%        SS     E      MM MM  A   A  P   P  H   H  O   O  R   R  E            %
%         SSS   EEE    M M M  AAAAA  PPPP   HHHHH  O   O  RRRR   EEE          %
%           SS  E      M   M  A   A  P      H   H  O   O  R R    E            %
%        SSSSS  EEEEE  M   M  A   A  P      H   H   OOO   R  R   EEEEE        %
%                                                                             %
%                                                                             %
%                        MagickCore Semaphore Methods                         %
%                                                                             %
%                              Software Design                                %
%                             William Radcliffe                               %
%                                John Cristy                                  %
%                                 June 2000                                   %
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
#include "magick/exception.h"
#include "magick/exception-private.h"
#include "magick/memory_.h"
#include "magick/semaphore.h"
#include "magick/string_.h"
#if defined(MAGICKCORE_HAVE_PTHREAD)
#include <pthread.h>
#endif
#if defined(__WINDOWS__)
#include <windows.h>
#endif

/*
  Struct declaractions.
*/
struct SemaphoreInfo
{
#if defined(MAGICKCORE_HAVE_PTHREAD)
  pthread_mutex_t
#elif defined(__WINDOWS__)
  CRITICAL_SECTION
#else
  unsigned long
#endif
    mutex;

#if defined(MAGICKCORE_DEBUG)
#if defined(MAGICKCORE_HAVE_PTHREAD)
  pthread_t
#elif defined(__WINDOWS__)
  DWORD
#else
  pid_t
#endif
    id;

  long
    reference_count;
#endif

  unsigned long
    signature;
};

/*
  Static declaractions.
*/
#if defined(MAGICKCORE_HAVE_PTHREAD)
static pthread_mutex_t
  semaphore_mutex = PTHREAD_MUTEX_INITIALIZER;
#elif defined(__WINDOWS__)
static LONG
  semaphore_mutex = 0;
#else
static long
  semaphore_mutex = 0;
#endif

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   A c q u i r e S e m a p h o r e I n f o                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  AcquireSemaphoreInfo() acquires a semaphore.
%
%  The format of the AcquireSemaphoreInfo method is:
%
%      void AcquireSemaphoreInfo(SemaphoreInfo **semaphore_info)
%
%  A description of each parameter follows:
%
%    o semaphore_info: Specifies a pointer to an SemaphoreInfo structure.
%
*/
MagickExport void AcquireSemaphoreInfo(SemaphoreInfo **semaphore_info)
{
  assert(semaphore_info != (SemaphoreInfo **) NULL);
#if defined(MAGICKCORE_HAVE_PTHREAD)
  if (pthread_mutex_lock(&semaphore_mutex) != 0)
    (void) fprintf(stderr,"pthread_mutex_lock failed %s\n",strerror(errno));
#elif defined(__WINDOWS__)
  while (InterlockedCompareExchange(&semaphore_mutex,1L,0L) != 0)
    Sleep(10);
#endif
  if (*semaphore_info == (SemaphoreInfo *) NULL)
    *semaphore_info=AllocateSemaphoreInfo();
#if defined(MAGICKCORE_HAVE_PTHREAD)
  if (pthread_mutex_unlock(&semaphore_mutex) != 0)
    (void) fprintf(stderr,"pthread_mutex_unlock failed %s\n",strerror(errno));
#elif defined(__WINDOWS__)
  InterlockedExchange(&semaphore_mutex,0L);
#endif
  (void) LockSemaphoreInfo(*semaphore_info);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   A l l o c a t e S e m a p h o r e I n f o                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  AllocateSemaphoreInfo() initializes the SemaphoreInfo structure.
%
%  The format of the AllocateSemaphoreInfo method is:
%
%      SemaphoreInfo *AllocateSemaphoreInfo(void)
%
*/
MagickExport SemaphoreInfo *AllocateSemaphoreInfo(void)
{
  SemaphoreInfo
    *semaphore_info;

  /*
    Allocate semaphore.
  */
  semaphore_info=(SemaphoreInfo *) malloc(sizeof(SemaphoreInfo));
  if (semaphore_info == (SemaphoreInfo *) NULL)
    ThrowFatalException(ResourceLimitFatalError,"MemoryAllocationFailed");
  (void) ResetMagickMemory(semaphore_info,0,sizeof(SemaphoreInfo));
  /*
    Initialize the semaphore.
  */
#if defined(MAGICKCORE_HAVE_PTHREAD)
  {
    int
      status;

    pthread_mutexattr_t
      mutex_info;

    status=pthread_mutexattr_init(&mutex_info);
    if (status != 0)
      {
        free(semaphore_info);
        ThrowFatalException(ResourceLimitFatalError,
          "UnableToInitializeSemaphore");
      }
    status=pthread_mutex_init(&semaphore_info->mutex,&mutex_info);
    (void) pthread_mutexattr_destroy(&mutex_info);
    if (status != 0)
      {
        free(semaphore_info);
        ThrowFatalException(ResourceLimitFatalError,
          "UnableToInitializeSemaphore");
      }
#if defined(MAGICKCORE_DEBUG)
    semaphore_info->id=pthread_self();
#endif
  }
#elif defined(__WINDOWS__)
  InitializeCriticalSection(&semaphore_info->mutex);
#if defined(MAGICKCORE_DEBUG)
  semaphore_info->id=GetCurrentThreadId();
#endif
#endif
#if defined(MAGICKCORE_DEBUG)
  semaphore_info->reference_count=0;
#endif
  semaphore_info->signature=MagickSignature;
  return(semaphore_info);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   D e s t r o y S e m a p h o r e                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DestroySemaphore() destroys the semaphore environment.
%
%  The format of the DestroySemaphore method is:
%
%      DestroySemaphore(void)
%
*/
MagickExport void DestroySemaphore(void)
{
#if defined(MAGICKCORE_HAVE_PTHREAD)
  if (pthread_mutex_destroy(&semaphore_mutex) != 0)
    (void) fprintf(stderr,"pthread_mutex_destroy failed %s\n",
      strerror(errno));
#endif
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   D e s t r o y S e m a p h o r e I n f o                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DestroySemaphoreInfo() destroys a semaphore.
%
%  The format of the DestroySemaphoreInfo method is:
%
%      void DestroySemaphoreInfo(SemaphoreInfo **semaphore_info)
%
%  A description of each parameter follows:
%
%    o semaphore_info: Specifies a pointer to an SemaphoreInfo structure.
%
*/
MagickExport void DestroySemaphoreInfo(SemaphoreInfo **semaphore_info)
{
  assert(semaphore_info != (SemaphoreInfo **) NULL);
  assert((*semaphore_info) != (SemaphoreInfo *) NULL);
  assert((*semaphore_info)->signature == MagickSignature);
#if defined(MAGICKCORE_HAVE_PTHREAD)
  (void) pthread_mutex_lock(&semaphore_mutex);
#elif defined(__WINDOWS__)
  while (InterlockedCompareExchange(&semaphore_mutex,1L,0L) != 0)
    Sleep(10);
#endif
#if defined(MAGICKCORE_HAVE_PTHREAD)
  (void) pthread_mutex_destroy(&(*semaphore_info)->mutex);
#elif defined(__WINDOWS__)
  DeleteCriticalSection(&(*semaphore_info)->mutex);
#endif
  (*semaphore_info)->signature=(~MagickSignature);
  free(*semaphore_info);
  *semaphore_info=(SemaphoreInfo *) NULL;
#if defined(MAGICKCORE_HAVE_PTHREAD)
  (void) pthread_mutex_unlock(&semaphore_mutex);
#elif defined(__WINDOWS__)
  InterlockedExchange(&semaphore_mutex,0L);
#endif
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   I n i t i a l i z e S e m a p h o r e                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  InitializeSemaphore() initializes the semaphore environment.
%
%  The format of the InitializeSemaphore method is:
%
%      InitializeSemaphore(void)
%
*/
MagickExport void InitializeSemaphore(void)
{
#if defined(MAGICKCORE_HAVE_PTHREAD)
  (void) pthread_mutex_init(&semaphore_mutex,
    (const pthread_mutexattr_t *) NULL);
#endif
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   L o c k S e m a p h o r e I n f o                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LockSemaphoreInfo() locks a semaphore.
%
%  The format of the LockSemaphoreInfo method is:
%
%      MagickBooleanType LockSemaphoreInfo(SemaphoreInfo *semaphore_info)
%
%  A description of each parameter follows:
%
%    o semaphore_info: Specifies a pointer to an SemaphoreInfo structure.
%
*/
MagickExport MagickBooleanType LockSemaphoreInfo(SemaphoreInfo *semaphore_info)
{
  assert(semaphore_info != (SemaphoreInfo *) NULL);
  assert(semaphore_info->signature == MagickSignature);
#if defined(MAGICKCORE_HAVE_PTHREAD)
  {
    int
      status;

    status=pthread_mutex_lock(&semaphore_info->mutex);
    if (status != 0)
      return(MagickFalse);
#if defined(MAGICKCORE_DEBUG)
    {
      pthread_t
        id;

      id=pthread_self();
      if ((semaphore_info->reference_count > 0) &&
          (pthread_equal(semaphore_info->id,id) != 0))
        {
          (void) fprintf(stderr,"Warning: unexpected recursive lock!\n");
          (void) fflush(stderr);
        }
      semaphore_info->id=id;
    }
#endif
  }
#elif defined(__WINDOWS__)
  {
    EnterCriticalSection(&semaphore_info->mutex);
#if defined(MAGICKCORE_DEBUG)
    {
      DWORD
        id;

      id=GetCurrentThreadId();
      if ((semaphore_info->reference_count > 0) &&
          (semaphore_info->id == id))
        {
          (void) fprintf(stderr,"Warning: unexpected recursive lock!\n");
          (void) fflush(stderr);
        }
      semaphore_info->id=id;
    }
#endif
  }
#endif
#if defined(MAGICKCORE_DEBUG)
  semaphore_info->reference_count++;
#endif
  return(MagickTrue);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   R e l i n g u i s h S e m a p h o r e I n f o                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  RelinquishSemaphoreInfo() relinquishes a semaphore.
%
%  The format of the RelinquishSemaphoreInfo method is:
%
%      RelinquishSemaphoreInfo(SemaphoreInfo *semaphore_info)
%
%  A description of each parameter follows:
%
%    o semaphore_info: Specifies a pointer to an SemaphoreInfo structure.
%
*/
MagickExport void RelinquishSemaphoreInfo(SemaphoreInfo *semaphore_info)
{
  assert(semaphore_info != (SemaphoreInfo *) NULL);
  assert(semaphore_info->signature == MagickSignature);
  (void) UnlockSemaphoreInfo(semaphore_info);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   U n l o c k S e m a p h o r e I n f o                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  UnlockSemaphoreInfo() unlocks a semaphore.
%
%  The format of the UnlockSemaphoreInfo method is:
%
%      MagickBooleanType UnlockSemaphoreInfo(SemaphoreInfo *semaphore_info)
%
%  A description of each parameter follows:
%
%    o semaphore_info: Specifies a pointer to an SemaphoreInfo structure.
%
*/
MagickExport MagickBooleanType UnlockSemaphoreInfo(
  SemaphoreInfo *semaphore_info)
{
  assert(semaphore_info != (SemaphoreInfo *) NULL);
  assert(semaphore_info->signature == MagickSignature);
#if defined(MAGICKCORE_DEBUG)
  if (semaphore_info->reference_count == 0)
    {
      (void) fprintf(stderr,"Warning: semaphore lock already unlocked!\n");
      (void) fflush(stderr);
      return(MagickFalse);
    }
  semaphore_info->reference_count--;
#endif
#if defined(MAGICKCORE_HAVE_PTHREAD)
  {
    int
      status;

#if defined(MAGICKCORE_DEBUG)
    assert(pthread_equal(semaphore_info->id,pthread_self()));
#endif
    status=pthread_mutex_unlock(&semaphore_info->mutex);
    if (status != 0)
      {
#if defined(MAGICKCORE_DEBUG)
        semaphore_info->reference_count++;
#endif
        return(MagickFalse);
      }
  }
#elif defined(__WINDOWS__)
#if defined(MAGICKCORE_DEBUG)
  assert(semaphore_info->id == GetCurrentThreadId());
#endif
  LeaveCriticalSection(&semaphore_info->mutex);
#endif
  return(MagickTrue);
}
