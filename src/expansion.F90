#include "param"
   Subroutine Expansion()
      Use partfield
#ifdef EXPANSION
      Use expansion_params
#endif
      Implicit None
#ifdef EXPANSION
#ifdef ACCELERATED_EXPANSION
!     if (t==0) then
!     read(73,*) ti,ai,daoai,wi,dwowi,usi
!     read(73,*) tf,af,daoaf,wf,dwowf,usf
!     dtwd=tf-ti
!     endif

      if (t .gt. ti+dtwd) then
      ti=tf; ai=af; daoai=daoaf
      wi=wf; dwowi=dwowf; usi=usf
      read(73,*) tf,af,daoaf,wf,dwowf,usf
      dtwd=tf-ti
      endif

      aa=(ai*(tf-t)+af*(t-ti))/dtwd
      daoa=(daoai*(tf-t)+daoaf*(t-ti))/(dtwd*wpiwci)
      ww=(wi*(tf-t)+wf*(t-ti))/dtwd
      dwow=(dwowi*(tf-t)+dwowf*(t-ti))/(dtwd*wpiwci)
      usw=(usi*(tf-t)+usf*(t-ti))/dtwd
#else
      ww=1.; dwow=0.
      aa = 1 + uu0*t
      daoa = uu0/aa
      usw = usw0
#endif
#else
     write(6,*) 'NO EXPANSION'
#endif
   End Subroutine Expansion
