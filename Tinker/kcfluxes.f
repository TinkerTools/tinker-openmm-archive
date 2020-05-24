c
c     ##################################################################
c     ##                                                              ##
c     ##  module cflux -- charge flux term forcefield parameters      ##
c     ##                                                              ##
c     ##################################################################
c
c
c     jb        charge flux over unit bond length
c     ja        charge flux over unit angle 
c
      module kcfluxes
      use sizes
      implicit none
      integer maxnbcf
      integer maxnacf
      parameter (maxnbcf=2000)
      parameter (maxnacf=2000)
      real*8 jbnd(maxnbcf)
      real*8 bp1l(maxnacf) 
      real*8 bp2l(maxnacf) 
      real*8 jbp1l(maxnacf)
      real*8 jbp2l(maxnacf)
      real*8 jtheta1l(maxnacf)
      real*8 jtheta2l(maxnacf)
      character*8 kcfb(maxnbcf)
      character*12 kcfa(maxnacf)
      save
      end
