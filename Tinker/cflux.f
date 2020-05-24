c
c     ##################################################################
c     ##                                                              ##
c     ##  module cflux -- charge flux term forcefield parameters      ##
c     ##                                                              ##
c     ##################################################################
c
c
c     jb        charge flux of unit bond length
c     ja        charge flux of unit angle
c     jbp1      charge flux of unit bond length for b'(1)
c     jbp2      charge flux of unit bond length for b'(2)
c     jtheta1   charge flux of unit angle for theta1 
c     jtheta2   charge flux of unit angle for theta2 
c     bp1       bond length in charge flux angle term  b'(1)
c     bp2       bond length in charge flux angle term  b'(2)
c     
c
      module cflux
      use sizes
      implicit none
      real*8, allocatable :: bp1(:),bp2(:)
      real*8, allocatable :: jbp1(:),jbp2(:)
      real*8, allocatable :: jtheta1(:),jtheta2(:)
      real*8, allocatable :: jb(:)
      real*8, allocatable :: pchrgflux(:)
      logical dobondcflux,doanglecflux 
      save
      end
