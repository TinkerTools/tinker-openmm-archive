c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################
c     ##                                            ## 
c     ##  charge flux subroutines for derivatives   ##   
c     ##                                            ##
c     ################################################

c==================================================================c
c    cfbondem calculate force due to charge flux-bond stretching   c
c==================================================================c
      subroutine cfbondem(dmppot,dcfemx,dcfemy,dcfemz)
      use sizes
      use atoms
      use atomid
      use bndstr
      use cflux
      implicit none
      real*8 ddqdx,ddqdy,ddqdz
      real*8 xa,ya,za
      real*8 xb,yb,zb
      real*8 xba,yba,zba
      real*8 rba,rba2
      real*8 pjb
      real*8 frcx,frcy,frcz
      integer ia,ib
      integer i
      real*8 dmppot(*)
      real*8 dcfemx(*)
      real*8 dcfemy(*)
      real*8 dcfemz(*)
c
c     loop over bond to calculate the force
c
      do i = 1, nbond
        ia = ibnd(1,i)
        ib = ibnd(2,i)
        pjb = jb(i)
        xa = x(ia) 
        ya = y(ia) 
        za = z(ia)
        xb = x(ib) 
        yb = y(ib) 
        zb = z(ib) 
        xba = xa - xb 
        yba = ya - yb 
        zba = za - zb
        rba2 = xba**2 + yba**2 + zba**2
        rba = sqrt(rba2)
        ddqdx = pjb*(xa-xb)/rba
        ddqdy = pjb*(ya-yb)/rba
        ddqdz = pjb*(za-zb)/rba
        frcx = -(dmppot(ia)-dmppot(ib))*ddqdx
        frcy = -(dmppot(ia)-dmppot(ib))*ddqdy
        frcz = -(dmppot(ia)-dmppot(ib))*ddqdz
        dcfemx(ia) = dcfemx(ia) + frcx
        dcfemy(ia) = dcfemy(ia) + frcy
        dcfemz(ia) = dcfemz(ia) + frcz
        dcfemx(ib) = dcfemx(ib) - frcx
        dcfemy(ib) = dcfemy(ib) - frcy
        dcfemz(ib) = dcfemz(ib) - frcz
      end do
      return
      end

c=========================================================c
c     cfangleem calculate force due to CF-angle bending   c
c=========================================================c

      subroutine cfangleem(dmppot,dcfemx,dcfemy,dcfemz)
      use sizes
      use atoms
      use atomid
      use angbnd 
      use cflux
      use math
      implicit none
      real*8 xa,ya,za
      real*8 xb,yb,zb
      real*8 xc,yc,zc
      real*8 xba,yba,zba
      real*8 xbc,ybc,zbc
      real*8 rba,rba2,rba3
      real*8 rbc,rbc2,rbc3
      real*8 frcxa1,frcya1,frcza1
      real*8 frcxb1,frcyb1,frczb1
      real*8 frcxc1,frcyc1,frczc1
      real*8 frcxa2,frcya2,frcza2
      real*8 frcxb2,frcyb2,frczb2
      real*8 frcxc2,frcyc2,frczc2
      real*8 pjbp1,pjbp2
      real*8 pjtheta1,pjtheta2
      real*8 term2xa,term2ya,term2za
      real*8 term2xc,term2yc,term2zc
      real*8 dot,term1
      integer ia,ib,ic
      integer i
      real*8 dmppot(*)
      real*8 dcfemx(*)
      real*8 dcfemy(*)
      real*8 dcfemz(*)

      do i = 1, nangle
        ia = iang(1,i)
        ib = iang(2,i)
        ic = iang(3,i)
        pjbp1 = jbp1(i)
        pjbp2 = jbp2(i)
        pjtheta1 = jtheta1(i) 
        pjtheta2 = jtheta2(i) 
        xa = x(ia) 
        ya = y(ia) 
        za = z(ia)
        xb = x(ib) 
        yb = y(ib) 
        zb = z(ib) 
        xc = x(ic) 
        yc = y(ic) 
        zc = z(ic) 
        xba = xa - xb 
        yba = ya - yb 
        zba = za - zb
        xbc = xc - xb 
        ybc = yc - yb 
        zbc = zc - zb
        rba2 = xba**2 + yba**2 + zba**2
        rba  = sqrt(rba2)
        rba3 = rba2*rba
        rbc2 = xbc**2 + ybc**2 + zbc**2
        rbc  = sqrt(rbc2)
        rbc3 = rbc2*rbc
c
c       terms due to bond stretching in angle

c       for site a
        frcxa1 = -(dmppot(ib)-dmppot(ic))*pjbp2*xba/rba
        frcya1 = -(dmppot(ib)-dmppot(ic))*pjbp2*yba/rba
        frcza1 = -(dmppot(ib)-dmppot(ic))*pjbp2*zba/rba
c       for site c
        frcxc1 = -(dmppot(ib)-dmppot(ia))*pjbp1*xbc/rbc
        frcyc1 = -(dmppot(ib)-dmppot(ia))*pjbp1*ybc/rbc
        frczc1 = -(dmppot(ib)-dmppot(ia))*pjbp1*zbc/rbc
c       for site b
        frcxb1 = -(frcxa1 + frcxc1) 
        frcyb1 = -(frcya1 + frcyc1) 
        frczb1 = -(frcza1 + frczc1) 
c
c       terms due to angle bending
c
        dot = xba*xbc + yba*ybc + zba*zbc
        term1 = -rba*rbc/sqrt(rba2*rbc2-dot**2)

        term2xa = xbc/(rba*rbc) -xba*dot/(rba3*rbc)
        term2xc = xba/(rba*rbc) -xbc*dot/(rba*rbc3)

        term2ya = ybc/(rba*rbc) -yba*dot/(rba3*rbc)
        term2yc = yba/(rba*rbc) -ybc*dot/(rba*rbc3)

        term2za = zbc/(rba*rbc) -zba*dot/(rba3*rbc)
        term2zc = zba/(rba*rbc) -zbc*dot/(rba*rbc3)

        frcxa2 = (dmppot(ia)*pjtheta1
     &           - dmppot(ib)*(pjtheta1+pjtheta2)
     &           + dmppot(ic)*pjtheta2)*term1*term2xa*radian
        frcya2 = (dmppot(ia)*pjtheta1
     &           - dmppot(ib)*(pjtheta1+pjtheta2)
     &           + dmppot(ic)*pjtheta2)*term1*term2ya*radian
        frcza2 = (dmppot(ia)*pjtheta1
     &           - dmppot(ib)*(pjtheta1+pjtheta2)
     &           + dmppot(ic)*pjtheta2)*term1*term2za*radian

        frcxc2 = (dmppot(ia)*pjtheta1
     &           - dmppot(ib)*(pjtheta1+pjtheta2)
     &           + dmppot(ic)*pjtheta2)*term1*term2xc*radian
        frcyc2 = (dmppot(ia)*pjtheta1
     &           - dmppot(ib)*(pjtheta1+pjtheta2)
     &           + dmppot(ic)*pjtheta2)*term1*term2yc*radian
        frczc2 = (dmppot(ia)*pjtheta1
     &           - dmppot(ib)*(pjtheta1+pjtheta2)
     &           + dmppot(ic)*pjtheta2)*term1*term2zc*radian

        frcxb2 = - (frcxa2 + frcxc2)
        frcyb2 = - (frcya2 + frcyc2)
        frczb2 = - (frcza2 + frczc2)

        dcfemx(ia) = dcfemx(ia) + frcxa1 + frcxa2 
        dcfemy(ia) = dcfemy(ia) + frcya1 + frcya2
        dcfemz(ia) = dcfemz(ia) + frcza1 + frcza2
        dcfemx(ib) = dcfemx(ib) + frcxb1 + frcxb2
        dcfemy(ib) = dcfemy(ib) + frcyb1 + frcyb2
        dcfemz(ib) = dcfemz(ib) + frczb1 + frczb2
        dcfemx(ic) = dcfemx(ic) + frcxc1 + frcxc2
        dcfemy(ic) = dcfemy(ic) + frcyc1 + frcyc2
        dcfemz(ic) = dcfemz(ic) + frczc1 + frczc2
      end do
      return
      end


c======================================================================c
c     compute the Ewald self-energy force term due to charge flux-bond c
c======================================================================c


      subroutine cfemselfb(ftrm,chrg,pchrgflx,dcfemx,dcfemy,dcfemz)
      use sizes
      use atoms
      use atomid
      use bndstr
      use cflux
      implicit none
      real*8 ddqdx,ddqdy,ddqdz
      real*8 xa,ya,za
      real*8 xb,yb,zb
      real*8 xba,yba,zba
      real*8 rba,rba2
      real*8 pjb
      real*8 ca,cb
      real*8 ftrm
      real*8 frcx,frcy,frcz
      integer ia,ib
      integer i
      real*8 dcfemx(*)
      real*8 dcfemy(*)
      real*8 dcfemz(*)
      real*8 pchrgflx(*)
      real*8 chrg(*)

      do i = 1, nbond
        ia = ibnd(1,i)
        ib = ibnd(2,i)
        pjb = jb(i)
        xa = x(ia) 
        ya = y(ia) 
        za = z(ia)
        xb = x(ib) 
        yb = y(ib) 
        zb = z(ib) 
        xba = xa - xb 
        yba = ya - yb 
        zba = za - zb
        rba2 = xba**2 + yba**2 + zba**2
        rba = sqrt(rba2)
        ca = chrg(ia) + pchrgflx(ia) 
        cb = chrg(ib) + pchrgflx(ib) 
        ddqdx = pjb*(xa-xb)/rba
        ddqdy = pjb*(ya-yb)/rba
        ddqdz = pjb*(za-zb)/rba
        frcx = -2.0d0*ftrm*(ca-cb)*ddqdx
        frcy = -2.0d0*ftrm*(ca-cb)*ddqdy
        frcz = -2.0d0*ftrm*(ca-cb)*ddqdz
        dcfemx(ia) = dcfemx(ia) + frcx
        dcfemy(ia) = dcfemy(ia) + frcy
        dcfemz(ia) = dcfemz(ia) + frcz
        dcfemx(ib) = dcfemx(ib) - frcx
        dcfemy(ib) = dcfemy(ib) - frcy
        dcfemz(ib) = dcfemz(ib) - frcz
      end do
      return
      end
c=====================================================================c
c compute the Ewald self-energy force term due to charge flux-angle   c
c=====================================================================c
      subroutine cfemselfa(ftrm,chrg,pchrgflx,dcfemx,dcfemy,dcfemz)
      use atoms
      use atomid
      use cflux
      use angbnd
      use math
      implicit none
      real*8 ca,cb,cc
      real*8 pjbp1,pjbp2
      real*8 pjtheta1,pjtheta2
      real*8 xa,xb,xc
      real*8 ya,yb,yc
      real*8 za,zb,zc
      real*8 xba,xbc
      real*8 yba,ybc
      real*8 zba,zbc
      real*8 rba,rbc
      real*8 rba2,rbc2
      real*8 rba3,rbc3
      real*8 frcxa1,frcya1,frcza1
      real*8 frcxb1,frcyb1,frczb1
      real*8 frcxc1,frcyc1,frczc1
      real*8 frcxa2,frcya2,frcza2
      real*8 frcxb2,frcyb2,frczb2
      real*8 frcxc2,frcyc2,frczc2
      real*8 ftrm
      real*8 term2xa,term2ya,term2za
      real*8 term2xc,term2yc,term2zc
      real*8 dot,term1,tmpterm
      integer i
      integer ia,ib,ic
      real*8 dcfemx(*)
      real*8 dcfemy(*)
      real*8 dcfemz(*)
      real*8 chrg(*)
      real*8 pchrgflx(*)
 
      do i = 1, nangle
        ia = iang(1,i)
        ib = iang(2,i)
        ic = iang(3,i)
        pjbp1 = jbp1(i)
        pjbp2 = jbp2(i)
        pjtheta1 = jtheta1(i) 
        pjtheta2 = jtheta2(i) 
        xa = x(ia) 
        ya = y(ia) 
        za = z(ia)
        xb = x(ib) 
        yb = y(ib) 
        zb = z(ib) 
        xc = x(ic) 
        yc = y(ic) 
        zc = z(ic) 
        ca = chrg(ia) + pchrgflx(ia) 
        cb = chrg(ib) + pchrgflx(ib) 
        cc = chrg(ic) + pchrgflx(ic) 
        xba = xa - xb 
        yba = ya - yb 
        zba = za - zb
        xbc = xc - xb 
        ybc = yc - yb 
        zbc = zc - zb

        rba2 = xba**2 + yba**2 + zba**2
        rba  = sqrt(rba2)
        rba3 = rba2*rba
        rbc2 = xbc**2 + ybc**2 + zbc**2
        rbc  = sqrt(rbc2)
        rbc3 = rbc2*rbc
c
c       terms due to bond stretching in angle
c
        frcxa1 = -2.0d0*ftrm*pjbp2*(cb-cc)*xba/rba
        frcya1 = -2.0d0*ftrm*pjbp2*(cb-cc)*yba/rba
        frcza1 = -2.0d0*ftrm*pjbp2*(cb-cc)*zba/rba
        frcxc1 = -2.0d0*ftrm*pjbp1*(cb-ca)*xbc/rbc
        frcyc1 = -2.0d0*ftrm*pjbp1*(cb-ca)*ybc/rbc
        frczc1 = -2.0d0*ftrm*pjbp1*(cb-ca)*zbc/rbc
        frcxb1 = -(frcxa1 + frcxc1) 
        frcyb1 = -(frcya1 + frcyc1) 
        frczb1 = -(frcza1 + frczc1) 

c
c       terms due to angle bending
c

        dot = xba*xbc + yba*ybc + zba*zbc
        term1 = -rba*rbc/sqrt(rba2*rbc2-dot**2)

        term2xa = xbc/(rba*rbc) -xba*dot/(rba3*rbc)
        term2xc = xba/(rba*rbc) -xbc*dot/(rba*rbc3)

        term2ya = ybc/(rba*rbc) -yba*dot/(rba3*rbc)
        term2yc = yba/(rba*rbc) -ybc*dot/(rba*rbc3)

        term2za = zbc/(rba*rbc) -zba*dot/(rba3*rbc)
        term2zc = zba/(rba*rbc) -zbc*dot/(rba*rbc3)

        tmpterm = (cb*(pjtheta1+pjtheta2)-ca*pjtheta1
     &           - cc*pjtheta2)*term1*radian

        frcxa2 = -2.0d0*ftrm*tmpterm*term2xa 
        frcya2 = -2.0d0*ftrm*tmpterm*term2ya 
        frcza2 = -2.0d0*ftrm*tmpterm*term2za 
        frcxc2 = -2.0d0*ftrm*tmpterm*term2xc 
        frcyc2 = -2.0d0*ftrm*tmpterm*term2yc 
        frczc2 = -2.0d0*ftrm*tmpterm*term2zc 
        frcxb2 = - (frcxa2 + frcxc2)
        frcyb2 = - (frcya2 + frcyc2)
        frczb2 = - (frcza2 + frczc2)

        dcfemx(ia) = dcfemx(ia) + frcxa1 + frcxa2
        dcfemy(ia) = dcfemy(ia) + frcya1 + frcya2
        dcfemz(ia) = dcfemz(ia) + frcza1 + frcza2
        dcfemx(ib) = dcfemx(ib) + frcxb1 + frcxb2
        dcfemy(ib) = dcfemy(ib) + frcyb1 + frcyb2
        dcfemz(ib) = dcfemz(ib) + frczb1 + frczb2
        dcfemx(ic) = dcfemx(ic) + frcxc1 + frcxc2
        dcfemy(ic) = dcfemy(ic) + frcyc1 + frcyc2
        dcfemz(ic) = dcfemz(ic) + frczc1 + frczc2
      end do
      return 
      end
