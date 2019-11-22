c     
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine chrgflux  --  charge flux calculation     ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "chrgflux" calculates the charge flux and return dQ
c
      subroutine chrgflux
      use cflux
      use atoms
      implicit none
      integer i
c
c     perform dynamic allocation of global array 
c
      if (allocated(pchrgflux)) deallocate (pchrgflux)

      allocate (pchrgflux(n))
c
c     zero out the charge flux at each site
c
      do i = 1, n
        pchrgflux(i) = 0.0d0
      end do
c
c     chose the method for summing over charge flux
c
      if (dobondcflux)  call chrgfluxb
      if (doanglecflux) call chrgfluxa
      return
      end
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine chrgfluxb  --  charge flux --bond         ##
c     ##                                                       ##
c     ###########################################################
c
c
c     charge flux due to bond stectching
c
      subroutine chrgfluxb
      use sizes
      use atomid
      use atoms
      use bndstr 
      use bound  
      use cflux
      use couple
      use group
      use potent
      use usage
      implicit none
      real*8 dq 
      real*8 xab,yab,zab,rab
      real*8 pjb,pb0 
      integer atoma,atomb
      integer ia,ib,i,j
      integer nha,nhb,n12a,n12b
      real*8 bigsign

c
c     calculate the bond stretching energy term
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         atoma = atomic(ia)
         atomb = atomic(ib)
         pjb = jb(i)
         pb0 = b0(i)
c
c     compute the value of the bond length deviation
c
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         if (use_polymer)  call image (xab,yab,zab)
         rab = sqrt(xab*xab + yab*yab + zab*zab)
         dq = pjb*(rab-pb0)
c
c     determine the BIGGER atom
c
         if (atoma .ne. atomb) then
            if (atoma .gt. atomb) then
              bigsign = -1.0d0
            else
              bigsign = 1.0d0
            end if
         else 
            n12a = n12(ia)
            n12b = n12(ib)
            if (n12a .ne. n12b) then
              if (n12a .gt. n12b) then
                bigsign = -1.0d0
              else
                bigsign = 1.0d0
              end if
            else
              nha = 0
              nhb = 0
              do j = 1, n12a
                if (atomic(i12(j,ia)) .eq. 1) then
                  nha = nha + 1
                end if
              end do
              do j = 1, n12b
                if (atomic(i12(j,ib)) .eq. 1) then
                  nhb = nhb + 1
                end if
              end do
              if (nha .ne. nhb) then
                if (nha .gt. nhb) then
                  bigsign = -1.0d0
                else
                  bigsign = 1.0d0
                end if
              else
                bigsign = 0.0d0
              end if
            end if
         end if
         pchrgflux(ia) = pchrgflux(ia) + dq*bigsign
         pchrgflux(ib) = pchrgflux(ib) - dq*bigsign
      end do
      return
      end

c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine chrgfluxa  --  charge flux --angle        ##
c     ##                                                       ##
c     ###########################################################
c
c
c    charge flux due to angle bending
c
      subroutine chrgfluxa
      use atomid
      use atoms 
      use cflux
      use sizes
      use angbnd
      use angpot
      use bound
      use energi
      use group
      use math
      use usage
      implicit none
      integer i,ia,ib,ic
      real*8 angle
      real*8 rab,rcb
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 dot,cosine
      real*8 ptheta0,pb10,pb20
      real*8 pjbp1,pjbp2
      real*8 pjtheta1,pjtheta2
      real*8 dq1,dq2

c
c     loop over all angle in the system 
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
c
c     assign the charge flux parameters  
c
         ptheta0 = theta0(i)
         pb10 = bond1(i)
         pb20 = bond2(i)
         pjbp1 = jbp1(i)
         pjbp2 = jbp2(i)
         pjtheta1 = jtheta1(i) 
         pjtheta2 = jtheta2(i)
c
c     calculate the bond length and angle 
c
         xia = x(ia)
         yia = y(ia)
         zia = z(ia)
         xib = x(ib)
         yib = y(ib)
         zib = z(ib)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xab = xia - xib
         yab = yia - yib
         zab = zia - zib
         xcb = xic - xib
         ycb = yic - yib
         zcb = zic - zib
         if (use_polymer) then
            call image (xab,yab,zab)
            call image (xcb,ycb,zcb)
         end if
         rab = sqrt(xab*xab + yab*yab + zab*zab)
         rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
         if (rab.ne.0.0d0 .and. rcb.ne.0.0d0) then
            dot = xab*xcb + yab*ycb + zab*zcb
            cosine = dot / (rab*rcb)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle = radian * acos(cosine)
         end if
c
c     charge flux increment for ange ia-ib-ic
c
         dq1 = pjbp1*(rcb-pb20)+pjtheta1*(angle-ptheta0)
         dq2 = pjbp2*(rab-pb10)+pjtheta2*(angle-ptheta0)
         pchrgflux(ia) = pchrgflux(ia) + dq1
         pchrgflux(ic) = pchrgflux(ic) + dq2
         pchrgflux(ib) = pchrgflux(ib) - (dq1+dq2)
      end do
      return
      end 
