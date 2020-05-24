c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kbond  --  bond stretch parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kbond" assigns a force constant and ideal bond length
c     to each bond in the structure and processes any new or
c     changed parameter values
c
c
      subroutine kcflux
      use sizes
      use atoms
      use angbnd
      use atomid
      use bndstr
      use cflux 
      use inform
      use iounit
      use kangs
      use kbonds
      use kcfluxes
      use keys
      use potent
      use usage
      implicit none
      integer i,j,k
      integer ia,ib,ic,ita,itb,itc
      integer na,nb
      integer size,next
      real*8 fc,bd,fj
      real*8 jtt1,jtt2,jb1,jb2
      character*4 pa,pb,pc
      character*8 blank2,pt2,pt,ptab,ptbc
      character*12 blank3,pt3
      character*20 keyword
      character*240 record
      character*240 string
c
c     process keywords containing bond stretch parameters
c
      blank2 = '        '
      blank3 = '            '
      size = 4
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=5,end=5)  ia,ib,fc,bd
   5       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnb
               if (kb(j).eq.blank2 .or. kb(j).eq.pt) then
                  blen(j) = bd
                  goto 8 
               end if
            end do
   8       continue
         else if (keyword(1:8) .eq. 'CFLUX-B ') then
            ia = 0
            ib = 0
            fj = 0.0d0
            dobondcflux = .true.
            string = record(next:240)
            read (string,*,err=10,end=10) ia,ib,fj
   10       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt2 = pa//pb
            else
               pt2 = pb//pa
            end if
            do j = 1, maxnbcf
               if (kcfb(j).eq.blank2 .or. kcfb(j).eq.pt2) then
                  kcfb(j) = pt2
                  jbnd(j) = fj
                  goto 15
               end if
            end do
   15       continue
         else if (keyword(1:8) .eq. 'CFLUX-A ') then
            ia = 0
            ib = 0
            ic = 0
            jtt1 = 0.0d0
            jtt2 = 0.0d0
            jb1 = 0.0d0
            jb2 = 0.0d0
            doanglecflux = .true.
            string = record(next:240)
            read (string,*,err=20,end=20) ia,ib,ic,jtt1,jtt2,jb1,jb2
  20        continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt3 = pa//pb//pc
            else
               pt3 = pc//pb//pa
            end if

            do j = 1, maxnacf
               if (kcfa(j).eq.blank3 .or. kcfa(j).eq.pt3) then
                  kcfa(j) = pt3
                  jtheta1l(j) = jtt1 
                  jtheta2l(j) = jtt2 
                  jbp1l(j) = jb1 
                  jbp2l(j) = jb2 
                  goto 25 
               end if
            end do
   25       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nb = maxnbcf
      do i = maxnbcf, 1, -1
         if (kcfb(i) .eq. blank2)  nb = i - 1
      end do
c
c     determine the total number of forcefield parameters
c
      na = maxnacf
      do i = maxnacf, 1, -1
         if (kcfa(i) .eq. blank3)  na = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(b0))  deallocate (b0)
      if (allocated(jb))  deallocate (jb)
      if (allocated(theta0))  deallocate (theta0)
      if (allocated(bp1))  deallocate (bp1)
      if (allocated(bp2))  deallocate (bp2)
      if (allocated(jbp1))  deallocate (jbp1)
      if (allocated(jbp2))  deallocate (jbp2)
      if (allocated(jtheta1))  deallocate (jtheta1)
      if (allocated(jtheta2))  deallocate (jtheta2)
      allocate (b0(nbond))
      allocate (jb(nbond))
      allocate (theta0(nangle))
      allocate (bp1(nangle))
      allocate (bp2(nangle))
      allocate (jbp1(nangle))
      allocate (jbp2(nangle))
      allocate (jtheta1(nangle))
      allocate (jtheta2(nangle))
c
c     assign bond stretching parameters for each bond
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt2 = pa//pb
         else
            pt2 = pb//pa
         end if
         b0(i) = 0.0d0
         jb(i) = 0.0d0
         do j = 1, nb
           if (kcfb(j) .eq. pt2) then
              jb(i) = jbnd(j)
           end if
         end do
      end do
c
c    assign ideal bond angle and force constant for each angle
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pt3 = pa//pb//pc
         else
            pt3 = pc//pb//pa
         end if

         if (ita .le. itb) then
            ptab = pa//pb
         else
            ptab = pb//pa
         end if

         if (itb .le. itc) then
            ptbc = pb//pc
         else
            ptbc = pc//pb
         end if

         bp1(i) = 0.0d0
         bp2(i) = 0.0d0
         jbp1(i) = 0.0d0
         jbp2(i) = 0.0d0
         jtheta1(i) = 0.0d0
         jtheta2(i) = 0.0d0
         do j = 1, na
           if (kcfa(j) .eq. pt3) then
             jtheta1(i) =jtheta1l(j)  
             jtheta2(i) =jtheta2l(j)  
             jbp1(i) = jbp1l(j) 
             jbp2(i) = jbp2l(j)
           end if
           do k = 1, nbond
              if (kb(k) .eq. ptab) then
                 bp1(i) = blen(k)
              end if
              if (kb(k) .eq. ptbc) then
                 bp2(i) = blen(k)
              end if
           end do
         end do
      end do
c
c     turn off the charge flux if bond and angle are not used 
c
      if (nb .eq. 0 .and. na .eq. 0)  use_cflux = .false.
      return
      end
