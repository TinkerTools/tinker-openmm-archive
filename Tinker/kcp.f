c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kcp  --  charge pene   parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kcp" assigns parameters to any additional user defined
c
c
      subroutine kcp
      use atoms
      use atomid
      use inform
      use iounit
      use keys
      use chgpen
      use kcps
      use potent
      use sizes
      implicit none
      integer i,k,next
      real*8 alphaf,coref 
      integer ncp
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(penalpha))  deallocate (penalpha)
      if (allocated(pencore))  deallocate (pencore)
      allocate (penalpha(n))
      allocate (pencore(n))
c
c     process keywords containing charge penetration 
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:3) .eq. 'CP ') then
            k = 0
            alphaf = -1.0d0 
            coref = -1.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  alphaf,coref 
   10       continue
            if ((k .gt. 0) .and. (k .le. maxtyp)) then
               apena(k) = alphaf
               apenc(k) = coref
            end if
         end if
      end do
c
c     find and store the charge penetration parameters 
c
      do i = 1, n
         penalpha(i) = apena(type(i))
         pencore(i) = apenc(type(i))
      end do
c
c     process keywords containing atom specific charge penetration 
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:3) .eq. 'CP ') then
            k = 0
            alphaf = 0.0d0
            coref = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:240)
               read (string,*,err=40,end=40) alphaf,coref 
   40          continue
               if (header) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Charge Penetration Parameters',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',8x,'Nuclear Charge'/)
               end if
               if (.not. silent) then
                  write (iout, 60)  k,alphaf,coref
   60             format (4x,i6,10x,f10.3,2x,f10.3)
               end if
               penalpha(k) = alphaf
               pencore(k) = coref
            end if
         end if
      end do
c
c     remove zero and undefined CP sites from the list
c
      ncp = 0
      if (use_mpole) then
         do i = 1, n
            if (penalpha(i).ne.0 .or. pencore(i).ne.0.0d0) then
               ncp = ncp + 1
               penalpha(ncp) = penalpha(i)
               pencore(ncp) = pencore(i)
            end if
         end do
      end if
      return
      end
