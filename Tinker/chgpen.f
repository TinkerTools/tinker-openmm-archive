c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  chgpen.f  --  specifics of charge penetration variables   ##
c     ##                                                            ##
c     ################################################################

c
c     penalpha        gordon or piquemal penetration parameter
c     pencore         effective nuclear charge 
c
      module chgpen
      use sizes
      implicit none
      real*8, allocatable :: penalpha(:) 
      real*8, allocatable :: pencore(:) 
      save
      end
