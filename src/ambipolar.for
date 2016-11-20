      module MOD_AMBIPOLAR
      real*8,parameter :: ak=1.38062259d-16  ! Boltzmann constant
      !*** VA Must be private for the final version
      real*8,allocatable,dimension(:) :: VA  ! Ambipolar Velocity I3.3, I4.7

      contains
      subroutine writeData(fname,data)
        character*(*) :: fname
        real*8 :: data(:)
        ifl=1012
        print *,'write file "dd-'//fname//'"'
        open(ifl,file='dd-'//fname,status='replace')
        write(ifl,'(e15.6)') data
        close(ifl)
      end subroutine
      
      pure function gradient(X,Y)
        real*8,dimension(:),intent(in) :: X,Y
        real*8,dimension(size(X)) :: gradient
        integer :: N
        
        N=size(X);
        if(N /= size(Y)) then 
          gradient=0./0.
        elseif (N==1) then
          gradient=0
        else
          gradient(1) = (Y(2)-Y(1))/(X(2)-X(1))
          gradient(2:N-1) = (Y(3:N)-Y(1:N-2))/(X(3:N)-X(1:N-2))
          gradient(N) = (Y(N)-Y(N-1))/(X(N)-X(N-1))
        endif
      end function
      !*********************************************************************!
      !* Calculate the ambipolar diffusion velocity                        *!
      !* See I: Fontenla, Avrett & Loeser "Energy Balance in the Solar     *!
      !* Transition Region. I. Hydrostatic Thermal Models with Ambipolar   *!
      !* Diffusion, Astro Jour. 355 p.700-718, 1-June-1990                 *!
      !*                                                                   *!
      !**********************************************************************!
      subroutine AMBIPOLAR(ND,N,T,ENTOT,RNE,LEVEL,RADIUS,POPNUM,RSTAR,timer)
      use UTILS,       only:assert 
      use MOD_GRADIFF, only:gradiff
      use MOD_TICTOC , only:writeToc
      implicit none;
      integer,intent(in) :: ND,N,  timer
      real*8,intent(in) :: RSTAR
      real*8,intent(in),dimension(ND)  :: T,ENTOT,RNE,RADIUS
      real*8,intent(in),dimension(ND,N):: POPNUM
      character,intent(in) :: LEVEL(:)*10 !* Name of Levels needed as a safeguard
      real*8,dimension(ND) :: XIO !* degree of ionization n_p/n_1 (I3.2)
      real*8,dimension(ND) :: press !* pressure due to hydrogen atoms, protons and electrons corresponding to this protons
      real*8,dimension(ND) :: Dx,DT !* see I3.1, I3.2
      integer,parameter :: LEVEL_H1 = 2   !** Index of the H1 Level
      integer,parameter :: LEVEL_H2 = 12  !** Index of the H2 Level
      real*8 :: DIFF(ND)
      
      integer :: ifl;
      
      if(.not.allocated(VA)) allocate(VA(ND))
      call assert(size(LEVEL)>=LEVEL_H2,'ambipolar:Level to small.')
      call assert(LEVEL(LEVEL_H1)=='H I......1','ambipolar:H I.')
      call assert(LEVEL(LEVEL_H2)=='H II......','ambipolar: H II')
      
      press= ENTOT*(1.D0+RNE)*ak*T
      Dx = 90.7*T**(1.76)/press  ! p 704, I3.1


      XIO  = POPNUM(:,LEVEL_H2)/POPNUM(:,LEVEL_H1)
      DT=64.1d0*(T**(1.76d0))/press*
     &      (XIO+2.57d0-4d3/(T*sqrt(XIO)))/(XIO+2d-2/XIO)  ! p. 704, I 3.2
      !*** Calculate VA, the ambipolar velocity = Dx*d(ln(XIO))/dz+DT*d(ln(T))/dz
      !call GRADIFF(ND,log(XIO),DIFF,RADIUS)
      VA = Dx*gradient(RADIUS,log(XIO))+DT*gradient(RADIUS,log(T))
      VA=VA/RSTAR  ! VA in cm/s

      !*** DEBUG from here on***
!      ifl=1012
!      open(ifl,file='VELO.AMBI',status='replace')
!      write(ifl,'(e12.3)') VA
!      close(ifl)
      print *,'AMBIV time: '//writeTOC(timer)
      print '("AMBIV: " 4e12.4e3)',VA
      end subroutine
      end module
