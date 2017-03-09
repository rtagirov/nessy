      module MOD_SABOLT
      contains
      !*** Calculate compression Factor depending on TYP
      !*** x must be in [0,1].
      !*** TYP must be one of JUMP, LINE, SPLI.
      !*** TYP = JUMP: Compression function = 1, i.e. full dilution
      !*** TYP = LINE: Linear compression from none to full dilution
      !*** TYP = SPLI(NE): spline compression s.t. cmpr(1)=1, cmpr(0)=0, cmpr'(1) = cmpr'(0)=0 ( _/" shape )
      function CMPRS(TYP,x)
      implicit none
      real*8 :: CMPRS
      character*4,intent(in) :: TYP
      real*8,intent(in) :: x
      selectcase(TYP(1:1))
        case('J')
          CMPRS=1  ! full dilution
        case('L')
          CMPRS=x
        case('S')
          CMPRS=3.*(x**2-2./3.*x**3)
      end select
      
      end function
      ! Calculate NPOT, TEFF, DIL and Tion_pot
      ! SAha BOLtzmann equation
      SUBROUTINE SABOLT (ENTOT,RNE,POPNUM,T,ND,
     $ N,NCHARG,WEIGHT,ELEVEL,EION,
     $ KODAT,NOM,MAXATOM,XTOT,HTOT,GTOT,
     $ RADIUS,MODHEAD,JOBNUM,
     $ Npot,Tion_pot,dil,teff,iTionsel)

      use MOD_DECF_SYN,only:CARDS
      use MOD_ERROR

      use auxfioss

      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(  out) :: npot
      real*8, intent(  out) :: TEFF,dil(nd),Tion_pot(ND,3)
      real*8, intent(in   ) :: ENTOT,RNE, POPNUM, T,WEIGHT,ELEVEL,EION
      real*8, intent(in   ) :: XTOT,HTOT,GTOT,RADIUS
      integer,intent(in   ) :: ND,N,NCHARG,KODAT,NOM,MAXATOM
      integer,intent(in   ) :: JOBNUM,iTionsel
      character,intent(in ) :: modhead*104
      character :: head*65, flname*12

      DIMENSION NCHARG(N)       ! The electron charge
      DIMENSION WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION nom(n)
      DIMENSION KODAT(MAXATOM)
      DIMENSION ENTOT(ND),RNE(ND),T(ND),POPNUM(ND,N)
     $ ,XTOT(ND),HTOT(ND),GTOT(ND)
     $ ,radius(nd)
      LOGICAL WRITEOUT_TION       !ms Writeout the Tion_{1,2,3}.mdisp files ?
      LOGICAL WRITEOUT_DILFAC     !ms Write out the Dilfac.mdisp file ?
C
c  change the "3" if the number of ionization ratios changes
      integer,parameter :: ionequ=3
      DIMENSION TION(3),Zeta(3)
      integer,dimension(ionequ) :: nreflw,nrefup,iact,nrepl

ccc      data ionequ,nreflw,nrefup/3,1,18,29,18,28,39/
      data zeta/-0.340,-0.44,-0.320/
      data sigr,pi/5.67d-5,3.14159/


      WRITEOUT_TION=.false.       ! No writeout
      WRITEOUT_DILFAC=.true.     ! No Writeout
c	do i=1,3 
c	print *,'SABOLT 1: zeta',zeta(i)
c	enddo

      ihelium = kodat(1)           ! KODAT(n) Number of elements/atom n
      ihydro  = kodat(2)           !n=1: helium, n=2: hydrogen, more see datom
      if (ihelium.gt.0) print *,' helium   included',ihelium
      if (ihydro .gt.0) print *,' hydrogen included',ihydro

      nreflw(1:ionequ)=0
      nrefup(1:ionequ)=0
      iact(1:ionequ)=0
      do i=1,n
        if (nom(i).eq.ihelium) then
          if (ncharg(i).eq.0) then
            !*** lower level He0/He+:
            if (nreflw(1).eq.0) then
              nreflw(1)=i
            endif
          elseif (ncharg(i).eq.1) then
              !***  upper level He0/He+ and lower level He+/He++
            if (nrefup(1).eq.0) nrefup(1)=i
            if (nreflw(2).eq.0) nreflw(2)=i
          else if (ncharg(i).eq.2) then
            if (nrefup(2).eq.0) nrefup(2)=i
          else
            print *,' unrecognized He stage',ncharg(i)
            stop ' SABOLT - level search He'
          endif
        else if (nom(i).eq.ihydro) then
          if (ncharg(i).eq.0) then
            !***   lower level H0/H+:
            if (nreflw(3).eq.0) nreflw(3)=i
          elseif (ncharg(i).eq.1) then
            !***     upper level H0/H+
            if (nrefup(3).eq.0) nrefup(3)=i
          else if (ncharg(i).ne.-1) then
              print *,' unrecognized H stage',ncharg(i)
              stop ' SABOLT - level search H'
          endif
        endif
      enddo
      npot=0
      do k=1,ionequ
        if (nreflw(k).ne.0 .and. nrefup(k).ne.0) then
          print '(A,2i3)',' ratio of levels ',nreflw(k),nrefup(k)
          iact(k)=1
          npot=npot+1
        else
          print *,' no levels for ratio number ',k
          iact(k)=0
        endif
      enddo
      !***  replace non-existing rations with another one
      !**   or replace all with one as selected 
      do k=1,ionequ
         if (iTionsel.eq.0) then
            if (iact(k).eq.0) then
               if (iact(2).eq.1) then
                  nrepl(k)=2
                  print *,' ratio',k,' replaced by 2'
               else if (iact(3).eq.1) then
                  nrepl(k)=3
                  print *,' ratio',k,' replaced by 3'
               else if (iact(1).eq.1) then
                  nrepl(k)=1
                  print *,' ratio',k,' replaced by 1'
               else
                  print *,' something wrong with tion replacement'
                  stop ' error SABOLT'
               endif
            else if (iact(k).eq.1) then
               nrepl(k)=k
            else
               print *,' iact has unknown value'
               stop 'error SABOLT'
            endif
         else if (iTionsel.le.3 .and. iTionsel.ge.1) then
            nrepl(k)=iTionsel
         else
            if (iTionsel.ne.-1) stop ' error SABOLT'
            print *,' Electron Temp selected for Tion'
            nrepl(k)=iTionsel
         endif
      enddo
      if (npot.lt.ionequ) npot=ionequ

      head=modhead(15:22)//modhead(35:91)

      print * ,' i   nreflw   nrefup   eion(nreflw)   zeta '
c***  check ground state:
      do i=1,ionequ
	print *,i,nreflw(i),nrefup(i),eion(nreflw(i)),zeta(i)
         if (iact(i).eq.1) then
         if (elevel(nreflw(i)).ne.0.) stop ' nreflw not ground stat'
         print *,i,nreflw(i),nrefup(i),eion(nreflw(i)),zeta(i)
	   zeta(i)=10.**zeta(i)
c	   print *,i,' SABOLT new zeta'
c         print *,i,nreflw(i),nrefup(i),eion(nreflw(i)),zeta(i)
         endif
      enddo

      !*** Calculate T_eff = (L/sigma)^1/4 *****************************
      teff4 = htot(1) * (4.d0*pi*radius(1)*radius(1)/sigr)
      if (teff4.le.0.) then
        print *,'sabolt: teff not defined'
        print '("sabolt: htot(1)   = ",f12.0)',htot(1)
        print '("sabolt: radius(1) = ",f12.0)',radius(1)
        print '("sabolt: sigr      = ",f12.0)',sigr

        call error(' Teff is not defined')
      endif
      teff  = teff4**(1./4.)
      print *,' Teff = ',teff

      !print *,'  not printed:  L,  log(entot),  W,      Te,   3 x Tion'
      !write (10,*) '   L,    log(entot),        Te,           3xTion'
      !*** Calculate the Dilution arry *********************************
      dil(1:ND)=1.
      NDTMIN=getLocalTMin(T,ND)
c      open (10,file='temp.cols',status='unknown')
c      open (11,file='pops.cols',status='unknown')
CMH	CACLULATION OF DILUTION ARRAY
CMH	CHECK WHETHER IT IS SET TO ONE IN LINOP !!
      DO 1515 L=ND,1,-1
        !*** Calculate diffusion factor
        w=1.
        
        if (xtot(l).gt.0.) w=1.d0-htot(l)/xtot(l)
        !!*** compress w
        dil(l) = w
        if(CARDS%COMPRESS /= '') then
          if(L<NDTMIN+CARDS%CMPRS_TMIN_OFFSET) then ! above temperature minimum
            if(L>=NDTMIN+CARDS%CMPRS_TMIN_OFFSET-CARDS%CMPRS_LENGTH)then
              dil(l)=w+(1.-w)*
     &           (1.-CMPRS(CARDS%COMPRESS,
     &                    dble(NDTMIN+CARDS%CMPRS_TMIN_OFFSET-L)
     &                   /dble(CARDS%CMPRS_LENGTH))
     &           )
            print *,l,dil(l)
            endif
          else
            dil(l)=1.
            w=1.
          endif
        endif
        if (l.eq.nd) then
          do k=1,3
            tion(k)=t(l)
            tion_pot(l,k)=tion(k)
          enddo
          tr=t(l)
        elseif (l.eq.1) then
          do k=1,3
            !*** last points are usually corrupted
            !**   use next inner point
            if (nrepl(k).ge.1 .and. nrepl(k).le.3) then
              tion_pot(l,k)=tion(nrepl(k))
            elseif (nrepl(k).eq.-1) then
              tion_pot(l,k)=t(l)
            else
              print *,' nrepl with wrong option'
              stop ' error SABOLT'
            endif
          enddo
          tr=t(l)
        else
        !*** CALCUATE A IONIZATION TEMPERATURE FOR A GIVEN ION RATIO 
        !**                                             (NREFLW- NREFUP)
        !**  ITERATE FOR THE IONIZATION TEMPERATURE
          DO 1616 i=1,ionequ
            IACT_EQ_1: if (iact(i).eq.1) then
              POPUP  = POPNUM(L,NREFUP(i))
              gup   = WEIGHT(NREFUP(i))
              POPLW  = POPNUM(L,nreflw(i))
              GLOW = WEIGHT(nreflw(i))
              RATIO = POPLW/POPUP
              DELERY = EION(nreflw(i))
              ! print *,i,nrefup(i),gup

              !-saha      FAC = RNE(L)*ENTOT(L)*2.07E-16*GLOW/gup
              !-pr               gratio=glow/gup
              !               print *,gratio
              !*** hier ist was ganz komisch wwwww
              FAC = RNE(L)*ENTOT(L)*2.07E-16*Glow/gup
              !-pr print *,' fac ' ,l, fac,rne(l),entot(l),glow,gup,gratio
              !-skip Zeta...     &              /(W* (zeta(i) + (1.-zeta(i))*w) )
     &              / W 
              !** FORUMLATION IN RSRLBV:
              !**        X3=CSAHA* (Z + ZXI(I)*(1.-Z)*WFAC(I) )
              RATIO = RATIO/FAC
              EXPON = DELERY*1.438/TION(i)
              IF (EXPON.LT.2.0) THEN
                TIONOL = (EXP(EXPON)/RATIO)**(2./3.)
              ELSE
                TIONOL = DELERY*1.438/LOG(RATIO*TION(i)*SQRT(TION(i)))
              ENDIF
              ITER=0
 111          ITER=ITER+1
              IF (TIONOL.LT.1000.) TIONOL=1000.

              RANEW = EXP(1.438*DELERY/TIONOL)/TIONOL/SQRT(T(l))
              DELRAT = RANEW-RATIO

              DRADT  = RANEW/TIONOL*(1.438*DELERY/TIONOL+1.d0)
              TIONOL = TIONOL+DELRAT/DRADT

              IF (ABS(DELRAT/DRADT).GT.1. .AND. ITER.LT.40) GOTO 111
              if (iter.ge.40) then
                print *,' ionisation iteration did not converge'
 !             ATTENTION!!!  
 !             stop 'error'
                endif
              IF (ITER.GT.30)PRINT *,' L, ITER, TION: ', L, ITER,TIONOL
              IF (TIONOL.LT.1000.) TIONOL=1000.
              TION(i) = TIONOL
            else IACT_EQ_1
              TION(i) = 0.
            endif IACT_EQ_1
 1616     enddo
          do i=1,ionequ
            if (nrepl(i).ge.1 .and. nrepl(i).le.3) then
              tion_pot(l,i)=tion(nrepl(i))
            else if (nrepl(i).eq.-1) then
              tion_pot(l,i)=t(l)
            else
              print *,' nrepl with wrong option'
              stop ' error SABOLT'
            endif
cc               tion_pot(l,i)=tion(nrepl(i))
          enddo

        !*** endif from if eq.nd else ...
        endif

        !c-pr PRINT 1517,L,LOG10(ENTOT(L)),W,T(L),(TION(I),I=1,ionequ),tr
        ! 1517    FORMAT(I5,2F10.3,10F8.0)

 1515 enddo

      if(WRITEOUT_TION) then
        do j=1,3
          write (flname,'(5HTion_,I1,6H.mdisp)') j
          open (10,file=flname,status='unknown')
          write (10,'(A,A65,I5)') '#TITLE =',head,jobnum
          write (10,'(A)'       ) '#XTITLE = log Density [cm !U-3!N]'
          write (10,'(A)'       ) '#YTITLE = log (population/density)'
          do i=1,nd
            write (10,*) log10(entot(i)),tion_pot(i,j)/1000.
          enddo
          close (10)
        enddo
      endif
      if(WRITEOUT_DILFAC) then
        flname='Dilfac.mdisp'
        open (10,file=flname,status='unknown')
        write (10,'(A,A65,I5)') '#TITLE =',head,jobnum
        write (10,'(A)'       ) '#XTITLE = log Density [cm !U-3!N]'
        write (10,'(A)'       ) '#YTITLE = log (population/density)'
        do i=1,nd
          write (10,*) log10(entot(i)),dil(i)
        enddo
        close (10)
      endif

      RETURN
      END subroutine
      end module
