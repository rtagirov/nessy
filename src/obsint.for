      module MOD_OBSINT

      contains

      SUBROUTINE OBSINT(LTOT,CORE,BCORE,DBDR,ip,
     $                  IRIND,RRAY,ZRAY,XCMF,
     $                  ND,NP,NVD,JP,NVOPA,VOPA0,DVOPA,
     $                  EMINT,XOBS0,DXOBS,NFOBR,XN,
     $                  ENTOT,RNE,SIGMAE,RSTAR,
     $                  XJ,DINT,RWLAE,DLAM)

!***  FORMAL INTEGRATION ALONG A RAY AT IMPACT PARAMETER JP
! ==> test version: see for "cpr" for opcity prints...

      use MOD_CUBIC
      use MOD_DECF_SYN
      use OPINT

      implicit real*8(a-h,o-z)
      real*8 :: RWLAE, DLAM(NFOBR)
      real*8 :: DINT(NVD, 2 * ND), XJ(ND, NVD)
      DIMENSION RRAY(LTOT),ZRAY(LTOT),XCMF(LTOT)
      DIMENSION EMINT(NFOBR)
      DIMENSION ENTOT(ND), RNE(ND)

      real*8, dimension(NP), intent(in) :: ip

      INTEGER   IRIND(LTOT)
      real*8 :: conf
      real*8,allocatable :: TAUSUM(:,:)  ! needed to print out the TAU=1 values
      real*8,allocatable :: ZFINE_SUM(:)  ! needed to print out the TAU=1 values
      integer ii
      real*8 ::mu_c, mu_l
      real*8,allocatable :: tau_r(:)
      real*8,allocatable :: contr_c(:)
      real*8,allocatable :: contr_l(:)
   


C***  THE FOLLOWING VECTORS ARE LOCAL IN ORDER TO SUPPORT AUTO-VECTORIZATION
cccc      DIMENSION 	EQUIVALENCE (TAU(1,1),OPAFINE(1,1))
      DIMENSION QIWI(NMAX),PIWI(NMAX)
      LOGICAL CORE

      IF (LTOT.LE.1) STOP ' OBSINT-1'

!      if(cards.PRINT_TAU) then

        allocate(TAUSUM(LTOT,NFOBR))
        allocate(ZFINE_SUM(LTOT))
        allocate(tau_r(NFOBR))
        allocate(contr_c(ND))
        allocate(contr_l(ND))

        ZFINE_SUM(:)=0.

!      endif

      !***  START AT FAR END OF GRID (from LTOT to 1)
      !***  (integration is from near to far within one grid step)
      LOOP_L: DO L=LTOT-1,1,-1

        !*** L  is the near side grid point
        !*** L1 is the far side point
        L1=L+1

        !***  DIFF IN VDOP UNITS  (0 for static atmospheres)
        DELX=XCMF(L1)-XCMF(L)

        !*** first integartion away from the core XCMF(L1) is still positive
        IF (CORE .AND. (L1 .EQ. ND)) DELX=-XCMF(L1)-XCMF(L)

        !***  N = NUMBER OF REFINED MESH POINTS - DELX POINTS PER VDOP
        N = (INT(DELX*XN)+1) +1   ! In the FIOSS-SUN this is always 2 because of static atmosphere

        !***  avoid detailed integration through the core
        IF (CORE .AND. (L1 .EQ. ND)) N=2
        DX=DELX/FLOAT(N-1)

        IF (N .GT. NMAX) THEN
          print *,' N, delx, dx: ',N, delx, dx
          STOP 'DIMENSION NMAX IS INSUFFICIENT - OBSINT'
        ENDIF

        if (n.gt.2) then
          !***        N>2 BRANCH ==> implies interpolation between grid points
          !
          !***  CUBIC INTERPOLATION OF Z AS FUNCTION OF XCMF
          !***  CALCULATE THE COEFFICIENTS P1 ... P4 FOR THE INTERVAL L,L1
          CALL CUBIC(L1,LTOT,ZRAY,XCMF,P1,P2,P3,P4)
          !***  find grid-z and correction
          zl  = sqrt((rray(l) -  ip(jp)) * (rray(l)  + ip(jp)))
          zll = sqrt((rray(ll) - ip(jp)) * (rray(ll) + ip(jp)))

          dizl  = zl  - abs(zray(l))
          dizl1 = zl1 - abs(zray(l1))

          !***  INTEGRATE FROM NEAR TO FAR
          LOOP_N:DO I=1,N
            XI=FLOAT(I-1)*DX+XCMF(L)
            XFINE(I)=XI
            DXA=XI-XCMF(L)
            DXB=XCMF(L1)-XI
            DXA3=DXA*DXA*DXA
            DXB3=DXB*DXB*DXB
            !***  CUBIC INTERPOLATION OF ZFINE
            ZFINE(I)=P1*DXA3+P2*DXA+P3*DXB3+P4*DXB
            !***  find the r value without using rr = zz+pp
            p=(zFINE(i)-zRAY(L1))/(zRAY(L)-zRAY(L1))
            q=1.-p
            diffz = Q*dizl1+P*dizl

            gridz = diffz+abs(zfine(I))

            rfine = SQRT(ip(jp) * ip(jp) + gridz * gridz)

            !***  now ... INTERPOLATED LINEARLY IN RADIUS R
            P=(RFINE-RRAY(L))/(RRAY(L1)-RRAY(L))
            Q=1.-P
            !***  INTERPOLATION WEIGHTS MODIFIED: LINEAR INTERPOLATION IN (R**2 *QUANTITY)
            ! **** correct??
            !-test this here might be really wrong...??!!!
            Q=Q*RRAY(L )*RRAY(L )/(RFINE*RFINE)
            P=P*RRAY(L1)*RRAY(L1)/(RFINE*RFINE)
            QIWI(I)=Q
            PIWI(I)=P
          ENDDO LOOP_N
        else
          !***  N=2 BRANCH
          QIWI(1)=1.
          QIWI(2)=0.
          PIWI(1)=0.
          PIWI(2)=1.
          XFINE(1)=Xcmf(l)
          XFINE(2)=Xcmf(l1)
          IF (CORE .AND. (L1 .EQ. ND)) XFINE(2)=-XCMF(L1)
          ZFINE(1)=zray(l)
          ZFINE(2)=zray(l1)

          !***  outward integration away from the core
          IF (CORE .AND. (L1 .EQ. ND)) ZFINE(2)=abs(zray(l1))
        endif

        if(cards.PRINT_TAU) ZFINE_SUM(L)=ZFINE(1)-ZFINE(N)+ZFINE_SUM(L+1)
        !***  loop over the number of inter-mesh grid points

        DO 3 I=1,N
          !***  loop over all observers-frame frequency points
          DO 13 K=1,NFOBR
            !***  XO is the observer's frame frequency in DOP units
            !***  K=1 is the bluemost but XOBS0 is a positive number
            !***  (dxobs is negative)
            !*** DOP units are positve for blue velocities!!!
            !*** (dxobs is a negative number and XOBS0 is positive)
            !*** K=1 is blue most frequency
            !*** XOBS0 is for index K=0
            XO=XOBS0+K*DXOBS
            !*** XFINE (= interpolated Xcmf) the doppler shift due to
            !*** the line of sight velocity
            !*** in DOP units (negative on near side)

            !***  XI is the comoving frame frequency that contributes to
            !*** XO at the given location in DOP units
            !*** note: on the front side, the element is moving AWAY from the
            !*** incident photon, thus the element sees the photon red-shifted,
            !*** Use a frequency red-shifted by V to absorb and emit at XO,
            !*** thus subtract the XFINE
            !*** HOWEVER:
            !***   near side: blue velocities ==> negative XCMF
            !*** THUS: it is correct to add XFINE
            XI=XO+XFINE(I)
            !***  VOPA0 is the frequency of the first element of the
            !***        opacity grid (this is in contrast to the
            !***        definition of XOBS0 which is the 0th)
            !***  KOPA  is the (nearest) index of the opacity grid
            !***        corresponding the XI get nearest index
            !***        ==> it seems that it is (sometimes?) important
            !***            to interpolate ..

            !*** version with interpolation
            nopa = (xi-vopa0)/dvopa + 1
            !test if (nopa.lt.1) stop 'nopa.lt.1'
            kopa = nopa+1
            !test if (kopa.gt.nvopa) stop ' kopa.gt.nvopa'
            !*** the stop must be commented out to allow vectorizing
            P = (XI-VOPA(NOPA))/(VOPA(KOPA)-VOPA(NOPA))
            Q = 1.-P
            OPATL  = Q*OPATOT(NOPA,IRIND(L))  + P*OPATOT(KOPA,IRIND(L))
            OPATL1 = Q*OPATOT(NOPA,IRIND(L1)) + P*OPATOT(KOPA,IRIND(L1))
            ETATL  = Q*ETATOT(NOPA,IRIND(L))  + P*ETATOT(KOPA,IRIND(L))
            ETATL1 = Q*ETATOT(NOPA,IRIND(L1)) + P*ETATOT(KOPA,IRIND(L1))

            !***pr
            ! if (abs(xi).lt.2.)
            !   $ print '(4I5,2F10.2,2e10.2)',L,L1,IRIND(L),IRIND(L1),
            !   $     XO,XI,opatl,opatl1
            !***  version with nearest index only
            !$$$ kopa = NINT((xi-vopa0)/dvopa + 1)
            !$$$ OPATL  = OPATOT(KOPA,IRIND(L ))
            !$$$ OPATL1 = OPATOT(KOPA,IRIND(L1))
            !$$$ ETATL  = ETATOT(KOPA,IRIND(L ))
            !$$$ ETATL1 = ETATOT(KOPA,IRIND(L1))
            IF (KOPA.LT.1 .OR. KOPA.GT.NVOPA) THEN
              PRINT *,' INDEX KOPA= ',KOPA,' OUT OF BOUND',KOPA,NVOPA
              stop ' ERROR OBSINT'
            ENDIF

            OPAFINE(K,I) = QIWI(I)*OPATL + PIWI(I)*OPATL1

            !*** Add the Thomson Emissivity to the Line Intensity
            !***
            !*** K is the observer's frame frequency index of the
            !*** resulting radiation KAPP is the apparent frequency
            !*** index seen by the scattering electron
            !*** note: on the front side, the electron is moving AWAY
            !***    from the incident photon, thus use the frequency
            !***    red-shifted by V to add to the Thomson Emissivity
            !*** XJ is the co-moving frame mean Intensity convolved with
            !***     the electron scattering profile
            kapp   = NINT((xi-vopa0)/dvopa + 1)
            etaelsc= entot(IRIND(L))*rne(IRIND(L))*SIGMAE*
     $               XJ(IRIND(L),KAPP)*rstar
            etatl  = etatl + etaelsc
            etaelsc= entot(IRIND(L1))*rne(IRIND(L1))*SIGMAE*
     $               XJ(IRIND(L1),KAPP)*rstar
            etatl1= etatl1+ etaelsc
            ETA   = QIWI(I)*ETATL + PIWI(I)*ETATL1
              !***  CALCULATE SOURCE FUNCTION
            SFINE(K,I)=ETA/OPAFINE(K,I)

 13       ENDDO
 3      ENDDO

        !***  ESTABLISH DTAU(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
        DO I=2,N
          DELTAZ=(ZFINE(I-1)-ZFINE(I))
          DO K=1,NFOBR
            DTAU(K,I)=0.5*(OPAFINE(K,I-1)+OPAFINE(K,I))*DELTAZ
            
          ENDDO
        ENDDO

        !***  ESTABLISH TAU(I) = OPTICAL DEPTH SCALE
        !***  note: i=1 is the near side (=L) and i=n the far (=L1)
        !*** Warning: TAU is the DTAU over the depth points L
        TAU(1:NFOBR,1)=0.
        DO I = 2, N

          TAU(1:NFOBR,I)=TAU(1:NFOBR,I-1)+DTAU(1:NFOBR,I)

        ENDDO

        if(CARDS.PRINT_TAU) TAUSUM(l,1:NFOBR)=TAU(1:NFOBR,N)

        if ((JP .eq. 1) .and. (L .le. ND-1)) then   ! do the job only for the central ray and for the near-side
        mu_c=1.  ! disc center
        mu_l=0.31225  ! R/R_sun=0.95
        tau_r(1:NFOBR)=0.    ! real optical depth

        do ii=LTOT-1, LTOT-l+1, -1  ! the far-side levels are already calculated and stored in tausum, so we use them here
        tau_r(1:NFOBR)=tau_r(1:NFOBR)+TAUSUM(ii, 1:NFOBR)   
        enddo

        contr_c(l)=0.
        contr_l(l)=0.

        do ii=1, NFOBR

           contr_c(l)=contr_c(l)+SFINE(ii,1)*OPAFINE(ii,1)*exp(-tau_r(ii)/mu_c)/mu_c
           contr_l(l)=contr_l(l)+SFINE(ii,1)*OPAFINE(ii,1)*exp(-tau_r(ii)/mu_l)/mu_l 

        enddo

        contr_c(l)=contr_c(l)/NFOBR
        contr_l(l)=contr_l(l)/NFOBR

        if (contr_c(l) .lt. 1.E-30) contr_c(l)=0.

        if (contr_l(l) .lt. 1.E-30) contr_l(l)=0.

        write(250,280), L, contr_c(l), contr_l(l) 

        endif

 280    format(I3, 2x, E12.4, E12.4)

        DO 53 K=1,NFOBR
          IF (TAU(K,N) .GT. 1.E-10) THEN
            EXPTAU=1.
            WTAU(K,1)=EXPTAU+ (EXP(-TAU(K,2))-EXPTAU)/DTAU(K,2)
          ENDIF
 53     ENDDO

        IF (N.GT.2) THEN
          DO 6 I=2,N-1
            DO 63 K=1,NFOBR
              IF (TAU(K,N) .GT. 1.E-10) THEN
                EXPTAU=EXP(-TAU(K,I))
                !***  WB FOR THE INTERVAL I-1, I
                WB=(EXP(-TAU(K,I-1))-EXPTAU)/DTAU(K,I)
                !***  WA FOR THE INTERVAL I, I+1
                WA=(EXP(-TAU(K,I+1))-EXPTAU)/DTAU(K,I+1)
                WTAU(K,I)=WA+WB
              ENDIF
63          ENDDO
 6        ENDDO
        ENDIF

        DO K=1,NFOBR
          IF (TAU(K,N) .GT. 1.E-10) THEN
            EXPTAU=EXP(-TAU(K,N))
            WTAU(K,N)=-EXPTAU-(EXPTAU-EXP(-TAU(K,N-1)))/DTAU(K,N)
          ENDIF
        ENDDO

        IF (CORE .AND. (L .EQ. ND)) THEN
          do 40 k=1,nfobr
            XO=XOBS0+K*DXOBS
            XI=XO+XCMF(ND)
            kopa = NINT((xi-vopa0)/dvopa + 1)
            XOPA = opaTOT(KOPA,IRIND(ND))
            PLUSI=BCORE+DBDR*ZRAY(ND)/XOPA
            EMINT(K)=PLUSI


            DINT(Kopa,L)=EMINT(K)
40        enddo
        ELSE
          !*** INTEGRATION SUM
          DO 92 K=1,NFOBR
            EXPTSC=EXP(-TAU(K,N))
            EMINT(K)=EMINT(K)*EXPTSC
92        ENDDO
          DO 93 I=1,N
            DO 83 K=1,NFOBR
              IF (TAU(K,N) .GT. 1.E-10) THEN

                EMINT(K)=EMINT(K)+SFINE(K,I)*WTAU(K,I)

              ENDIF
 83         ENDDO
 93       ENDDO

          !*** the formal integration is finished
          !*** now store the intensity in the array DINT
          !*** DINT is the co-moving frame intensity

          DO KOPA=1,NVOPA
            XI=vopa0+(kopa-1)*dvopa
            XO=XI-Xcmf(L)
            K=nint((XO-Xobs0)/dXobs)
            K=max(1,K)
            K=min(NFOBR,K)
            DINT(Kopa,L)=EMINT(K)
          ENDDO
        ENDIF

      ENDDO LOOP_L

      !***  END OF LOOP OVER GRID POINTS
      !*** sum over the depthpoints and when > 1 print out

      IF(CARDS.PRINT_TAU) call PRINTTAU(NFOBR,TAUSUM,ZFINE_SUM,ND,RSTAR,RWLAE,DLAM,JP)
      if(allocated(ZFINE_SUM)) deallocate(ZFINE_SUM)
      if(allocated(TAUSUM)   ) deallocate(TAUSUM)

      RETURN

      END subroutine

      subroutine PRINTTAU(NFOBR,TAUSUM,ZFINE_SUM,ND,RSTAR,RWLAE,DLAM,JP)

      use interpolation
      use utils

      implicit none

      integer,intent(in):: NFOBR,ND,JP  !# of Freq/Depth points, Impact Parameter
      real*8, intent(in):: TAUSUM(:,:),RSTAR ! L,K
      real*8, intent(in):: RWLAE,DLAM(:)     ! central frequency, delta to RWLAE
      real*8, intent(in):: ZFINE_SUM(:) ! fine radial grid(LTOT)
      real*8 :: Z(ND)           ! Fine corrected grid
      real*8 :: TAUSSUM(ND)     ! accumulated tau
      real*8 :: LAMBDA,y,dy     ! Freq point, tau=1,tau=1 error
      real*8 :: coeff
      real*8 :: FH(NFOBR), FHA(40)
      real*8 :: ZZ1,ZZ2, tau1, tau2
      integer :: K,L, iFROM,iTO ! loop variables, index
      character*(*),parameter :: FORMAT = '(3(e18.10e4,x),i1,x,i3)'
      if(JP/=1) return
      if(size(ZFINE_SUM)<ND) then
c*** call PRINTTAU only for the central ray
        print *,'obsint: printTAU: warning: ZFINE to small' 
        return
      endif


      !*** 1) calculate the depth-grid.
      !***  !!! This only works for the central ray !!!
      Z(1:ND)=(ZFINE_SUM(1:ND)-ZFINE_SUM(ND))*RSTAR/1e5  ! RSTAR in cm to km

      !*** 2) Calculate the tau=1 for every frequency point
      DO K=1,NFOBR
        !*** 1) add up the opacities  taus(l) = \int_depth(l)^\inf tau(l)
        TAUSSUM(1)=0.
        do L=2,ND
          TAUSSUM(L)=TAUSSUM(L-1)+TAUSUM(L,K)
          if(TAUSSUM(L)==TAUSSUM(L-1)) then
            print *,'obsint: warning: TAUSUM(L,K) < eps'
            return
          endif
        enddo
        !*** 2) find the point where TAU switches to 1
        do L=1,ND  
          if(TAUSSUM(L)>1.) exit
        enddo
        ZZ1=Z(L-1)
        ZZ2=Z(L)
        tau1=taussum(L-1)
        tau2=taussum(L)
        FH(K)=ZZ1-(tau1-1.)/(tau1-tau2)*(ZZ1-ZZ2)


         if ((JP .eq. 1) .and. (K .eq. 1000)) then

         open(270, file='atmosinfo.txt',status='unknown') 
         do L=1, ND         

         write(270,*), L, Z(L), TAUSSUM(L)

         enddo
         close(270)
      


         endif


        !*** 3) write the TAU=1 to file
        LAMBDA = RWLAE+DLAM(K) ! The frequency point
        !*** The file format is Frequency, Formation height,
        !***      error-estimate, status (1 is highest depth point,
        !***                              2 lowest, 0 in between)
        !*** i.e. if status==1  then the atmosphere did not go high up
        !*** enough if status == 2 then the atmosphere is transparent
        !*** and otherwise at some point in the atmosphere there is
        !*** absorption.
        IF(L<=2)THEN
          write(9998,FORMAT)  LAMBDA,Z(L),0. , 1,L
        ELSEIF(L==ND)THEN
          write(9998,FORMAT)  LAMBDA,Z(L),0. , 2,L
        ELSE
          
          coeff=((TAUSSUM(L))-1.)/(1.-TAUSSUM(L-1.))
            
   !       iFROM=max(L-3,1);  iTO=min(L+1,ND)
   !       call assert(TAUSSUM(iFROM)/= TAUSSUM(iTo))
   !       call polint(TAUSSUM(iFROM:iTO), Z(iFROM:iTO),1.,y,dy)
           dy=0.
           y=(Z(L)+Z(L-1)*coeff)/(1.+coeff)

         write(9998,FORMAT)  LAMBDA,    y   ,dy , 0,L
        ENDIF
      ENDDO


      open(unit=200, file='../form.height', access='append')
   
      do k=1, 100

      write(200, *), RWLAE+DLAM(k*20-10), sum(FH(20*(k-1)+1:20*k))/20.

      enddo

      close(200)


      end subroutine PRINTTAU
      end module
