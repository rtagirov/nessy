      MODULE MOD_LINOP_MS

      PUBLIC :: LINOP_MS

      CONTAINS
      !*****************************************************************
      ! Calculate the line absorption and emission for every frequency point
      ! (ABLIN, EMLIN) based on the line profile, the depthpoint ID

      ! Called by inibl
      SUBROUTINE LINOP_MS(ID,ABLIN,EMLIN)

      use MOD_SYNSUBM, only:PROFIL,PHE1,PHE2
      use MOD_TICTOC,  only:WRITETOC
      use MOD_ERROR,   only:ERROR
      use MOD_DECF_SYN
      use MOD_chemeq
      use constants

      USE COMMON_BLOCK
      USE PHYS

C     TOTAL LINE OPACITY (ABLIN) AND EMISSIVITY (EMLIN)
      implicit none
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/LINDAT.FOR'
      integer,intent(in) :: ID
      logical :: LPR,ltrad

!RT:  VARIABLES FOR THE CALCULATION OF NLTE LINES

!=======================================================================

!      LOGICAL :: CONDITION

      INTEGER :: NT

      REAL*8 ::  NLTE_ABS, NLTE_EMI

      REAL*8 ::  WAV, FRE, BLU, DCR, IEF, NF1, NF2

!=======================================================================

      real*8,dimension(MFREQ) :: ABLIN,EMLIN, abs_her
      real*8,allocatable,dimension(:) :: ABLINN,XFR, xfr1
      real*8,allocatable :: ABL_A(:) ! ABL Array Version
   
      real*8 ::  DOPA1,VDWC    ! commons PRFQUA
      integer :: nltoff,iemoff ! commons linoff
      real*8 ::  velmax        ! commons

      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
      common/linoff/velmax,nltoff,iemoff

      real*8, PARAMETER::sqrp=1.77245385

      real*8, PARAMETER::DP0=3.33564E-11
      real*8, PARAMETER::CLm=299792458.

      integer :: il,l, ind(1)      ! outer,inner Loop Counter
      integer :: innlt,ISP,ION,IAT ! Helper Variable for array elements

      real*8 :: TEM1,PLANW         ! Helper Variables
      real*8 :: AGAM, agam1, DOP1  ! AGAM=profil, Helper Variable
                                   ! FREQ0 is the threshold frequency from
                                   ! level I of ion KI. Threshold cross-sections will be of the order

      real*8, allocatable :: AB0(:), opacw(:) ! AB0: Absorption at the line center
      integer :: AB0_COUNTER
      real*8 :: SL0                           ! Source function at the line center
      real*8 :: BI,BJ,STIMB,ABL
      real*8 :: WDID, frw, freq1, freqN
      real*8 :: freqb1, freqb2, xb
      real*8 :: opres, opres1, Delta, AB0mean, adc
      real*8 :: Ffudge, freqt, lambdat
      real*8 :: linefreq, lambdaAn
      real*8 :: freqBN1, freqBN2
      real*8 :: xtest(1000)
      real*8 :: opacwm, ABLINmean, ABLINmin
      integer :: iD1, iD2, Nc
      integer :: iL1, iL2, ILmin, ILmax
      integer :: mn, j, ii, i

      real*8, allocatable :: ABLIN_old(:)
      real*8 :: testint

      NF1 = LIGHT_SPEED**(2.0D0) / PLANCK / 2.0D0 ! FIRST NORMALIZATION FACTOR FOR NLTE LINE TREATMENT

      AB0_COUNTER = 0
      print '("linop_ms: NLIN0 = ",i7,",  ID=",i3,'//
     *      '", time=",A," C(AB0) = ",$)',
     *        NLIN0,ID,writeTOC()

      if(MFREQ < NFREQ) call ERROR('LINOP_MS: MFREQ < NFREQ ')

      ABLIN(:NFREQ) = 0.
      EMLIN(:NFREQ) = 0.

      IF (NLIN0 .EQ. 0) THEN

          PRINT '(A)', '-0'

          RETURN ! No lines to compute

      ENDIF

      NDPMIN = getLocalTMin(TEMP,ND)
      allocate(ABLINN(NFREQ))
      allocate(ABL_A(NFREQ))
      allocate(XFR(NFREQ))
      allocate(XFR1(NFREQ))
      allocate(opacw(NFREQ))
      
      ABLINN(:)= 0.
      !******************************************
      !* overall loop over contributing lines
      !------------------------------------------
      ! ATTENTION !!!!!!

      if(CARDS%COMPRESS=='') wdil(id)=1.

      !lvi=Velw(id).gt.velmax.and.iemoff.eq.0                      ! should be always false (Velw < velmax)
      !lne=velw(id).gt.velmax.and.nltoff.gt.0.and.iemoff.gt.0      ! should be always false  -- micha
      ltrad=trad(1,id).eq.trad(2,id).and.trad(1,id).eq.trad(3,id)
      ltrad=ltrad.and.trad(1,id).eq.temp(id)
      planw=xjcon(id)

      !*** Loop over all lines, calculate for every line the Absorption and Emission
       
      freq1=FREQ(1)
      freqN=FREQ(NFREQ)

      CALL calcmolopac(NFREQ, freq, temp(ID), ID)

!     **********   Molecular loop       *************************                             
      DO mn = 1, Nmol

         IF (imax(mn) .GT. imin(mn)) THEN

             DO j = imin(mn), imax(mn)

                linefreq=clight_cgs*Mollines(mn)%Etran(j) 
!               now the trick converting it to "air" :)) frequency

                lambdaAn=(clight_cgs/linefreq)*1.d8

                linefreq=1.d8*clight_cgs/(airlambda(lambdaAn))

                dop1=1./(dp0*vturbmol(mn,ID)*linefreq)
         
                XFR(:)=ABS(FREQ(1:NFREQ)-linefreq)*DOP1

                freqb1=linefreq-4./dop1
                freqb2=linefreq+4./dop1

                iD1 = min(max(int((freq1*freqN/freqb1-freqN)*(Nfreq-1)
     *                /(freq1-freqN)),1),NFREQ)

                iD2 = min(int((freq1*freqN/freqb2-freqN)*(Nfreq-1)/
     *                (freq1-freqN))+1, Nfreq)

                iD2 = max(iD2,iD1)

                IF (iD1 .LE. NFREQ) THEN

                    IF (iD2 .GT. NFREQ) iD2 = NFREQ

                    ABLIN(iD1:iD2)=ABLIN(iD1:iD2)+dop1*alpha0(mn,j)*
     *              VOIGTK_MS(0., XFR(ID1:iD2), iD2-iD1+1)/sqrp

                ENDIF

             ENDDO

         ENDIF

      ENDDO

      ILmin = -1
      opres =  0.
      opres1 = 0.

      PRINT*, 'NLIN0 = ', NLIN0

      ALLOCATE(AB0(NLIN0))

      DO IL = 1, NLIN0

         agam =  profil(il, indat(il) / 100, id)
         INNLT = INDNLT(IL)
         tem1 =  1d0 / trad(ipotl(il), id)

         IAT =  INDAT(IL) / 100
         ION =  MOD(INDAT(IL), 100)
         DOP1 = DOPA1(IAT, ID)

         freqBN1 = CLm*freq1/(CLm+freq1*1.d-10)
   
         freqBN2 = CLm*freqN/(CLm-freqN*1.d-10)

         adc = 0.

         if (abs(freq0(IL)-1.072400126201562E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.069653987843886E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.051119263314454E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.040375995374756E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-7.621203392359742E+014) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-7.554360585908571E+014) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.548147717044358E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.551507183793536E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.809239653640779E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.808025216495057E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.808258593515800E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.808834660026463E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.193976511764430E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.193554397271716E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.191490437830097E+015) .lt. 10. ) adc=1.
         if (abs(freq0(IL)-1.190102205277549E+015) .lt. 10. ) adc=1.

         IF (((FREQ0(IL) .gt. freqBN1) .and. (FREQ0(IL) .lt. freqBN2)) .or. (adc .eq. 1.)) THEN

             IF ((ILmin==-1) .and. (adc .eq. 0.)) ILmin = IL

             IF (adc .eq. 0.) ILmax = IL

             IF (INNLT == 0) THEN

                 IF (excl0(il) < 2000.) THEN ! Attention: Currently WDID is always 1 (see wdil(id) = 1. above)

                     wdid = 1.

                 ELSE

                     wdid = wdil(id)

                 ENDIF
 
                 AB0(IL)=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*DOP1*STIM(ID)*wdid

                 IF (.not. ltrad) sl0 = planw ! = xjcon(id)

             ELSEIF (INNLT > 0) THEN

                 AB0(IL) = ABCENT(INNLT, ID)
                 SL0 =     SLIN(INNLT, ID)
 
             ELSE ! INNLT < 0

                 BI = 1d0
                 BJ = 1d0

                 IF (ILOWN(IL)> 0) BI = BFAC(ILOWN(IL), ID)
                 IF (IUPN(IL) > 0) BJ = BFAC(IUPN(IL), ID)

                 STIMB=BI*(1.-EXHK(ID)*BJ/BI)

                 AB0(IL)=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*DOP1*STIMB

                 SL0=PLAN(ID)*(1d0-EXHK(ID))/(BI/BJ-EXHK(ID))

             ENDIF    

         ENDIF

         IF (adc .eq. 1.) THEN

             agam=profil(il,indat(il)/100,id)

             XFR(:)=ABS(FREQ(1:NFREQ)-FREQ0(IL))*DOP1

             ABLIN(1:NFREQ)=ABLIN(1:NFREQ)+ AB0(IL)*VOIGTK_MS(AGAM,XFR(1:NFREQ),NFREQ)

             ab0(il)=0.

         ENDIF

      ENDDO

      AB0mean=sum(AB0(ILmin:ILmax))/(ILmax-ILmin+1)

      Nc=0

        LINE_LOOP: DO IL=ILmin, ILmax

        agam=profil(il,indat(il)/100,id)
        INNLT=INDNLT(IL)
        tem1=1d0/trad(ipotl(il),id)

        IAT=INDAT(IL)/100
        ION=MOD(INDAT(IL),100)
        DOP1=DOPA1(IAT,ID)




        ISP=ISPRF(IL)
        IF (ISP.GE.6) cycle LINE_LOOP     ! If    ISP > 5 goto end of loop (why the LPR >=5 above?)
        agam=profil(il,indat(il)/100,id)
    
        if(IAT<=0) print *,IAT
        LPR=.NOT.(ISP.GT.1.AND.ISP.LE.5)  ! LPR = ISP not in (1,5] 
        
     

        freqb1=FREQ0(IL)-4./dop1

        freqb2=FREQ0(IL)+4./dop1

        xb=4.



        XFR(:)=ABS(FREQ(1:NFREQ)-FREQ0(IL))*DOP1

        XFR1(:)=(FREQ(1:NFREQ)-FREQ0(IL))*DOP1

        Delta=(XFR1(NFREQ)-XFR1(1))/(NFREQ-1)



        if (AB0(IL) .gt. 10.*AB0mean) then
        Nc=Nc+1

       IF(INNLT.EQ.0.and.ltrad) THEN


               IF(LPR) THEN


            ABLIN(1:NFREQ)=ABLIN(1:NFREQ)+ AB0(IL)*VOIGTK_MS(AGAM,XFR(1:NFREQ),NFREQ)

          ELSE

            do l=1,nfreq
              ABLIN(l)=ABLIN(l)+AB0(IL)*PHE1(ID,FREQ(l),ISP-1)
            enddo
          END IF




          ELSE
          IF_LPR1: if(LPR) then
     
!            ABL_A = AB0(IL)*VOIGTK_MS(AGAM,XFR,NFREQ)
!            ABLINN=ABLINN+ABL_A
!                      EMLIN(1:NFREQ)=EMLIN(1:NFREQ)+ABL_A*SL0
            
            
          else IF_LPR1
           
 

 !             do l=1,nfreq
      

 !             ABL=AB0(IL)*PHE1(ID,FREQ(l),ISP-1)
 !             ABLINN(l)=ABLINN(l)+ABL
            
!              EMLIN(l)=EMLIN(l)+ABL*SL0
            
!            enddo
          endif IF_LPR1     
        ENDIF





        cycle LINE_LOOP
        endif
        



!        if (((freq1 .gt. freqb2) .or. (freqN .lt. freqb1)) .and. 
!     *  (INNLT.EQ.0.and.ltrad)) cycle LINE_LOOP
         

        iD1=max(int((freq1*freqN/freqb1-freqN)*(Nfreq-1)
     *  /(freq1-freqN)),1)
        iD2=min(int((freq1*freqN/freqb2-freqN)*(Nfreq-1)/
     *  (freq1-freqN))+1, Nfreq)
  
         iD2=max(iD2,iD1)
      

 
        if(IAT<=0) print *,IAT
        LPR=.NOT.(ISP.GT.1.AND.ISP.LE.5)  ! LPR = ISP not in (1,5]

 
        if(AB0(IL)==0.) cycle LINE_LOOP !** only do all the work if its worth it

    
        AB0_COUNTER = AB0_COUNTER + 1

        IF (INNLT .EQ. 0 .and. ltrad) THEN

          !*************
          !* LTE lines *
          !*************

!       test - no background

            if (iD1 .le. NFREQ) then

                if (iD2 .gt. NFREQ) iD2 = NFREQ

                ABLIN(iD1 : iD2) = ABLIN(iD1:iD2)+AB0(IL)*VOIGTK_MS(AGAM, XFR(ID1:iD2), iD2-iD1+1)

            endif


        if (AGAM .gt. 0.01) then

        freqb1=FREQ0(IL)-400*agam/dop1

        freqb2=FREQ0(IL)+400.*agam/dop1

        xb=400.*agam

        iL1=int((freq1*freqN/freqb1-freqN)*(Nfreq-1)/(freq1-freqN))+1   
        iL2=int((freq1*freqN/freqb2-freqN)*(Nfreq-1)/(freq1-freqN)) 

        iL1=max(1,iL1)
        iL2=min(iL2,NFREQ)

        if ((iL1 .lt. iD1) .and. (iL2 .gt. iD2)) then

!           test - no background
            ABLIN(iL1:iD1)=ABLIN(iL1:iD1)+AB0(IL)*VOIGTK_MS(AGAM, XFR(IL1:iD1), iD1-iL1+1)
      
            ABLIN(iD2:iL2)=ABLIN(iD2:iL2)+AB0(IL)*VOIGTK_MS(AGAM, XFR(ID2:iL2), iL2-iD2+1)

        endif
     
        endif



!       test - no background


!          opres=opres+2.*AB0(IL)*agam/(Pi*xb*Delta)
 
!        CHECK!!!
         opres=opres+2.*AB0(IL)*agam/(sqrp*xb*Delta)
!        CHECK!!!         
         

    !     if (AB0(IL) .lt. 100.*AB0mean) then
         opres1=opres1+2.*AB0(IL)/(Delta)
    !     endif


        !**************
        !* NLTE LINES *
        !**************
        ELSE
          IF_LPR: if(LPR) then  !{
            ! ABLINN(:)=ABLINN(:)+AB0*VOIGTK_MS(AGAM,XFR(1:NFREQ),NFREQ)
            ABL_A = AB0(IL)*VOIGTK_MS(AGAM,XFR,NFREQ)
            ABLINN=ABLINN+ABL_A
            ! if(.not. lne) then
            EMLIN(1:NFREQ)=EMLIN(1:NFREQ)+ABL_A*SL0
            !endif
            !again, special expressions for 4 selected He I lines
          else IF_LPR
            do l=1,nfreq
              ABL=AB0(IL)*PHE1(ID,FREQ(l),ISP-1)
              ABLINN(l)=ABLINN(l)+ABL
              ! if(.not. lne) then
              EMLIN(l)=EMLIN(l)+ABL*SL0
              ! endif
            enddo
          endif IF_LPR      !}
        ENDIF
      enddo LINE_LOOP ! Main loop over all lines

      Ffudge=0.
      ABLIN(1:NFREQ)=ABLIN(1:NFREQ)+(opres)/(NFREQ-1)

      ABLINmean=sum(ABLIN)/NFREQ
      ABLINmin=minval(ABLIN)

      print '(i7)', AB0_COUNTER

C     LTE line emissivity

CMH   IF THE DEPTPOINT IS OUTWARD OF THE TEMPERATURE MIMIMUM
CMH   THE PLANCK-FUNCION (OR TEMPERATURE) IS SET TO THE
CMH   VALUE AT THE TEMPERATURE MINIMUM (ID = 55)
C     print *,ID,'LINOP: USE THE PLANCK FUNKTION OF ID=55 (FOR MODEL C)'

!     always true

!      CONDITION = velw(id) .le. velmax

!      PRINT*, 'SUBSTITUTING NDPMIN:', MAX(NDPMIN, ID), WDIL(ID), CONDITION

      if (velw(id) .le. velmax) EMLIN(1 : NFREQ) = EMLIN(1:NFREQ)+ABLIN(1:NFREQ)*PLAN(max(NDPMIN,id))*WDIL(ID)
!      if (velw(id) .le. velmax) EMLIN(1 : NFREQ) = EMLIN(1:NFREQ)+ABLIN(1:NFREQ)*PLAN(60)*WDIL(ID)

      deallocate(opacw)

!RT   HERE WE ADD A LOOP OVER ALL NLTE TRANSITIONS CALCULATED IN HMINUS (EXCEPT FOR HYDROGEN) BECAUSE WE DID NOT MANAGE TO
!RT   FIGURE OUT HOW TO TREAT NLTE LINES IN FIOSS IN...NLTE.

!      PRINT*, 'LINOP NTC:', NTC

!      IF (NTC .GE. 1) THEN

!          PRINT*, 'ENTERED THE LOOP'

!          DO NT = 1, NTC

!RT          FOR THE REFERENCE SEE RUTTEN, RADIATIVE TRANSFER IN STELLAR ATMOSPHERES, P.23

!             PRINT*, 'ELID :', ELID(NTI(NT))

!             IF (ELID(NTI(NT)) .EQ. 1) CYCLE ! HYDROGEN IS ALREADY TREATED IN NLTE ELSEWHERE IN THE CODE (NO IDEA WHERE THOUGH)
 
!             WAV = WAVTRA(NTI(NT)) * 1D-8 ! WAVELENGTH SHOULD BE IN CM

!             PRINT*, 'WAVTRA LINOP:', WAVTRA(NTI(NT))
 
!             FRE = LIGHT_SPEED / WAV
 
!             XFR(1 : NFREQ) = ABS(FREQ(1 : NFREQ) - FRE) * DOPA1(ELID(NTI(NT)), ID)! * 10.0D0 ! FREQUENCY IN THE DIMENSIONLESS UNITS (X)
 
!             BLU = NF1 * (SWU(NTI(NT)) / SWL(NTI(NT))) * FRE**(-3.0D0) * AUL(NTI(NT))
 
!             DCR = NTDU(ID, NT) / NTDL(ID, NT) ! RATIO OF THE NLTE DEPARTURE COEFFICIENTS OF THE UPPER AND LOWER TRANSITION LEVELS
 
!             IEF = 1.0D0 - DCR * DEXP(-PLANCK * FRE / BOLTZ / TEMP(ID)) ! INDUCED EMISSION FACTOR
 
!             NF2 = PLANCK * FRE / 4.0D0 / PAI ! SECOND NORMALIZATION FACTOR FOR NLTE LINE TREATMENT
 
!             NLTE_ABS = NF2 * NTPL(ID, NT) * BLU * IEF * 0.0D0
 
!             NLTE_EMI = NF2 * NTPU(ID, NT) * AUL(NTI(NT)) * 0.0D0
 
!             ABLIN(1 : NFREQ) = ABLIN(1 : NFREQ) + NLTE_ABS * VOIGTK_MS(0.0D0, XFR(1 : NFREQ), NFREQ)! / SQRT(PAI)
!             EMLIN(1 : NFREQ) = EMLIN(1 : NFREQ) + NLTE_EMI * VOIGTK_MS(0.0D0, XFR(1 : NFREQ), NFREQ)! / SQRT(PAI)

!             WRITE(*, '(A,2x,I2,2x,I2,4(2x,E15.7))') 'SF:', 
!     $       NT, ID, TEMP(ID), NLTE_EMI / NLTE_ABS / PLAN(ID), NTDU(ID, NT), NTDL(ID, NT)

!          ENDDO

!      ENDIF

!     add opacity of NLTE lines (RT: this is a bit from the original NLTE treatment,
!     the one the functioning of which we didn't quite understand, it doesn't pertain to what we implemented in the loop above)

      ABLIN(1:NFREQ)=ABLIN(1:NFREQ)+ABLINN(:)

!     special routine for selected He II lines

      IF(NSP.NE.0) THEN
        DO l=1,NSP
          ISP=ISP0(l)
          if(velw(id).gt.velmax.and.nltoff.gt.0) cycle
          IF(ISP.GE.6.AND.ISP.LE.24) CALL PHE2(ISP,ID,ABLIN,EMLIN)
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE LINOP_MS

      PURE FUNCTION VOIGTK_MS(A,V,NFREQ)
      !====================
      !*****************************************************************
      !** Voigt function after Kurucz (in Computational Astrophysics)
      !** also see Mihalas 2nd Edition, p. 279
      implicit none
      integer,parameter :: MVOI=2001
      real*8, parameter :: SQRT2  =  1.4142135623730950488_8
      real*8, parameter :: C11283 =  1.12838_8
      !real*8, parameter :: C32    =  3.2_8
      real*8, parameter :: C05642 =   .5642_8
      real*8, parameter :: C79788 =   .79788_8
      real*8, parameter :: C14    =  1.4_8
      real*8, parameter :: C37613 =   .37613_8
      real*8, parameter :: C23    =  2._8/3._8
      real*8, parameter :: CV1    = - .122727278_8
      real*8, parameter :: CV2    =   .532770573_8
      real*8, parameter :: CV3    = - .96284325_8
      real*8, parameter :: CV4    =   .979895032_8
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      integer,intent(in) ::NFREQ
      real*8,dimension(NFREQ),intent(in) ::V
      real*8 :: VOIGTK_MS(NFREQ)
      real*8 :: H0TAB,H1TAB,H2TAB,AA, H0TABi,H2TABi,Vi
      real*8,intent(in) :: A
      real*8  :: VV,HH1,HH2,HH3,HH4, U,UU
      integer :: IV
      integer :: i
      ! Method 1: A < 0.2
      IF(A.LT. 0.2) THEN
        do i=1,NFREQ  !** DO Loops are faster than where/elsewhere
          VI=V(i)
          if(VI<=10d0) then
            IV=V(i)*200d0+1.5_8
            VOIGTK_MS(i)=(H2TAB(IV)*A+H1TAB(IV))*A+H0TAB(IV)
          else
            VOIGTK_MS(i)=C05642*A/(VI**2)
          endif
        enddo
      else
        do i=1,NFREQ
          Vi = V(i)
          VV=Vi**2
          if(A<=C14 .and. A+Vi <= 3.2d0) then
            IV=Vi*200._4+1.5_4
            H0TABI = H0TAB(IV)
            H2TABI = H2TAB(IV)
            IV=Vi*200d0+1.5_8

            HH1=H1TAB(IV)+H0TABi*C11283
            HH2=H2TABi+HH1*C11283-H0TABi
            HH3=(1d0-H2TABi)*C37613-HH1*C23*VV+HH2*C11283
            HH4=(3d0*HH3-HH1)*C37613+H0TABi*C23*VV*VV
            VOIGTK_MS(i)=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0TABi)*
     *           (((CV1*A+CV2)*A+CV3)*A+CV4)
          else
            AA=A**2
            U=(AA+VV)*SQRT2
            UU=U**2
            VOIGTK_MS(i)=((((AA-10.*VV)*AA*3.+15.*VV**2)/UU+3d0*VV-AA)
     *                    /UU+1d0)*A*C79788/U
          endif
        enddo
      endif
      END FUNCTION VOIGTK_MS


      function getLocalTMin(T,ND)
      implicit none
      integer :: getLocalTMin
      integer,intent(in) ::ND
      real*8,intent(in) :: T(ND)
      integer :: NDPMIN,L
      IF(ND==1) THEN
        getLocalTMin=1
        RETURN
      ELSE
      NDPMIN = 0                    ! Sanity check 1
      DO L=2,ND-1                   ! Find local minimum in T
        if ((T(L) .LT. T(L-1)) .AND. (T(L) .LT. T(L+1))) THEN
          NDPMIN = L
        endif
      ENDDO

!      PRINT*, 'TEMP MIN', NDPMIN, T(NDPMIN)

      if (NDPMIN .eq. 0) then       ! Sanity check 1 - finish
        print *,'linop_ms: something wrong! NDPMIN = ',NDPMIN,
     $   ' setting NDPMIN to 1'
        NDPMIN = 1              ! continue, dont abort
      endif
      if ((T(NDPMIN) .gt. 5000d0) .and. (ND .gt. 1)) then !Sanity check 2
        print *,'something wrong with Temperature minimum!',NDPMIN, T(NDPMIN)
        write(6,*)'something wrong with Temperature minimum!',
     $                        NDPMIN,T(NDPMIN)
        pause
      endif
      ENDIF
    !! ATTENTION
     
 !      NDPMIN=1

    !! ATTENTION

      getLocalTMin = NDPMIN
      end function getLocalTMin

      end module MOD_LINOP_MS
