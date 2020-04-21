!  This module is a common block similar in all AFT Model programs and is
!  written in FORTRAN 90.
!                     J G Li   26 Oct 2000
!! Adapted for multiple cell 2D advection tests using UNO schemes.
!!                    J G Li    8 Aug 2007
!! Adapted for global multiple cell grid   J G Li   12 Nov 2007
!! Modified for SMC extended over Arctic Ocean  J G Li   26 Nov 2008
!! Adapted for multi-resolution UK3 to Global 25km grid  J G Li   8 Feb 2010 
!! Modified to add minimum y-size in V-flux array.       J G Li  16 Feb 2010 
!! Adapted for 6-25km global ocean SMC grid.  J G Li   22 Feb 2010
!! Modified to use new rules on second cell selection.  J G Li   26 Feb 2010 
!! Modified for SMC625 grid global part only.   J G Li   12 Dec 2011 
!! Restore second cell selection by one face end equal.   JGLi28Feb2012 
!! Modified to use new cell count line for SMC625 model.  JGLi01Oct2012 
!! Adapted for SMC6125 grid with refined UK waters.   JGLi08Jan2013 
!!

MODULE Constants
   IMPLICIT NONE

! Parameters fixed in the program
   INTEGER,PARAMETER::NCL=2250000, NFC=2360000

! Variables to be used for data storage
   INTEGER:: NU, NV, NC, N9, N8, N4, N2, N1, NLon, NLONMAX, NCMAX
   INTEGER:: ICE(4,NCL), KG(NCL)
   INTEGER, DIMENSION(7,NFC)::  ISD
   INTEGER, DIMENSION(8,NFC)::  JSD
   INTEGER:: I,II,IJ,IJK,J,JJ,JK,JKL,K,KK,KL,KLM,L,LL,LM,LMN,M,MM,MN,N,NN
   INTEGER:: SBD, ABD
   LOGICAL:: ARCTICG

END MODULE Constants

!!
!! This program generates the one level 2D spherical multiple cell grid
!!  face arrays.
!! v1.0 Adapted for multiple cell 2D advection tests using UNO schemes.
!!       J G Li   26 Jul 2007
!! v2.0 Preliminary work to generalise for any SMC grid.
!!       A Saulter 30 Apr 2020
!!

 PROGRAM AdapGrid 
   USE Constants
   IMPLICIT NONE

   REAL:: CNST, CNST1, CNST2, CNST3, BX, BY

   REAL:: BXSB, BYSB, BX0, BY0
   INTEGER:: NBLAT, NBLON, NLEVS, NLat, NPol
   LOGICAL:: ARCTIC
   CHARACTER*80:: FNAME 
   NAMELIST /GRID_NML/ NLEVS, NBLAT, NBLON, BXSB, BYSB, BX0, BY0, ARCTIC, FNAME 

!  Read the grid NAMELIST
   OPEN(UNIT=11, FILE='smcSides.nml',STATUS='OLD',IOSTAT=nn,ACTION='READ')
   READ (UNIT=11, NML=GRID_NML) 
   CLOSE(11)

!  Calculate cell and dx,dy values for highest tier
   NLon = NBLON * 2**(NLEVS-1)
   NLat = NBLAT * 2**(NLEVS-1)
   NPol = NLat
   BX = BXSB / 2.0**(NLEVS-1)
   BY = BYSB / 2.0**(NLEVS-1)
   NLONMAX = NINT(360.0 / BX)
   ARCTICG = ARCTIC
   
!  Grid corners set to lower left from read in centre value
   BX0 = BX0 - BXSB/2.0
   BY0 = BY0 - BYSB/2.0

!  Read Global Multiple-Cell info

   OPEN(UNIT=9, FILE=TRIM(FNAME),STATUS='OLD',IOSTAT=nn,ACTION='READ')
   IF(nn /= 0) PRINT*,' Cells .dat file was not opened! '
   IF (ARCTIC) THEN
      READ (9,*) NC, SBD, ABD
   ELSE
      IF (NLEVS .EQ. 1) THEN
         READ (9,*) NC, N1 
      ELSE IF (NLEVS .EQ. 2) THEN
         READ (9,*) NC, N1, N2
      ELSE IF (NLEVS .EQ. 3) THEN
         READ (9,*) NC, N1, N2, N4 
      ELSE IF (NLEVS .EQ. 4) THEN
         READ (9,*) NC, N1, N2, N4, N8 
      ENDIF
   ENDIF
   DO J=1,NC
      READ (9,*) ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), KG(J)
   END DO
   CLOSE(9)
   PRINT*, ' Cells .dat read done  NC=', NC

!  Output a few to check input values
   DO J=NC/10,NC,NC/10
      WRITE(6,'(i8,2i6,2i4,i6)') J, ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), KG(J)
   END DO

!  Call subroutines to generate flux faces

   CALL CellSide

!  Open files to store writups
   IF (ARCTIC) THEN
      OPEN(UNIT=16,FILE='ww3ASide.txt',STATUS='UNKNOWN',IOSTAT=nn,ACTION='WRITE')
   ELSE
      OPEN(UNIT=16,FILE='ww3GSide.txt',STATUS='UNKNOWN',IOSTAT=nn,ACTION='WRITE')
   ENDIF
   IF(nn /= 0) PRINT*,' Output file ww3Side.txt was not opened! '

!  Header messages and configuration information 
   WRITE(UNIT=16,FMT='(1x/   &
     &  "  Multiple Cell 2D Grid Generation Output " /)' )

   WRITE(UNIT=16,FMT='(1x," Cell Units  BX  BY deg = ",2f14.11)')  BX, BY
   WRITE(UNIT=16,FMT='(1x," Cell Units XLon0 YLat0 = ",2f14.8 )')  BX0, BY0
   WRITE(UNIT=16,FMT='(1x," Horizontal cell number = ",i8)' )  NC
   IF (ARCTIC) THEN
      WRITE(UNIT=16,FMT='(1x," Arctic boundary number = ",3i8)')  SBD, ABD
      WRITE(UNIT=16,FMT='(1x," Lon/Lat/NPol grid No.s = ",3i8)')  NLon, NLat, NPol
   ELSE
      IF (NLEVS .EQ. 1) THEN
         WRITE(UNIT=16,FMT='(1x," Size 1 cell number = ",i8)')  N1
      ELSE IF (NLEVS .EQ. 2) THEN
         WRITE(UNIT=16,FMT='(1x," Size 1 2 cell number = ",2i8)')  N1, N2
      ELSE IF (NLEVS .EQ. 3) THEN
         WRITE(UNIT=16,FMT='(1x," Size 1 2 4 cell number = ",3i8)')  N1, N2, N4
      ELSE IF (NLEVS .EQ. 4) THEN
         WRITE(UNIT=16,FMT='(1x," Size 1 2 4 8 cell number = ",4i8)')  N1, N2, N4, N8
      ENDIF
      WRITE(UNIT=16,FMT='(1x," Lon/Lat grid No.s = ",3i8)')  NLon, NLat
   ENDIF
   WRITE(UNIT=16,FMT='(1x," Max Lon grid No.s = ",i8)')  NLONMAX
   WRITE(UNIT=16,FMT='(1x," Total number of U-face = ",i8)' )  NU
   WRITE(UNIT=16,FMT='(1x," Total number of V-face = ",i8)' )  NV

 3912 FORMAT(1x,i4,3F9.1,ES12.3)
 638  FORMAT(5(i6,f8.2))

!  Close all files
   CLOSE(16)

 9999  PRINT*, ' AdapGrid completed '

 END PROGRAM AdapGrid 
!  End of main program


! Subroutine that generates the cell side information
 SUBROUTINE CellSide
   USE Constants
   IMPLICIT NONE
   REAL:: CNST, CNST1, CNST2, CNST3
   LOGICAL:: VERBOSE=.FALSE.

!!    Test integer division for boundary cell numbers
!!    Size 2**n cell should bounded by -n cells
      DO II=0, 8
         JJ=2**II
         CNST=FLOAT(JJ)
         K=-INT( LOG(CNST)/LOG(2.) + 0.1)
      Write(6,*) " Cell size and boundary cell index =", JJ, K
      ENDDO

!     Generate i & j inter-sides for all cells
      WRITE(6,*) " Start creating inner face ..."

!     Generate i & j inter-sides for all cells at the highest level

      II=0
      JJ=0
      IF (ARCTICG) THEN
         ! Exclude last cell, the North Polar cell.
         NCMAX = NC-1
      ELSE
         NCMAX = NC
      ENDIF

      DO L=1, NCMAX

         IF(MOD(L, 5000) .eq. 0) WRITE(6,*) " Done L=", L, II, JJ

!!  Cyclic boundary for i-side at L-cell east side
         LM=ICE(1,L)+ICE(3,L)
         IF(LM .ge. NLONMAX) THEN
            IF(VERBOSE) print*, "Wrapping lon" , LM, NLONMAX
            LM=LM-NLONMAX
         ENDIF

!!  Cell height size is different from width size sometimes
         KL=ICE(2,L)+ICE(4,L)

         DO M=1, NCMAX

!!  Cyclic boundary for i-side at M-cell east side
            MN=ICE(1,M) + ICE(3,M)
            IF(MN .ge. NLONMAX) THEN 
               IF(VERBOSE) print*, "Wrapping lon" , MN, NLONMAX
               MN=MN-NLONMAX
            ENDIF

!!  U-faces
            IF(( ICE(1,M) .eq. LM ) .AND.          &
     &         ( ICE(2,M)+ICE(4,M) .eq. KL .OR.    &
     &           ICE(2,M) .eq. ICE(2,L) ))  THEN  
               II=II+1
               ISD(1,II)=ICE(1,M)
               ISD(2,II)=MAX(ICE(2,M), ICE(2,L)) 
               ISD(3,II)=MIN(ICE(4,M), ICE(4,L)) 
               ISD(5,II)=L
               ISD(6,II)=M 
            ENDIF

!!   V-faces
            IF(( ICE(2,M) .eq. KL ) .AND.        &
               ( ICE(1,M) .eq. ICE(1,L) .OR. MN .eq. LM ))  THEN 
               JJ=JJ+1
               JSD(1,JJ)=MAX(ICE(1,M), ICE(1,L)) 
               JSD(2,JJ)=ICE(2,M)
               JSD(3,JJ)=MIN(ICE(3,M), ICE(3,L)) 
               JSD(5,JJ)=L
               JSD(6,JJ)=M 
!!  Minimum Y-size of the two bounding cells will be used to sort 
!!  cell sizes for multi-step implementation.
               JSD(8,JJ)=MIN(ICE(4,M), ICE(4,L)) 
            ENDIF
          END DO
      END DO
 
      IJK=II
      LMN=JJ

!     Set boundary u faces
      WRITE(6,*) " Start creating u boundary face II JJ=", II, JJ

      DO 111 L=1, NCMAX
         I=0
         J=0
         IJ=0
         K=0
         N=0
         KK=0
         
         IF(MOD(L, 10000) .eq. 0) WRITE(6,*) " Done L II=", L, II

!!    Cyclic boundary need to be taken into account
         LM=ICE(1,L)+ICE(3,L)
         IF(LM .ge. NLONMAX) THEN
            IF(VERBOSE) print*, "Wrapping lon" , LM, NLONMAX
            LM=LM-NLONMAX
         ENDIF

!!    Loop through all inner faces 
         DO M=1, IJK

!!    See if the L cell west face is covered
            IF( ISD(1,M) .eq. ICE(1,L) ) THEN
                IF( ISD(2,M) .eq. ICE(2,L) ) THEN
                    I=1 
                    IJ=IJ+ISD(3,M)
                ELSEIF( ISD(2,M)+ISD(3,M) .eq. ICE(2,L)+ICE(4,L) ) THEN
                    J=1
                    IJ=IJ+ISD(3,M)
                ENDIF
            ENDIF

!!          and see if the L cell east face is covered
            IF( ISD(1,M) .eq. LM ) THEN
                IF( ISD(2,M) .eq. ICE(2,L) )  THEN 
                    K=1
                    KK=KK+ISD(3,M)
                ELSEIF( ISD(2,M)+ISD(3,M) .eq. ICE(2,L)+ICE(4,L) ) THEN
                    N=1
                    KK=KK+ISD(3,M)
                ENDIF
            ENDIF

!!  End of inner face M=1,IJK loop
         END DO

         IF(KK+IJ .gt. 2*ICE(4,L) )  WRITE(6,*) "Over done i-side for cell L,IJ,KK=", L, IJ, KK
         IF(KK+IJ .ge. 2*ICE(4,L) )  GOTO  111

          IF(IJ .eq. 0)  THEN
!!  Full boundary cell for west side
               II=II+1
               ISD(1,II)=ICE(1,L)
               ISD(2,II)=ICE(2,L)
               ISD(3,II)=ICE(4,L)
!!  New boundary cells proportional to cell y-sizes 
!!  Updated for any 2**n sizes
               ISD(5,II)=-INT( LOG(FLOAT(ISD(3,II)))/LOG(2.) + 0.01 )
               ISD(6,II)=L
          ENDIF
          IF(KK .eq. 0)  THEN
!!  Full boundary cell for east side
               II=II+1
               ISD(1,II)=LM
               ISD(2,II)=ICE(2,L)
               ISD(3,II)=ICE(4,L)
               ISD(5,II)=L
!!  Updated for any 2**n sizes
               ISD(6,II)=-INT( LOG(FLOAT(ISD(3,II)))/LOG(2.) + 0.01 )
          ENDIF

!!  Half cell size west boundary faces
          IF(IJ .gt. 0  .and. IJ .lt. ICE(4,L) )  THEN
             IF( I .eq. 0 )  THEN
!!  lower half west cell face
               II=II+1
               ISD(1,II)=ICE(1,L)
               ISD(2,II)=ICE(2,L)
               ISD(3,II)=ICE(4,L)/2
!!  Updated for any 2**n sizes
               ISD(5,II)=-INT( LOG(FLOAT(ISD(3,II)))/LOG(2.) + 0.01 )
!!  Size 1 for cell 0, size 2 uses cell -1 and size 4 uses cell -2
               ISD(6,II)=L
             ENDIF
             IF( J .eq. 0 )  THEN
!!  Upper half west cell face
               II=II+1
               ISD(1,II)=ICE(1,L)
               ISD(2,II)=ICE(2,L)+ICE(4,L)/2
               ISD(3,II)=ICE(4,L)/2
!!  Updated for any 2**n sizes
               ISD(5,II)=-INT( LOG(FLOAT(ISD(3,II)))/LOG(2.) + 0.01 )
               ISD(6,II)=L
             ENDIF
          ENDIF

!!  Half cell size east boundary faces
          IF(KK .gt. 0  .and. KK .lt. ICE(4,L) )  THEN
             IF( K .eq. 0 )  THEN
!!  lower half east cell face
               II=II+1
               ISD(1,II)=LM
               ISD(2,II)=ICE(2,L)
               ISD(3,II)=ICE(4,L)/2
!!  Size 1 for cell 0, size 2 uses cell -1 and size 4 uses cell -2
               ISD(5,II)=L
!!  Updated for any 2**n sizes
               ISD(6,II)=-INT( LOG(FLOAT(ISD(3,II)))/LOG(2.) + 0.01 )
             ENDIF
             IF( N .eq. 0 )  THEN
!!  Upper half west cell face
               II=II+1
               ISD(1,II)=LM
               ISD(2,II)=ICE(2,L)+ICE(4,L)/2
               ISD(3,II)=ICE(4,L)/2
!!  Size 1 for cell 0, size 2 uses cell -1 and size 4 uses cell -2
               ISD(5,II)=L
!!  Updated for any 2**n sizes
               ISD(6,II)=-INT( LOG(FLOAT(ISD(3,II)))/LOG(2.) + 0.01 )
             ENDIF
          ENDIF

 111  CONTINUE


!     Set boundary v faces
      WRITE(6,*) " Start creating v boundary face II JJ=", II, JJ

      DO 222 L=1, NCMAX
         I=0
         J=0
         IJ=0
         K=0
         N=0
         NN=0
         
         IF(MOD(L, 10000) .eq. 0) WRITE(6,*) " Done L JJ=", L, JJ

!!    Loop through all V faces already set 
         DO M=1, LMN

!!    See if the L cell south face is covered
            IF( JSD(2,M) .eq. ICE(2,L) ) THEN
               IF( JSD(1,M) .eq. ICE(1,L) ) THEN
                  I=1
                  IJ=IJ+JSD(3,M)
               ELSEIF( JSD(1,M)+JSD(3,M) .eq. ICE(1,L)+ICE(3,L) ) THEN
                  J=1
                  IJ=IJ+JSD(3,M)
               ENDIF
            ENDIF
!!          and see if the L cell north face is covered
            IF( JSD(2,M) .eq. ICE(2,L) + ICE(4,L) )  THEN 
               IF( JSD(1,M) .eq. ICE(1,L) ) THEN 
                  K=1
                  NN=NN+JSD(3,M)
               ELSEIF( JSD(1,M)+JSD(3,M) .eq. ICE(1,L)+ICE(3,L) )  THEN
                  N=1
                  NN=NN+JSD(3,M)
               ENDIF
            ENDIF

!!   End M=1, LMN V-side loop
         END DO

         IF(NN+IJ .gt. 2*ICE(3,L) )  WRITE(6,*)  "Over done j-side for L, IJ, NN=", L, IJ, NN
         IF(NN+IJ .ge. 2*ICE(3,L) )  GOTO  222

          IF(IJ .eq. 0)  THEN
!!  Full boundary cell for south side
               JJ=JJ+1
               JSD(1,JJ)=ICE(1,L)
               JSD(2,JJ)=ICE(2,L)
               JSD(3,JJ)=ICE(3,L)
!!  New boundary cells proportional to cell sizes 
!!  Updated for any 2**n sizes
               JSD(5,JJ)=-INT( LOG(FLOAT(ICE(3,L)))/LOG(2.) + 0.01 )
               JSD(6,JJ)=L
               JSD(8,JJ)=ICE(4,L)
!!  No cells over Antarctic land so there is no S Polar cell.
          ENDIF
          IF(NN .eq. 0)  THEN
!!  Full boundary cell for north side
               JJ=JJ+1
               JSD(1,JJ)=ICE(1,L)
               JSD(2,JJ)=ICE(2,L)+ICE(4,L)
               JSD(3,JJ)=ICE(3,L)
               JSD(5,JJ)=L
               IF (ARCTICG) THEN
                  !  North polar cell takes the whole last 4 rows above JSD=ICE(2,NC).
                  !  Note ICE(2,L) represents lower-side of the cell.  
                  !  Polar cell is the last cell NC.
                  IF( ICE(2,L)+ICE(4,L) .eq. ICE(2,NC) ) THEN
                     JSD(6,JJ)=NC
                     WRITE(6,*) "Set north pole v face for cell L", L
                  ELSE
                     !  Updated for any 2**n sizes
                     JSD(6,JJ)=-INT( LOG(FLOAT(ICE(3,L)))/LOG(2.) + 0.01 )
                  ENDIF
               ELSE
                  !  Updated for any 2**n sizes
                  JSD(6,JJ)=-INT( LOG(FLOAT(ICE(3,L)))/LOG(2.) + 0.01 )
               ENDIF
               JSD(8,JJ)=ICE(4,L)
          ENDIF

!!  Half cell size south boundary faces
          IF(IJ .gt. 0  .and. IJ .lt. ICE(3,L) )  THEN
             IF( I .eq. 0 )  THEN
!!  left half cell face
               JJ=JJ+1
               JSD(1,JJ)=ICE(1,L)
               JSD(2,JJ)=ICE(2,L)
               JSD(3,JJ)=ICE(3,L)/2
!!  New boundary cells proportional to cell sizes 
!!  Updated for any 2**n sizes
               JSD(5,JJ)=-INT( LOG(FLOAT(JSD(3,JJ)))/LOG(2.) + 0.01 )
               JSD(6,JJ)=L
               JSD(8,JJ)=ICE(4,L)
             ENDIF
             IF( J .eq. 0 )  THEN
!!  right half cell face
               JJ=JJ+1
               JSD(1,JJ)=ICE(1,L)+ICE(3,L)/2
               JSD(2,JJ)=ICE(2,L)
               JSD(3,JJ)=ICE(3,L)/2
!!  New boundary cells proportional to cell sizes 
!!  Updated for any 2**n sizes
               JSD(5,JJ)=-INT( LOG(FLOAT(JSD(3,JJ)))/LOG(2.) + 0.01 )
               JSD(6,JJ)=L
               JSD(8,JJ)=ICE(4,L)
             ENDIF
          ENDIF

!!  Half cell size north boundary faces
          IF(NN .gt. 0  .and. NN .lt. ICE(3,L) )  THEN
             IF( K .eq. 0 )  THEN
!!  left half north cell face
               JJ=JJ+1
               JSD(1,JJ)=ICE(1,L)
               JSD(2,JJ)=ICE(2,L)+ICE(4,L)
               JSD(3,JJ)=ICE(3,L)/2
               JSD(5,JJ)=L
!!  New boundary cells proportional to cell sizes 
!!  Updated for any 2**n sizes
               JSD(6,JJ)=-INT( LOG(FLOAT(JSD(3,JJ)))/LOG(2.) + 0.01 )
               JSD(8,JJ)=ICE(4,L)
             ENDIF
             IF( N .eq. 0 )  THEN
!!  right half north cell face
               JJ=JJ+1
               JSD(1,JJ)=ICE(1,L)+ICE(3,L)/2
               JSD(2,JJ)=ICE(2,L)+ICE(4,L)
               JSD(3,JJ)=ICE(3,L)/2
               JSD(5,JJ)=L
!!  New boundary cells proportional to cell sizes 
!!  Updated for any 2**n sizes
               JSD(6,JJ)=-INT( LOG(FLOAT(JSD(3,JJ)))/LOG(2.) + 0.01 )
               JSD(8,JJ)=ICE(4,L)
             ENDIF
          ENDIF

 222  CONTINUE

!   Store top level U V side numbers in NU NV 
      NU=II
      NV=JJ

!!  Loop over all u faces to find the second cells next to the L and M cells
!!  Boundary cells will be duplicated for second cells
      WRITE(6,*) " Find extra second u cell II JJ=", II, JJ
      DO I=1, NU
         L=ISD(5,I)
         M=ISD(6,I)
         KK=0
         NN=0

!!  Boundary L cell just duplicate it as LL cell
         IF(L .LE. 0) THEN
            ISD(4,i)=L
            KK=1
         ELSE
!!  Find the second LL cell by loop over all faces again
!!  The two U faces have to share at least one y-end.
!!  The second U face should be no less than this face.
!!  Restore one face end equal and suspend no less requirement. JGLi28Feb2012
            DO K=1, NU
!              IF( (L .EQ. ISD(6,K)) .AND. (ISD(3,I) .LE. ISD(3,K)) .AND.  &
               IF( (L .EQ. ISD(6,K)) .AND.   &
     &             ( (ISD(2,I)+ISD(3,I) .eq. ISD(2,K)+ISD(3,K)) .or.       &
     &               (ISD(2,I) .eq. ISD(2,K)) ) )  THEN
                   ISD(4,I)=ISD(5,K)
                   KK=1
               ENDIF
            ENDDO
         ENDIF

!!  Boundary M cell just duplicate it as MM cell
         IF(M .LE. 0) THEN
            ISD(7,I)=M
            NN=1
         ELSE
!!  Find the second MM cell by loop over all faces again
!!  The two U faces have to share at least one y-end.
            DO N=1, NU
!              IF( (M .EQ. ISD(5,N)) .AND. (ISD(3,I) .LE. ISD(3,N)) .AND.   &
               IF( (M .EQ. ISD(5,N)) .AND.                                  &
     &             ( (ISD(2,I)+ISD(3,I) .eq. ISD(2,N)+ISD(3,N)) .or.        &
     &               (ISD(2,I) .eq. ISD(2,N)) ) )  THEN
                   ISD(7,I)=ISD(6,N)
                   NN=1
               ENDIF
            ENDDO
         ENDIF

!!  Duplicate central cell if upstream cells are not selected.
         IF( KK .eq. 0) THEN
            WRITE(6,*) " L cell duplicated for U face I=", I
            ISD(4,I)=L
         ENDIF
         IF( NN .eq. 0) THEN
            WRITE(6,*) " M cell duplicated for U face I=", I
            ISD(7,I)=M
         ENDIF

         IF(MOD(I, 10000) .eq. 0) WRITE(6,*) " Done U face I=", I

!!  End of u face loop
      ENDDO

!!  Loop over all v faces to find the second cells next to the L and M cells
!!  Boundary cells will be duplicated for second cells
      WRITE(6,*) " Find extra second v cell II JJ=", II, JJ

      DO J=1, NV
         L=JSD(5,J)
         M=JSD(6,J)
         KK=0
         NN=0

!!  Boundary L cell just duplicate it as LL cell
         IF(L .LE. 0) THEN
            JSD(4,J)=L
            KK=1
         ELSE
!!  Find the second LL cell by loop over all faces again
!!  The two V faces have to share at least one x-end.
!!  Restore one face end equal and suspend no less requirement. JGLi28Feb2012
            DO K=1, NV
!              IF( (L .EQ. JSD(6,K)) .AND. (JSD(3,J) .LE. JSD(3,K)) .AND.  &
               IF( (L .EQ. JSD(6,K)) .AND.                                 &
     &             ( (JSD(1,J)+JSD(3,J) .eq. JSD(1,K)+JSD(3,K)) .or.       &
     &               (JSD(1,J) .eq. JSD(1,K)) ) )  THEN 
                   JSD(4,J)=JSD(5,K)
                   KK=1
               ENDIF
            ENDDO
         ENDIF

!!  Boundary M cell just duplicate it as MM cell
!!  Duplicate north polar cell as well
         IF((ARCTICG) .AND. (M .EQ. NC)) THEN
            JSD(7,J)=M
            NN=1
         ELSEIF(M .LE. 0 ) THEN
            JSD(7,J)=M
            NN=1
         ELSE
!!  Find the second LL cell by loop over all faces again
!!  The two V faces have to share at least one x-end.
            DO N=1, JJ
!              IF( (M .EQ. JSD(5,N)) .AND. (JSD(3,J) .LE. JSD(3,N)) .AND.  &
               IF( (M .EQ. JSD(5,N)) .AND.                                 &
     &             ( (JSD(1,J)+JSD(3,J) .eq. JSD(1,N)+JSD(3,N)) .or.       &
     &               (JSD(1,J) .eq. JSD(1,N)) ) )  THEN 
                   JSD(7,J)=JSD(6,N)
                   NN=1
               ENDIF
            ENDDO
         ENDIF

         IF( KK .eq. 0) THEN
             WRITE(6,*) " L cell duplicated for V face J=", J
             JSD(4,J)=L
         ENDIF
         IF( NN .eq. 0) THEN
             WRITE(6,*) " M cell duplicated for V face J=", J
             JSD(7,J)=M
         ENDIF

         IF(MOD(J, 10000) .eq. 0) WRITE(6,*) " Done V face J=", J

!!  End of v face loop
      ENDDO

!!  Check whether any overlaping exists
      WRITE(6,*) " Check any overlaping NU NV=", NU, NV

      IJ=0
      DO I=1, NU-1
         L=ISD(1,I)
         M=ISD(2,I)
            DO K=I+1, NU
               IF( L .EQ. ISD(1,K) .AND. M .EQ. ISD(2,K) )  THEN
                 IJ=IJ+1
                 WRITE(6,*) IJ, ' Overlaping U face K, I, J, L, MM, M, N, NN' 
                 WRITE(6,333) I, (ISD(N,I), N=1,7)
                 WRITE(6,333) K, (ISD(N,K), N=1,7)
               ENDIF
            ENDDO
         IF(MOD(I, 10000) .eq. 0) WRITE(6,*) " Checked U face I=", I
      ENDDO

 333  FORMAT(7I8)

      IJ=0
      DO J=1, NV-1
         L=JSD(1,J)
         M=JSD(2,J)
            DO K=J+1, NV
               IF( L .EQ. JSD(1,K) .AND. M .EQ. JSD(2,K) )  THEN
                 IJ=IJ+1
                 WRITE(6,*) IJ, ' Overlaping V face K, I, J, L, MM, M, N, NN'
                 WRITE(6,333) J, (JSD(N,J), N=1,7)
                 WRITE(6,333) K, (JSD(N,K), N=1,7)
               ENDIF
            ENDDO
         IF(MOD(J, 10000) .eq. 0) WRITE(6,*) " Checked V face J=", J
      ENDDO


!!  Output ISD JSD variables for later use
      WRITE(6,*) " Storing face array info NU,NV=", NU, NV

   IF (ARCTICG) THEN
      OPEN(UNIT=10,FILE='ww3AISide.d',STATUS='UNKNOWN',IOSTAT=nn)
   ELSE
      OPEN(UNIT=10,FILE='ww3GISide.d',STATUS='UNKNOWN',IOSTAT=nn)
   ENDIF
   IF(nn /= 0) PRINT*,' File ISide.d was not opened! '
      DO I=1,NU
         WRITE(10,FMT='(2i6,i4,4i8)')  (ISD(N,I), N=1,7)
      END DO
   CLOSE(10)

   IF (ARCTICG) THEN
      OPEN(UNIT=10,FILE='ww3AJSide.d',STATUS='UNKNOWN',IOSTAT=nn)
   ELSE
      OPEN(UNIT=10,FILE='ww3GJSide.d',STATUS='UNKNOWN',IOSTAT=nn)
   ENDIF
   IF(nn /= 0) PRINT*,' File JSide.d was not opened! '
      DO J=1,NV
         WRITE(10,FMT='(2i6,i4,4i8,i4)')  (JSD(N,J), N=1,8)
      END DO
   CLOSE(10)

   PRINT*, ' I J-Sides output done '

 999  PRINT*, ' Sub CellSide ended.'

      RETURN

 END SUBROUTINE CellSide


