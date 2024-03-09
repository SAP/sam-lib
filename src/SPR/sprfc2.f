C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRFC2   ( MSPAR , MTREES, MSIFA , SM    , TOL   ,
     $                        IWORK , RWORK , LPU   , IERR  , MVMULT,
     $                        MMMULT                                  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           LPU , IERR
      INTEGER           MSPAR(*)              , MTREES(*)             ,
     $                  MSIFA(*)              , IWORK(*)

      DOUBLE PRECISION  SM(*)                 , TOL(3)                ,
     $                  RWORK(*)

C ======================================================================
C  S A M  library routine :  SPRFC2                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To perform the Cholesky decomposition CC'.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Dec. 10, 1996 (acd)
C                 Corrected errors in the length for IWORK.
C                 Applies only to the test in this routine since all
C                 other routines work on array of fixed length.
*                 Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
*                 Feb. 13, 2003 (kmo)
*                 IERR = -3 only on numerical errors (singular matrix).
*                 Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*                 Mar. 19, 2003 (kmo)
*                 Suppress error messages if LPU < 0.
*                 Mar. 26, 2003 (kmo)
*                 Added more error output when the matrix is singular.
*                 Aug. 14, 2003 (kmo)
*                 Added MSPAR(57) in call to SPRF21.
C
C     MSPAR - INTEGER(*)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed, except for MSPAR(1) that is set to 6 in
C              case of sucessful exit, and to -6 else.
C     MTREES - INTEGER(NTREES)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     SM - DOUBLE PRECISION(NSM)
C      Entry : As output from SPRADM.
C      Exit  : The L and D factors, A is destroyed.
C     TOL - DOUBLE PRECISION(3)
C      Entry : TOL(1) is set to 'numerical zero', may be 0.0.
C      Exit  : TOL(2) = max(a(i)/d(i))=min(d(i)/a(i))
C              TOL(3) = trace(a)/min(d(i))
C     LPU - INTEGER
C      Entry : Output unit.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call,
C              = -1 Error from SPRTRS is not cleared.
C              = -2 MTREES, MSIFA, SM, IWORK and/or RWORK
C                   are not great enough.
C              = -3 Singular matrix.
C              = -4 Other error from SPRF21.
C
C     Working arrays
C     --------------
C
C     IWORK - INTEGER(NIWORK)
C      Entry : Not defined. Partitioned to internal work array,
C              see below.
C      Exit  : Need not be saved.
C     RWORK - DOUBLE PRECISION(NRWORK)
C      Entry : Not defined.
C              Used to store the multifrontal stack.
C      Exit  : Need not be saved.
C
C     Procedures
C     ----------
C
C     SPRF21
C     SPRER1
C     MVMULT
C     MMMULT
C
C     Intrinsic
C     ---------
C
C     None
C
C     Include Blocks
C     --------------
C
C     None
C
C     Common Blocks
C     -------------
C
C     None
C
C     ------------------------------------------------------------------

      INTEGER
     $  NIWORK, NMSIFA, NRWORK, NSM   , NTREES

      INTEGER
     $  DIAG  , ELSEQ , FIRST , IWUSED, LFIRST, LIERR , LINDX , LNZ   ,
     $  MAXSUP, MFRONT, MLUSED, MTUSED, NB    , NBLOCK, NELACT, NEQ   ,
     $  NOFSUB, NSUPER, NUMNEG, NZEROL, RWUSED, SIBL  , SMUSED, SPLIT ,
     $  STACK , SUPEQT, SUPSUP, WSTACK, XBLOCK, XELSEQ, XLINDX, XLNZ  ,
     $  XSIBL , XSPLIT, XSUPER

      EXTERNAL          SPRF21, SPRER1, MVMULT, MMMULT

C     ==================================================================
#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRFC2'
#endif

C     -------------------
C     GET INFO FROM MSPAR
C     -------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -6
         IERR = -1
         CALL SPRER1 ( 30, 'SPRFC2', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         MSPAR(1) = 6
         IERR = 0
      ENDIF

      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NB     = MSPAR(12)
      MAXSUP = MSPAR(13)
      MFRONT = MSPAR(14)
      NOFSUB = MSPAR(15)
      NZEROL = MSPAR(16)
      WSTACK = MSPAR(17)
      NELACT = MSPAR(19)
      NBLOCK = MSPAR(22)

      NMSIFA = MSPAR( 3)
      NTREES = MSPAR( 4)
      NSM    = NEQ + NZEROL
      NIWORK = MFRONT + 1
      NRWORK = WSTACK

C     ---------------------------------
C     DIVIDE IWORK AND CHECK ITS LENGTH
C     ---------------------------------
      FIRST  = 1
      LFIRST = MFRONT + 1
      IWUSED = FIRST  + LFIRST - 1

C     -------------------------------------------
C     SET POINTERS TO MTREES AND CHECK ITS LENGTH
C     -------------------------------------------
      XSUPER = MSPAR(45)
      XSIBL  = XSUPER + NSUPER + 1
      SIBL   = XSIBL  + NSUPER + 1
      XLINDX = MSPAR(46)

* === SUPSUP-acd
      SUPSUP = XLINDX + NSUPER + 1
      ELSEQ  = SUPSUP + NSUPER
      XELSEQ = ELSEQ  + NELACT
      XBLOCK = XELSEQ + NB     + 1
      MTUSED = XBLOCK + NB

C     ------------------------------------------
C     SET POINTERS TO MSIFA AND CHECK ITS LENGTH
C     ------------------------------------------
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)
      XSPLIT = XLNZ   + NSUPER + 1
      SPLIT  = XSPLIT + NSUPER
      MLUSED = SPLIT  + NBLOCK

C     ---------------------------------------
C     SET POINTERS TO SM AND CHECK ITS LENGTH
C     ---------------------------------------
      DIAG   = 1
      LNZ    = DIAG   + NEQ
      SMUSED = LNZ    + NZEROL - 1

C     ------------------------------------------
C     SET POINTERS TO RWORK AND CHECK ITS LENGTH
C     ------------------------------------------
      STACK  = 1
      RWUSED = STACK  + WSTACK - 1

      IF ( IWUSED .GT. NIWORK .OR. MTUSED .GT. NTREES .OR.
     $     MLUSED .GT. NMSIFA .OR. SMUSED .GT. NSM    .OR.
     $     RWUSED .GT. NRWORK                              ) THEN
         MSPAR(1) = -2
         IF ( LPU .GE. 0 ) THEN
            WRITE(LPU,'(/A)') ' *** ERROR from SPRFC2'
            IF ( IWUSED .GT. NIWORK ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'IWORK too small',NIWORK,IWUSED
            ENDIF
            IF ( MTUSED .GT. NTREES ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'MTREES too small',NTREES,MTUSED
            ENDIF
            IF ( MLUSED .GT. NMSIFA ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'MSIFA too small',NMSIFA,MLUSED
            ENDIF
            IF ( SMUSED .GT. NSM ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'SM too small',NSM,SMUSED
            ENDIF
            IF ( RWUSED .GT. NRWORK ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'RWORK too small',NRWORK,RWUSED
            ENDIF
         ENDIF
         IERR = -2
         CALL SPRER1 ( 36, 'SPRFC2', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

C     ==================================================================

C     ----------------------------------------------------
C     COMPUTE D AND C BY A CALL TO SPRF21 AND STORE IN SM.
C     ----------------------------------------------------
      CALL SPRF21 ( NSUPER, NEQ, NOFSUB, LFIRST, WSTACK, NZEROL,
     $              NBLOCK, MSPAR(57), NUMNEG, MTREES(XSUPER),
     $              MTREES(XSIBL), MTREES(SIBL), MTREES(XLINDX),
     $              MTREES(SUPSUP), MSIFA(LINDX), MSIFA(XSPLIT),
     $              MSIFA(SPLIT), MSIFA(XLNZ), IWORK(FIRST), SM(DIAG),
     $              SM(LNZ), TOL, RWORK(STACK), LIERR, SUPEQT,
     $              MVMULT, MMMULT )

      IF ( LIERR .NE. 0 ) THEN
         IF ( LIERR/10 .EQ. -3 ) THEN
            IERR = -3
         ELSE
            IERR = -4
         ENDIF
         CALL SPRER1 ( 40, 'SPRFC2', NUMNEG, SUPEQT, LIERR, LPU, IERR )
         MSPAR(27) = SUPEQT
         IF ( LIERR .EQ. -30 .AND. LPU.GE.0) THEN
            WRITE(LPU,600) TOL, TOL(2)/ABS(TOL(3))
         ENDIF
      ENDIF

      MSPAR(26) = NUMNEG

C     ==================================================================
#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRFC2'
#endif
      RETURN

  600 FORMAT(5X,'ZERO FACTOR THRESHOLD,   TOL(1) = ',1PE12.5,
     $     / 5X,'DIAGONAL ENTRY (PIVOT),   DIAG  = ',1PE12.5,
     $     / 5X,'PREVIOUS PIVOT,          ADIAG  = ',1PE12.5,
     $     / 5X,'DIAGONAL DECAY, DIAG/ABS(ADIAG) = ',1PE12.5,
     $          ' < TOL(1)' )

C     ------------------------------------------------------------------
C     END OF MODULE SPRFC2
      END


      SUBROUTINE   SPRF21   ( NSUPER, NEQ   , NOFSUB, LFIRST, WSTACK,
     $                        NZEROL, NBLOCK, LPAR57, NUMZRO, XSUPER,
     $                        XSIBL , SIBL  , XLINDX, SUPSUP, LINDX ,
     $                        XSPLIT, SPLIT , XLNZ  , FIRST , DIAG  ,
     $                        LNZ   , TOL   , STACK , ERROR , SUPEQT,
     $                        MVMULT, MMMULT  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , NOFSUB, LFIRST, WSTACK, NZEROL,
     $                  NBLOCK, LPAR57, NUMZRO, ERROR , SUPEQT
      INTEGER           XSUPER(NSUPER+1)      , XSIBL(NSUPER+1)       ,
     $                  SIBL(NSUPER)          , XLINDX(NSUPER+1)      ,
     $                  SUPSUP(NSUPER)        , LINDX(NOFSUB)         ,
     $                  XSPLIT(*)             , SPLIT(*)              ,
     $                  XLNZ(NSUPER+1)        , FIRST(LFIRST)

      DOUBLE PRECISION  LNZ(NZEROL)           , STACK(WSTACK)         ,
     $                  DIAG(NEQ)             , TOL(3)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRF21
C
C --- TO PERFORM A SUPERNODE MULTIFRONT COMPLETE or partial (dependent
*     on the contents of SUPSUP) FACTORIZATION STEP. THE NODES MUST BE
*     IN POSTORDER FOR THE INTERNAL STACK MANAGEMENT TO WORK.
*     THIS VERSION HAS NOT IMPLEMENTED THE LAST CHILD OVERLAP TECHNIQUE.
*     THE ROUTINE ASSUMES THAT THE MATRIX COEFFICIENTS ARE ASSEMBLED
*     INTO LNZ ON ENTRY TO THE ROUTINE.
C
C     THIS VERSION INCORPORATES THE USE OF BLAS3 ROUTINES IN THE
C     ELIMINATION STEPS, AND IT ALSO WORKS ON A COMPLETE FRONT MATRIX ON
C     TOP OF THE STACK.
C
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : DEC. 10, 1996 (ACD)
C                 CORRECTIONS OF LFIRST TEXT.
C                 APR. 09, 2003 (KMO)
C                 ADDED CHECK FOR ZERO PIVOTS BEFORE THE FACTORIZATION.
C                 AUG. 14, 2003 (KMO)
C                 ADDED INPUT PARAMETER LPAR57.
C                 JAN. 02, 2008 (KMO)
C                 DON'T CHECK FOR ZERO PIVOTS AMONG THE RETAINED DOFS.
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     NOFSUB : INTEGER
C              NUMBER OF INDICES USED TO REPRESENT THE STRUCTURE OF L.
C     LFIRST : INTEGER
C              LENGTH OF WORK ARRAY FIRST, MFRONT+1.
C              MFRONT IS THE ORDER OF THE GREATEST FRONT MATRIX.
C     WSTACK : INTEGER
C              WORKING STACK STORAGE REQUIREMENT FOR THE FACTORIZATION
C              STEP.
C     NZEROL : INTEGER
C              SIZE OF THE FACTOR MATRIX INCLUDING THE DIAGONAL.
C     NBLOCK : INTEGER
C              NUMBER OF PANELS FOR BETTER STACK USE.
C     LPAR57 : INTEGER
C              IF (LPAR57.LE.0) DO A FULL SOLVE EVEN IF THE EQUATION
C              SYSTEM IS PARTITIONED INTO INTERNAL AND RETAINED DOFS.
C     XSUPER : INTEGER(NSUPER+1)
C              THE FUNDAMENTAL SUPERNODE PARTITION.
C     (XSIBL,
C     SIBL)  : INTEGER(NSUPER+1), INTEGER(NSUPER-1)
C              THE CHILDREN OF EACH SUPERNODE.
C     (XLINDX,
C     LINDX) : INTEGER(NSUPER+1), INTEGER(NOFSUB)
C              THE NODE BLOCK NON-ZERO STRUCTURE OF L.
C     (XSPLIT,
C     SPLIT) : INTEGER(NSUPER), INTEGER(NBLOCK+1)
C              INFORMATION ABOUT THE PANELS FOR BETTER USE OF THE
C              CACHE.
C     XLNZ   : INTEGER(NSUPER+1)
C              POINTER TO START OF SUBMATRICES IN LNZ FOR EACH
C              SUPERNODE.
C
C     ON EXIT
C     -------
C
C     NUMZRO : INTEGER
C              THE NUMBER OF ZERO PIVOTS FOUND.
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              THE CHOLESKY SUPERNODE SUBMATRICES OF C STORED CON-
C              SEQUTIVELY ON DENSE FORMAT INCLUDED THE DIAGONAL.
C     TOL    : DOUBLE PRECISION(3)
C              TOL(1): FACTOR ZERO THRESHOLD.
C              TOL(2): MAX(A(I,I)/D(I))=MIN(D(I)/A(I,I))
C                      SET TO ZERO IF A(I,I) = ZERO.
C              TOL(3): TRACE(A)/MIN(D(I))
C                      SET TO ZERO IF D(I) = ZERO.
C     ERROR  : INTEGER
C              ERROR FLAG.
C              =   0 : NORMAL RETURN
C              = -10 : INSUFFICIENT STACK SPACE. OCCURRED FOR
C                      THE SUPERNODE WITH FIRST EQUATION SUPEQT.
C              = -20 : INSUFFICIENT SPACE FOR C. OCCURRED FOR
C                      THE SUPERNODE WITH FIRST EQUATION SUPEQT.
C              = -32 : PIVOT REDUCED TO .LE. ZERO IN EQUATION SUPEQT.
C              = -35 : ZERO PIVOTS IN NUMZRO EQUATIONS, FIRST ONE SUPEQT.
C              = -40 : TRIES TO MOVE THE UPDATE MATRIX BEYOND THE LOWER
C                      LIMIT FOR SUPERNODE WITH FIRST EQUATION SUPEQT.
C     SUPEQT : INTEGER
C              EQUATION FOR WHICH ERROR OCCURED.
C
C     WORKING ARRAYS
C     --------------
C
C     FIRST  : INTEGER(MFRONT+1)
C              LOCAL INDEX SETS FOR EACH CHILD, LENGTH MFRONT.
C              START FOR COLUMNS IN FRONT, LENGTH MFRONT+1.
C     STACK  : DOUBLE PRECISION(WSTACK)
C              STORAGE FOR THE MULTIFRONTAL STACK.
C
C     PROCEDURES
C     ----------
C
C     SPRF12
C     SPRF13
C     SPRF14
C     SPRF25
C     SPRF16
C     SPRF18
C     MVMULT
C     MMMULT
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           CHDLEN, CHDSZE, CHILD , CSKIP , CSTRT , FSIZE ,
     $                  FSTVAR, I     , IPJS  , IPKS  , IPLNZ , IPSTCK,
     $                  ISTART, ISTOP , J     , JS    , LENVAR, PSTRT ,
     $                  SUPSZE, SUPVAR, UPDSZE

      DOUBLE PRECISION  MDIAG , ONE   , ZERO
      PARAMETER       ( ONE = 1.0D0, ZERO = 0.0D0 )

      LOGICAL           PANEL

      EXTERNAL          SPRF12, SPRF13, SPRF14, SPRF25, SPRF16, SPRF18,
     $                  MVMULT, MMMULT

C     ==================================================================

      ERROR  = 0
      SUPEQT = 0
      NUMZRO = 0
      IPSTCK = WSTACK + 1
      IPLNZ  = 1
      PANEL  = NBLOCK .GT. 0
      TOL(2) = ZERO
      TOL(3) = ZERO

C     ==================================================================

C     ------------------------------------------
C     INITIALIZE DIAG AND CHECK FOR ZERO PIVOTS.
C     ------------------------------------------
C     KMO 02/01/08 -------------------------------
      IF (LPAR57.GT.0) CALL SPRF18 (NEQ,DIAG,ZERO)
C     --------------------------------------------
      DO 200 JS = 1, NSUPER
C        KMO 02/01/08 ----------------------------------
         IF (SUPSUP(JS).GT.0 .AND. LPAR57.GT.0) GOTO 200
C        -----------------------------------------------
         CHDLEN = XLINDX(JS+1) - XLINDX(JS)
         IPJS = XLNZ(JS)
         DO 100 J = XSUPER(JS), XSUPER(JS+1)-1
            DIAG(J) = LNZ(IPJS)
C           KMO 09/04/03 ------------
            IF (DIAG(J).EQ.ZERO) THEN
               NUMZRO = NUMZRO + 1
               IF (NUMZRO.EQ.1) THEN
                  SUPEQT = J
                  ERROR = -35
               ENDIF
            ENDIF
C           -------------------------
            TOL(3) = TOL(3) + DIAG(J)
            IPJS = IPJS + CHDLEN
            CHDLEN = CHDLEN - 1
  100    CONTINUE
  200 CONTINUE
      IF (ERROR .LT. 0) RETURN

C     ----------------------
C     FOR EACH SUPERNODE ...
C     ----------------------
      DO 800 JS = 1, NSUPER

C        -----------------------
C        SET UP CHARACTERISTICS.
C        -----------------------
         FSTVAR = XSUPER(JS)
         SUPVAR = XSUPER(JS+1) - FSTVAR
         LENVAR = XLINDX(JS+1) - XLINDX(JS)
         FSIZE  = (LENVAR*(LENVAR+1))/2
         UPDSZE = ((LENVAR-SUPVAR)*(LENVAR-SUPVAR+1))/2
         SUPSZE = FSIZE - UPDSZE

C        ---------------------
C        ALLOCATE STACK SPACE.
C        ---------------------
         IPJS = IPSTCK - FSIZE
         IF ( IPJS .LE. 0 ) THEN
            SUPEQT = FSTVAR
            ERROR  = -10
            WSTACK = WSTACK - IPJS + 1
            RETURN
         ENDIF

C        ------------------------------------------
C        COPY MATRIX COEFFICIENTS INTO FULLY SUMMED
C        PART OF THE FRONT MATRIX, AND INITIALIZE
C        THE UPDATE PART TO ZERO.
C        ------------------------------------------
         DO 300 I = 1, SUPSZE
            STACK(IPJS+I-1) = LNZ(IPLNZ+I-1)
  300    CONTINUE
         DO 400 I = 1, UPDSZE
            STACK(IPJS+SUPSZE+I-1) = ZERO
  400    CONTINUE

C        ----------------------------------
C        POINTERS TO CHILDREN LIST IN SIBL.
C        ----------------------------------
         ISTART = XSIBL(JS)
         ISTOP  = XSIBL(JS+1) - 1
         IF ( ISTART .LE. ISTOP ) THEN

C           ------------------
C           FOR EACH CHILD ...
C           ------------------
            PSTRT  = XLINDX(JS)
            DO 500 I = ISTOP, ISTART, -1

C              -------------------------------------
C              SET UP CHARACTERISTICS FOR THE CHILD.
C              -------------------------------------
               CHILD  = SIBL(I)
               IPKS   = IPSTCK
               CSKIP  = XSUPER(CHILD+1) - XSUPER(CHILD)
               CSTRT  = XLINDX(CHILD)
               CHDLEN = XLINDX(CHILD+1) - CSTRT - CSKIP
               CSTRT  = CSTRT + CSKIP
               CHDSZE = (CHDLEN*(CHDLEN+1))/2

C              ---------------------------
C              ASSEMBLE THE UPDATE MATRIX.
C              ---------------------------
               IF ( LENVAR .GT. CHDLEN ) THEN

C                 ----------------------------
C                 NOT A TRIVIAL ASSEMBLY.
C                 COMPUTE THE LOCAL INDEX SET.
C                 ----------------------------
                  CALL SPRF12 ( LENVAR, CHDLEN, LINDX(PSTRT),
     $                          LINDX(CSTRT), FIRST )
                  IF ( LENVAR .GT. SUPVAR ) THEN

C                    -------------------------------
C                    ASSEMBLE INTO BOTH FULLY SUMMED
C                    PART AND UPDATE PART OF MATRIX.
C                    -------------------------------
                     CALL SPRF13 ( LENVAR, SUPVAR, CHDLEN,
     $                             FIRST, STACK(IPKS), STACK(IPJS),
     $                             STACK(IPJS+SUPSZE) )
                  ELSE

C                    -----------------------
C                    ONLY FULLY SUMMED PART.
C                    -----------------------
                     CALL SPRF13 ( SUPVAR, SUPVAR, CHDLEN,
     $                             FIRST, STACK(IPKS), STACK(IPJS),
     $                             STACK(IPJS) )
                  ENDIF
               ELSE

C                 --------------------------------------------
C                 THE SIZE OF THE UPDATE MATRIX IS THE SAME
C                 AS THE PARENTS FRONT MATRIX, TRIVIAL UPDATE.
C                 N O T E THIS MUST BE TWO CALLS
C                         WHEN PARTITIONED MATRIX IS USED.
C                 --------------------------------------------
                  CALL SPRF14 ( FSIZE, STACK(IPKS), STACK(IPJS) )
               ENDIF

C              -------------------------------
C              DECREMENT STACK POINTER IPSTCK.
C              -------------------------------
               IPSTCK = IPKS + CHDSZE
  500       CONTINUE

         ENDIF

         IF ( SUPSUP(JS).LE.0 .OR. LPAR57.LE.0 ) THEN

C           -----------------------------------------------------
C           PERFORM SUPVAR STEPS OF ELIMINATION since this is not
*           a retained supernode.
C           -----------------------------------------------------
            IF ( PANEL ) THEN

               J = XSPLIT(JS)
               CALL SPRF25 ( LENVAR, SUPVAR, PANEL, SPLIT(J),
     $                       STACK(IPJS), DIAG(FSTVAR), TOL, ERROR,
     $                       FIRST, MVMULT, MMMULT )

            ELSE

               CALL SPRF25 ( LENVAR, SUPVAR, PANEL, FIRST,
     $                       STACK(IPJS), DIAG(FSTVAR), TOL, ERROR,
     $                       FIRST, MVMULT, MMMULT )

            ENDIF

            IF ( ERROR .LT. 0 ) THEN
               SUPEQT = FSTVAR - ERROR - 1
               ERROR  = -32
               RETURN
            ENDIF

         ENDIF

         IF ( SUPSZE .GT. 0 ) THEN

C           -----------------------------
C           ALLOCATE THE FACTOR SPACE ...
C           -----------------------------
            IF ( IPLNZ+SUPSZE-1 .GT. NZEROL ) THEN
               SUPEQT = FSTVAR
               ERROR  = -20
               RETURN
            ENDIF

C           ----------------------------
C           ... STORE THE FACTORED PART.
C           ----------------------------
            DO 600 I = 1, SUPSZE
               LNZ(IPLNZ+I-1) = STACK(IPJS+I-1)
  600       CONTINUE

C           --------------------------------
C           ... AND UPDATE POINTER INTO LNZ.
C           --------------------------------
            IPLNZ = IPLNZ + SUPSZE
         ENDIF

         IF ( UPDSZE .GT. 0 ) THEN

C           ---------------------------
C           NODE HAS UPDATE MATRIX ...
C              ... SHIFT IT TO THE TOP.
C           ---------------------------
            IF ( (IPJS+FSIZE-1) .LT. (IPSTCK-UPDSZE) ) THEN

C              ----------------------------------
C              ... NO OVERLAP, CALL COPY ROUTINE.
C              ----------------------------------
               CALL SPRF16 ( UPDSZE, STACK(IPJS+SUPSZE),
     $                       STACK(IPSTCK-UPDSZE) )
            ELSEIF ( (IPJS+FSIZE-1) .GE. (IPSTCK-UPDSZE) .AND.
     $               (IPJS+FSIZE)   .LE. IPSTCK                ) THEN

C              ------------------------------
C              ... OVERLAP, USE IN-LINE LOOP.
C              ------------------------------
               IF ( (IPJS+FSIZE) .EQ. IPSTCK ) THEN

C                 ----------------------------------
C                 NO ACTION SINCE THE UPDATE MATRIX
C                 IS ALREADY PLACED WHERE IT SHOULD.
C                 ----------------------------------
               ELSE

C                 -----------------------------
C                 MOVE UPDATE MATRIX BACKWARDS.
C                 -----------------------------
                  DO 700 I = UPDSZE, 1, -1
                     STACK(IPSTCK-UPDSZE+I-1) = STACK(IPJS+SUPSZE+I-1)
  700             CONTINUE
               ENDIF

            ELSE

C              --------------------------------
C              THIS IS AN ERROR, CANNOT HAPPEN.
C              DEBUGG EXIT.
C              --------------------------------
               SUPEQT = FSTVAR
               ERROR  = -40
               RETURN
            ENDIF

C           -----------------------------
C           ... UPDATE THE STACK POINTER.
C           -----------------------------
            IPSTCK = IPSTCK - UPDSZE
         ENDIF

  800 CONTINUE

C     -------------
C     FINALIZE TOL.
C     -------------
C     KMO 02/01/08 -------------------------------------------------
      MDIAG = ZERO
      DO 900 I = 1, NEQ
         IF (DIAG(I).GT.ZERO) THEN
            IF (MDIAG.EQ.ZERO .OR. MDIAG.GT.DIAG(I)) MDIAG = DIAG(I)
         ENDIF
  900 CONTINUE
C     --------------------------------------------------------------

      IF ( TOL(2) .NE. ZERO ) THEN
         TOL(2) = ONE/TOL(2)
      ELSE
         TOL(2) = ZERO
      ENDIF
      IF ( MDIAG .NE. ZERO ) THEN
         TOL(3) = TOL(3)/(MDIAG*MDIAG)
      ELSE
         TOL(3) = ZERO
      ENDIF

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRF21
      END


      SUBROUTINE   SPRF25   ( M     , N     , PANEL , SPLIT , X     ,
     $                        DIAG  , TOL   , ERROR , XPNT  , MVMULT,
     $                        MMMULT                                  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           M     , N     , ERROR
      INTEGER           SPLIT(*)              , XPNT(M+1)

      DOUBLE PRECISION  X(M*(M+1)/2)          , DIAG(N)               ,
     $                  TOL(3)

      LOGICAL           PANEL

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRF25
C
C --- TO ELIMINATE THE N COLUMNS OF THE CURRENT FRONT MATRIX AND TO UP-
C     DATE THE REMAINING M-N COLUMNS OF THE MATRIX.
C     THE VERSION IS A BLOCK-BLOCK UPDATE ROUTINE BASED ON BLAS3 METHOD.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     M      : INTEGER
C              ORDER OF THE CURRENT FRONT MATRIX.
C     N      : INTEGER
C              NUMBER OF FULLY SUMMED COLUMNS OF THE FRONT MATRIX, I.E.
C              THE NUMBER OF COLUMNS TO BE ELIMINATED IN THE CURRENT
C              STEP.
C     PANEL  : LOGICAL
C              OPTION FOR OPTIMIZED CACHE USE.
C              .EQ. .TRUE.  : ARRAY SPLIT CONTAINS THE BLOCK INFO FOR
C                             THE CURRENT SUPERNODE.
C              .EQ. .FALSE. : NO PANELLING AND THE ARRAY SPLIT IS NOT
C                             REFERENCED IN THE ROUTINE.
C     SPLIT  : INTEGER(*)
C              INFORMATION ABOUT PANELS FOR BETTER USE OF THE CACHE.
C
C     UPDATED PARAMETERS
C     ------------------
C
C     X      : DOUBLE PRECISION(M*(M+1)/2)
C              ON INPUT IT CONTAINS THE NEWLY FORMED FRONT MATRIX,
C              AND ON RETURN IT CONTAINS IN THE FIRST N COLUMNS THE
C              NEW FACTORS AND IN ITS LAST M-N COLUMNS THE UPDATE
C              MATRIX TO ITS PARENT.
C     DIAG   : DOUBLE PRECISION(N)
C              THE DIAGONAL ENTRIES OF THE LDL' FACTOR MATRIX.
C     TOL    : DOUBLE PRECISION(3)
C              TOL(2) UPDATED.
C
C     ON EXIT
C     -------
C
C     ERROR  : INTEGER
C              ERROR FLAG.
C              = 0: NORMAL RETURN.
C              < 0: PIVOT REDUCED TO .LE. ZERO IN COLUMN ABS(ERROR).
C
C     WORKING ARRAYS
C     --------------
C
C     XPNT   : INTEGER(M+1)
C              POINTERS TO THE FIRST ENTRY IN EACH COLUMN OF THE FRONT
C              MATRIX.
C
C     PROCEDURES
C     ----------
C
C     MVMULT : MATRIX-VECTOR MULTIPLY.
C     MMMULT : MATRIX-MATRIX MULTIPLY.
C     SPRF27 : PARTIAL CC' FACTORIZATION.
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           FSTBLK, FSTCOL, I     , J     , JBLK  , JPNT  ,
     $                  MM    , NN    , NXTCOL, OFFSET, Q

      EXTERNAL          MVMULT, MMMULT, SPRF27

      OFFSET(I,J) = ( I*J ) - ( J*( J-1 )/2 )

C     ==================================================================

C     -------------------------
C     QUICK RETURN IF POSSIBLE.
C     -------------------------
      IF ( M .LE. 0 .OR. N .LE. 0 ) RETURN

C     -------------
C     COMPUTE XPNT.
C     -------------
CCDIR$ IVDEP
      DO 100 I = 1, M
         XPNT(I) = OFFSET(M,I-1) + 1
  100 CONTINUE
      XPNT(M+1) = XPNT(M) + 1
C *** The reason why IWORK in this routine set must have length
C *** MFRONT+1.

C     ==================================================================

C     ------------------------
C     START THE FACTORIZATION.
C     ------------------------
      JBLK   = 1
      FSTCOL = 1
      MM     = M
      JPNT   = XPNT(FSTCOL)
      ERROR  = 0
      IF ( PANEL ) THEN

C        ------------------------
C        ARRAY SPLIT GIVES THE
C        BLOCK PARTITION FOR
C        BETTER USE OF THE CACHE.
C        ------------------------
         FSTBLK = SPLIT(JBLK)
  200    IF ( FSTCOL .LE. N ) THEN
            JBLK   = JBLK   + 1
            NN     = SPLIT(JBLK) - FSTBLK
            FSTBLK = SPLIT(JBLK)

C           ---------------------------------------------
C           FACTORIZE THE FIRST NN COLUMNS OF THE MATRIX.
C           ---------------------------------------------
            CALL SPRF27 ( MM, NN, XPNT(FSTCOL), X, M, DIAG(FSTCOL),
     $                    TOL, ERROR, MVMULT )
            IF ( ERROR .NE. 0 ) RETURN
            NXTCOL = FSTCOL + NN
            Q      = M      - NXTCOL + 1
            MM     = MM     - NN
            JPNT   = XPNT(NXTCOL)
            IF ( Q .GT. 0 ) THEN

C              -------------------------------------
C              THEN UPDATE THE TRAILING SUBMATRIX
C              BY THE RANK NN UPDATE COMPUTED ABOVE.
C              -------------------------------------
               CALL MMMULT ( MM, NN, Q, XPNT(FSTCOL), X,
     $                       X(JPNT), MM )
            ENDIF

C           -----------------
C           UPDATE FSTCOL AND
C           GO TO NEXT PANEL.
C           -----------------
            FSTCOL = NXTCOL
            GOTO 200
         ENDIF
      ELSE

C        --------------------------------------------
C        FACTORIZE THE FIRST N COLUMNS OF THE MATRIX.
C        --------------------------------------------
         CALL SPRF27 ( M, N, XPNT, X, M, DIAG, TOL, ERROR, MVMULT )
         IF ( ERROR .NE. 0 ) RETURN
         NXTCOL = FSTCOL + N
         Q      = M      - NXTCOL + 1
         MM     = MM     - N
         JPNT   = XPNT(NXTCOL)
         IF ( Q .GT. 0 ) THEN

C           ------------------------------------
C           THEN UPDATE THE TRAILING SUBMATRIX
C           BY THE RANK N UPDATE COMPUTED ABOVE.
C           ------------------------------------
            CALL MMMULT ( MM, N, Q, XPNT, X, X(JPNT), MM )
         ENDIF

      ENDIF

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRF25
      END


      SUBROUTINE   SPRF27   ( M, N, XPNT, X, LDX, LZDIAG, TOL, IERR,
     $                        MVMULT )

      IMPLICIT NONE

      INTEGER             M, N, XPNT(N+1), LDX, IERR

      DOUBLE PRECISION    X(LDX*(LDX+1)/2), LZDIAG(N), TOL(3)

C***********************************************************************
C
C     CREATED:   JAN. 11, 1994 (ACD)
C     REVISIONS: MAR. 26, 2003 (KMO)
C                RETURN DIAG AND ADIAG THROUGH TOL WHEN SINGULAR MATRIX.
C
C     PURPOSE - THIS ROUTINE PERFORMS CCT
C               FACTORIZATION ON THE COLUMNS OF A SUPERNODE
C               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS
C               EXTERNAL TO THE SUPERNODE.
C
C     INPUT PARAMETERS -
C        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN).
C        N      - NUMBER OF COLUMNS IN THE SUPERNODE.
C        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END
C                 OF THE J-TH COLUMN OF THE SUPERNODE.
C        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO
C                 BE FACTORED.
C        LDX    - THE LEADING DIMENSION OF THE FRONT MATRIX
C                 OF WHICH THE LOWER TRIANGLE INCLUDING THE
C                 THE DIAGONAL IS STORED.
C
C     OUTPUT PARAMETERS -
C        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF
C                 THE SUPERNODE.
C        IERR   - UNCHANGED IF THERE IS NO ERROR.
C                 = -EQUATION WHERE DIAGONAL ENTRY .LE. 0 IS
C                    ENCOUNTERED.
C
C        MVMULT - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY.
C
C***********************************************************************

      INTEGER             I     , JPNT  , JCOL  , MM

      DOUBLE PRECISION    ADIAG , DIAG  , ONE
      PARAMETER         ( ONE = 1.0D0 )

      INTRINSIC           ABS   , MAX   , SQRT

      EXTERNAL            MVMULT

C***********************************************************************

      MM = M
      JPNT = XPNT(1)
      DO 200 JCOL = 1, N

C        ----------------------------------
C        UPDATE JCOL WITH PREVIOUS COLUMNS.
C        ----------------------------------
         IF ( JCOL .GT. 1 ) THEN
            CALL MVMULT ( MM, JCOL-1, X(JPNT), XPNT, X )
         ENDIF

C        -----------------------------------
C        CHECK AND STORE THE DIAGONAL ENTRY.
C        -----------------------------------
         ADIAG = LZDIAG(JCOL)
         DIAG  = X(JPNT)
         IF ( DIAG .LE. TOL(1)*ABS(ADIAG) ) THEN
            IERR   = -JCOL
            TOL(2) = DIAG
            TOL(3) = ADIAG
            RETURN
         ENDIF

C        --------------
C        UPDATE TOL(2).
C        --------------
         TOL(2) = MAX(TOL(2),ADIAG/DIAG)

C        ------------------------------------------
C        TAKE THE SQUARE-ROOT OF THE DIAGONAL ENTRY
C        AND UPDATE LZDIAG AND X WITH THE RESULT.
C        ------------------------------------------
         DIAG = SQRT(DIAG)
         LZDIAG(JCOL) = DIAG
         X(JPNT) = DIAG

C        ----------------------------------------------------
C        SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY.
C        ----------------------------------------------------
         DIAG = ONE / DIAG
         MM   = MM - 1
         JPNT = JPNT + 1
         DO 100 I = 1, MM
            X(JPNT) = X(JPNT) * DIAG
            JPNT = JPNT + 1
  100    CONTINUE

  200 CONTINUE

      RETURN
      END
