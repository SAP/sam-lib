C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C **********************************************************************
C This module contains slightly modified versions of some subroutines
C in the ASM module, not (yet) part of the official SAM library.

      SUBROUTINE ADDEM1 (EM,TTCC,MPAR,MADOF,MEQN,
     +                   MPMNPC,MMNPC,MPMCEQ,MMCEQ,MSKY,
     +                   IEL,NEDOF,LPU,NSV,
     +                   SM,SV,MEEN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDEM1                   GROUP 7 / PUBLIC
C
C     Identical to ADDEM, except that it uses ELMEQ2 instead of ELMEQ
C     when MPAR(18).GT.0
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ,ELMEQ2,NOPRED,ASMERR   (SAM-7)
C                                   ABS                (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-01-14 / 1.0   (ADDEM)
C                       02-07-20 / 2.0   K.M.Okstad
C                       07-09-05 / 2.1   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF,NSV
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*),MSKY(*)
      DOUBLE PRECISION  TTCC(*),EM(NEDOF,NEDOF),SM(*),SV(*)
C
      INTEGER           I,ICEQ,IEQ,II,IP,J,JCEQ,JEQ,JJ,JP,JV,K,KEQ,KP,
     +                  NECEQ,NELDOF,NENOD,NEPRD,NESLV,NMSTI,NMSTJ,
     +                  nndof,NV,NOPRED
      DOUBLE PRECISION  C0,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR,ELMEQ,ELMEQ2,NOPRED
C ----------------------------------------------------------------------
      IERR = 0
      IF (IEL.GT.0 .AND. IEL.LE.MPAR(2)) GO TO 10
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
C
   10 NV = ABS(NSV)
      IF (NV.LT.2)             GO TO 20
      IF (NV.EQ.MPAR(9))       GO TO 20
      CALL ASMERR (17,IEL,NSV,LPU,IERR)
      GO TO 1000
C
   20 NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
      IP    = MPMNPC(IEL)
      IF (MPAR(18).GT.0)       GO TO 25
C     Assumes that NNDOF = MADOF(I+I) - MADOF(I) for node I
      CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +            MEEN,NELDOF,NESLV,NEPRD)
      GO TO 30
   25 nndof = NEDOF/nenod
C     Assumes that NEDOF = NNDOF*NENOD, NNDOF is constant
      CALL ELMEQ2 (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,nndof,
     +             MEEN,NELDOF,NESLV,NEPRD)
   30 IF (NELDOF.EQ.NEDOF)     GO TO 40
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
   40 NECEQ = NESLV+NEPRD
      IF (NSV.LT.0)            GO TO 220
C ----------------------------------------------------------------------
C  ADD ELEMENTS CORRESPONDING TO FREE  D O F S  IN  EM  INTO  SM
C ----------------------------------------------------------------------
      DO 200 J=1,NEDOF
         JEQ = MEEN(J)
         IF (JEQ.LT.1)         GO TO 200
         DO 150 I=1,J
            IEQ = MEEN(I)
            IF (IEQ.LT.1)      GO TO 150
            K   = MSKY(IEQ)
            IF (IEQ.LT.JEQ)    K = MSKY(JEQ)-JEQ+IEQ
            IF (IEQ.GT.JEQ)    K = K-IEQ+JEQ
            SM(K) = SM(K) + EM(I,J)
  150    CONTINUE
  200 CONTINUE
C
  220 IF (NECEQ.EQ.0)                  GO TO 1000
      IF (NSV.LT.-1 .AND. NEPRD.EQ.0)  GO TO 1000
C ----------------------------------------------------------------------
C  ADD (APPROPRIATELY WEIGHTED) ELEMENTS CORRESPONDING TO CONSTRAINED
C  (DEPENDENT AND PRESCRIBED)  D O F S  IN  EM  INTO  SM  AND/OR  SV
C ----------------------------------------------------------------------
      DO 600 J=1,NEDOF
         IF (MEEN(J).GE.0)     GO TO 600
         JCEQ =-MEEN(J)
         JP   = MPMCEQ(JCEQ)
         IF (NSV.EQ.0)         GO TO 320
         C0   = TTCC(JP)
         IF (C0.EQ.ZERO)       GO TO 310
         KP   = 0
         IF (NV.EQ.1)          GO TO 250
         JV   = NOPRED(MPMCEQ,JCEQ)
         IF (JV.EQ.0)          GO TO 310
         KP   = (JV-1)*MPAR(11)
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SV  (R-H SIDE)
  250    DO 300 I=1,NEDOF
            IEQ = MEEN(I)
            IF (IEQ.EQ.0)      GO TO 300
            IF (IEQ.LT.0)      GO TO 260
            KEQ = KP+IEQ
            SV(KEQ) = SV(KEQ) - C0*EM(I,J)
            GO TO 300
  260       ICEQ =-IEQ
            nmsti = mpmceq(iceq+1) - mpmceq(iceq) - 1
            IF (NMSTI.EQ.0)    GO TO 300
            IP = MPMCEQ(ICEQ)
            DO 280 II=1,NMSTI
               IP  = IP+1
               if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) go to 280
               keq = kp+meqn(mmceq(ip))
               SV(KEQ) = SV(KEQ) - C0*TTCC(IP)*EM(I,J)
  280       CONTINUE
  300    CONTINUE
C
  310    IF (NSV.LT.0)         GO TO 600
  320    nmstj = mpmceq(jceq+1) - mpmceq(jceq) - 1
         IF (NMSTJ.EQ.0)       GO TO 600
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SM
         DO 500 JJ=1,NMSTJ
            JP  = JP+1
            if (mmceq(jp).le.0 .or. ttcc(jp).eq.zero) go to 500
            jeq = meqn(mmceq(jp))
            DO 400 I=1,NEDOF
               IEQ = MEEN(I)
               IF (IEQ.EQ.0)   GO TO 400
               IF (IEQ.LT.0)   GO TO 360
               IF (IEQ.LT.JEQ) GO TO 330
               IF (IEQ.GT.JEQ) GO TO 340
               K     = MSKY(IEQ)
               SM(K) = SM(K) + (TTCC(JP)+TTCC(JP))*EM(I,J)
               GO TO 400
  330          K     = MSKY(JEQ)-JEQ+IEQ
               GO TO 350
  340          K     = MSKY(IEQ)-IEQ+JEQ
  350          SM(K) = SM(K) + TTCC(JP)*EM(I,J)
               GO TO 400
C
  360          ICEQ  =-IEQ
               nmsti = mpmceq(iceq+1) - mpmceq(iceq) - 1
               IF (NMSTI.EQ.0) GO TO 400
               IP    = MPMCEQ(ICEQ)
               DO 380 II=1,NMSTI
                  IP  = IP+1
                  if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) go to 380
                  ieq = meqn(mmceq(ip))
                  IF (IEQ.GT.JEQ) GO TO 380
                  K   = MSKY(JEQ)-JEQ+IEQ
                  SM(K) = SM(K) + TTCC(IP)*TTCC(JP)*EM(I,J)
  380          CONTINUE
  400       CONTINUE
  500    CONTINUE
  600 CONTINUE
C ----------------------------------------------------------------------
 1000 RETURN
      END
      SUBROUTINE ADDEM2 (EM,TTCC,MPAR,MADOF,MEQN,
     +                   MPMNPC,MMNPC,MPMCEQ,MMCEQ,
     +                   IEL,NEDOF,NEQ,LPU,NSV,
     +                   SM,SV,MEEN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDEM2                   GROUP 7 / PUBLIC
C
C     Identical to ADDEM1, except that the system matrix is dense.
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ,ELMEQ2,NOPRED,ASMERR   (SAM-7)
C                                   ABS                (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-01-14 / 1.0   (ADDEM)
C                       04-01-15 / 2.0   K.M.Okstad
C                       07-09-05 / 2.1   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF,NEQ,NSV
      INTEGER           MADOF(*),MMCEQ(*),MEEN(NEDOF),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*)
      DOUBLE PRECISION  TTCC(*),EM(NEDOF,NEDOF),SM(NEQ,NEQ),SV(NEQ)
C
      INTEGER           I,ICEQ,IEQ,II,IP,J,JCEQ,JEQ,JJ,JP,JV,KEQ,KP,
     +                  NECEQ,NELDOF,NENOD,NEPRD,NESLV,NMSTI,NMSTJ,
     +                  nndof,NV,NOPRED
      DOUBLE PRECISION  C0,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR,ELMEQ,ELMEQ2,NOPRED
C ----------------------------------------------------------------------
      IERR = 0
      IF (IEL.GT.0 .AND. IEL.LE.MPAR(2))  GO TO 10
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
C
   10 NV = ABS(NSV)
      IF (NV.LT.2)             GO TO 20
      IF (NV.EQ.MPAR(9))       GO TO 20
      CALL ASMERR (17,IEL,NSV,LPU,IERR)
      GO TO 1000
C
   20 NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
      IP    = MPMNPC(IEL)
      IF (MPAR(18).GT.0)       GO TO 25
C     Assumes that NNDOF = MADOF(I+I) - MADOF(I) for node I
      CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +            MEEN,NELDOF,NESLV,NEPRD)
      GO TO 30
   25 nndof = NEDOF/nenod
C     Assumes that NEDOF = NNDOF*NENOD, NNDOF is constant
      CALL ELMEQ2 (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,nndof,
     +             MEEN,NELDOF,NESLV,NEPRD)
   30 IF (NELDOF.EQ.NEDOF)     GO TO 40
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
   40 NECEQ = NESLV+NEPRD
      IF (NSV.LT.0)            GO TO 220
C ----------------------------------------------------------------------
C  ADD ELEMENTS CORRESPONDING TO FREE  D O F S  IN  EM  INTO  SM
C ----------------------------------------------------------------------
      DO 200 J=1,NEDOF
         JEQ = MEEN(J)
         IF (JEQ.LT.1)         GO TO 200
         DO 150 I=1,J
            IEQ = MEEN(I)
            IF (IEQ.LT.1)      GO TO 150
            SM(ieq,jeq) = SM(ieq,jeq) + EM(I,J)
            IF (IEQ.EQ.JEQ)    GO TO 150
            SM(jeq,ieq) = SM(jeq,ieq) + EM(J,I)
  150    CONTINUE
  200 CONTINUE
C
  220 IF (NECEQ.EQ.0)                  GO TO 1000
      IF (NSV.LT.-1 .AND. NEPRD.EQ.0)  GO TO 1000
C ----------------------------------------------------------------------
C  ADD (APPROPRIATELY WEIGHTED) ELEMENTS CORRESPONDING TO CONSTRAINED
C  (DEPENDENT AND PRESCRIBED)  D O F S  IN  EM  INTO  SM  AND/OR  SV
C ----------------------------------------------------------------------
      DO 600 J=1,NEDOF
         IF (MEEN(J).GE.0)     GO TO 600
         JCEQ =-MEEN(J)
         JP   = MPMCEQ(JCEQ)
         IF (NSV.EQ.0)         GO TO 320
         C0   = TTCC(JP)
         IF (C0.EQ.ZERO)       GO TO 310
         KP   = 0
         IF (NV.EQ.1)          GO TO 250
         JV   = NOPRED(MPMCEQ,JCEQ)
         IF (JV.EQ.0)          GO TO 310
         KP   = (JV-1)*MPAR(11)
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SV  (R-H SIDE)
  250    DO 300 I=1,NEDOF
            IEQ = MEEN(I)
            IF (IEQ.EQ.0)      GO TO 300
            IF (IEQ.LT.0)      GO TO 260
            KEQ = KP+IEQ
            SV(KEQ) = SV(KEQ) - C0*EM(I,J)
            GO TO 300
  260       ICEQ =-IEQ
            nmsti = mpmceq(iceq+1) - mpmceq(iceq) - 1
            IF (NMSTI.EQ.0)    GO TO 300
            IP = MPMCEQ(ICEQ)
            DO 280 II=1,NMSTI
               IP  = IP+1
               if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) goto 280
               KEQ = KP+meqn(mmceq(ip))
               SV(KEQ) = SV(KEQ) - C0*TTCC(IP)*EM(I,J)
  280       CONTINUE
  300    CONTINUE
C
  310    IF (NSV.LT.0)         GO TO 600
  320    nmstj = mpmceq(jceq+1) - mpmceq(jceq) - 1
         IF (NMSTJ.EQ.0)       GO TO 600
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SM
         DO 500 JJ=1,NMSTJ
            JP  = JP+1
            if (mmceq(jp).le.0 .or. ttcc(jp).eq.zero) go to 500
            jeq = meqn(mmceq(jp))
            DO 400 I=1,NEDOF
               IEQ = MEEN(I)
               IF (IEQ.EQ.0)   GO TO 400
               IF (IEQ.LT.0)   GO TO 360
               SM(ieq,jeq) = SM(ieq,jeq) + TTCC(JP)*EM(I,J)
               SM(jeq,ieq) = SM(jeq,ieq) + TTCC(JP)*EM(J,I)
               GO TO 400
C
  360          ICEQ  =-IEQ
               nmsti = mpmceq(iceq+1) - mpmceq(iceq) - 1
               IF (NMSTI.EQ.0) GO TO 400
               IP    = MPMCEQ(ICEQ)
               DO 380 II=1,NMSTI
                  IP  = IP+1
                  if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) go to 380
                  ieq = meqn(mmceq(ip))
                  SM(ieq,jeq) = SM(ieq,jeq) + TTCC(IP)*TTCC(JP)*EM(I,J)
  380          CONTINUE
  400       CONTINUE
  500    CONTINUE
  600 CONTINUE
C ----------------------------------------------------------------------
 1000 RETURN
      END
      SUBROUTINE ADDED (EM,TTCC,MPAR,MADOF,MEQN,
     +                  MPMNPC,MMNPC,MPMCEQ,MMCEQ,
     +                  IEL,NEDOF,LPU,NSV,
     +                  SM,SV,MEEN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDED                    GROUP 7 / PUBLIC
C
C     T A S K :  TO ADD/SUBTRACT THE ELEMENTS CORRESPONDING TO FREE
C     AND CONSTRAINED  D O F S  OF THE ELEMENT diagonal MATRIX  EM
C     (OF ELEMENT NO. IEL) TO/FROM THE CURRENT CONTENT OF THE
C     APPROPRIATE ELEMENTS OF
C       MATRIX     SM  IF  NSV.GE.0  (ADDITION WITHOUT/WITH WEIGHTING)
C       VECTOR(S)  SV  IF  NSV.NE.0  (SUBTRACTION WITH WEIGHTING)
C     SM IS A diagonal SYSTEM MATRIX STORED IN a vector, AND
C     SV IS A SYSTEM VECTOR (ABS(NSV)=1) OR A SET OF SYSTEM VECTORS
C     (ABS(NSV)=NPDOF) CONTAINING THE CONTRIBUTIONS TO THE RIGHT-HAND
C     SIDE(S) FROM PRESCRIBED AND (POSSIBLY) DEPENDENT D O F S.
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ,ELMEQ2,NOPRED,ASMERR   (SAM-7)
C                                   ABS                (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Knut Morten Okstad (based on ADDEM)
C     DATE/VERSION  :   01-08-27 / 1.0
C                       02-07-20 / 1.1
C                       07-09-05 / 1.2
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF,NSV
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*)
      DOUBLE PRECISION  TTCC(*),EM(NEDOF),SM(*),SV(*)
C
      INTEGER           II,IP,J,JCEQ,JDOF,JEQ,JJ,JP,JV,KP,
     +                  NECEQ,NELDOF,NENOD,NEPRD,NESLV,NMSTI,NMSTJ,
     +                  nndof,NV,NOPRED
      DOUBLE PRECISION  C0,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR,ELMEQ,ELMEQ2,NOPRED
C ----------------------------------------------------------------------
      IERR = 0
      IF (IEL.GT.0 .AND. IEL.LE.MPAR(2))  GO TO 10
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
C
   10 NV = ABS(NSV)
      IF (NV.LT.2)             GO TO 20
      IF (NV.EQ.MPAR(9))       GO TO 20
      CALL ASMERR (17,IEL,NSV,LPU,IERR)
      GO TO 1000
C
   20 NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
      IP    = MPMNPC(IEL)
      IF (MPAR(18).GT.0)       GO TO 25
C     Assumes that NNDOF = MADOF(I+I) - MADOF(I) for node I
      CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +            MEEN,NELDOF,NESLV,NEPRD)
      GO TO 30
   25 nndof = NEDOF/nenod
C     Assumes that NEDOF = NNDOF*NENOD, NNDOF is constant
      CALL ELMEQ2 (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,nndof,
     +             MEEN,NELDOF,NESLV,NEPRD)
   30 IF (NELDOF.EQ.NEDOF)     GO TO 40
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
   40 NECEQ = NESLV+NEPRD
      IF (NSV.LT.0)            GO TO 220
C ----------------------------------------------------------------------
C  ADD ELEMENTS CORRESPONDING TO FREE  D O F S  IN  EM  INTO  SM
C ----------------------------------------------------------------------
      DO 200 J=1,NEDOF
         JEQ = MEEN(J)
         IF (JEQ.GT.0)         SM(JEQ) = SM(JEQ) + EM(J)
  200 CONTINUE
C
  220 IF (NECEQ.EQ.0)                  GO TO 1000
      IF (NSV.LT.-1 .AND. NEPRD.EQ.0)  GO TO 1000
C ----------------------------------------------------------------------
C  ADD (APPROPRIATELY WEIGHTED) ELEMENTS CORRESPONDING TO CONSTRAINED
C  (DEPENDENT AND PRESCRIBED)  D O F S  IN  EM  INTO  SM  AND/OR  SV
C ----------------------------------------------------------------------
      DO 600 J=1,NEDOF
         JCEQ =-MEEN(J)
         IF (JCEQ.LE.0)        GO TO 600
         JP   = MPMCEQ(JCEQ)
         IF (NSV.EQ.0)         GO TO 320
         C0   = TTCC(JP)
         IF (C0.EQ.ZERO)       GO TO 310
         KP   = 0
         IF (NV.EQ.1)          GO TO 250
         JV   = NOPRED(MPMCEQ,JCEQ)
         IF (JV.EQ.0)          GO TO 310
         KP   = (JV-1)*MPAR(11)
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SV  (R-H SIDE)
  250    nmsti = mpmceq(jceq+1) - mpmceq(jceq) - 1
         IF (NMSTI.EQ.0)       GO TO 310
         DO 280 JJ=1,NMSTI
            JP  = JP+1
            if (mmceq(jp).le.0 .or. ttcc(jp).eq.zero) go to 280
            JEQ = KP+meqn(mmceq(jp))
            SV(JEQ) = SV(JEQ) - C0*TTCC(JP)*EM(J)
  280    CONTINUE
C
  310    IF (NSV.LT.0)         GO TO 600
  320    nmstj = mpmceq(jceq+1) - mpmceq(jceq) - 1
         IF (NMSTJ.EQ.0)       GO TO 600
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SM
         JP = MPMCEQ(JCEQ)
         DO 500 JJ=1,NMSTJ
            JP   = JP+1
            jdof = mmceq(jp)
            if (jdof.le.0 .or. ttcc(jp).eq.zero) go to 500
            jeq  = meqn(jdof)
C                                               ** Also add off-diagonal
C                                                  terms to the diagonal
            IP = MPMCEQ(JCEQ)
            DO 380 II=1,NMSTJ
               IP = IP+1
               if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) go to 380
               SM(JEQ) = SM(JEQ) + TTCC(IP)*TTCC(JP)*EM(J)
  380       CONTINUE
  500    CONTINUE
  600 CONTINUE
C ----------------------------------------------------------------------
 1000 RETURN
      END
      SUBROUTINE ADDEV1 (EV,TTCC,MPAR,MADOF,MEQN,
     +                   MPMNPC,MMNPC,MPMCEQ,MMCEQ,MPREAC,
     +                   IEL,NELDOF,IADDRF,LPU,SV,REAC,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDEV1                    GROUP 7 / PUBLIC
C
C     T A S K :  TO ADD THE ELEMENTS CORRESPONDING TO FREE, SUPPRESSED
C     AND DEPENDENT  D O F S  OF THE ELEMENT VECTOR  EV  (ELEMENT NO.
C     IEL) TO THE CURRENT CONTENT OF THE APPROPRIATE ELEMENTS OF THE
C     SYSTEM VECTOR  SV  AND REACTION FORCE VECTOR  REAC.
C     IF  IADDRF  IS NEGATIVE, THE REACTION FORCE CONTRIBUTIONS ARE
C     SUBTRACTED FROM  REAC  INSTEAD.
C
C     ROUTINES CALLED/REFERENCED :  ADDNV1 AND ASMERR  (SAM-7)
C
C     PROGRAMMED BY : Knut Morten Okstad
C     DATE/VERSION  : 05-10-11 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IADDRF,IEL,IERR,LPU,NELDOF
      INTEGER           MADOF(*),MMCEQ(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*),MPREAC(*)
      DOUBLE PRECISION  TTCC(*),EV(NELDOF),SV(*),REAC(*)
C
      INTEGER           I,INODE,NEDOF,NENOD,NNDOF
C
      EXTERNAL          ADDNV1, ASMERR
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (IEL .LT. 1 .OR. IEL .GT. MPAR(2)) THEN
         CALL ASMERR (12,IEL,IEL,LPU,IERR)
         RETURN
      END IF
C
      NEDOF = 0
      NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
C
      DO 100 I = 1, NENOD
         INODE = MMNPC(MPMNPC(IEL)+I-1)
C
         IF (NEDOF .LT. NELDOF) then
            CALL ADDNV1 (EV(NEDOF+1),TTCC,MPAR,MADOF,
     +                   MEQN,MPMCEQ,MMCEQ,MPREAC,
     +                   INODE,IADDRF,LPU,SV,REAC,NNDOF)
         ELSE IF (INODE .GE. 1 .AND. INODE .LE. MPAR(1)) THEN
            NNDOF = MADOF(INODE+1) - MADOF(INODE)
         ELSE
            CALL ASMERR (18,INODE,INODE,LPU,NNDOF)
         END IF
         IF (NNDOF .LT. 0) THEN
            IERR = NNDOF
            RETURN
         END IF
C
         NEDOF = NEDOF + NNDOF
  100 CONTINUE
C
      IF (NELDOF .NE. NEDOF) THEN
         CALL ASMERR (12,IEL,IEL,LPU,IERR)
         WRITE(LPU,"(5X,'NELDOF, NEDOF =',2I6/)") NELDOF, NEDOF
      END IF
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE ADDEV2 (EV,TTCC,MPAR,MADOF,MEQN,
     +                   MPMNPC,MMNPC,MPMCEQ,MMCEQ,
     +                   IEL,NEDOF,LPU,SV,MEEN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDEV2                   GROUP 7 / PUBLIC
C
C     Identical to ADDEV, except that it uses ELMEQ2 instead of ELMEQ,
C     when MPAR(18).GT.0
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ,ELMEQ2 AND ASMERR  (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-04-01 / 1.0   (ADDEV)
C                       02-20-07 / 2.0   K.M.Okstad
C                       07-09-05 / 2.1   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*)
      DOUBLE PRECISION  TTCC(*),EV(NEDOF),SV(*)
C
      INTEGER           I,ICEQ,IEQ,II,IP,NELDOF,NENOD,NEPRD,NESLV,NMST,
     +                  nndof
C
      EXTERNAL          ASMERR,ELMEQ,ELMEQ2
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (IEL.GT.0 .AND. IEL.LE.MPAR(2))  GO TO 10
      CALL ASMERR (12,IEL,IEL,LPU,IERR)
      GO TO 100
C
   10 NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
      IP    = MPMNPC(IEL)
      IF (MPAR(18).GT.0)    GO TO 15
C     Assumes that NNDOF = MADOF(I+I) - MADOF(I) for node I
      CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +            MEEN,NELDOF,NESLV,NEPRD)
      GO TO 18
   15 nndof = NEDOF/nenod
C     Assumes that NEDOF = NNDOF*NENOD, NNDOF is constant
      CALL ELMEQ2 (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,nndof,
     +             MEEN,NELDOF,NESLV,NEPRD)
   18 IF (NELDOF.EQ.NEDOF)  GO TO 20
      CALL ASMERR (12,IEL,IEL,LPU,IERR)
      GO TO 100
C
   20 DO 50 I=1,NEDOF
         IEQ = MEEN(I)
         IF (IEQ.EQ.0)      GO TO 50
         IF (IEQ.LT.0)      GO TO 30
         SV(IEQ) = SV(IEQ) + EV(I)
         GO TO 50
C
   30    ICEQ =-IEQ
         nmst = mpmceq(iceq+1) - mpmceq(iceq) - 1
         IF (NMST.EQ.0)     GO TO 50
         IP   = MPMCEQ(ICEQ)
         DO 40 II=1,NMST
            IP  = IP+1
            if (mmceq(ip).le.0) go to 40
            IEQ = meqn(mmceq(ip))
            SV(IEQ) = SV(IEQ) + TTCC(IP)*EV(I)
   40    CONTINUE
   50 CONTINUE
C ----------------------------------------------------------------------
  100 RETURN
      END
      SUBROUTINE ADDNV1 (VNOD,TTCC,MPAR,MADOF,
     +                   MEQN,MPMCEQ,MMCEQ,MPREAC,
     +                   INOD,IADDRF,LPU,SV,REAC,NNDOF)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDNV1                    GROUP 7 / PUBLIC
C
C     TASK :  To add the elements corresponding to free, suppressed and
C             dependent  d o f s  of the vector VNOD, which contains
C             one element per  d o f  of node INOD, to the current
C             content of the appropriate elements of the system vector
C             SV and/or the reaction force vector REAC.
C             The elements of VNOD are assumed to be expressed in the
C             same coordinates as the corresponding elements of the
C             system vector SV and reaction force vector REAC.
C             If IADDRF is negative, the reaction force contributions
C             are subtracted from REAC instead.
C
C     ROUTINES CALLED/REFERENCED :  ASMERR  (SAM-7)
C
C     PROGRAMMED BY : Knut Morten Okstad (based on ADDNV)
C     DATE/VERSION  : 05-10-11 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           INOD,IADDRF,LPU,NNDOF
      INTEGER           MADOF(*),MMCEQ(*),MEQN(*),MPAR(*)
      INTEGER           MPMCEQ(*),MPREAC(*)
      DOUBLE PRECISION  TTCC(*),VNOD(*),SV(*),REAC(*)
C
      INTEGER           I,ICEQ,IDOFM,IEQ,IP,IREAC
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------

      NNDOF = 0
      IF (INOD .LT. 1 .OR. INOD .GT. MPAR(1)) THEN
         CALL ASMERR (18,INOD,INOD,LPU,NNDOF)
         RETURN
      END IF
C
      DO 200 I = MADOF(INOD), MADOF(INOD+1)-1
         NNDOF = NNDOF + 1
         IEQ   = MEQN(I)
         IREAC = MPREAC(I)
         IF (IEQ .GT. 0) THEN
C
C     ..... Free DOF, add to system RHS vector
            SV(IEQ) = SV(IEQ) + VNOD(NNDOF)
C
         ELSE IF (IREAC .GT. 0) THEN
C
C     ..... Fixed or prescribed DOF, add to reaction forces
            IF (IADDRF .GE. 0) THEN
               REAC(IREAC) = REAC(IREAC) + VNOD(NNDOF)
            ELSE
               REAC(IREAC) = REAC(IREAC) - VNOD(NNDOF)
            END IF
C
         ELSE IF (IEQ .LT. 0) THEN
C
C     ..... Slave DOF, add to the associated master DOFs
            ICEQ = -IEQ
            DO 100 IP = MPMCEQ(ICEQ)+1, MPMCEQ(ICEQ+1)-1
               IDOFM = MMCEQ(IP)
               IF (IDOFM .GT. 0) THEN
C
C     ........... Free master DOF
                  IEQ = MEQN(IDOFM)
                  SV(IEQ) = SV(IEQ) + TTCC(IP)*VNOD(NNDOF)
C
               ELSE IF (IDOFM .LT. 0) THEN
C
C     ........... Fixed or prescribed master DOF
                  IREAC = MPREAC(-IDOFM)
                  IF (IADDRF .GE. 0) THEN
                     REAC(IREAC) = REAC(IREAC) + TTCC(IP)*VNOD(NNDOF)
                  ELSE
                     REAC(IREAC) = REAC(IREAC) - TTCC(IP)*VNOD(NNDOF)
                  END IF
C
               END IF
  100       CONTINUE
C
         END IF
  200 CONTINUE
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE EXTEV2 (SV,MADOF,MPMNPC,MMNPC,IEL,nndof,EV,NEDOF)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EXTEV2                   GROUP 7 / PUBLIC
C
C     Same as EXTEV. However, this version allows that the number of
C     nodal dofs for the element (NNDOF) is less than given by MADOF.
C
C     ROUTINES CALLED/REFERENCED :  MIN  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-29 / 1.0   (EXTEV)
C                       02-07-20 / 2.0   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IEL,NEDOF,nndof,MADOF(*),MMNPC(*),MPMNPC(*)
      DOUBLE PRECISION  EV(*),SV(*)
C
      INTEGER           I,IE,INOD,IS,J,JE,JS
C ----------------------------------------------------------------------
      NEDOF = 0
      IS    = MPMNPC(IEL)
      IE    = MPMNPC(IEL+1) - 1
      DO 100 I=IS,IE
         INOD = MMNPC(I)
         JS   = MADOF(INOD)
         JE   = MIN(JS+nndof,MADOF(INOD+1)) - 1
         DO 50 J=JS,JE
            NEDOF     = NEDOF+1
            EV(NEDOF) = SV(J)
   50    CONTINUE
  100 CONTINUE
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE ELMEQ2 (MADOF,MNPC,MPMCEQ,MEQN,NENOD,nndof,
     +                   MEEN,NEDOF,NESLV,NEPRD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ELMEQ2                   GROUP 7 / PRIVATE
C
C     Same as ELMEQ. However, this version allows that the number of
C     nodal dofs for the element (NNDOF) is less than given by MADOF.
C
C     ROUTINES CALLED/REFERENCED :  MIN  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-25 / 1.0
C                       02-07-20 / 2.0   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           NEDOF,NENOD,NEPRD,NESLV,nndof
      INTEGER           MADOF(*),MEEN(*),MEQN(*),MNPC(NENOD),MPMCEQ(*)
C
      INTEGER           I,ICEQ,IE,INOD,IS,NMST,NN
C ----------------------------------------------------------------------
      NEDOF = 0
      NESLV = 0
      NEPRD = 0
      DO 100 INOD=1,NENOD
         NN = MNPC(INOD)
         IS = MADOF(NN)
         IE = MIN(IS+nndof,MADOF(NN+1)) - 1
         DO 50 I=IS,IE
            NEDOF       = NEDOF+1
            MEEN(NEDOF) = MEQN(I)
            IF (MEQN(I).GE.0)  GO TO 50
C
            ICEQ =-MEQN(I)
            NMST = MPMCEQ(ICEQ+1) - MPMCEQ(ICEQ) - 1
C                                               ** PRESCRIBED  D O F
            IF (NMST.EQ.0)     NEPRD = NEPRD+1
C                                               ** DEPENDENT  D O F
            IF (NMST.GT.0)     NESLV = NESLV+1
   50    CONTINUE
  100 CONTINUE
C ----------------------------------------------------------------------
      RETURN
      END
