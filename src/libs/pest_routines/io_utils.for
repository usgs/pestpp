c     ******************************************************************
      SUBROUTINE PARNAM(IFAIL,J1,J2,TPAR,CLINE)
      implicit none
c     PEST routine by J. Doherty
C --  SUBROUTINE PARNAM EXTRACTS A PARAMETER NAME FROM A STRING
c     ******************************************************************

      INTEGER IFAIL
      INTEGER J1,J2,I,J
      CHARACTER*50 TPAR
      CHARACTER*(*) CLINE

      IFAIL=0
      TPAR=' '
      IF(J2-J1.LE.1) THEN
        IFAIL=1
        RETURN
      END IF
      DO 10 I=J1+1,J2-1
      IF(CLINE(I:I).EQ.' ') GO TO 10
      GO TO 30
10    CONTINUE
      IFAIL=2
      RETURN
30    J=MIN(50,J2-I)
      TPAR(1:J)=CLINE(I:I+J-1)
      RETURN
#ifdef PESTMOD
      END SUBROUTINE PARNAM
#else
      END
#endif

c     ******************************************************************
      SUBROUTINE WHICH1(IFAIL,NPAR,IPAR,APAR,TPAR)
      implicit none
c     PEST routine by J. Doherty
C --  SUBROUTINE WHICH1 LOCATES A STRING IN AN ARRAY
c     ******************************************************************

      INTEGER NPAR,IPAR,I
      INTEGER IFAIL
      CHARACTER*(*) TPAR
      CHARACTER*(*) APAR(NPAR)

      IFAIL=0
      IF((IPAR.LT.1).OR.(IPAR.GT.NPAR)) IPAR=1
      CALL LOWCAS(TPAR)
      IF(TPAR.EQ.APAR(IPAR)) RETURN
      IF(IPAR.NE.NPAR)THEN
        DO 20 I=IPAR+1,NPAR
        IF(TPAR.EQ.APAR(I))THEN
          IPAR=I
          RETURN
        END IF
20      CONTINUE
      END IF
      IF(IPAR.NE.1)THEN
        DO 40 I=IPAR-1,1,-1
        IF(TPAR.EQ.APAR(I)) THEN
          IPAR=I
          RETURN
        END IF
40      CONTINUE
      END IF
      IFAIL=1
      RETURN
#ifdef PESTMOD
      END SUBROUTINE WHICH1
#else        
      END
#endif


c     ******************************************************************
      SUBROUTINE WRTSIG(IFAIL,VAL,WORD,NW,PRECIS,TVAL,NOPNT)
      implicit none
c     PEST routine by J. Doherty
C --
C --  SUBROUTINE WRTSIG WRITES A NUMBER INTO A CONFINED SPACE WITH MAXIMUM
C --  PRECISION
C --

c       failure criteria:
c           ifail= 1 ...... number too large or small for single precision type
c           ifail= 2 ...... number too large or small for double precision type
c           ifail= 3 ...... field width too small to represent number
c           ifail=-1 ...... internal error type 1
c           ifail=-2 ...... internal error type 2
c           ifail=-3 ...... internal error type 3
c     ******************************************************************

      INTEGER PRECIS,LW,POS,INC,D,P,W,J,JJ,K,JEXP,N,JFAIL,NW,
     +        EPOS,PP,NOPNT,KEXP,IFLAG,LEXP
      INTEGER IFAIL
      DOUBLE PRECISION VAL,TVAL
      CHARACTER*29 TWORD,TTWORD,FMT*14
      CHARACTER*50 WORD

C     The following line overcomes what appears to be a bug in the LF90
C     compiler

#ifdef LAHEY
      if(abs(val).lt.1.0d-300) val=0.0d0
#endif

      LEXP=0
      IFLAG=0
      WORD=' '
      POS=1
      IF(VAL.LT.0.0D0)POS=0
#ifdef USE_D_FORMAT
      WRITE(TWORD,'(1PD23.15D3)') VAL
#else
      WRITE(TWORD,'(1PE23.15E3)') VAL
#endif
      READ(TWORD(20:23),'(I4)') JEXP
      EPOS=1
      IF(JEXP.LT.0)EPOS=0
      JFAIL=0
      IFAIL=0
      IF(PRECIS.EQ.0)THEN
        LW=MIN(15,NW)
      ELSE
        LW=MIN(23,NW)
      END IF

      N=0
      IF(NOPNT.EQ.1)N=N+1
      IF(POS.EQ.1)N=N+1
      IF(PRECIS.EQ.0)THEN
        IF(ABS(JEXP).GT.38)THEN
          IFAIL=1
c          WORD=trim(adjustl(WORD))//CHAR(0)
          RETURN
        END IF
        IF(POS.EQ.1) THEN
          IF(LW.GE.13) THEN
            WRITE(WORD,'(1PE13.7)',ERR=80) VAL
            GO TO 200
          END IF
        ELSE
          IF(LW.GE.14)THEN
            WRITE(WORD,'(1PE14.7)',ERR=80) VAL
            GO TO 200
          END IF
        END IF
        IF(LW.GE.14-N) THEN
          LW=14-N
          GO TO 80
        END IF
      ELSE
        IF(ABS(JEXP).GT.275)THEN
          IFAIL=2
c          WORD=trim(adjustl(WORD))//CHAR(0)
          RETURN
        END IF
        IF(POS.EQ.1) THEN
          IF(LW.GE.22) THEN
#ifdef USE_D_FORMAT
            WRITE(WORD,'(1PD22.15D3)',ERR=80) VAL
#else
            WRITE(WORD,'(1PE22.15E3)',ERR=80) VAL
#endif
            GO TO 200
          END IF
        ELSE
          IF(LW.GE.23) THEN
#ifdef USE_D_FORMAT
            WRITE(WORD,'(1PD23.15D3)',ERR=80) VAL
#else
            WRITE(WORD,'(1PE23.15E3)',ERR=80) VAL
#endif
            GO TO 200
          END IF
        END IF
        IF(LW.GE.23-N)THEN
          LW=23-N
          GO TO 80
        END IF
      END IF

      IF(NOPNT.EQ.1)THEN
        IF((JEXP.EQ.LW-2+POS).OR.(JEXP.EQ.LW-3+POS))THEN
          WRITE(FMT,15)LW+1
15        FORMAT('(F',I2,'.0)')
          WRITE(WORD,FMT,ERR=19) VAL
          IF(INDEX(WORD,'*').NE.0) GO TO 19
          IF(WORD(1:1).EQ.' ') GO TO 19
          WORD(LW+1:LW+1)=' '
          GO TO 200
        END IF
      END IF
19    D=MIN(LW-2+POS,LW-JEXP-3+POS)
20    IF(D.LT.0) GO TO 80
      WRITE(FMT,30) LW,D
30    FORMAT('(F',I2,'.',I2,')')
      WRITE(WORD,FMT,ERR=80) VAL
      IF(INDEX(WORD,'*').NE.0) THEN
        D=D-1
        GO TO 20
      END IF
      K=INDEX(WORD,'.')
      IF(K.EQ.0)THEN
        IFAIL=-1
c        WORD=trim(adjustl(WORD))//CHAR(0)
        RETURN
      END IF
      IF((K.EQ.1).OR.((POS.EQ.0).AND.(K.EQ.2)))THEN
        DO 70 J=1,3
        IF(K+J.GT.LW) GO TO 75
        IF(WORD(K+J:K+J).NE.'0') GO TO 200
70      CONTINUE
        GO TO 80
75      IFAIL=3
c        WORD=trim(adjustl(WORD))//CHAR(0)
        RETURN
      END IF
      GO TO 200

80    WORD=' '
      IF(NOPNT.EQ.0)THEN
        D=LW-7
        IF(POS.EQ.1) D=D+1
        IF(EPOS.EQ.1) D=D+1
        IF(ABS(JEXP).LT.100) D=D+1
        IF(ABS(JEXP).LT.10) D=D+1
        IF((JEXP.GE.100).AND.(JEXP-(D-1).LT.100))THEN
          P=1+(JEXP-99)
          D=D+1
          LEXP=99
        ELSE IF((JEXP.GE.10).AND.(JEXP-(D-1).LT.10))THEN
          P=1+(JEXP-9)
          D=D+1
          LEXP=9
        ELSE IF((JEXP.EQ.-10).OR.(JEXP.EQ.-100)) THEN
          IFLAG=1
          D=D+1
        ELSE
          P=1
        END IF
        INC=0
85      IF(D.LE.0) GO TO 300
        IF(IFLAG.EQ.0)THEN
          WRITE(FMT,100,ERR=300) P,D+7,D-1
        ELSE
          WRITE(FMT,100,ERR=300) 0,D+8,D
        END IF
        WRITE(TWORD,FMT) VAL
        IF(IFLAG.EQ.1) GO TO 87
        READ(TWORD(D+4:D+7),'(I4)',ERR=500) KEXP
        IF(((KEXP.EQ.10).AND.((JEXP.EQ.9).OR.(LEXP.EQ.9))).OR.
     +    ((KEXP.EQ.100).AND.((JEXP.EQ.99).OR.LEXP.EQ.99))) THEN
          IF(INC.EQ.0)THEN
            IF(LEXP.EQ.0)THEN
              IF(D-1.EQ.0) THEN
                D=D-1
              ELSE
                P=P+1
              END IF
            ELSE IF(LEXP.EQ.9)THEN
              IF(JEXP-(D-2).LT.10) THEN
                P=P+1
              ELSE
                D=D-1
              END IF
            ELSE IF(LEXP.EQ.99)THEN
              IF(JEXP-(D-2).LT.100)THEN
                P=P+1
              ELSE
                D=D-1
              END IF
            END IF
            INC=INC+1
            GO TO 85
          END IF
        END IF
#ifdef USE_D_FORMAT
87      J=INDEX(TWORD,'D')
#else
87      J=INDEX(TWORD,'E')
#endif
        GO TO 151
      END IF
      INC=0
      P=LW-2
      PP=JEXP-(P-1)
      IF(PP.GE.10)THEN
        P=P-1
        IF(PP.GE.100)P=P-1
      ELSE IF(PP.LT.0)THEN
        P=P-1
        IF(PP.LE.-10)THEN
          P=P-1
          IF(PP.LE.-100)P=P-1
        END IF
      END IF
      IF(POS.EQ.0)P=P-1
90    CONTINUE
      D=P-1
      W=D+8
      WRITE(FMT,100) P,W,D
      IF(D.LT.0)THEN
        IF(JFAIL.EQ.1) GO TO 300
        JFAIL=1
        P=P+1
        GO TO 90
      END IF
#ifdef USE_D_FORMAT
100   FORMAT('(',I2,'pD',I2,'.',I2,'D3)')
#else
100   FORMAT('(',I2,'pE',I2,'.',I2,'E3)')
#endif
      WRITE(TWORD,FMT) VAL
#ifdef USE_D_FORMAT
      J=INDEX(TWORD,'D')
#else
      J=INDEX(TWORD,'E')
#endif
      IF(TWORD(J-1:J-1).NE.'.')THEN
        IFAIL=-1
c        WORD=trim(adjustl(WORD))//CHAR(0)
        RETURN
      END IF
      N=1
      IF(TWORD(J+1:J+1).EQ.'-') N=N+1
      IF(TWORD(J+2:J+2).NE.'0') THEN
        N=N+2
        GO TO 120
      END IF
      IF(TWORD(J+3:J+3).NE.'0') N=N+1
120   N=N+1
      IF(J+N-2-POS.LT.LW)THEN
        IF(INC.EQ.-1) GO TO 150
        TTWORD=TWORD
        P=P+1
        INC=1
        GO TO 90
      ELSE IF(J+N-2-POS.EQ.LW) THEN
        GO TO 150
      ELSE
        IF(INC.EQ.1)THEN
          TWORD=TTWORD
          GO TO 150
        END IF
        IF(JFAIL.EQ.1) GO TO 300
        P=P-1
        INC=-1
        GO TO 90
      END IF

150   J=INDEX(TWORD,'.')
151   IF(POS.EQ.0)THEN
        K=1
      ELSE
       K=2
      END IF
      WORD(1:J-K)=TWORD(K:J-1)
      JJ=J
      J=J-K+1
      IF(PRECIS.EQ.0)THEN
        WORD(J:J)='E'
      ELSE
#ifdef USE_D_FORMAT
        WORD(J:J)='D'
#else
        WORD(J:J)='E'
#endif
      END IF
      JJ=JJ+2
      IF(NOPNT.EQ.0) JJ=JJ-1
      IF(TWORD(JJ:JJ).EQ.'-')THEN
        J=J+1
        WORD(J:J)='-'
      END IF
      IF(TWORD(JJ+1:JJ+1).NE.'0')THEN
        J=J+2
        WORD(J-1:J)=TWORD(JJ+1:JJ+2)
        GO TO 180
      END IF
      IF(TWORD(JJ+2:JJ+2).NE.'0')THEN
        J=J+1
        WORD(J:J)=TWORD(JJ+2:JJ+2)
      END IF
180   J=J+1
      WORD(J:J)=TWORD(JJ+3:JJ+3)
      IF(IFLAG.EQ.1)THEN
        IF(POS.EQ.1)THEN
          JJ=1
        ELSE
          JJ=2
        END IF
        N=LEN_TRIM(WORD)
        DO 190 J=JJ,N-1
190     WORD(J:J)=WORD(J+1:J+1)
        WORD(N:N)=' '
      END IF

200   IF(LEN_TRIM(WORD).GT.LW)THEN
        IFAIL=-2
c        WORD=trim(adjustl(WORD))//CHAR(0)
        RETURN
      END IF
      WRITE(FMT,30) LW,0
      READ(WORD,FMT,ERR=400) TVAL
c      WORD=trim(adjustl(WORD))//CHAR(0)
      RETURN
300   IFAIL=3
c      WORD=trim(adjustl(WORD))//CHAR(0)
      RETURN
400   IFAIL=-3
c      WORD=trim(adjustl(WORD))//CHAR(0)
      RETURN
500   IFAIL=-2
c      WORD=trim(adjustl(WORD))//CHAR(0)
      RETURN
#ifdef PESTMOD
      END SUBROUTINE WRTSIG
#else        
      END
#endif


c     ******************************************************************
      SUBROUTINE LOWCAS(ASTRNG)
      implicit none
c     PEST routine by J. Doherty
C --  SUBROUTINE LOWCAS CONVERTS A STRING TO LOWER CASE
c     ******************************************************************

      INTEGER I,J
      CHARACTER*(*) ASTRNG

      DO 10 I=1,LEN_TRIM(ASTRNG)
      J=ICHAR(ASTRNG(I:I))
      IF((J.GE.65).AND.(J.LE.90)) ASTRNG(I:I)=CHAR(J+32)
10    CONTINUE
      RETURN
#ifdef PESTMOD
      END SUBROUTINE LOWCAS
#else        
      END
#endif

c     ******************************************************************
      SUBROUTINE TABREM(CLINE)
      implicit none
c     PEST routine by J. Doherty
C --  SUBROUTINE TABREM REMOVES TABS FROM A STRING
c     ******************************************************************

      INTEGER I
      CHARACTER*(*) CLINE

      DO 10 I=1,LEN(CLINE)
10    IF(ICHAR(CLINE(I:I)).EQ.9) CLINE(I:I)=' '

      RETURN
#ifdef PESTMOD
      END SUBROUTINE TABREM
#else
      END
#endif

c     ******************************************************************
      SUBROUTINE CMPRSS(CLINE)
      implicit none
c     PEST routine by J. Doherty
C --  SUBROUTINE CMPRSS COMPRESSES AN INSTRUCTION LINE BY REMOVING EXCESS
C --  BLANK CHARACTERS
c     ******************************************************************

      INTEGER NBLC,J
      CHARACTER*(*) CLINE

      IF(CLINE.EQ.' ') RETURN
10    NBLC=LEN_TRIM(CLINE)
      J=INDEX(CLINE(1:NBLC),'  ')
      IF(J.NE.0) THEN
c        CALL SHIFTL(CLINE(J+1:NBLC))
        CLINE(J+1:NBLC)=ADJUSTL(CLINE(J+1:NBLC))
        GO TO 10
      END IF
      RETURN

#ifdef PESTMOD
      END SUBROUTINE CMPRSS
#else
      END
#endif

c     ******************************************************************
      SUBROUTINE OUTRD(JFAIL,NINSTR,NOUFL,ASIZE,NUML,NOBS,NBLBMX,
     +                 LCINS,LL,OBSN1,OBSN2,IIOBS,OBS,AOBS,A,MRKDEL,
     +                 CLINE,BUF)
c     PEST routine by J. Doherty
c     ******************************************************************
      implicit none
C -- SUBROUTINE OUTRD READS MODEL OUTPUT FILES USING INTRUCTIONS

ctm        USE PESTDATA, ONLY : INST

        LOGICAL LOPENED
        INTEGER MCALL,CIL,IFAIL,ASIZE,NINSTR,NCALL,JFAIL
        INTEGER INS,NBLB,NBLC,I,J,N1,N2,NOL,MRKTYP,J1,J2,inst,
     +  NUML,JOBS,IOBS,NOBS,N3,NUM1,NUM2,IFILE,IL,NOUFL,NBLBMX,
     +  INSNUM,ALMARK,BEGINS,INSFLE,DUMFLG
        INTEGER LCINS(NINSTR)
        INTEGER LL(NUML),OBSN1(NOBS),OBSN2(NOBS),IIOBS(NOBS)
        DOUBLE PRECISION OBS(NOBS),RTEMP
        CHARACTER*15 FMT,OBSNAM*50,MKRDEL*1,AA*1
        CHARACTER MRKDEL(NOUFL)
        CHARACTER A(ASIZE)
        CHARACTER*(*) AOBS(NOBS)
        CHARACTER*(*) CLINE,BUF
        CHARACTER*200 FLENME

ctm        COMMON /FLENME/FLENME
ctm        COMMON /MODCAL/NCALL,MCALL

        JFAIL=0
        IFILE=0
        IL=0
        JOBS=0
        MKRDEL=MRKDEL(1)
        CIL=0
        IOBS=1
        BEGINS=0
        BUF=' '
        
ctm     ------------------        
        inst=35
        mcall=1
ctm     ------------------        

        INS=1
10      IF(INS.LT.NINSTR)THEN
          NBLB=LCINS(INS+1)-LCINS(INS)
        ELSE
          NBLB=ASIZE-LCINS(INS)+1
        END IF
c        BUF(1:NBLBMX)=' '
        BUF(1:200)=' '
        DO 20 I=1,NBLB
20      BUF(I:I)=A(LCINS(INS)+I-1)
25      N2=0
        INSNUM=0

50      CALL GETINT(IFAIL,BUF,N1,N2,NBLB,MKRDEL)
        IF(IFAIL.NE.0) THEN
ctm           CALL STPERR(68,5,BUF,0,' ',CLINE)
           GO TO 9891
        END IF
51      IF(N1.EQ.0) GO TO 1000
        INSNUM=INSNUM+1
        IF(INSNUM.EQ.1)THEN
          IF(BUF(N1:N1).NE.'&') THEN
            MRKTYP=0
            ALMARK=1
            BEGINS=0
          ELSE
            IF(INS.EQ.INSFLE+1) THEN
ctm              CALL STPERR(73,5,BUF,0,' ',CLINE)
              GO TO 9891
            END IF
            IF(BEGINS.EQ.1)THEN
              INS=INS-1
              GO TO 10
            END IF
          END IF
        END IF
        IF(ICHAR(BUF(N1:N1)).EQ.2)THEN
          IF(IFILE.NE.0) CLOSE(UNIT=INST)
          DO 60 I=N1+1,NBLB
          IF(BUF(I:I).NE.' ') GO TO 70
60        CONTINUE
70        FLENME=BUF(I:NBLB)
          open(inst,file=flenme)
ct,          CALL FFOPEN(JFAIL,INST,'r',' ',65,CLINE)
          IF(JFAIL.NE.0) GO TO 9891
          IFILE=IFILE+1
          CIL=0
          MKRDEL=MRKDEL(IFILE)
          INSFLE=INS
          GO TO 1000
        ELSE IF((BUF(N1:N1).EQ.'l').OR.(BUF(N1:N1).EQ.'L'))THEN
          ALMARK=0
          IL=IL+1
          IF(MCALL.EQ.1)THEN
            WRITE(FMT,150) N2-N1
150         FORMAT('(I',I4,')')
            READ(BUF(N1+1:N2),FMT,ERR=9050) NOL
            LL(IL)=NOL
          ELSE
            NOL=LL(IL)
          END IF
          IF(NOL.GT.1) THEN
            DO 160 I=1,NOL-1
            READ(INST,*,END=9100)
160         CIL=CIL+1
          END IF
          READ(INST,22,END=9100) CLINE
22        FORMAT(A)
          IF(INDEX(CLINE,CHAR(9)).NE.0) CALL TABREP(CLINE)
          CIL=CIL+1
          NBLC=LEN_TRIM(CLINE)
          MRKTYP=1
          J1=0
        ELSE IF(BUF(N1:N1).EQ.MKRDEL)THEN
          IF(MRKTYP.EQ.0)THEN
200         READ(INST,22,END=9100) CLINE
            IF(INDEX(CLINE,CHAR(9)).NE.0) CALL TABREP(CLINE)
            CIL=CIL+1
            J1=INDEX(CLINE,BUF(N1+1:N2-1))
            IF(J1.EQ.0) GO TO 200
            NBLC=LEN_TRIM(CLINE)
            J1=J1+N2-N1-2
            MRKTYP=1
          ELSE
            IF(J1.GE.NBLC) THEN
              IF(ALMARK.EQ.1) THEN
                BEGINS=1
                GO TO 25
              END IF
              GO TO 9200
            END IF
            J2=INDEX(CLINE(J1+1:NBLC),BUF(N1+1:N2-1))
            IF(J2.EQ.0) THEN
              IF(ALMARK.EQ.1) THEN
                BEGINS=1
                GO TO 25
              END IF
              GO TO 9200
            END IF
            J1=J1+J2
            J1=J1+N2-N1-2
          END IF
        ELSE IF(BUF(N1:N1).EQ.'&')THEN
          IF(INSNUM.NE.1) THEN
ctm            CALL STPERR(72,5,BUF,0,' ',CLINE)
            GO TO 9891
          END IF
        ELSE IF((BUF(N1:N1).EQ.'w').OR.(BUF(N1:N1).EQ.'W'))THEN
          ALMARK=0
          IF(J1.GE.NBLC) GO TO 9400
          J2=INDEX(CLINE(J1+1:NBLC),' ')
          IF(J2.EQ.0) GO TO 9400
          J1=J1+J2
          DO 210 I=J1,NBLC
          IF(CLINE(I:I).NE.' ') GO TO 220
210       CONTINUE
          I=NBLC+1
220       J1=I-1
        ELSE IF((BUF(N1:N1).EQ.'t').OR.(BUF(N1:N1).EQ.'T'))THEN
          ALMARK=0
          WRITE(FMT,150) N2-N1
          READ(BUF(N1+1:N2),FMT,ERR=9000) J2
          IF(J2.LT.J1) THEN
ctm            CALL STPERR(81,4,BUF,CIL,' ',CLINE)
            GO TO 9891
          END IF
          J1=J2
          IF(J1.GT.NBLC) THEN
ctm            CALL STPERR(70,4,BUF,CIL,' ',CLINE)
            GO TO 9891
          END IF
        ELSE IF((BUF(N1:N1).EQ.'[').OR.(BUF(N1:N1).EQ.'('))THEN
          ALMARK=0
          AA=BUF(N1:N1)
          JOBS=JOBS+1
          IF(MCALL.EQ.1)THEN
            IF(AA.EQ.'[')THEN
              N3=INDEX(BUF(N1:N2),']')
            ELSE
              N3=INDEX(BUF(N1:N2),')')
            END IF
            N3=N3+N1-1
            OBSNAM=BUF(N1+1:N3-1)
            CALL WHICH1(IFAIL,NOBS,IOBS,AOBS,OBSNAM)
            IF(IFAIL.NE.0) GO TO 9700
            CALL GETNUM(IFAIL,BUF,N3,N2,NUM1,NUM2,FMT)
            IF(IFAIL.NE.0) THEN
ctm              CALL STPERR(64,5,BUF,0,' ',CLINE)
              GO TO 9891
            END IF
            IF(NUM1.LE.0)THEN
ctm              CALL STPERR(64,5,BUF,0,' ',CLINE)
              GO TO 9891
            END IF
            OBSN1(JOBS)=NUM1
            OBSN2(JOBS)=NUM2
            IIOBS(JOBS)=IOBS
          ELSE
            NUM1=OBSN1(JOBS)
            NUM2=OBSN2(JOBS)
            IOBS=IIOBS(JOBS)
          END IF
          IF(AA.EQ.'(') THEN
            CALL GETTOT(IFAIL,CLINE,NUM1,NUM2,NBLC)
            IF(IFAIL.NE.0) THEN
ctm              CALL STPERR(88,3,AOBS(IOBS)(:LEN_TRIM(AOBS(IOBS))),
ctm     +        CIL,' ',CLINE)
              GO TO 9891
            END IF
          ELSE
           IF(NUM1.GT.NBLC) THEN
ctm             CALL STPERR(88,3,
ctm     +       AOBS(IOBS)(:LEN_TRIM(AOBS(IOBS))),CIL,' ',CLINE)
             GO TO 9891
           END IF
           IF(NUM2.GT.NBLC) NUM2=NBLC
           IF(CLINE(NUM1:NUM2).EQ.' ') THEN
ctm             CALL STPERR(88,3,
ctm     +       AOBS(IOBS)(:LEN_TRIM(AOBS(IOBS))),CIL,' ',CLINE)
             GO TO 9891
           END IF
          END IF
          WRITE(FMT,250) NUM2-NUM1+1
250       FORMAT('(F',I4,'.0)')
          READ(CLINE(NUM1:NUM2),FMT,ERR=260) OBS(IOBS)
          J1=NUM2
          GO TO 50
ctm260       CALL STPERR(82,3,AOBS(IOBS)(:LEN_TRIM(AOBS(IOBS))),
ctm     +      CIL,' ',CLINE)
260            GO TO 9891
        ELSE IF(BUF(N1:N1).EQ.'!') THEN
          ALMARK=0
          CALL LOWCAS(BUF(N1+1:N2-1))
          IF((N2-N1.NE.4).OR.(BUF(N1+1:N2-1).NE.'dum'))THEN
            JOBS=JOBS+1
            IF(MCALL.EQ.1) THEN
              OBSNAM=BUF(N1+1:N2-1)
              CALL WHICH1(IFAIL,NOBS,IOBS,AOBS,OBSNAM)
              IF(IFAIL.NE.0) GO TO 9700
              IIOBS(JOBS)=IOBS
            ELSE
              IOBS=IIOBS(JOBS)
            END IF
            DUMFLG=0
          ELSE
            DUMFLG=1
          END IF
          CALL GETNXT(IFAIL,CLINE,J1,NUM1,NUM2,NBLC)
          IF(IFAIL.NE.0) THEN
            IF(DUMFLG.EQ.0) THEN
ctm              CALL STPERR(88,3,AOBS(IOBS)(:LEN_TRIM(AOBS(IOBS))),
ctm     +          CIL,' ',CLINE)
                GO TO 9891
            ELSE
ctm              CALL STPERR(88,3,'dum',CIL,' ',CLINE)
              GO TO 9891
            END IF
          END IF
          WRITE(FMT,250) NUM2-NUM1+1
          READ(CLINE(NUM1:NUM2),FMT,ERR=270) RTEMP
          IF(DUMFLG.EQ.0) OBS(IOBS)=RTEMP
          J1=NUM2
          GO TO 50
270       CALL GETINT(IFAIL,BUF,N1,N2,NBLB,MKRDEL)
          IF(IFAIL.NE.0) THEN
ctm            CALL STPERR(68,5,BUF,0,' ',CLINE)
            GO TO 9891
          END IF
          IF(N1.EQ.0)THEN
            IF(DUMFLG.EQ.1) GO TO 9900
            GO TO 9800
          END IF
          IF(BUF(N1:N1).NE.MKRDEL) THEN
            IF(DUMFLG.EQ.1) GO TO 9900
            GO TO 9800
          END IF
          J2=INDEX(CLINE(J1+1:NBLC),BUF(N1+1:N2-1))
          IF(J2.EQ.0) THEN
            IF(DUMFLG.EQ.1) GO TO 9900
            GO TO 9800
          END IF
          NUM2=J1+J2-1
          IF(NUM2.LT.NUM1)THEN
            IF(DUMFLG.EQ.1) GO TO 9900
            GO TO 9800
          END IF
          WRITE(FMT,250) NUM2-NUM1+1
          IF(DUMFLG.EQ.1)THEN
            READ(CLINE(NUM1:NUM2),FMT,ERR=9900) RTEMP
          ELSE
            READ(CLINE(NUM1:NUM2),FMT,ERR=9800) OBS(IOBS)
          END IF
          J1=NUM2
          GO TO 51
        ELSE
ctm          CALL STPERR(64,5,BUF,0,' ',CLINE)
          GO TO 9891
        END IF
        GO TO 50
1000    INS=INS+1
        IF(INS.LE.NINSTR) GO TO 10

        IF(MCALL.EQ.1)THEN
          DO 1100 I=1,NOBS
          DO 1050 J=1,JOBS
          IF(IIOBS(J).EQ.I) GO TO 1100
1050      CONTINUE
ctm          CALL STPERR(45,1,AOBS(I)(:LEN_TRIM(AOBS(I))),0,' ',CLINE)
          GO TO 9891
1100      CONTINUE
        END IF

        CLOSE(UNIT=INST)
        RETURN

ctm9000    CALL STPERR(74,5,BUF,0,' ',CLINE)
9000    GO TO 9891
ctm9050    CALL STPERR(75,5,BUF,0,' ',CLINE)
9050    GO TO 9891
ctm9100    CALL STPERR(66,6,BUF,0,' ',CLINE)
9100    GO TO 9891
ctm9200    CALL STPERR(67,4,BUF,CIL,' ',CLINE)
9200    GO TO 9891
ctm9400    CALL STPERR(69,4,BUF,CIL,' ',CLINE)
9400    GO TO 9891
ctm9700    CALL STPERR(84,6,BUF,0,OBSNAM(:LEN_TRIM(OBSNAM)),CLINE)
9700    GO TO 9891
ctm9800    CALL STPERR(82,3,AOBS(IOBS)(:LEN_TRIM(AOBS(IOBS))),
ctm     +    CIL,' ',CLINE)
9800    GO TO 9891
ctm9900    CALL STPERR(82,3,'dum',CIL,' ',CLINE)
9900    GO TO 9891

9891    JFAIL=1
        INQUIRE(UNIT=INST,OPENED=LOPENED)
        IF(LOPENED)CLOSE(UNIT=INST)
        RETURN

#ifdef PESTMOD
        END SUBROUTINE OUTRD
#else
        END
#endif



        SUBROUTINE GETINT(IFAIL,BUF,N1,N2,NBLB,MRKDEL)
      implicit none

C -- SUBROUTINE GETINT GETS THE NEXT STORED INSTRUCTION FOR PROCESSING

        INTEGER N1,N2,NBLB,I,II
        INTEGER IFAIL
        CHARACTER MRKDEL
        CHARACTER*(*) BUF

        IFAIL=0
        IF(N2.GE.NBLB) THEN
          N1=0
          RETURN
        END IF
        DO 10 I=N2+1,NBLB
        IF((BUF(I:I).NE.' ').AND.(ICHAR(BUF(I:I)).NE.9)) GO TO 50
10      CONTINUE
        N1=0
        RETURN
50      N1=I
        IF(BUF(N1:N1).NE.MRKDEL)THEN
          I=INDEX(BUF(N1:NBLB),' ')
          II=INDEX(BUF(N1:NBLB),CHAR(9))
          IF((I.EQ.0).AND.(II.EQ.0))THEN
            I=0
          ELSE IF(I.EQ.0)THEN
            I=II
          ELSE IF(II.EQ.0) THEN
            I=I
          ELSE
            I=MIN(I,II)
          END IF
          IF(I.NE.0) THEN
            N2=N1+I-2
          ELSE
            N2=NBLB
          END IF
        ELSE
          IF(N1.EQ.NBLB)THEN
            IFAIL=1
            RETURN
          END IF
          I=INDEX(BUF(N1+1:NBLB),MRKDEL)
          IF(I.EQ.0) THEN
            IFAIL=1
            RETURN
          END IF
          N2=N1+I
        END IF

        RETURN
#ifdef PESTMOD
        END SUBROUTINE GETINT
#else        
        END
#endif        


        SUBROUTINE TABREP(CLINE)
      implicit none

C -- SUBROUTINE TABREP REPLACES A TAB BY BLANK SPACE(S)

        INTEGER LLEN,I,J,K,NBLC
        CHARACTER*(*) CLINE

        LLEN=LEN(CLINE)
        DO 10 I=LLEN,1,-1
        IF(CLINE(I:I).NE.' ') GO TO 20
10      CONTINUE
        RETURN
20      NBLC=I

        I=0
30      I=I+1
        IF(I.GT.NBLC)RETURN
        IF(ICHAR(CLINE(I:I)).NE.9) GO TO 30
        J=((I-1)/8+1)*8-I
        IF(J.EQ.0) THEN
          CLINE(I:I)=' '
        ELSE
          CLINE(I:I)=' '
          NBLC=NBLC+J
          IF(NBLC.GT.LLEN) NBLC=LLEN
          DO 50 K=NBLC,((I-1)/8+1)*8,-1
          CLINE(K:K)=CLINE(K-J:K-J)
50        CONTINUE
          DO 60 K=I+1,MIN(NBLC,I+J)
          CLINE(K:K)=' '
60        CONTINUE
          I=I+J
        END IF
        GO TO 30

#ifdef PESTMOD
        END SUBROUTINE TABREP
#else
        END
#endif


        SUBROUTINE GETNUM(IFAIL,BUF,N3,N2,NUM1,NUM2,FMT)
      implicit none

C -- SUBROUTINE GETNUM RETRIEVES CHARACTER POSITIONS FROM FIXED AND
C -- SEMI-FIXED OBSERVATION INSTRUCTIONS

        INTEGER N3,NUM1,NUM2,I,N2
        INTEGER IFAIL
        CHARACTER*(*) BUF
        CHARACTER*(*) FMT

        IFAIL=0
        I=INDEX(BUF(N3+1:N2),':')
        IF(I.EQ.0) GO TO 100
        WRITE(FMT,20) I-1
20      FORMAT('(I',I3,')')
        READ(BUF(N3+1:N3+I-1),FMT,ERR=100) NUM1
        N3=N3+I
        I=N2-N3
        IF(I.LT.1) GO TO 100
        WRITE(FMT,20) I
        READ(BUF(N3+1:N2),FMT,ERR=100) NUM2
        RETURN
100     IFAIL=1
        RETURN
#ifdef PESTMOD
        END SUBROUTINE GETNUM
#else        
        END
#endif


        SUBROUTINE GETTOT(IFAIL,CLINE,J1,J2,NBLC)
      implicit none

C -- SUBROUTINE GETTOT DETERMINES THE EXACT POSITION OCCUPIED BY A NUMBER

        INTEGER IFAIL
        INTEGER J1,J2,NBLC,I
        CHARACTER*(*) CLINE

        IFAIL=0
        IF(J1.GT.NBLC)THEN
          IFAIL=1
          RETURN
        END IF
        IF(J2.GT.NBLC)J2=NBLC
        IF(CLINE(J2:J2).EQ.' ') THEN
          DO 10 I=J2,J1,-1
          IF(CLINE(I:I).NE.' ')THEN
            J2=I
            GO TO 100
          END IF
10        CONTINUE
          IFAIL=1
          RETURN
        ELSE
          IF(J2.EQ.NBLC) GO TO 100
          DO 20 I=J2,NBLC
          IF(CLINE(I:I).EQ.' ') THEN
            J2=I-1
            GO TO 100
          END IF
20        CONTINUE
          J2=NBLC
        END IF
100     IF(J1.EQ.1) GO TO 200
        DO 120 I=J1,1,-1
        IF(CLINE(I:I).EQ.' ') THEN
          J1=I+1
          GO TO 200
        END IF
120     CONTINUE
        J1=1
200     RETURN

#ifdef PESTMOD
        END SUBROUTINE GETTOT
#else
        END
#endif

        SUBROUTINE GETNXT(IFAIL,CLINE,J1,NUM1,NUM2,NBLC)
      implicit none

C -- SUBROUTINE GETNXT GETS THE NEXT SPACE-DELIMITED WORD

        INTEGER IFAIL
        INTEGER J1,NUM1,NUM2,NBLC,I
        CHARACTER*(*) CLINE

        IFAIL=0
        DO 20 I=J1+1,NBLC
        IF(CLINE(I:I).NE.' ') GO TO 50
20      CONTINUE
        IFAIL=1
        RETURN
50      NUM1=I
        I=INDEX(CLINE(NUM1:NBLC),' ')
        IF(I.EQ.0) THEN
          NUM2=NBLC
        ELSE
          NUM2=NUM1+I-2
        END IF

        RETURN
#ifdef PESTMOD
        END SUBROUTINE GETNXT
#else
        END
#endif
