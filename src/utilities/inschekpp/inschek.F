

        PROGRAM INSCHK

C -- Number of observations
        INTEGER ZNOBS
        PARAMETER (ZNOBS = 500000)

        INTEGER I,J,NBLC,MAXWID,IERR,MCASE,LENIN,LENOUT,CHKTYP,II,
     +  J1,EXTLEN,LENMOD,MMLINE,NOBS,MAXOBS,IS,IOBS,N1,
     +  N2,ITEMP,N3,NUM1,NUM2,BEGINS,NBLB,MMMLIN,INSNUM,MRKTYP,
     +  ALLMRK,J2,DUMFLG
        INTEGER JLINE,IFAIL,CIL,INS,NOL
        DOUBLE PRECISION RTEMP
        DOUBLE PRECISION OBSVAL(:)
        CHARACTER*10 AMXWID,ATEMP*30,TPAR*20,AMXOBS,PARDEL*1,FILADD*1
        CHARACTER*12 AVERSION
        CHARACTER*130 COMLIN,INSFLE*80,OUTFLE*80,MODFLE*80,
     +  FMT*20,AA*1
        CHARACTER*200 AALINE,AAALIN
        CHARACTER*20 AOBS(:)
        CHARACTER*2120 ERRMSG
C -- File size
        CHARACTER*2010 CLINE,BUF

        ALLOCATABLE:: AOBS,OBSVAL

        INTEGER NBLNK

        COMMON /IERR/IERR,MCASE
        COMMON /MMLINE/MMLINE
        COMMON /ERRMSG/ERRMSG
        COMMON /AALINE/AALINE

        J1=0
        FILADD=' '
        NOBS=0
        IERR=0
C -- File size
        MAXWID=2000
        CALL WRTINT(AMXWID,MAXWID)
        MAXOBS=ZNOBS
        CALL WRTINT(AMXOBS,MAXOBS)

C        INCLUDE 'version.inc'
C        WRITE(6,1,ERR=9000) TRIM(AVERSION)
1       FORMAT(' INSCHEK Version ',A,'. ',
     +  'Watermark Numerical Computing.',/)

        ALLOCATE(AOBS(MAXOBS),STAT=IFAIL)
        IF(IFAIL.NE.0) THEN
          WRITE(ERRMSG,5)
5         FORMAT(' Cannot allocate sufficient memory to run ',
     +    'INSCHEK.')
          CALL WRTERR
          GO TO 9800
        END IF

C -- THE COMMAND LINE IS PARSED

#ifdef NO_CMDLINE
3       WRITE(*,'(A)') ' Enter command line arguments > '
        READ(*,'(A)',ERR=3) COMLIN
#else
        CALL PGETCL(COMLIN)
#endif
        CALL SPACESUB(COMLIN)
        NBLC=NBLNK(COMLIN)
        IF(NBLC.EQ.0) THEN
         GO TO 9900
        END IF
        DO 20 I=NBLC,NBLC-2,-1
        IF(COMLIN(I:I).EQ.'/')THEN
          FILADD=COMLIN(I+1:I+1)
          COMLIN(I:NBLC)=' '
          NBLC=NBLNK(COMLIN)
          GO TO 25
        END IF
20      CONTINUE
25      CONTINUE
        DO 50 I=1,NBLC
        IF(COMLIN(I:I).NE.' ') GO TO 60
50      CONTINUE
        GO TO 9900
60      DO 70 J=I,NBLC
        IF(COMLIN(J:J).EQ.' ') GO TO 90
70      CONTINUE
        CHKTYP=1
        INSFLE=COMLIN(I:NBLC)
#ifdef CAPFILE
        CALL UPCAS(INSFLE)
#endif
        LENIN=NBLC-I+1
        GO TO 200
90      J=J-1
        INSFLE=COMLIN(I:J)
#ifdef CAPFILE
        CALL UPCAS(INSFLE)
#endif
        LENIN=J-I+1
        DO 100 I=J+1,NBLC
        IF(COMLIN(I:I).NE.' ') GO TO 110
100     CONTINUE
110     DO 120 J=I,NBLC
        IF(COMLIN(J:J).EQ.' ') GO TO 9900
120     CONTINUE
        CHKTYP=2
        MODFLE=COMLIN(I:NBLC)
#ifdef CAPFILE
        CALL UPCAS(MODFLE)
#endif
        LENMOD=NBLC-I+1

200     CONTINUE

        CALL REMCHAR(INSFLE,CHAR(211))
c        CALL SHIFTL(INSFLE)
        CALL REMCHAR(MODFLE,CHAR(211))
c        CALL SHIFTL(MODFLE)

        CALL EXTCHK(EXTLEN,LENIN,INSFLE)
        IF(EXTLEN.LT.LENIN)THEN
#ifdef CAPFILE
          IF(INSFLE(EXTLEN:LENIN).EQ.'.OBF')THEN
            WRITE(ERRMSG,202)
202         FORMAT(' Instruction file must not have an extension of ',
     +      '".OBF".')
#else
          IF(INSFLE(EXTLEN:LENIN).EQ.'.obf')THEN
            WRITE(ERRMSG,202)
202         FORMAT(' Instruction file must not have an extension of ',
     +      '".obf".')
#endif
            CALL WRTERR
            GO TO 9800
          END IF
        END IF
c        IF(LENIN.GT.51) THEN
c          WRITE(ERRMSG,220)
c220       FORMAT(' Name of instruction file must not exceed 51 ',
c     +    'characters.')
c          CALL WRTERR
c          GO TO 9800
c        END IF

C -- THE INSTRUCTION FILE IS OPENED

#ifdef UNICOS
        OPEN(UNIT=40,FILE=INSFLE,STATUS='OLD',RECL=1000,ERR=230)
#else
        OPEN(UNIT=40,FILE=INSFLE,STATUS='OLD',ERR=230)
#endif
        GO TO 250
230     WRITE(ERRMSG,240) INSFLE(:NBLNK(INSFLE))
240     FORMAT(' Cannot open instruction file ',A,'.')
        CALL WRTERR
        GO TO 9800
250     OUTFLE=INSFLE
        DO 310 I=EXTLEN,LEN(INSFLE)
        OUTFLE(I:I)=' '
310     CONTINUE
        OUTFLE(EXTLEN:EXTLEN+3)='.obf'
        IF(FILADD.NE.' ') OUTFLE='pest'//FILADD//'.tmp'
#ifdef CAPFILE
        CALL UPCAS(OUTFLE)
#endif
        LENOUT=NBLNK(OUTFLE)
c        IF(LENOUT.GT.51) THEN
c          WRITE(ERRMSG,330) OUTFLE(:NBLNK(OUTFLE))
c330       FORMAT(' Cannot write file ',A,': name length limited ',
c     +    'to 51 characters.')
c          CALL WRTERR
c          GO TO 9800
c        END IF

C -- IF REQUESTED, THE MODEL OUTPUT FILE IS OPENED

        IF(CHKTYP.EQ.2)THEN
          IF(MODFLE(:NBLNK(MODFLE)).EQ.INSFLE(:NBLNK(INSFLE)))THEN
            WRITE(ERRMSG,212)
212         FORMAT(' Model output file cannot have the same name as ',
     +      'the instruction file.')
            CALL WRTERR
            GO TO 9800
          END IF
c          IF(LENMOD.GT.51)THEN
c            WRITE(ERRMSG,206)
c206         FORMAT(' Name of model output file may not exceed ',
c     +      '51 characters.')
c            CALL WRTERR
c            GO TO 9800
c          END IF
#ifdef UNICOS
          OPEN(UNIT=35,FILE=MODFLE,STATUS='OLD',RECL=3000,
     +    ERR=207)
#else
          OPEN(UNIT=35,FILE=MODFLE,STATUS='old',ERR=207)
#endif
          GO TO 209
207       WRITE(ERRMSG,208) MODFLE(:NBLNK(MODFLE))
208       FORMAT(' Cannot open model output file ',A,'.')
          CALL WRTERR
          GO TO 9800
        END IF
209     CONTINUE
        WRITE(6,210,ERR=9000) INSFLE(:NBLNK(INSFLE))
210     FORMAT(' Errors in file ',A,' ----->')

C -- THE INSTRUCTION FILE IS READ AND CHECKED

2640    JLINE=1
        CALL WRTINT(ATEMP,JLINE)
        AALINE=' Line '//ATEMP(1:INDEX(ATEMP,' ')-1)//' of instruction '
     +  //'file '//INSFLE(1:LENIN)//': '
        MMLINE=NBLNK(AALINE)+1
        READ(40,'(A)',ERR=2900,END=2950) CLINE
        CALL TABREM(CLINE)
        CALL LOWCAS(CLINE)
        IF((CLINE(1:3).EQ.'ptf').or.(CLINE(1:3).eq.'jtf'))THEN
          WRITE(ERRMSG,2645) INSFLE(1:LENIN)
2645      FORMAT(' File ',A,' is a template file: an instruction ',
     +    'file is expected.')
          CALL WRTERR
          GO TO 9800
        END IF
        IF(((CLINE(1:4).NE.'pif ').AND.
     +      (CLINE(1:4).NE.'jif ')).OR.(NBLNK(CLINE).LT.5))THEN
          WRITE(ERRMSG,2650) AALINE(1:MMLINE)
2650      FORMAT(A,'first line must be "pif" or "jif" + ',
     +    'marker delimiter.')
          CALL WRTERR
          GO TO 9800
        END IF
        PARDEL=CLINE(5:5)
        II=ICHAR(PARDEL)
        IF(PARDEL.EQ.' ')THEN
          WRITE(ERRMSG,2660) AALINE(1:MMLINE)
2660      FORMAT(A,'marker delimiter not found in fifth character ',
     +    'position.')
          CALL WRTERR
          GO TO 9800
        ELSE IF((PARDEL.EQ.'(').OR.(PARDEL.EQ.')').OR.(PARDEL.EQ.'[')
     +    .OR.(PARDEL.EQ.']').OR.(PARDEL.EQ.'!').OR.
     +    ((II.GE.48).AND.(II.LE.57)).OR.
     +    ((II.GE.97).AND.(II.LE.122)).OR.
     +    (PARDEL.EQ.'&').OR.(PARDEL.EQ.':'))THEN
          WRITE(ERRMSG,2670) AALINE(1:MMLINE)
2670      FORMAT(A,'illegal marker delimiter.')
          CALL WRTERR
          PARDEL=' '
          GO TO 9800
        END IF
        IS=0
2680    JLINE=JLINE+1
        CALL WRTINT(ATEMP,JLINE)
        AALINE=' Line '//ATEMP(1:INDEX(ATEMP,' ')-1)//
     +  ' of instruction file '//INSFLE(1:LENIN)//': '
        MMLINE=NBLNK(AALINE)+1
        READ(40,'(A)',ERR=2900,END=2950) CLINE
        IF(INDEX(CLINE,CHAR(9)).NE.0)THEN
          IF(INDEX(CLINE,PARDEL).EQ.0) CALL TABREM(CLINE)
        ENDIF
        IF(CLINE.EQ.' ') GO TO 2680
        NBLC=NBLNK(CLINE)
        IF(NBLC.GT.MAXWID)THEN
          WRITE(ERRMSG,2705) AALINE(1:MMLINE),AMXWID(:NBLNK(AMXWID))
2705      FORMAT(A,'line exceeds ',A,' characters.')
          CALL WRTERR
          GO TO 2680
        END IF
        CALL SHIFTL(CLINE)
        IF((IS.EQ.0).AND.(CLINE(1:1).EQ.'&'))THEN
          WRITE(ERRMSG,2706) AALINE(1:MMLINE)
2706      FORMAT(A,'first instruction line in file cannot begin ',
     +    'with continuation character.')
          CALL WRTERR
        END IF
        IS=1
        N2=0
        IF((CLINE(1:1).NE.'l').AND.(CLINE(1:1).NE.'L').AND.
     +  (CLINE(1:1).NE.'&').AND.(CLINE(1:1).NE.PARDEL))THEN
          WRITE(ERRMSG,2685) AALINE(1:MMLINE)
2685      FORMAT(A,'instruction line must begin with "l", a ',
     +    'marker delimiter or "&" (for continuation).')
          CALL WRTERR
        END IF
        IF(CLINE(1:1).NE.'&') J1=0
2690    NBLC=NBLNK(CLINE)
2700    CALL GETINT(IFAIL,CLINE,N1,N2,NBLC,PARDEL)
        IF(IFAIL.NE.0)THEN
          WRITE(ERRMSG,2710) AALINE(1:MMLINE)
2710      FORMAT(A,'unbalanced marker delimiters.')
          CALL WRTERR
          GO TO 2680
        END IF
        IF(N1.EQ.0) GO TO 2680
        IF((CLINE(N1:N1).EQ.'l').OR.(CLINE(N1:N1).EQ.'L'))THEN
          IF(N1.EQ.N2)THEN
            CALL ILINT1(CLINE(N1:N2))
            GO TO 2700
          END IF
          WRITE(FMT,2730) N2-N1
2730      FORMAT('(I',I4,')')
          READ(CLINE(N1+1:N2),FMT,ERR=2740) NOL
          GO TO 2750
2740      CALL ILINT1(CLINE(N1:N2))
          GO TO 2700
2750      IF(N1.NE.1)THEN
            WRITE(ERRMSG,2720) AALINE(1:MMLINE)
2720        FORMAT(A,'line advance instruction can only occur at the ',
     +      'beginning of an instruction line.')
            CALL WRTERR
          END IF
          IF(NOL.LE.0) THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2755)
2755        FORMAT('  (Integer following "l" must be positive.)')
            CALL WRTERR
          END IF
        ELSE IF(CLINE(N1:N1).EQ.PARDEL)THEN
          IF(INDEX(CLINE(N1:N2),CHAR(9)).NE.0)THEN
            WRITE(ERRMSG,2686) AALINE(1:MMLINE)
2686        FORMAT(A,'tab character embedded in marker; space ',
     +      'characters must be used.')
            CALL WRTERR
          END IF
          IF(N2-N1.EQ.1)THEN
            WRITE(ERRMSG,2760)AALINE(1:MMLINE)
2760        FORMAT(A,'marker has zero length.')
            CALL WRTERR
          ELSE
            J1=J1+N2-N1-1
            IF(J1.GT.MAXWID)THEN
              WRITE(ERRMSG,2765) AALINE(1:MMLINE),
     +        AMXWID(:NBLNK(AMXWID))
2765          FORMAT(A,'model output file width cannot exceed ',A,
     +        ' characters.')
              CALL WRTERR
              GO TO 2680
            END IF
          ENDIF
        ELSE IF(CLINE(N1:N1).EQ.'&')THEN
          IF(N1.NE.N2)THEN
            CALL ILINT1(CLINE(N1:N2))
          ELSE
            IF(N1.GT.1)THEN
              WRITE(ERRMSG,2770) AALINE(1:MMLINE)
2770          FORMAT(A,'"&" is a continuation character and can only ',
     +        'occur at the beginning of an instruction line.')
              CALL WRTERR
            END IF
          END IF
        ELSE IF((CLINE(N1:N1).EQ.'w').OR.(CLINE(N1:N1).EQ.'W'))THEN
          IF(N1.NE.N2)THEN
            CALL ILINT1(CLINE(N1:N2))
          ELSE
            J1=J1+1
            IF(J1.GT.MAXWID) THEN
              WRITE(ERRMSG,2765) AALINE(1:MMLINE),
     +        AMXWID(:NBLNK(AMXWID))
              CALL WRTERR
              GO TO 2680
            END IF
          ENDIF
        ELSE IF((CLINE(N1:N1).EQ.'t').OR.(CLINE(N1:N1).EQ.'T'))THEN
          IF(N1.EQ.N2)THEN
            CALL ILINT1(CLINE(N1:N2))
            GO TO 2700
          END IF
          WRITE(FMT,2730)N2-N1
          READ(CLINE(N1+1:N2),FMT,ERR=2780) ITEMP
          GO TO 2790
2780      CALL ILINT1(CLINE(N1:N2))
          GO TO 2700
2790      IF(ITEMP.LE.0)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2848)
2848        FORMAT('  (Number following "t" must be greater than ',
     +      'zero.)')
            CALL WRTERR
          ELSE IF(ITEMP.LT.J1)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2800)
2800        FORMAT('  (Tab position results in backward move - model ',
     +      'output file must be read from left to right.)')
            CALL WRTERR
          ELSE IF(ITEMP.GT.MAXWID)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2847) AMXWID(:NBLNK(AMXWID))
2847        FORMAT('  (Model output file width cannot exceed ',A,
     +      ' characters.)')
            CALL WRTERR
            GO TO 2680
          ELSE
            J1=ITEMP
          END IF
        ELSE IF((CLINE(N1:N1).EQ.'[').OR.(CLINE(N1:N1).EQ.'('))THEN
          AA=CLINE(N1:N1)
          IF(AA.EQ.'[')THEN
            N3=INDEX(CLINE(N1:N2),']')
          ELSE
            N3=INDEX(CLINE(N1:N2),')')
          END IF
          IF(N3.EQ.0)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2810)
2810        FORMAT('  (Instruction contains unbalanced parentheses.)')
            CALL WRTERR
            GO TO 2700
          END IF
          N3=N3+N1-1
          IF((N3-1.GT.N1).AND.(N3-1-(N1+1)+1.LE.20))TPAR=
     +    CLINE(N1+1:N3-1)
          CALL LOWCAS(TPAR)
          IF(N3-1.EQ.N1)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2815)
2815        FORMAT('  (Instruction missing observation name.)')
            CALL WRTERR
          ELSE IF(N3-1-(N1+1)+1.GT.20)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2816)
2816        FORMAT('  (Observation name greater than 20 characters ',
     +      'long.)')
            CALL WRTERR
          ELSE IF((N3-1-(N1+1)+1.EQ.3).AND.(TPAR.EQ.'dum'))
     +      THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2817)
2817        FORMAT('  (Dummy observation name "dum" only allowed for ',
     +      'non-fixed observations.)')
            CALL WRTERR
          ELSE
            IF(NOBS.EQ.0) THEN
              NOBS=1
              AOBS(1)=TPAR
              IOBS=1
            ELSE
              CALL WHICH1(IFAIL,NOBS,IOBS,AOBS,TPAR)
              IF(IFAIL.NE.0)THEN
                NOBS=NOBS+1
                IF(NOBS.GT.MAXOBS)THEN
                  WRITE(ERRMSG,2825) AMXOBS(:NBLNK(AMXOBS)),
     +            INSFLE(1:LENIN)
2825              FORMAT(' Program INSCHEK can handle only ',A,
     +            ' different observations: more than this are cited ',
     +            'in file ',A,'.')
                  CALL WRTERR
                  GO TO 9800
                END IF
                AOBS(NOBS)=TPAR
                IOBS=NOBS
              ELSE
                WRITE(ERRMSG,2830)AALINE(1:MMLINE),TPAR(:NBLNK(TPAR))
2830            FORMAT(A,'observation "',A,'" already cited.')
                CALL WRTERR
              END IF
            END IF
            IF(N2.EQ.N3)THEN
              CALL ILINT1(CLINE(N1:N2))
              WRITE(ERRMSG,2840)
              CALL WRTERR
              GO TO 2700
            END IF
            CALL GETNUM(IFAIL,CLINE,N3,N2,NUM1,NUM2,FMT)
            IF(IFAIL.NE.0)THEN
              CALL ILINT1(CLINE(N1:N2))
              WRITE(ERRMSG,2840)
2840          FORMAT('  (Cannot read numbers in format "N1:N2" after ',
     +        'observation name.)')
              CALL WRTERR
            ELSE IF(NUM1.EQ.0)THEN
              CALL ILINT1(CLINE(N1:N2))
              WRITE(ERRMSG,2841)
2841          FORMAT('  (First Number in format "N1:N2" following ',
     +        'observation name must not be zero.)')
              CALL WRTERR
            ELSE IF(NUM2.LE.NUM1)THEN
              CALL ILINT1(CLINE(N1:N2))
              WRITE(ERRMSG,2845)
2845          FORMAT('  (Numbers in format "N1:N2" following ',
     +        'observation name must be in increasing order.)')
              CALL WRTERR
            ELSE IF (NUM2.GT.MAXWID) THEN
              CALL ILINT1(CLINE(N1:N2))
              WRITE(ERRMSG,2847) AMXWID(:NBLNK(AMXWID))
              CALL WRTERR
              GO TO 2680
            ELSE
              IF(NUM1.LT.J1)THEN
                WRITE(ERRMSG,2850) AALINE(1:MMLINE)
2850            FORMAT(A,
     +          'model output file must be read from left to right.')
                CALL WRTERR
              END IF
              J1=NUM2
            END IF
          END IF
        ELSE IF(CLINE(N1:N1).EQ.'!')THEN
          IF(N1.EQ.N2)THEN
            CALL ILINT1(CLINE(N1:N2))
            GO TO 2700
          END IF
          ITEMP=INDEX(CLINE(N1+1:N2),'!')
          IF(ITEMP.EQ.0)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2860)
2860        FORMAT('  ("!"s not balanced.)')
            CALL WRTERR
            GO TO 2700
          END IF
          IF(CLINE(N2:N2).NE.'!')THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2870)
2870        FORMAT('  (Instruction must end in "!".)')
            CALL WRTERR
            GO TO 2700
          END IF
          IF(N1+1.EQ.N2)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2815)
            CALL WRTERR
          ELSE IF(N2-1-(N1+1)+1.GT.20)THEN
            CALL ILINT1(CLINE(N1:N2))
            WRITE(ERRMSG,2816)
            CALL WRTERR
          ELSE
            TPAR=CLINE(N1+1:N2-1)
            CALL LOWCAS(TPAR)
            IF(TPAR.EQ.'dum') GO TO 2700
            IF(NOBS.EQ.0) THEN
              NOBS=1
              AOBS(1)=TPAR
              IOBS=1
            ELSE
              CALL WHICH1(IFAIL,NOBS,IOBS,AOBS,TPAR)
              IF(IFAIL.NE.0)THEN
                NOBS=NOBS+1
                IF(NOBS.GT.MAXOBS)THEN
                  WRITE(ERRMSG,2825) AMXOBS(:NBLNK(AMXOBS)),
     +            INSFLE(1:LENIN)
                  CALL WRTERR
                  GO TO 9800
                END IF
                AOBS(NOBS)=TPAR
                IOBS=NOBS
              ELSE
                WRITE(ERRMSG,2830) AALINE(1:MMLINE),TPAR(:NBLNK(TPAR))
                CALL WRTERR
              END IF
            END IF
          END IF
        ELSE
          CALL ILINT1(CLINE(N1:N2))
        END IF
        GO TO 2700

2900    WRITE(ERRMSG,2910) AALINE(1:MMLINE)
2910    FORMAT(A,'cannot read line.')
        CALL WRTERR
        GO TO 9800
2950    IF(JLINE.EQ.1)THEN
          WRITE(ERRMSG,2960) AALINE(1:MMLINE)
2960      FORMAT(A,'unexpected end to file.')
          CALL WRTERR
          GO TO 9800
        END IF
        IF(IERR.NE.0) GO TO 9800

        IF(NOBS.EQ.0)THEN
          WRITE(ERRMSG,295) INSFLE(:NBLNK(INSFLE))
295       FORMAT(' No observations identified in file ',A,'.')
          CALL WRTERR
        END IF
        IF(IERR.NE.0) GO TO 9800
        WRITE(6,300,ERR=9000)
300     FORMAT(' No errors encountered.')
        IF(CHKTYP.GE.2) GO TO 500

C -- A FILE LISTING ALL OBSERVATION NAMES IS WRITTEN

        OPEN(UNIT=50,FILE=OUTFLE,STATUS='UNKNOWN',ERR=335)
        GO TO 350
335     WRITE(ERRMSG,340) OUTFLE(:NBLNK(OUTFLE)),
     +  INSFLE(:NBLNK(INSFLE))
340     FORMAT(' Cannot open file ',A,' to list observations ',
     +  'identified in instruction file ',A,'.')
        CALL WRTERR
        GO TO 9800
350     DO 370 I=1,NOBS
        WRITE(50,'(1X,A)',ERR=9700) AOBS(I)
370     CONTINUE
        CLOSE(UNIT=50)
        WRITE(6,*,ERR=9000)
        CALL WRTINT(ATEMP,NOBS)
        WRITE(ERRMSG,380) ATEMP(:NBLNK(ATEMP)),
     +  INSFLE(:NBLNK(INSFLE)),
     +  OUTFLE(:NBLNK(OUTFLE))
380     FORMAT(' ',A,' observations identified in file ',A,
     +  ': these are listed in file ',A,'.')
        CALL WRTERR
        WRITE(6,*,ERR=9000)
        IERR=IERR-1
        GO TO 9999

500     WRITE(6,*,ERR=9000)

        ALLOCATE(OBSVAL(NOBS),STAT=IFAIL)
        IF(IFAIL.NE.0)THEN
          WRITE(ERRMSG,502)
502       FORMAT(' Insufficient memory available to continue ',
     +    'INSCHEK execution.')
          CALL WRTERR
          GO TO 9820
        END IF

        DO 503 I=1,NOBS
503     OBSVAL(I)=-1.0D300

C -- IF REQUESTED, THE MODEL OUTPUT FILE IS READ

        WRITE(ERRMSG,540) MODFLE(1:LENMOD)
540     FORMAT(' Reading model output file ',A,' ----->')
        CALL WRTERR
        IERR=IERR-1
        REWIND(UNIT=40,ERR=9600)
        BEGINS=0
        CIL=0
        INS=0
550     INS=INS+1
        READ(40,'(A)',ERR=9500,END=800) BUF
        IF(INDEX(BUF,CHAR(9)).NE.0)CALL TABREM(BUF)
        NBLB=NBLNK(BUF)
        IF(NBLB.EQ.0)THEN
          IF(BEGINS.EQ.1) THEN
            BACKSPACE(UNIT=40,ERR=9550)
            BACKSPACE(UNIT=40,ERR=9550)
            INS=INS-2
          END IF
          GO TO 550
        END IF
        IF(INS.EQ.1) GO TO 550
        INSNUM=0
        CALL WRTINT(ATEMP,INS)
        AALINE=' (Instruction on line '
     +  //ATEMP(1:INDEX(ATEMP,' ')-1)//' of file '//
     +  INSFLE(1:LENIN)//'.)'
        MMLINE=NBLNK(AALINE)
        N2=0

650     CALL GETINT(IFAIL,BUF,N1,N2,NBLB,PARDEL)
        IF(IFAIL.NE.0) GO TO 9500
651     IF(N1.EQ.0) GO TO 550
        INSNUM=INSNUM+1
        IF(INSNUM.EQ.1)THEN
          IF(BUF(N1:N1).NE.'&') THEN
            MRKTYP=0
            ALLMRK=1
            BEGINS=0
          ELSE
            IF(BEGINS.EQ.1)THEN
              INS=INS-2
              BACKSPACE(UNIT=40,ERR=9550)
              BACKSPACE(UNIT=40,ERR=9550)
              GO TO 550
            END IF
          END IF
        END IF
        IF((BUF(N1:N1).EQ.'l').OR.(BUF(N1:N1).EQ.'L'))THEN
          ALLMRK=0
          WRITE(FMT,649) N2-N1
649       FORMAT('(I',I4,')')
          READ(BUF(N1+1:N2),FMT,ERR=9500) NOL
          IF(NOL.GT.1) THEN
            DO 660 I=1,NOL-1
            READ(35,*,ERR=9070,END=9100)
660         CIL=CIL+1
          END IF
          READ(35,'(A)',ERR=9070,END=9100) CLINE
          IF(INDEX(CLINE,CHAR(9)).NE.0) CALL TABREP(CLINE)
          CIL=CIL+1
          CALL WRTINT(ATEMP,CIL)
          AAALIN=' Line '//ATEMP(1:INDEX(ATEMP,' ')-1)//
     +    ' of file '//MODFLE(1:LENMOD)//': '
          MMMLIN=NBLNK(AAALIN)+1
          NBLC=NBLNK(CLINE)
          IF(NBLC.GT.MAXWID) GO TO 9050
          MRKTYP=1
          J1=0
        ELSE IF(BUF(N1:N1).EQ.PARDEL)THEN
          IF(MRKTYP.EQ.0)THEN
700         READ(35,'(A)',ERR=9070,END=9150) CLINE
            IF(INDEX(CLINE,CHAR(9)).NE.0) CALL TABREP(CLINE)
            CIL=CIL+1
            NBLC=NBLNK(CLINE)
            IF(NBLC.GT.MAXWID) GO TO 9050
            J1=INDEX(CLINE,BUF(N1+1:N2-1))
            IF(J1.EQ.0) GO TO 700
            CALL WRTINT(ATEMP,CIL)
            AAALIN=' Line '//ATEMP(1:INDEX(ATEMP,' ')-1)//
     +      ' of file '//MODFLE(1:LENMOD)//': '
            MMMLIN=NBLNK(AAALIN)+1
            J1=J1+N2-N1-2
            MRKTYP=1
          ELSE
            IF(J1.GE.NBLC) THEN
              IF(ALLMRK.EQ.1) THEN
                BEGINS=1
                N2=0
                INSNUM=0
                GO TO 650
              END IF
              GO TO 9200
            END IF
            J2=INDEX(CLINE(J1+1:NBLC),BUF(N1+1:N2-1))
            IF(J2.EQ.0) THEN
              IF(ALLMRK.EQ.1) THEN
                BEGINS=1
                N2=0
                INSNUM=0
                GO TO 650
              END IF
              GO TO 9200
            END IF
            J1=J1+J2
            J1=J1+N2-N1-2
          END IF
        ELSE IF(BUF(N1:N1).EQ.'&') THEN
          IF(INSNUM.NE.1) GO TO 9500
        ELSE IF((BUF(N1:N1).EQ.'w').OR.(BUF(N1:N1).EQ.'W'))THEN
          ALLMRK=0
          IF(J1.GE.NBLC) GO TO 9250
          J2=INDEX(CLINE(J1+1:NBLC),' ')
          IF(J2.EQ.0) GO TO 9250
          J1=J1+J2
          DO 710 I=J1,NBLC
          IF(CLINE(I:I).NE.' ') GO TO 720
710       CONTINUE
          I=NBLC+1
720       J1=I-1
        ELSE IF((BUF(N1:N1).EQ.'t').OR.(BUF(N1:N1).EQ.'T'))THEN
          ALLMRK=0
          WRITE(FMT,649) N2-N1
          READ(BUF(N1+1:N2),FMT,ERR=9500) J2
          IF(J2.LT.J1) GO TO 9500
          J1=J2
          IF(J1.GT.NBLC) THEN
            WRITE(ERRMSG,725) AAALIN(1:MMMLIN),AALINE(1:MMLINE)
725         FORMAT(A,'tab position beyond end of line.',A)
            CALL WRTERR
            GO TO 9820
          END IF
        ELSE IF((BUF(N1:N1).EQ.'[').OR.(BUF(N1:N1).EQ.'('))THEN
          ALLMRK=0
          AA=BUF(N1:N1)
          IF(AA.EQ.'[')THEN
            N3=INDEX(BUF(N1:N2),']')
          ELSE
            N3=INDEX(BUF(N1:N2),')')
          END IF
          N3=N3+N1-1
          TPAR=BUF(N1+1:N3-1)
          CALL WHICH1(IFAIL,NOBS,IOBS,AOBS,TPAR)
          IF(IFAIL.NE.0) GO TO 9500
          CALL GETNUM(IFAIL,BUF,N3,N2,NUM1,NUM2,FMT)
          IF(IFAIL.NE.0) GO TO 9500
          IF(AA.EQ.'(') THEN
           CALL GETTOT(IFAIL,CLINE,NUM1,NUM2,NBLC)
           IF(IFAIL.NE.0) GO TO 9300
          ELSE
           IF(NUM1.GT.NBLC) GO TO 9300
           IF(NUM2.GT.NBLC) NUM2=NBLC
           IF(CLINE(NUM1:NUM2).EQ.' ') GO TO 9300
          END IF
          WRITE(FMT,750) NUM2-NUM1+1
750       FORMAT('(F',I4,'.0)')
          READ(CLINE(NUM1:NUM2),FMT,ERR=9350) OBSVAL(IOBS)
          J1=NUM2
        ELSE IF(BUF(N1:N1).EQ.'!') THEN
          ALLMRK=0
          CALL LOWCAS(BUF(N1+1:N2-1))
          IF((N2-N1.NE.4).OR.(BUF(N1+1:N2-1).NE.'dum'))THEN
            TPAR=BUF(N1+1:N2-1)
            CALL WHICH1(IFAIL,NOBS,IOBS,AOBS,TPAR)
            IF(IFAIL.NE.0) GO TO 9500
            DUMFLG=0
          ELSE
            TPAR='dum'
            DUMFLG=1
          END IF
          CALL GETNXT(IFAIL,CLINE,J1,NUM1,NUM2,NBLC)
          IF(IFAIL.NE.0) GO TO 9300
          IF(CLINE(NUM2:NUM2).EQ.',') GO TO 770  !to get around a bug in lf90
          WRITE(FMT,750) NUM2-NUM1+1
          READ(CLINE(NUM1:NUM2),FMT,ERR=770) RTEMP
          IF(DUMFLG.EQ.0) OBSVAL(IOBS)=RTEMP
          J1=NUM2
          GO TO 650
770       CALL GETINT(IFAIL,BUF,N1,N2,NBLB,PARDEL)
          IF(IFAIL.NE.0) GO TO 9500
          IF(N1.EQ.0) GO TO 9350
          IF(BUF(N1:N1).NE.PARDEL) GO TO 9350
          J2=INDEX(CLINE(J1+1:NBLC),BUF(N1+1:N2-1))
          IF(J2.EQ.0) GO TO 9350
          NUM2=J1+J2-1
          IF(NUM2.LT.NUM1) GO TO 9350
          WRITE(FMT,750) NUM2-NUM1+1
          READ(CLINE(NUM1:NUM2),FMT,ERR=9350) RTEMP
          IF(DUMFLG.NE.1) OBSVAL(IOBS)=RTEMP
          J1=NUM2
          GO TO 651
        ELSE
          GO TO 9500
        END IF
        GO TO 650


800     CLOSE(UNIT=35)
        DO 805 I=1,NOBS
        IF(OBSVAL(I).LT.-1.0D299)GO TO 9500
805     CONTINUE

C -- THE OBSERVATION VALUES ARE WRITTEN TO THE INSCHEK OUTPUT FILE

        WRITE(6,300,ERR=9000)
        OPEN(UNIT=50,FILE=OUTFLE,STATUS='UNKNOWN',ERR=810)
        GO TO 830
810     WRITE(ERRMSG,811) OUTFLE(1:LENOUT),MODFLE(1:LENMOD)
811     FORMAT(' Cannot open file ',A,' to list observation values ',
     +  'read from file ',A,'.')
        CALL WRTERR
        GO TO 9820
830     DO 850 I=1,NOBS
        WRITE(50,860,ERR=9700) AOBS(I),OBSVAL(I)
860     FORMAT(1X,A,2X,1PG15.7)
850     CONTINUE
        CLOSE(UNIT=50)
        WRITE(6,*,ERR=9000)
        CALL WRTINT(ATEMP,NOBS)
        WRITE(ERRMSG,870) ATEMP(:NBLNK(ATEMP)),INSFLE(1:LENIN),
     +  OUTFLE(1:LENOUT),MODFLE(1:LENMOD)
870     FORMAT(' ',A,' observations identified in file ',A,
     +  ': these are listed in file ',A,' together with their values ',
     +  'as read from file ',A,'.')
        CALL WRTERR
        IERR=IERR-1
        WRITE(6,*,ERR=9000)
        GO TO 9999

C -- SOME ERROR MESSAGES .....

9000    CALL EXIT(100)
9050    WRITE(ERRMSG,9060) AMXWID(:NBLNK(AMXWID))
9060    FORMAT(' Model output file width exceeds ',A,' characters.')
        CALL WRTERR
        GO TO 9820
9070    WRITE(ERRMSG,9080) AALINE(1:MMLINE)
9080    FORMAT(' Inexplicable error while reading model output file ',
     +  A,'.')
        GO TO 9820
9100    WRITE(ERRMSG,9110) AALINE(1:MMLINE)
9110    FORMAT(' Line advance moves cursor beyond end of file.',A)
        CALL WRTERR
        GO TO 9820
9150    WRITE(ERRMSG,9160) AALINE(1:MMLINE)
9160    FORMAT(' End of file encountered while searching for ',
     +  'marker.',A)
        CALL WRTERR
        GO TO 9820
9200    WRITE(ERRMSG,9210) AAALIN(1:MMMLIN),AALINE(1:MMLINE)
9210    FORMAT(A,'unable to find secondary marker.',A)
        CALL WRTERR
        GO TO 9820
9250    WRITE(ERRMSG,9260) AAALIN(1:MMMLIN),AALINE(1:MMLINE)
9260    FORMAT(A,'unable to find requested whitespace or ',
     +  'whitespace preceeds end of line.',A)
        CALL WRTERR
        GO TO 9820
9300    WRITE(ERRMSG,9310) AAALIN(1:MMMLIN),TPAR(:NBLNK(TPAR)),
     +  AALINE(1:MMLINE)
9310    FORMAT(A,'cannot locate observation "',A,'".',A)
        CALL WRTERR
        GO TO 9820
9350    WRITE(ERRMSG,9360) AAALIN(1:MMMLIN),TPAR(:NBLNK(TPAR)),
     +  AALINE(1:MMLINE)
9360    FORMAT(A,'error reading observation "',A,'".',A)
        CALL WRTERR
        GO TO 9820
9500    WRITE(ERRMSG,9510) INSFLE(:NBLNK(INSFLE))
9510    FORMAT(' Error re-reading instruction file ',A,'.')
        CALL WRTERR
        GO TO 9820
9550    WRITE(ERRMSG,9560) INSFLE(:NBLNK(INSFLE))
9560    FORMAT(' Unable to backspace in file ',A,'.')
        CALL WRTERR
        GO TO 9820
9600    WRITE(ERRMSG,9610) INSFLE(:NBLNK(INSFLE))
9610    FORMAT(' Cannot rewind instruction file ',A,'.')
        CALL WRTERR
        GO TO 9820
9700    WRITE(ERRMSG,9710) OUTFLE(:NBLNK(OUTFLE))
9710    FORMAT(' Cannot write to file ',A)
        CALL WRTERR
        GO TO 9800
9800    WRITE(6,*,ERR=9000)
C        OPEN(UNIT=36,FILE='pest.tmp',STATUS='UNKNOWN',ERR=9000)
C        WRITE(36,'(a)',ERR=9000)'1'
C        CLOSE(UNIT=36,ERR=9000)
9820    CALL EXIT(50)
9900    WRITE(6,9910,ERR=9000)
9910    FORMAT(' INSCHEK is run using the command:',/)
        WRITE(6,9920,ERR=9000)
9920    FORMAT('     INSCHEK insfile [modfile]',/,/,
     +  ' where',/)
        WRITE(6,9930,ERR=9000)
9930    FORMAT('     "insfile" is a PEST instruction file, and')
        WRITE(6,9940,ERR=9000)
9940    FORMAT('     "modfile" is an [optional] model output file ',
     +  'to be read by INSCHEK.')
        GO TO 9800

9999    CONTINUE
C        OPEN(UNIT=36,FILE='pest.tmp',STATUS='UNKNOWN',ERR=9000)
C        WRITE(36,'(a)',ERR=9000)'0'
C        CLOSE(UNIT=36)
9998    CALL EXIT(0)
        END


        SUBROUTINE EXTCHK(EXTLEN,LENIN,INFLE)

C -- SUBROUTINE EXTCHK DETERMINES THE LOCATION OF THE FILENAME
C -- EXTENSION

        INTEGER EXTLEN,LENIN,I
        CHARACTER*(*) INFLE

        DO 10 I=LENIN,1,-1
        IF(INFLE(I:I).EQ.'.') GO TO 30
#ifdef UNIX
        IF(INFLE(I:I).EQ.'/') GO TO 50
#else
        IF(INFLE(I:I).EQ.CHAR(92)) GO TO 50
#endif
10      CONTINUE
        EXTLEN=LENIN+1
        RETURN
30      EXTLEN=I
        RETURN
50      EXTLEN=LENIN+1
        RETURN
        END


        SUBROUTINE WRTERR

C -- SUBROUTINE WRTERR FORMATS AND WRITES AN ERROR MESSAGE

        INTEGER NBLNK
        INTEGER IERR,J,JEND,I,NBLC,ITAKE,MCASE
        CHARACTER*2120 ERRMSG

        COMMON /IERR/IERR,MCASE
        COMMON /ERRMSG/ERRMSG

        IERR=IERR+1
        IF(IERR.GT.100)THEN
          WRITE(6,*,ERR=200)
C          OPEN(UNIT=36,FILE='pest.tmp',STATUS='unknown',ERR=200)
C          WRITE(36,'(a)',ERR=200)'1'
C          CLOSE(UNIT=36,ERR=200)
          CALL EXIT(50)
        END IF
        ITAKE=0
        J=0
        NBLC=NBLNK(ERRMSG)
5       JEND=J+78-ITAKE
        IF(JEND.GE.NBLC) GO TO 100
        DO 10 I=JEND,J+1,-1
        IF(ERRMSG(I:I).EQ.' ') THEN
          IF(ITAKE.EQ.0) THEN
            WRITE(6,'(a)',ERR=200) ERRMSG(J+1:I)
            ITAKE=3
          ELSE
            WRITE(6,'(a)',ERR=200) '   '//ERRMSG(J+1:I)
          END IF
          J=I
          GO TO 5
        END IF
10      CONTINUE
        IF(ITAKE.EQ.0)THEN
          WRITE(6,'(a)',ERR=200) ERRMSG(J+1:JEND)
          ITAKE=3
        ELSE
          WRITE(6,'(a)',ERR=200) '   '//ERRMSG(J+1:JEND)
        END IF
        J=JEND
        GO TO 5
100     JEND=NBLC
        IF(ITAKE.EQ.0)THEN
          WRITE(6,'(a)',ERR=200) ERRMSG(J+1:JEND)
        ELSE
          WRITE(6,'(a)',ERR=200) '   '//ERRMSG(J+1:JEND)
        END IF
        RETURN

200     CALL EXIT(100)
        END



        SUBROUTINE WHICH1(IFAIL,NPAR,IPAR,APAR,TPAR)

C -- SUBROUTINE WHICH1 LOCATES A PARTICULAR STRING IN AN ARRAY OF STRINGS

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
20        CONTINUE
        END IF
        IF(IPAR.NE.1)THEN
          DO 40 I=IPAR-1,1,-1
          IF(TPAR.EQ.APAR(I)) THEN
            IPAR=I
            RETURN
          END IF
40        CONTINUE
        END IF
        IFAIL=1
        RETURN
        END



        SUBROUTINE LOWCAS(ASTRNG)

C -- SUBROUTINE LOWCAS CONVERTS A STRING TO LOWER CASE

        INTEGER NBLNK
        INTEGER I,J
        CHARACTER*(*) ASTRNG

        DO 10 I=1,NBLNK(ASTRNG)
        J=ICHAR(ASTRNG(I:I))
        IF((J.GE.65).AND.(J.LE.90)) ASTRNG(I:I)=CHAR(J+32)
10      CONTINUE
        RETURN
        END


        SUBROUTINE WRTINT(ATEMP,IVAL)

C -- SUBROUTINE WRITINT WRITES AN INTEGER TO A STRING

        INTEGER IVAL
        CHARACTER*(*) ATEMP

        WRITE(ATEMP,'(i10)')IVAL
        CALL SHIFTL(ATEMP)
        RETURN
        END



        SUBROUTINE GETNUM(IFAIL,BUF,N3,N2,NUM1,NUM2,FMT)

C -- SUBROUTINE GETNUM DETERMINES THE CHARACTER POSITIONS OF FIXED
C -- AND SEMI-FIXED OBSERVATIONS

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
        END


        SUBROUTINE ILINT1(AINSTR)

C -- SUBROUTINE ILLINT1 ASSISTS IN WRITING ERROR MESSAGES

        INTEGER MMLINE
        CHARACTER*200 ERRMSG*2120,AALINE
        CHARACTER*(*) AINSTR

        COMMON /ERRMSG/ ERRMSG
        COMMON /AALINE/ AALINE
        COMMON /MMLINE/MMLINE

        WRITE(ERRMSG,10) AALINE(1:MMLINE),AINSTR
10      FORMAT(A,'illegal instruction - "',A,'".')
        CALL WRTERR

        RETURN
        END



        SUBROUTINE GETNXT(IFAIL,CLINE,J1,NUM1,NUM2,NBLC)

C -- SUBROUTINE GETNXT RETRIEVES THE NEXT SPACE-DELIMITED SUBSTRING

        INTEGER IFAIL
        INTEGER J1,NUM1,NUM2,NBLC,I,II
        CHARACTER*(*) CLINE

        IFAIL=0
        DO 20 I=J1+1,NBLC
        IF(CLINE(I:I).NE.' ') GO TO 50
20      CONTINUE
        IFAIL=1
        RETURN
50      NUM1=I
        I=INDEX(CLINE(NUM1:NBLC),' ')
        II=INDEX(CLINE(NUM1:NBLC),',')
        IF(II.NE.0)THEN
          IF(I.EQ.0)THEN
            I=II
          ELSE
            I=MIN(I,II)
          END IF
        END IF
        IF(I.EQ.0) THEN
          NUM2=NBLC
        ELSE
          NUM2=NUM1+I-2
        END IF

        RETURN
        END


        SUBROUTINE GETTOT(IFAIL,CLINE,J1,J2,NBLC)

C -- SUBROUTINE GETTOT DETERMINES THE SPACE OCCUPIED BY A NUMBER ON
C -- A MODEL OUTPUT FILE

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

        END


        SUBROUTINE SHIFTL(AA)

C -- SUBROUTINE SHIFTL REMOVES LEADING BLANK CHARACTERS FROM A
C -- STRING

        INTEGER L,I,J,II
        CHARACTER*(*) AA

        L=LEN(AA)
        DO 10 I=1,L
        IF((AA(I:I).NE.' ').AND.(ICHAR(AA(I:I)).NE.9)) GO TO 50
10      CONTINUE
        RETURN
50      IF(I.EQ.1) RETURN
        II=I-1
        DO 100 J=I,L
100     AA(J-II:J-II)=AA(J:J)
        DO 110 J=1,II
110     AA(L-J+1:L-J+1)=' '
        RETURN
        END


        SUBROUTINE GETINT(IFAIL,BUF,N1,N2,NBLB,MRKDEL)

C -- SUBROUTINE GETINT RETRIEVES THE NEXT INSTRUCTION

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
        END


        SUBROUTINE TABREM(CLINE)

C -- SUBROUTINE TABREM REMOVES TABS FROM A STRING

        INTEGER I
        CHARACTER*(*) CLINE

        DO 10 I=1,LEN(CLINE)
10      IF(ICHAR(CLINE(I:I)).EQ.9) CLINE(I:I)=' '

        RETURN
        END


        SUBROUTINE TABREP(CLINE)

C -- SUBROUTINE TABREP REPLACES A TAB WITH THE CORRECT NUMBER OF
C -- BLANK SPACES

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

        END

#ifdef CAPFILE
        SUBROUTINE UPCAS(ASTRNG)

C -- SUBROUTINE UPCAS CONVERTS A STRING TO UPPER CASE

        INTEGER NBLNK
        INTEGER I,J
        CHARACTER*(*) ASTRNG

        DO 10 I=1,NBLNK(ASTRNG)
        J=ICHAR(ASTRNG(I:I))
        IF((J.GE.97).AND.(J.LE.122)) ASTRNG(I:I)=CHAR(J-32)
10      CONTINUE
        RETURN
        END
#endif

