


#ifndef NO_CMDLINE

#ifdef LAHEY
        SUBROUTINE PGETCL(COMLIN)

        CHARACTER COMLIN*(*)

        CALL GETCL(COMLIN)
        RETURN
        END
#else
        SUBROUTINE PGETCL(COMLIN)
        CHARACTER*(*) COMLIN
        INTEGER IARGC
        INTEGER LLEN,NARG,IB,I,NB,IBB
        CHARACTER*120 ARG(10)

        LLEN=LEN(COMLIN)
        NARG=IARGC()
        COMLIN=' '
        IF(NARG.EQ.0) RETURN
        IB=0
        DO 100 I=1,MIN(NARG,10)
          CALL GETARG(I,ARG(I))
          NB=LEN_TRIM(ARG(I))
          IBB=MIN(IB+NB+1,LLEN)
          COMLIN(IB+1:IBB)= ARG(I)(1:NB)
          IB=IBB
          IF(IB.GE.LLEN)RETURN
100     CONTINUE
        RETURN
        END
#endif

#else
        SUBROUTINE DUMMY
        RETURN
        END
#endif
C-- EOF

