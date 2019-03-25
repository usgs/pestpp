      SUBROUTINE SPACESUB(ASTRING)

      INTEGER I,J,K,N
      CHARACTER*1 BB
      CHARACTER*(*) ASTRING

      BB=CHAR(211)
      N=LEN_TRIM(ASTRING)
      K=1
10    CONTINUE
      IF(K.GT.N) GO TO 100
      DO I=K,N
        IF((ASTRING(I:I).EQ.'''').OR.(ASTRING(I:I).EQ.'"'))THEN
          ASTRING(I:I)=' '
          DO J=I+1,N
            IF((ASTRING(J:J).EQ.'''').OR.(ASTRING(J:J).EQ.'"'))THEN
              ASTRING(J:J)=' '
              K=J+1
              GO TO 10
            END IF
            IF(ASTRING(J:J).EQ.' ')ASTRING(J:J)=BB
          END DO
          GO TO 100
        END IF
      END DO

100   CONTINUE
      RETURN
#ifdef PESTMOD
      END SUBROUTINE SPACESUB
#else
      END
#endif


      SUBROUTINE REMCHAR(ASTRING,AA)

      INTEGER I,N
      CHARACTER*(*) ASTRING,AA

      N=LEN_TRIM(ASTRING)
      DO I=1,N
        IF(ASTRING(I:I).EQ.AA)ASTRING(I:I)=' '
      END DO

      RETURN
#ifdef PESTMOD
      END SUBROUTINE REMCHAR
#else
      END
#endif

