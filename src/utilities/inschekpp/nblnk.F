
        INTEGER FUNCTION NBLNK(ASTRNG)

        CHARACTER*(*) ASTRNG 

        NBLNK=LEN_TRIM(ASTRNG)
        RETURN
#ifdef PESTMOD
        END FUNCTION NBLNK
#else        
        END
#endif        


