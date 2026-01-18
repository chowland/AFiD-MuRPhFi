        !COMPILER-GENERATED INTERFACE MODULE: Wed Feb 19 10:14:01 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PRIMEFACTORS__genmod
          INTERFACE 
            SUBROUTINE PRIMEFACTORS(NUM,FACTORS,NFACT)
              INTEGER(KIND=4), INTENT(IN) :: NUM
              INTEGER(KIND=4), INTENT(OUT) :: FACTORS(*)
              INTEGER(KIND=4), INTENT(INOUT) :: NFACT
            END SUBROUTINE PRIMEFACTORS
          END INTERFACE 
        END MODULE PRIMEFACTORS__genmod
