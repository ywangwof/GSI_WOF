        !COMPILER-GENERATED INTERFACE MODULE: Mon Jul 19 23:01:08 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FFTPACK_RADF3__genmod
          INTERFACE 
            SUBROUTINE FFTPACK_RADF3(IDO,L1,CC,CH,WA1,WA2)
              INTEGER(KIND=4) :: L1
              INTEGER(KIND=4) :: IDO
              REAL(KIND=4) :: CC(IDO,L1,3)
              REAL(KIND=4) :: CH(IDO,3,L1)
              REAL(KIND=4) :: WA1(IDO)
              REAL(KIND=4) :: WA2(IDO)
            END SUBROUTINE FFTPACK_RADF3
          END INTERFACE 
        END MODULE FFTPACK_RADF3__genmod
