        !COMPILER-GENERATED INTERFACE MODULE: Mon Jul 19 23:01:09 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FFTPACK_RFFTB__genmod
          INTERFACE 
            SUBROUTINE FFTPACK_RFFTB(N,R,WSAVE)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: R(N)
              REAL(KIND=4) :: WSAVE(2*N+15)
            END SUBROUTINE FFTPACK_RFFTB
          END INTERFACE 
        END MODULE FFTPACK_RFFTB__genmod
