        !COMPILER-GENERATED INTERFACE MODULE: Sat Oct 21 14:25:28 2017
        MODULE RK4__genmod
          INTERFACE 
            SUBROUTINE RK4(N,X,H,Y,RP)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: X
              REAL(KIND=4) :: H
              REAL(KIND=4) :: Y(5)
              EXTERNAL RP
            END SUBROUTINE RK4
          END INTERFACE 
        END MODULE RK4__genmod
