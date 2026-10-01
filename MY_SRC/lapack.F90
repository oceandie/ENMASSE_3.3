MODULE lapack
   !! Explicit interfaces to the LAPACK routines called in splines.F90
   IMPLICIT NONE
   PUBLIC

   INTERFACE

      SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE par_kind, ONLY: wp
         INTEGER            INFO, LDA, LDB, N, NRHS
         INTEGER            IPIV( * )
         REAL(wp)           A( LDA, * ), B( LDB, * )
      END SUBROUTINE DGESV

      SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
         USE par_kind, ONLY: wp
         INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
         INTEGER            IPIV( * )
         REAL(wp)           AB( LDAB, * ), B( LDB, * )
      END SUBROUTINE DGBSV

   END INTERFACE

END MODULE lapack
