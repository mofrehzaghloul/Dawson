    ! Module set_rk provides the kind type parameter needed to define
    ! the precision of a complete package along with values
    ! for commonly used standard precisions

    Module set_rk
      Use iso_fortran_env, Only: real32, real64, real128
      Implicit None
      Public :: rk, sp, dp, qp

      ! Precisions
      Integer, Parameter :: sp = real32  ! Single Precision
      Integer, Parameter :: dp = real64  ! Double Precision
      Integer, Parameter :: qp = real128 ! Quadruple Precision

      ! Set the precision for the whole package
      Integer, Parameter :: rk = qp

      ! To change the default package-precision change
      ! the parameter assigned to rk above to the required
      ! precision and recompile the complete package.

    End Module
