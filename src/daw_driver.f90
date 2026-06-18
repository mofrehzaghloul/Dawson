  ! -----------------------------------------------------------------------
  ! Daw_Driver is a driving code for thr elementary function Daw_rk(x)
  ! included in the module elemental_daw_mod_rk found in the file
  ! elemental_daw_mod_rk.f90.
  ! --
  ! For compilation:
  ! gfortran -O3 set_rk.f90 elemental_daw_mod_rk.f90
  ! Daw_Driver.f90 -o Daw_Driver.exe
  ! !  --warn-all -Wmissing-declarations -std=c90 -C
  ! ! or
  ! ifort -O3 set_rk.f90 elemental_daw_mod_rk.f90
  ! Daw_Driver.f90 -o Daw_Driver.exe
  ! ! or
  ! use the accompanying "makefile"
  ! -----------------------------------------------------------------------
  Program daw_driver
    Use set_rk, Only: rk, qp
    Use elemental_daw_mod_rk, Only: daw_rk
    Implicit None

    Integer, Parameter :: mmax = 400001, nrep = 10
    Real (rk), Dimension (mmax) :: xref, yref
    Real (rk), Dimension (mmax) :: yp
    Real (rk) :: a, b, time_begin, time_end
    Integer :: i, n
    Character (Len=45) :: filename


    ! Accuracy & Efficiency Comparison
    filename = 'daw_ref_values.txt'
    Open (Unit=27, File=filename, Status='unknown')

    Do i = 1, 400001
      Read (27, *) xref(i), yref(i)
    End Do

    Write (*, *) ' '
    Write (*, *) 'PRECISION =', rk
    Write (*, *) '**********************************************************'
    Write (*, '(A38,ES15.7E3,A3,ES15.7E3,A8)') ' !!**** Accuracy check using data in [', xref(1), &
      ' , ', xref(400001), '] ****!!'
    Write (*, *) '**********************************************************'


    ! :-Accuracy check against externally generated reference data
    Call cpu_time(time_begin)

    Do n = 1, nrep
      yp = daw_rk(xref)
    End Do

    Call cpu_time(time_end)
    Write (*, *) 'Present    : Daw_rk'
    Write (*, '(A37,ES15.7E3)') 'Max Re. Error |y-yref|/|yref|:  = ', maxval(abs(yp-yref)/abs(yref))

    ! :-
    ! Efficiency Comparison and Timing Experiments
    ! :-
    a = 1.0E-6_rk                  ! x_min
    b = 1.0E+6_rk                  ! x_max
    Call check_daw(a, b)

  Contains
    ! -----
    Subroutine check_daw(xmin, xmax)
      ! -----
      ! This subroutine receives the smallest(xmin) and largest(xmax) values
      ! of the real input and returns the average CPU time consumed in
      ! evaluating the Dawson integral for an array of mmax values of
      ! the input "xref" using the present code 
      ! The integer "mmax" is defined in the main program.
      ! The values of the input variable are equally spaced on the
      ! logarithmic scale
      ! -----

      Real (rk), Intent (In) :: xmin, xmax
      Real (rk) :: lg10_xmin, lg10_xmax
      Integer :: k, n

      Write (*, *) ' '
      Write (*, *) ' '
      Write (*, *) ' '
      Write (*, *) '*********************************************************'
      Write (*, '(A34,ES15.7E3,A3,ES15.7E3,A8)') ' !!**** Efficiency Comparison in [', xmin, ' , ', &
        xmax, '] ****!!'
      Write (*, *) '*********************************************************'

      If (xmin>0.0_rk .And. xmax>0.0_rk) Then
        lg10_xmin = log10(xmin)
        lg10_xmax = log10(xmax)
        Do k = 1, mmax
          xref(k) = 10.0_rk**(lg10_xmin+real((k-1),kind=rk)*(lg10_xmax-lg10_xmin)/real((mmax- &
            1),kind=rk))
        End Do
      Else
        Do k = 1, mmax
          xref(k) = xmin + real((k-1), kind=rk)*(xmax-xmin)/real((mmax-1), kind=rk)
        End Do
      End If

      Write (*, *) 'nrep=', nrep

      Write (*, *) '    '
      Call cpu_time(time_begin)
      Do n = 1, nrep
        yp = daw_rk(xref)
      End Do
      Call cpu_time(time_end)
      Write (*, '(A37,ES15.7E3,A10)') ' Elapsed time: Present Daw_rk      =   ', &
        (time_end-time_begin)/nrep, ' seconds'

    End Subroutine
  End Program
