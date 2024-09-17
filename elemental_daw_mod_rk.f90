    Module elemental_daw_mod_rk
      ! :----------------
      ! This module provides generic interfaces for the Dawson function of a
      ! real
      ! argument,  Dawson(x)=EXP(-x^2)*INT(0 <= t <= x) EXP(t^2 ) dt, in any
      ! of the
      ! standard precision arithmetic (single, double or quad) based on
      ! the choice of an integer "rk" in the subsidiary module "set_rk" in the
      ! file "set_rk.f90.
      ! :------------------

      Use set_rk
      Implicit None

      Include 'cheb_t100_daw_parameters_sdq.f90'

      Public :: daw_rk
    Contains
      ! -------------------------
      Elemental Function daw_rk(x)
        ! Daw_rk evaluates the Dawson integral of a real argument,x;
        ! Dawson(x)=EXP(-x^2)*INT(0 <= t <= x) EXP(t^2 ) dt
        ! THE CALLING SEQUENCE FOR THIS FUNCTION IS
        ! y = Daw_rk(x)
        ! The function can be evaluated using single, double and quadruple
        ! precision depending on the choice of the integer "rk" in the
        ! associated module "set_rk"
        ! :-
        ! Developed in March 2023
        ! by Mofreh R. Zaghloul
        ! United Arab Emirates University
        ! AlAin, Abu Dhabi, UAE
        ! !---------------

        Implicit None
        Real (rk), Intent (In) :: x
        Real (rk) :: ax, daw_rk
        Real (rk) :: y100, t
        Integer :: j, ycase, jmin


        ax = abs(x)

        If (big_border_12<ax) Then
          daw_rk = half/x
        Else If (ax<=xsmall) Then
          daw_rk = x
        Else If (ax<=small_border_1) Then
          daw_rk = priv_ncont_frac0(1, x)
        Else If (ax<=small_border_2) Then
          daw_rk = priv_ncont_frac0(2, x)
        Else If (ax<=small_border_3) Then
          daw_rk = priv_ncont_frac0(3, x)
        Else If (ax<=small_border_4) Then
          daw_rk = priv_ncont_frac0(4, x)
        Else If (ax<=small_border_5) Then
          daw_rk = priv_ncont_frac0(5, x)
        Else If (ax<=small_border_6) Then
          daw_rk = priv_ncont_frac0(6, x)
        Else If (ax<=small_border_7) Then
          daw_rk = priv_ncont_frac0(7, x)
          
           ! :- Chebyshev subinterval polynomial approximation
        Else If (ax<=cheb_ax) Then
          y100 = one80/(ax+one_pt_8)
          ycase = y100
          t = two*y100 - (two*ycase+one)
          jmin = ycase*np_pls_1 + 1
          daw_rk = t*cffs(jmin)
         
          Do j = jmin + 1, jmin + np_min_1
            daw_rk = t*(daw_rk+(cffs(j)))
          End Do
          daw_rk = (daw_rk+cffs(jmin+np))*((1-c_ax)*ax+c_ax) & ! c_ax=(rk/qp)+(sp/rk))
            *sign(one, x)

        Else If (ax<=big_border_1) Then
          daw_rk = priv_ncont_frac(12, x)
        Else If (ax<=big_border_2) Then
          daw_rk = priv_ncont_frac(11, x)
        Else If (ax<=big_border_3) Then
          daw_rk = priv_ncont_frac(10, x)
        Else If (ax<=big_border_4) Then
          daw_rk = priv_ncont_frac(9, x)
        Else If (ax<=big_border_5) Then
          daw_rk = priv_ncont_frac(8, x)
        Else If (ax<=big_border_6) Then
          daw_rk = priv_ncont_frac(7, x)
        Else If (ax<=big_border_7) Then
          daw_rk = priv_ncont_frac(6, x)
        Else If (ax<=big_border_8) Then
          daw_rk = priv_ncont_frac(5, x)
        Else If (ax<=big_border_9) Then
          daw_rk = priv_ncont_frac(4, x)
        Else If (ax<=big_border_10) Then
          daw_rk = priv_ncont_frac(3, x)
        Else If (ax<=big_border_11) Then
          daw_rk = priv_ncont_frac(2, x)
        Else If (ax<=big_border_12) Then
          daw_rk = x/(two*x*x-one)
        End If

      End Function

      
      !:- Laplace continued fraction for large values of x
      !
      Elemental Function priv_ncont_frac(m, x) Result (y)
        Implicit None
        Real (rk), Intent (In) :: x
        Integer, Intent (In) :: m
        Real (rk) :: y
        Integer :: k

        y = real(m, kind=rk)/x
        Do k = m - 1, 1, -1
          y = real(k, kind=rk)/(x-half*y)
        End Do
        y = one/(two*x-y)

        Return
      End Function

      
      !:- Continued fraction expansion for |x| near the origin
      ! 
      Elemental Function priv_ncont_frac0(m, x) Result (y)
        Implicit None
        Real (rk), Intent (In) :: x
        Integer, Intent (In) :: m
        Real (rk) :: y, x_sqr
        Integer :: k

        x_sqr = x*x
        y = lc_0_cff(m)*x_sqr
        Do k = m - 1, 1, -1
          y = lc_0_cff(k)*x_sqr/(one+y)
        End Do
        y = x/(one+y)
        Return

      End Function
    End Module
