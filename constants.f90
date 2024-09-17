  ! Constants
Real(rk),Parameter ::                                              &
  two_by_sqrt_pi   = 1.128379167095512573896158903121545171688_rk, & !2/sqrt(pi)
  one_by_sqrt_pi   = 0.564189583547756286948079451560772585844_rk, & !1/sqrt(pi)
  one    = 1.0_rk,                                                 &
  two    = 2.0_rk,                                                 &
  zero   = 0.0_rk,                                                 &
  half   = 0.5_rk,                                                 &
  two_p_one = 2.1_rk,                                              &
  eps    = EPSILON( one ),                                         &
  xhuge  = HUGE( one ),                                            &
  sqrt_lg_xhuge  = SQRT(LOG(xhuge)),                               &
  zeps12  = eps**(one/12.0_rk)
