program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  ! variables:
  ! BE: beta*energy,
  ! dE: change in energy between new and initial config
  ! S: array containing Spins indexed as row, column

  real(dp), allocatable :: BE(:), c_ss(:), r(:), c_ss_fit(:)
  real(dp)              :: BJ, nu, err_nu, chi_s, chi, err_chi, Mag, &
                           err_Mag, Cv, err_Cv
  integer, allocatable  :: S(:,:)
  integer               :: runtime, L, method, r_max, n_corr
  logical               :: calc_css
  
  call get_usr_args(method,calc_css) 
  call user_in(BJ,L,r_max,n_corr)
  allocate(S(L,L),BE(n_meas),c_ss(r_max),c_ss_fit(r_max),r(r_max))
  call init_random_seed()
  call init_lattice(S,L)

  call run_sim(S,method,r_max,n_corr,BE,BJ,r,Mag,err_Mag,runtime,&
    calc_css,c_ss,c_ss_fit,nu,err_nu,chi_s,chi,err_chi,Cv,err_Cv)
  
  call results_out(BE,BJ,L,r,runtime,calc_css,c_ss,c_ss_fit,nu,err_nu,chi_s, &
    chi,err_chi,Mag,err_Mag,Cv,err_Cv)
  deallocate(S,r,BE,c_ss,c_ss_fit)
end program
