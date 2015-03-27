program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  ! variables:
  ! BE: beta*energy,
  ! dE: change in energy between new and initial config
  ! h: external field
  ! S: array containing Spins indexed as row, column

  real(dp), allocatable :: BE(:), c_ss(:), r(:), c_ss_fit(:)
  real(dp)              :: BJ, nu, chi, Mag, Cv, h = 0._dp
  integer, allocatable  :: S(:,:), t(:)
  integer               :: runtime, L, method, r_max, n_corr
  
  call user_in(BJ,L,method,r_max,n_corr)
  allocate(S(L,L),t(n_meas),BE(n_meas),c_ss(r_max),c_ss_fit(r_max),r(r_max))

  call init_random_seed()
  call init_lattice(S,L)
  call run_sim(S,L,method,r_max,n_corr,BE,BJ,h,t,r,Mag,runtime,c_ss,&
    c_ss_fit,nu,chi,Cv)
  call results_out(BE,BJ,t,r,h,runtime,c_ss,c_ss_fit,nu,chi,Mag,Cv)
  
  deallocate(S,t,r,BE,c_ss,c_ss_fit)
end program
