program main
  use constants
  use initialize
  use markov
  use plotroutines
  use io
  implicit none
  
  ! variables:
  ! BE: beta*energy,
  ! dE: change in energy between new and initial config
  ! h: external field
  ! S: array containing Spins indexed as row, column

  real(dp), allocatable :: BE(:), c_ss(:), r(:), c_ss_fit(:)
  real(dp)              :: BJ, nu, chi, Mag, Cv, h = 0._dp
  integer, allocatable  :: S(:,:), m(:), t(:)
  integer               :: runtime, L, N, r_max, n_corr
  
  call user_in(BJ,L,N,r_max,n_corr)
  allocate(S(L,L),m(n_meas),t(n_meas),BE(n_meas),c_ss(r_max),&
    c_ss_fit(r_max),r(r_max))

  call init_random_seed()
  call init_lattice(S,L)
  call animate_lattice('')
  
  call run_sim(S,L,r_max,n_corr,BE,BJ,h,t,r,m,Mag,runtime,c_ss,c_ss_fit,nu,&
    chi, Cv)
  
  call close_lattice_plot()
  call results_out(BJ,h,runtime,nu,chi,Mag,Cv)
  call line_plot(real(t,dp),BE,'t','energy','','',1)
  call line_plot(real(t,dp),real(m,dp),'t','magnetization','','',2)
  call line_plot(r,c_ss,'r','corr','corr','',3,c_ss_fit,'fit')
  
  deallocate(S,m,t,r,BE,c_ss,c_ss_fit)
end program
