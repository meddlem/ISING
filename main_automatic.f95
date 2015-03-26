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
  real(dp)              :: BJ, alpha, chi, Cv, h = 0._dp
  integer, allocatable  :: S(:,:), m(:), t(:)
  integer               :: runtime, L, N, r_max, n_corr, i, j

  integer, parameter  :: lattice_lengths(1:7) = (/8, 16, 32, 64, 128, 256, 512/)
  integer, parameter  :: N_BJ_STEPS = 10
  real(dp) :: spin_coupling_energies(N_BJ_STEPS) 
  call linspace(spin_coupling_energies, 0d0,1d0,N_BJ_STEPS)

  do i = 1, size(lattice_lengths)
  do j = 1, N_BJ_STEPS

    L = lattice_lengths(i)
    BJ = spin_coupling_energies(j)

    call init_constant_variables(L,N,r_max,n_corr) ! variables scaling with lattice size and energy

    allocate(S(L,L),m(n_meas),t(n_meas),BE(n_meas),c_ss(r_max),&
      c_ss_fit(r_max),r(r_max))

    call init_random_seed()
    call init_lattice(S,L)
    call animate_lattice('')
    
    call run_sim(S,L,r_max,n_corr,BE,BJ,h,t,r,m,runtime,c_ss,c_ss_fit,alpha,&
      chi, Cv)
    
    call close_lattice_plot()
    call results_out(BJ,BE(n_meas),h,runtime,alpha, chi, Cv)
    call line_plot(real(t,dp),BE,'t','energy','','',1)
    call line_plot(real(t,dp),real(m,dp),'t','magnetization','','',2)
    call line_plot(r,c_ss,'r','corr','corr','',3,c_ss_fit,'fit')
    
    deallocate(S,m,t,r,BE,c_ss,c_ss_fit)

  end do
  end do

contains

  subroutine init_constant_variables(L,N,r_max,n_corr) 
    integer, intent(in)  :: L
    integer, intent(out)  :: r_max, n_corr, N
    N = L**2
    n_corr = L/3 ! number of spins used to calculate correlation
    r_max = L/4 ! distances over which to calc correlation function
  end subroutine

  subroutine linspace(output_var, start_var, stop_var, N_steps)
      integer, intent(in)   :: N_steps
      real(dp), intent(in)  :: start_var, stop_var
      integer               :: i
      real(dp), dimension(N_steps), intent(out) :: output_var
      forall (i=1:N_steps) output_var(i) = (stop_var - start_var)/N_steps*i
  end subroutine

end program
