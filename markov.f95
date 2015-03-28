module markov 
  use constants
  use process_data 
  use Wolff
  use Swendsen_Wang
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(S,L,method,r_max,n_corr,BE,BJ,h,t,r,Mag,runtime, &
      calc_css,c_ss,c_ss_fit,nu,chi,Cv)
    integer, intent(inout)  :: S(:,:)
    real(dp), intent(inout) :: BE(:), BJ
    integer, intent(in)     :: L, method, r_max, n_corr
    logical, intent(in)     :: calc_css
    real(dp), intent(in)    :: h
    integer, intent(out)    :: t(:), runtime
    real(dp), intent(out)   :: c_ss(:), r(:), c_ss_fit(:), Mag, nu, chi, Cv

    integer(lng), allocatable :: N_SW(:), m(:)
    real(dp), allocatable     :: g(:,:)
    integer(lng) :: start_time, m_tmp, N_SW_tmp, end_time
    integer      :: i, j 
    real(dp)     :: p
    
    allocate(g(n_meas,r_max),N_SW(n_meas),m(n_meas))
    ! initialize needed variables
    j = 0
    t = (/(i,i=0,n_meas-1)/)
    r = real((/(i,i=1,r_max)/),dp)
    p = 1 - exp(-2._dp*BJ)
    N_SW_tmp = 0
    m_tmp = 0

    call animate_lattice()
    
    call system_clock(start_time)
    do i=1,steps
      call gen_config(S,L,m_tmp,N_SW_tmp,p,method)

      if (i > meas_start) then
        j = j+1
        m(j) = m_tmp ! record magnetization
        N_SW(j) = N_SW_tmp ! record clustersize

        if (calc_css) call s_corr(g(j,:),S,L,r_max,n_corr)
        call calc_energy(BE(j),S,L,BJ,h)
      endif

      if (mod(i,plot_interval) == 0) call write_lattice(S,L) ! pipe
    enddo    
    call system_clock(end_time)
    
    call close_lattice_plot()
    call proc_sim_output(L,N_SW,m,start_time,end_time,g,r,BE,&
      calc_css,c_ss_fit,c_ss,nu,Mag,Cv,runtime,Chi)
    deallocate(g,N_SW,m)
  end subroutine

  subroutine gen_config(S,L,m,N_SW,p,method)
    integer, intent(inout)    :: S(:,:)
    real(dp), intent(in)      :: p
    integer, intent(in)       :: L, method
    integer(lng), intent(out) :: m, N_SW ! fix dit nog 
    
    if (method == 1) then
      call SwWa_clust(S,L,m,N_SW,p)
    elseif (method == 2) then
      call Wolff_clust(S,L,m,N_SW,p)
    endif
  end subroutine
  
  pure subroutine calc_energy(BE,S,L,BJ,h)
    integer, intent(in)   :: S(:,:), L
    real(dp), intent(in)  :: h, BJ
    real(dp), intent(out) :: BE
    
    integer :: i, j, k, nn(4,2)

    BE = 0._dp ! initialze energy 

    do i = 1,L
      do j = 1,L
        nn = nn_idx([i,j],L) ! get nearest neighbors of spin i,j
        do k = 1,4
          BE = BE - BJ*S(i,j)*S(nn(k,1),nn(k,2))
        enddo
      enddo
    enddo

    BE = 0.5_dp*BE ! account for double counting of pairs
    BE = BE - h*sum(S) ! add external field
  end subroutine

  pure subroutine s_corr(g,S,L,r_max,n_corr)
    integer, intent(in)   :: S(:,:), L, r_max, n_corr
    real(dp), intent(out) :: g(:)
    
    real(dp)              :: g_tmp(n_corr,r_max)
    integer               :: i, r_0, r_1
                          
    r_0 = (L-n_corr)/2
    r_1 = r_0 + 1
    do i=1,n_corr
      g_tmp(i,:) = S(i+r_0,i+r_0)*S(i+r_0,i+r_1:i+r_0+r_max) 
    enddo

    g = sum(g_tmp,1)/n_corr 
  end subroutine
end module
