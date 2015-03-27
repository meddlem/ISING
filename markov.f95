module markov 
  use constants
  use main_routines
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(S,L,r_max,n_corr,BE,BJ,h,t,r,m,Mag,runtime,c_ss,c_ss_fit,&
      nu,chi,Cv)
    integer, intent(inout)  :: S(:,:)
    real(dp), intent(inout) :: BE(:), BJ
    integer, intent(in)     :: L, r_max, n_corr
    real(dp), intent(in)    :: h
    integer, intent(out)    :: t(:), m(:), runtime
    real(dp), intent(out)   :: c_ss(:), r(:), c_ss_fit(:), nu, chi, Cv

    integer, allocatable    :: N_SWC(:)
    real(dp), allocatable   :: g(:,:)
    integer  :: i, j, start_time, m_tmp, N_SWC_tmp, end_time
    real(dp) :: p, Mag
    
    allocate(g(n_meas,r_max),N_SWC(n_meas))
    ! initialize needed variables
    j = 0
    t = (/(i,i=0,n_meas-1)/)
    r = real((/(i,i=1,r_max)/),dp)
    p = 1 - exp(-2._dp*BJ)

    call system_clock(start_time)
    do i=1,steps
      call gen_config(S,L,m_tmp,N_SWC_tmp,p)

      if ((mod(i,meas_step) == 0) .and. (i > meas_start)) then
        j = j+1
        m(j) = m_tmp ! record magnetization
        N_SWC(j) = N_SWC_tmp ! record clustersize

        call s_corr(g(j,:),S,L,r_max,n_corr)
        call calc_energy(BE(j),S,L,BJ,h)
      endif

      if (mod(i,plot_interval) == 0) call write_lattice(S,L) ! pipe
    enddo    
    call system_clock(end_time)
    call sim_run_output(L,N_SWC,m,start_time,end_time,g,r,BE,&
      c_ss_fit,c_ss,nu,Mag,Cv,runtime,Chi)

    deallocate(g,N_SWC)
  end subroutine

  subroutine gen_config(S,L,m,N_SWC,p)
    ! generates wolff cluster
    integer, intent(inout) :: S(:,:)
    real(dp), intent(in)   :: p
    integer, intent(in)    :: L
    integer, intent(out)   :: m, N_SWC 

    integer, allocatable :: C(:,:)
    integer :: i, j, S_init, x(2), nn(4,2)

    allocate(C(L**2,2))
    ! initialize variables 
    i = 1 ! labels spin in cluster
    N_SWC = 1 ! number of spins in cluster
    C = 0 ! init array that holds indices of all spins in cluster
    call random_spin(x,L) ! start cluster by choosing 1 spin

    S_init = S(x(1),x(2)) ! save state of chosen spin
    C(1,:) = x ! add chosen spin to cluster     
    S(x(1),x(2)) = -S_init ! flip initial spin
    
    do while (i<=N_SWC)
      x = C(i,:) ! pick a spin x in the cluster
      nn = nn_idx(x,L) ! get nearest neighbors of spin x
      
      do j = 1,4 ! iterate over neighbors of x
        call try_add(S,C,N_SWC,S_init,nn(j,:),p)
      enddo
      i = i+1 ! move to next spin in cluster
    enddo

    m = sum(S) ! calculate instantaneous magnetization
    deallocate(C)
  end subroutine

  subroutine try_add(S,C,N_SWC,S_init,s_idx,p)
    integer, intent(inout) :: S(:,:), N_SWC, C(:,:)
    integer, intent(in)    :: S_init, s_idx(:)
    real(dp), intent(in)   :: p
    
    real(dp) :: r

    if (S(s_idx(1),s_idx(2)) == S_init) then 
      call random_number(r)

      if (r<p) then ! add spin to cluster with probability p
        N_SWC = N_SWC+1

        C(N_SWC,:) = s_idx 
        S(s_idx(1),s_idx(2)) = -S_init ! flip spin
      endif
    endif
  end subroutine
  
  pure function nn_idx(x, L)
    ! returns indices of nearest neighbors of x_ij, accounting for PBC
    integer, intent(in) :: x(2), L
    
    integer :: nn_idx(4,2)

    nn_idx(1,:) = merge(x + [1,0], [1,x(2)], x(1) /= L)
    nn_idx(2,:) = merge(x + [0,1], [x(1),1], x(2) /= L) 
    nn_idx(3,:) = merge(x - [1,0], [L,x(2)], x(1) /= 1) 
    nn_idx(4,:) = merge(x - [0,1], [x(1),L], x(2) /= 1) 
  end function
  
  subroutine random_spin(x, L)
    ! returns index of randomly picked spin
    integer, intent(in)  :: L
    integer, intent(out) :: x(:)
    
    real(dp) :: u(2)

    call random_number(u)
    u = L*u + 0.5_dp
    x = nint(u) ! index of spin to flip
  end subroutine
  
  pure subroutine calc_energy(BE,S,L,BJ,h)
    integer, intent(in)   :: S(:,:), L
    real(dp), intent(in)  :: h, BJ
    real(dp), intent(out) :: BE
    
    integer :: i, j, k, nn(4,2)

    if (size(S,1) < 2) return !check
    
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
