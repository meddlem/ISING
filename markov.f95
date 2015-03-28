module markov 
  use constants
  use process_data 
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

    call animate_lattice()
    
    call system_clock(start_time)
    do i=1,steps
      call gen_config(S,L,m_tmp,N_SW_tmp,p,method)

      if ((mod(i,meas_step) == 0) .and. (i > meas_start)) then
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
    call sim_proc_output(L,N_SW,m,start_time,end_time,g,r,BE,&
      calc_css,c_ss_fit,c_ss,nu,Mag,Cv,runtime,Chi)
    deallocate(g,N_SW,m)
  end subroutine

  subroutine gen_config(S,L,m,N_SW,p,method)
    integer, intent(inout)    :: S(:,:)
    real(dp), intent(in)      :: p
    integer, intent(in)       :: L, method
    integer(lng), intent(out) :: m, N_SW ! fix dit nog 

    logical, allocatable :: Bond(:,:,:), Mrkd(:,:)
    integer, allocatable :: N_SW_rec(:), C(:,:)
    integer(lng)  :: k, N
    integer       :: i, j 
    
    N = int(L,lng)**2
    allocate(Bond(2,L,L),C(2,8*N),Mrkd(L,L),N_SW_rec(N))

    ! initialize variables 
    Bond = .false. ! init array that holds bonds in x,y dirs
    Mrkd = .false. ! init marked by grow_cluster routine
    N_SW_rec = 0
    k = 0

    ! scan lattice and form bonds
    call freeze_bonds(S,Bond,L,p)
    
    ! build clusters recursively
    do i=1,L
      do j=1,L
        N_SW = 0 ! init cluster size
        call grow_cluster(i,j,S,L,Bond,Mrkd,C,N_SW)
        
        if (N_SW > 0) then
          k = k+1
          N_SW_rec(k) = N_SW ! record cluster size
        endif
      enddo
    enddo

    ! remember to return N_SW_rec average or something
    m = sum(S) ! calculate instantaneous magnetization
    deallocate(Bond,Mrkd,N_SW_rec)
  end subroutine

  subroutine freeze_bonds(S,Bond,L,p)
    ! create bonds between neighboring spins
    logical, intent(inout)  :: Bond(:,:,:)
    integer, intent(in)     :: S(:,:), L 
    real(dp), intent(in)    :: p
    
    integer   :: i, j
    real(dp)  :: r
    
    do i=1,L
      do j=1,L
        if (S(i,j)==S(modulo(i,L)+1,j)) then
          call random_number(r)
          if (r<p) Bond(1,i,j) = .true.
        endif
        
        if (S(i,j)==(S(i,modulo(j,L)+1))) then
          call random_number(r)
          if (r<p) Bond(2,i,j) = .true.
        endif
      enddo
    enddo
  end subroutine

  subroutine grow_cluster(i_init,j_init,S,L,Bond,Mrkd,C,N_SW)
    ! try to form cluster around spin i,j
    integer, intent(inout)  :: S(:,:), N_SW, C(:,:)
    logical, intent(inout)  :: Mrkd(:,:)
    integer, intent(in)     :: i_init, j_init, L
    logical, intent(in)     :: Bond(:,:,:) 

    integer   :: i, j, k, N_stack
    real(dp)  :: r
    logical   :: flsp
    
    if (Mrkd(i_init,j_init)) then
      return
    else
      ! decide if cluster spin will be flipped
      call random_number(r)
      flsp = .false. 
      if (r<0.5_dp) flsp = .true.

      ! init stack holding sites to be scanned
      k = 1
      N_stack = 1
      C(:,1) = [i_init,j_init]
      
      do while (k<=N_stack)
        ! pick spin from stack
        i = C(1,k) 
        j = C(2,k) 

        if (.not. Mrkd(i,j)) then
          Mrkd(i,j) = .true. ! mark site as visited
          N_SW = N_SW + 1 ! increase cluster size
          if (flsp) S(i,j) = -S(i,j) ! flip spin
          
          ! scan bonds with neighbors 
          call grow_stack(i,j,L,Bond,C,N_stack)
        endif
        k = k+1
      enddo
    endif
  end subroutine

  pure subroutine grow_stack(i,j,L,Bond,C,N_stack)
    ! checks bonds with neighbors, adds to stack if bonded
    integer, intent(inout)  :: C(:,:), N_stack
    logical, intent(in)     :: Bond(:,:,:)
    integer, intent(in)     :: i, j, L
    
    if (Bond(1,i,j)) then
      N_stack = N_stack+1
      C(:,N_stack) = [modulo(i,L)+1,j] ! add to stack
    endif
    
    if (Bond(1,modulo(i-2,L)+1,j)) then
      N_stack = N_stack+1
      C(:,N_stack) = [modulo(i-2,L)+1,j] 
    endif
    
    if (Bond(2,i,j)) then 
      N_stack = N_stack+1
      C(:,N_stack) = [i,modulo(j,L)+1] 
    endif

    if (Bond(2,i,modulo(j-2,L)+1)) then 
      N_stack = N_stack+1
      C(:,N_stack) = [i,modulo(j-2,L)+1] 
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
