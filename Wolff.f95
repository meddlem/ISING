module Wolff
  use constants
  implicit none
  private
  public :: Wolff_clust, nn_idx

contains
  subroutine Wolff_clust(S,L,m,N_SWC,p)
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
end module
