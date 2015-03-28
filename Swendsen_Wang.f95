module Swendsen_Wang
  use constants
  implicit none
  private 
  public :: SwWa_clust

contains
  subroutine SwWa_clust(S,L,m,N_SW,p)
    integer, intent(inout)    :: S(:,:)
    real(dp), intent(in)      :: p
    integer, intent(in)       :: L
    integer(lng), intent(out) :: m, N_SW 
    
    logical, allocatable :: Bond(:,:,:), Mrkd(:,:)
    integer, allocatable :: N_SW_rec(:), C(:,:)
    integer(lng)  :: k, N
    integer       :: i, j 
    
    ! allocate needed arrays
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
        ! make x,y bonds with probability p when spins are equal
        if (S(i,j) == S(modulo(i,L)+1 ,j)) then
          call random_number(r)
          if (r < p) Bond(1,i,j) = .true.
        endif
        
        if (S(i,j) == (S(i, modulo(j,L)+1))) then
          call random_number(r)
          if (r < p) Bond(2,i,j) = .true.
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

      ! initialize stack holding sites to be scanned
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
          call add_stack(i,j,L,Bond,C,N_stack)
        endif
        ! move to next spin in stack
        k = k+1 
      enddo
    endif
  end subroutine

  pure subroutine add_stack(i,j,L,Bond,C,N_stack)
    ! checks bonds with neighbors, adds to stack if bonded
    integer, intent(inout)  :: C(:,:), N_stack
    logical, intent(in)     :: Bond(:,:,:)
    integer, intent(in)     :: i, j, L
    
    ! if a bond exists add it to the stack
    if (Bond(1,i,j)) then
      N_stack = N_stack+1
      C(:,N_stack) = [modulo(i,L)+1, j] 
    endif
    
    if (Bond(1, modulo(i-2,L)+1, j)) then
      N_stack = N_stack+1
      C(:,N_stack) = [modulo(i-2,L)+1, j] 
    endif
    
    if (Bond(2,i,j)) then 
      N_stack = N_stack+1
      C(:,N_stack) = [i, modulo(j,L)+1] 
    endif

    if (Bond(2, i, modulo(j-2,L)+1)) then 
      N_stack = N_stack+1
      C(:,N_stack) = [i, modulo(j-2,L)+1] 
    endif
  end subroutine
end module
