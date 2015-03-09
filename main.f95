program main
  use constants
  use initialize
  use markov
  use plotroutines
  use io
  implicit none
  
  ! variables:
  ! BE = beta*energy,
  ! dE = change in energy between new and initial config
  ! h = external field
  ! array containing Spins indexed as row,column

  real(dp) :: BJ, BE, BE_branch(n_br), BE_init, dE(n_br), h
  integer, allocatable :: S(:,:,:)
  integer :: i, j, start_time, end_time, runtime
  
  allocate(S(n_br,L+2,L+2))
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S(1,:,:))
  call init_energy(BE,S(1,:,:),BJ,h)
  
  BE_init = BE
  BE_branch = BE
  call sync_branches(S,BE_branch)

  call gnu_lattice_plot(S(1,:,:),1,'initial state')
 
  call system_clock(start_time)
  do i = 1,100*N
    do j = 1,n_br 
      call gen_config(S(j,:,:),dE(j),BJ,h)
    enddo
    
    BE_branch = BE_branch + dE
    if (mod(i,5*N)==0) call sync_branches(S,BE_branch)
  enddo
  call system_clock(end_time)
  
  call sync_branches(S,BE_branch)
  BE = BE_branch(1)

  runtime = (end_time - start_time)/1000
  call results_out(BJ,BE,BE_init,h,runtime)
  call gnu_lattice_plot(S(1,:,:),2,'final state')
end program
