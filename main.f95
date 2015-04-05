program main
  use constants
  use initialize
  use markov
  use io
  use plotroutines
  implicit none
  
  ! variables:
  ! BE: beta*energy,
  ! dE: change in energy between new and initial config
  ! S: array containing Spins indexed as row, column

  real(dp), allocatable :: BJ(:), chi_s(:,:), chi(:,:), err_chi(:,:), &
                           Mag(:,:), err_Mag(:,:), Cv(:,:), err_Cv(:,:)
  integer, allocatable  :: L(:) 
  integer               :: method, L_s, T_s
  logical               :: calc_css
  
  ! initialize variables
  L_s = 7
  T_s = 26
  
  allocate(L(L_s),BJ(T_s),chi_s(L_s,T_s),chi(L_s,T_s),err_chi(L_s,T_s),&
    Mag(L_s,T_s),err_Mag(L_s,T_s),Cv(L_s,T_s),err_Cv(L_s,T_s))

  call get_usr_args(method,calc_css) 
  call init_random_seed()
  call autorun_sim(method,L,BJ,Cv,err_Cv,Mag,err_Mag,chi_s,chi,err_chi)
  call auto_results(L,BJ,Mag,chi_s,Cv)

  deallocate(L,BJ,chi_s,chi,err_chi,Mag,err_Mag,Cv,err_Cv)
end program
