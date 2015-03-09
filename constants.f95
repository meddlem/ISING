module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)

  integer, parameter :: L = 100 ! side of lattice 
  integer, parameter :: N = L**2 ! number of spins
  integer, parameter :: n_br = 3 ! number of branches
  integer, parameter :: steps = 20*N ! number of iterations 
  integer, parameter :: s_intvl = 5*N
end module
