module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)

  integer, parameter :: steps = 20000 ! number of iterations
  integer, parameter :: meas_start = 10000 ! start measurement after .. steps 
  integer, parameter :: n_meas = steps - meas_start 
  integer, parameter :: plot_interval = 99 ! plot every .. steps
end module
