module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer
  ! steps: number of iterations to perform for fixed lattice size, temp
  ! meas_start: start of measurment after .. steps
  ! n_meas: total number of measurements
  ! plot_interval: plot every .. iterations

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)

  integer, parameter :: steps = 20000 
  integer, parameter :: meas_start = 12000
  integer, parameter :: n_meas = steps - meas_start 
  integer, parameter :: plot_interval = 99 
end module
