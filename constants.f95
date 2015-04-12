module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer
  ! steps: number of iterations to perform for fixed lattice size, temp
  ! n_avg: number of steps to take block average over
  ! n_blocks: number of data blocks
  ! meas_start: start of measurment after .. steps
  ! n_meas: total number of measurements
  ! plot_interval: plot every .. iterations
  ! BJ_st : number of temperature values for autorun
  ! BJ_c : central point used for temp interval for autorun
  ! BJ_intv : interval between temp measurements for autorun

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter  :: dp = selected_real_kind(15,307)
  integer, parameter  :: lng = selected_int_kind(8)

  integer, parameter :: steps = 20000 
  integer, parameter :: n_avg = 100
  integer, parameter :: meas_start = 4000
  integer, parameter :: n_meas = steps - meas_start 
  integer, parameter :: n_blocks = n_meas/n_avg 
  integer, parameter :: plot_interval = 99 

  integer, parameter  :: BJ_st = 51
  real(dp), parameter :: BJ_c = 0.44_dp 
  real(dp), parameter :: BJ_intv = 0.002_dp 
end module
