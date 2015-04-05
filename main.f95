program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  integer :: method
  logical :: calc_css, auto
  
  call get_usr_args(method,calc_css,auto) 
  call init_random_seed()

  if (auto) then 
    call autorun(method,auto)
  else 
    call singlerun(method,calc_css,auto)
  endif
end program
