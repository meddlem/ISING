module io
  use constants
  use plotroutines
  implicit none
  private
  public :: get_usr_args, user_in, results_out
contains

  subroutine get_usr_args(method,calc_css)
    integer, intent(out) :: method
    logical, intent(out) :: calc_css 
    
    character(10) :: arg
    integer       :: i

    ! defaults
    method = 1
    calc_css = .false.

    ! check command line arguments
    do i=1,iargc()
      call getarg(i,arg)
      if (trim(arg) == '-S') then
        method = 1
      elseif (trim(arg) == '-W') then
        method = 2
      elseif (trim(arg) == '-M') then
        method = 3
      elseif (trim(arg) == '-c') then
        calc_css = .true.
      endif
    enddo
  end subroutine

  subroutine user_in(method,BJ,L,h,r_max,n_corr)
    integer, intent(in)   :: method
    real(dp), intent(out) :: BJ, h
    integer, intent(out)  :: L, r_max, n_corr
    
    
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Beta*J = " 
    read(*,*) BJ
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L
    if (method == 3) then
      write(*,'(A)',advance='no') "h = " 
      read(*,*) h
    endif
    write(*,'(A)') "Running simulation..."
    
    ! set variables
    n_corr = L/3 ! number of spins used to calculate correlation
    r_max = L/4 ! distances over which to calc correlation function
  end subroutine

  subroutine results_out(BE,BJ,t,r,h,runtime,calc_css,c_ss,c_ss_fit,nu, &
      chi,Mag,Cv) 
    real(dp), intent(in) :: BE(:), BJ, r(:), h, c_ss(:), c_ss_fit(:), &
      nu, chi, Mag, Cv
    logical, intent(in)  :: calc_css
    integer, intent(in)  :: t(:), runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "Beta*J :", BJ
      write(12,*) "field :", h 
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "specific heat", Cv
      write(12,*) "Magnetization", Mag
      write(12,*) "(unsubtracted) susceptibility", chi
      if (calc_css) write(12,*) "nu: ", nu
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    ! plot results
    call line_plot(real(t,dp),BE,'t','energy','','',1)
    
    if (calc_css) then
      call line_plot(r,c_ss,'r','corr','corr','',3,c_ss_fit,'fit')
    endif
    
    call system('cat output.txt')
  end subroutine
end module
