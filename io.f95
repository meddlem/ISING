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
    method = 2 ! use wolff by default
    calc_css = .false.

    ! check command line arguments
    do i=1,iargc()
      call getarg(i,arg)
      if ((trim(arg) == '--SW') .or. (trim(arg) == '-S')) then
        method = 1
      elseif ((trim(arg) == '--Wolff') .or. (trim(arg) == '-W')) then
        method = 2
      elseif (trim(arg) == '-c') then
        calc_css = .true.
      endif
    enddo
  end subroutine

  subroutine user_in(BJ,L,r_max,n_corr)
    real(dp), intent(out) :: BJ
    integer, intent(out)  :: L, r_max, n_corr
    
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Beta*J = " 
    read(*,*) BJ
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L
    write(*,'(A)') "Running simulation..."
    
    ! set variables
    n_corr = L/3 ! number of spins used to calculate correlation
    r_max = L/4 ! distances over which to calc correlation function
  end subroutine

  subroutine results_out(BE,BJ,r,runtime,calc_css,c_ss,c_ss_fit,nu,err_nu, &
      chi_s,chi,err_chi,Mag,err_Mag,Cv,err_Cv) 
    real(dp), intent(in) :: BE(:), BJ, r(:), c_ss(:), c_ss_fit(:), &
      nu, err_nu, chi_s, chi, err_chi, Mag, err_Mag, Cv, err_Cv
    logical, intent(in)  :: calc_css
    integer, intent(in)  :: runtime

    real(dp), allocatable :: t(:)
    integer               :: i
    character(30)         :: output_fmt

    allocate(t(n_meas))
    ! init
    forall(i=0:n_meas-1) t(i+1) = real(i,dp)
    output_fmt = '(A,T25,F8.4,A,F8.4)'

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "Beta*J :", BJ
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,output_fmt) "specific heat", Cv, " ± ", err_Cv 
      write(12,output_fmt) "Magnetization", Mag, " ± ", err_Mag
      write(12,output_fmt) "unsubtracted susceptibility", chi, " ± ", err_chi
      write(12,output_fmt) "susceptibility", chi_s, " ± ", err_chi 
      if (calc_css) write(12,output_fmt) "nu: ", nu, " ± ", err_nu 
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    ! plot results
    call line_plot(real(t,dp),BE,'t','energy','','',1)
    
    if (calc_css) then
      call line_plot(r,c_ss,'r','corr','corr','',3,c_ss_fit,'fit')
    endif
    
    call system('cat output.txt')
    deallocate(t)
  end subroutine
end module
