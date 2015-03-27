module io
  use constants
  use plotroutines
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(BJ,L,N,r_max,n_corr)
    real(dp), intent(out) :: BJ
    integer, intent(out)  :: L, N, r_max, n_corr
    real(dp)              :: L_tmp
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Beta*J = " 
    read(*,*) BJ
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L_tmp
    write(*,'(A)') "Running simulation..."
    
    L = int(L_tmp)
    N = L**2
    n_corr = L/3 ! number of spins used to calculate correlation
    r_max = L/4 ! distances over which to calc correlation function
  end subroutine

  subroutine results_out(BE,BJ,t,r,h,runtime,c_ss,c_ss_fit,nu,chi,Mag,Cv) 
    real(dp), intent(in) :: BE(:), BJ, r(:), h, c_ss(:), c_ss_fit(:), &
      nu, chi, Mag, Cv
    integer, intent(in)  :: t(:), runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "Beta*J :", BJ
      write(12,*) "field :", h 
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "nu: ", nu
      write(12,*) "specific heat", Cv
      write(12,*) "Magnetization", Mag
      write(12,*) "(unsubtracted) susceptibility", chi
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    ! plot results
    call line_plot(real(t,dp),BE,'t','energy','','',1)
    call line_plot(r,c_ss,'r','corr','corr','',3,c_ss_fit,'fit')
    
    call system('cat output.txt')
  end subroutine
end module
