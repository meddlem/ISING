module io
  use constants
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(BJ,L,N)
    real(dp), intent(out) :: BJ
    integer, intent(out)  :: L, N
    real(dp)              :: L_tmp
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Beta*J = " 
    read(*,*) BJ
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L_tmp
    L = int(L_tmp)
    N = L**2
    print *, BJ, L
    write(*,'(A)') "Running simulation..."
  end subroutine

  subroutine results_out(BJ,BE,h,runtime,alpha, chi, Cv) 
    real(dp), intent(in) :: BJ, BE, h, alpha, chi, Cv
    integer, intent(in) :: runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "Beta*J :", BJ
      write(12,*) "field :", h 
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Final Energy", BE
      write(12,*) "alpha: ", alpha
      write(12,*) "specific heat", Cv
      write(12,*) "(unsubtracted) susceptibility", chi
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    call system('cat output.txt')
  end subroutine
end module
