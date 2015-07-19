module io
  use constants
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(BJ,h)
    real(dp), intent(out) :: BJ, h
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Beta*J = " 
    read(*,*) BJ
    write(*,'(A)',advance='no') "h = " 
    read(*,*) h
    write(*,'(A)') "Running simulation..."
  end subroutine


  subroutine results_out(BJ,BE,h,runtime, Cv, alpha) 
    real(dp), intent(in) :: BJ, BE, h, Cv, alpha
    integer, intent(in) :: runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "Beta*J :", BJ
      write(12,*) "field :", h 
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Final Energy", BE
      write(12,*) "Specific Heat", Cv
      write(12,*) "alpha: ", alpha
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    call system('cat output.txt')
  end subroutine
end module
