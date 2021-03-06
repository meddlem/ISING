program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  integer :: method
  logical :: calc_css, auto
  
  call init_random_seed()
  call get_usr_args(method,calc_css,auto) 

  if (auto) then 
    call autorun(method,auto)
  else 
    call singlerun(method,calc_css,auto)
  endif

contains
  subroutine singlerun(method,calc_css,auto)
    integer, intent(in) :: method
    logical, intent(in) :: calc_css, auto

    real(dp), allocatable :: c_ss(:), r(:), c_ss_fit(:)
    real(dp)              :: BJ, Q, eta, err_eta, chi_s, chi_s_err, chi, err_chi, Mag, &
                             err_Mag, Cv, err_Cv
    integer, allocatable  :: S(:,:)
    integer               :: runtime, L, r_max, n_corr

    call user_in(auto,L,r_max,n_corr,BJ)
    allocate(S(L,L),c_ss(r_max),c_ss_fit(r_max),r(r_max))
    call init_lattice(S,L)

    call markov_chain(S,method,auto,r_max,n_corr,BJ,r,Q,Mag,err_Mag,runtime,&
      calc_css,c_ss,c_ss_fit,eta,err_eta,chi_s,chi_s_err,chi,err_chi,Cv,err_Cv)
    
    call results_out(BJ,L,r,runtime,calc_css,c_ss,c_ss_fit,eta,err_eta,&
      chi_s,chi_s_err,chi,err_chi,Mag,err_Mag,Cv,err_Cv)
    deallocate(S,r,c_ss,c_ss_fit) 
  end subroutine

  subroutine autorun(method,auto)
    integer, intent(in) :: method
    logical, intent(in) :: auto
    
    integer, allocatable  :: S(:,:)
    real(dp), allocatable :: BJ(:), chi_s(:), chi_s_err(:), chi(:), &
      err_chi(:), Mag(:), err_Mag(:), Cv(:), err_Cv(:), Q(:), r(:), c_ss(:), &
      c_ss_fit(:)
    integer  :: i, L, runtime, r_max, n_corr 
    real(dp) :: eta, err_eta
    logical  :: calc_css = .false.

    ! initialize
    call user_in(auto,L,r_max,n_corr)
    allocate(S(L,L),c_ss(r_max),c_ss_fit(r_max),r(r_max),BJ(BJ_st),&
      chi_s(BJ_st),chi_s_err(BJ_st),chi(BJ_st),err_chi(BJ_st),Mag(BJ_st),&
      err_Mag(BJ_st),Cv(BJ_st),err_Cv(BJ_st),Q(BJ_st))
    call init_lattice(S,L)
    call init_BJ_vals(BJ)
    
    do i = 1,BJ_st
      write(*,'(A,F6.3)') 'BJ= ', BJ(i)

      call markov_chain(S,method,auto,r_max,n_corr,BJ(i),r,Q(i),Mag(i), &
        err_Mag(i),runtime,calc_css,c_ss,c_ss_fit,eta,err_eta,chi_s(i), &
        chi_s_err(i),chi(i),err_chi(i),Cv(i),err_Cv(i))
    enddo
    
    call auto_results(L,BJ,Q,Mag,chi_s,Cv)
    deallocate(S,c_ss,c_ss_fit,r,BJ,chi_s,chi_s_err,chi,err_chi,Mag,err_Mag,&
      Cv,err_Cv,Q)
  end subroutine
end program
