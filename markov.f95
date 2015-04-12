module markov 
  use constants
  use output_processing 
  use Wolff
  use Swendsen_Wang
  use plotroutines
  implicit none
  private
  public :: markov_chain 

contains
  subroutine markov_chain(S,method,auto,r_max,n_corr,BJ,r,Q,Mag,err_Mag, &
      runtime,calc_css,c_ss,c_ss_fit,eta,err_eta,chi_s,chi_s_err,chi, &
      err_chi,Cv,err_Cv)
    integer, intent(inout)  :: S(:,:)
    real(dp), intent(inout) :: BJ
    integer, intent(in)     :: method, r_max, n_corr
    logical, intent(in)     :: calc_css, auto
    integer, intent(out)    :: runtime
    real(dp), intent(out)   :: c_ss(:), r(:), c_ss_fit(:), Q, Mag, err_Mag, eta, &
      err_eta, chi_s, chi_s_err, chi, err_chi, Cv, err_Cv

    integer(lng), allocatable :: N_SW(:), N_SW_2(:), m(:)
    real(dp), allocatable     :: g(:,:), BE(:)
    integer(lng) :: start_time, m_tmp, N_SW_tmp, N_SW_2_tmp, end_time
    integer      :: i, j, L
    real(dp)     :: p
    
    ! initialize needed variables
    L = size(S,1)
    j = 0
    N_SW_tmp = 0
    N_SW_2_tmp = 0
    m_tmp = 0
    r = real((/(i,i=1,r_max)/),dp)
    p = 1 - exp(-2._dp*BJ)
    
    allocate(g(n_meas,r_max),BE(n_meas),N_SW(n_meas),N_SW_2(n_meas),m(n_meas))
    
    if (.not. auto) call animate_lattice()
    call system_clock(start_time)
    
    ! run markov chain
    do i=1,steps
      call gen_config(S,L,m_tmp,N_SW_tmp,N_SW_2_tmp,p,method)

      if (i > meas_start) then
        j = j+1
        m(j) = m_tmp ! record magnetization
        N_SW(j) = N_SW_tmp ! record clustersize
        N_SW_2(j) = N_SW_2_tmp ! squared cluster size

        if (calc_css) call s_corr(g(j,:),S,L,r_max,n_corr)
        call calc_energy(BE(j),S,L,BJ)
      endif

      if ((mod(i,plot_interval) == 0) .and. (.not. auto)) then
        call write_lattice(S,L) ! write lattice to pipe
      endif
    enddo    
    
    call system_clock(end_time)
    if (.not. auto) call close_lattice_plot()
    
    ! calculate ensemble averages
    call calc_M_chi(L,N_SW,N_SW_2,m,Q,Mag,err_Mag,chi_s,chi_s_err,chi,err_chi,&
      method)
    call calc_spec_heat(BE,L,Cv,err_Cv)
    if (calc_css) call calc_corr_function(g,r,c_ss_fit,c_ss,eta,err_eta)
    
    runtime = (end_time - start_time)/1000

    deallocate(g,BE,N_SW,N_SW_2,m)
  end subroutine

  subroutine gen_config(S,L,m,N_SW,N_SW_2,p,method)
    integer, intent(inout)    :: S(:,:)
    real(dp), intent(in)      :: p
    integer, intent(in)       :: L, method
    integer(lng), intent(out) :: m, N_SW, N_SW_2 
    
    if (method == 1) then
      call SW_clust(S,L,m,N_SW,N_SW_2,p)
    elseif (method == 2) then
      call Wolff_clust(S,L,m,N_SW,p)
      N_SW_2 = N_SW**2
    endif
  end subroutine
  
  pure subroutine calc_energy(BE,S,L,BJ)
    integer, intent(in)   :: S(:,:), L
    real(dp), intent(in)  :: BJ
    real(dp), intent(out) :: BE
    
    integer :: i, j, k, nn(4,2)

    BE = 0._dp ! initialze energy 

    do i = 1,L
      do j = 1,L
        nn = nn_idx([i,j],L) ! get nearest neighbors of spin i,j
        do k = 1,4
          BE = BE - BJ*S(i,j)*S(nn(k,1),nn(k,2))
        enddo
      enddo
    enddo

    BE = 0.5_dp*BE ! account for double counting of pairs
    !BE = BE - h*sum(S) ! add external field
  end subroutine

  pure subroutine s_corr(g,S,L,r_max,n_corr)
    ! calculate spin-spin correlation function
    integer, intent(in)   :: S(:,:), L, r_max, n_corr
    real(dp), intent(out) :: g(:)
    
    real(dp) :: g_tmp(n_corr,r_max)
    integer  :: i, r_0
    ! intialize starting position             
    r_0 = (L-n_corr)/2
    
    ! iterate over n_corr spins on diagonal, get correlation up to r_max
    do i=1,n_corr
      g_tmp(i,:) = S(i+r_0, i+r_0)*S(i+r_0, i+r_0+1:i+r_0+r_max) 
    enddo

    g = sum(g_tmp,1)/n_corr 
  end subroutine
end module
