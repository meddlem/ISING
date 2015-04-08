module output_processing 
  use constants
  implicit none
  private 
  public :: calc_M_chi, calc_corr_function, calc_spec_heat 

contains
  
  pure subroutine calc_M_chi(L,N_SW,N_SW_2,m,Q,Mag,err_Mag,chi_s,chi,err_chi,&
      method)
    ! calculates magnetization and susceptibility
    integer, intent(in)      :: L, method
    integer(lng), intent(in) :: m(:), N_SW(:), N_SW_2(:)
    real(dp), intent(out)    :: Q, Mag, err_Mag, chi, chi_s, err_chi

    real(dp) :: N_SW_mean, N
    real(dp), allocatable :: m_r(:)

    allocate(m_r(n_meas))

    ! initialize variables
    N = real(L,dp)**2
    N_SW_mean = sum(real(N_SW,dp))/n_meas
    m_r = real(m,dp)
    
    ! calculate magnetization
    Mag = sum(abs(m_r))/(n_meas*N)

    ! calculate susceptibility
    if (method==1) then
      chi = sum(real(N_SW_2,dp)/N**2)/n_meas
      err_chi = std_err(real(N_SW_2,dp)/N**2)
    elseif (method==2) then
      chi = N_SW_mean/N 
      err_chi = std_err(real(N_SW,dp)/N)
    endif

    Q = chi_s/(N*Mag**2) + 1._dp
    chi_s = (sum(abs(m_r)**2)/n_meas - sum(abs(m_r)/n_meas)**2)/N
    ! also calculate the error for chi_s..
    err_Mag = std_err(m_r/N)
    deallocate(m_r)
  end subroutine

  pure subroutine calc_spec_heat(BE,L,Cv,err_Cv)
    ! calculates specific heat, per particle
    real(dp), intent(in)     :: BE(:)
    integer(lng), intent(in) :: L
    real(dp), intent(out)    :: Cv, err_Cv
    
    real(dp)     :: N, mu_BE, mu_BE_2, mu_BE_4

    N = real(L,dp)**2
    mu_BE = sum(BE)/n_meas
    mu_BE_2 = sum((BE-mu_BE)**2)/n_meas
    mu_BE_4 = sum(block_avg(BE-mu_BE)**4)/n_blocks**2

    Cv = mu_BE_2/N
    err_Cv = sqrt(mu_BE_4)/N
  end subroutine

  pure subroutine calc_corr_function(g,r,c_ss_fit,c_ss,nu,err_nu)
    ! calculate correlation function 
    real(dp), intent(in)  :: g(:,:), r(:)
    real(dp), intent(out) :: c_ss(:), c_ss_fit(:), nu, err_nu
    
    real(dp) :: offset

    c_ss = sum(g,1)/n_meas 
    call lin_fit(nu,err_nu,offset,-log(c_ss),log(r))
    c_ss_fit = exp(-offset)*r**(-nu)
  end subroutine 
    
  pure function std_err(A)
    ! calculates std error of A from blocked data
    real(dp), intent(in) :: A(:)
    real(dp) :: std_err, sigma_blocks_2, Avg(n_blocks)

    Avg = block_avg(A)
  
    sigma_blocks_2 = sum((Avg - sum(Avg)/n_blocks)**2)/(n_blocks-1)
    std_err = sqrt(sigma_blocks_2/n_blocks)
  end function

  pure function block_avg(A)
    ! returns array containing block average of A
    real(dp), intent(in) :: A(:)
    real(dp) :: block_avg(n_blocks)
    integer :: j

    do j = 0,(n_blocks-1)
      block_avg(j+1) = sum(A(n_avg*j+1:n_avg*(j+1)))/n_avg
    enddo
  end function
  
  pure subroutine lin_fit(slope,err_slope,offset,y,x)
    real(dp), intent(out) :: slope, err_slope, offset
    real(dp), intent(in)  :: y(:), x(:)
    real(dp) :: mu_x, mu_y, ss_yy, ss_xx, ss_yx, s
    ! linear regression
    ! see also: http://mathworld.wolfram.com/LeastSquaresFitting.html

    mu_x = sum(x)/size(x)
    mu_y = sum(y)/size(y)
    ss_yy = sum((y - mu_y)**2)
    ss_xx = sum((x - mu_x)**2)
    ss_yx = sum((x - mu_x)*(y - mu_y))
    
    slope = ss_yx/ss_xx
    offset = mu_y - slope*mu_x
    s = sqrt((ss_yy - slope*ss_yx)/(size(x)-2))

    err_slope = s/(sqrt(ss_xx)*6._dp) 
  end subroutine
end module
