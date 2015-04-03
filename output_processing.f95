module output_processing 
  use constants
  implicit none
  private 
  public :: calc_chi, calc_corr_function, calc_spec_heat 

contains
  pure subroutine calc_chi(L,N_SW,N_SW_2,m,Mag,err_Mag,chi,err_chi,method)
    ! calculates magnetization and susceptibility
    integer, intent(in)      :: L, method
    integer(lng), intent(in) :: m(:), N_SW(:), N_SW_2(:)
    real(dp), intent(out) :: Mag, err_Mag, chi, err_chi

    real(dp)     :: N_SW_mean
    integer(lng) :: N

    ! initialize variables
    N = L**2
    N_SW_mean = sum(real(N_SW,dp))/n_meas
    
    ! calculate magnetization for low temp
    if (N_SW_mean > N/2) then
      Mag = sum(real(abs(m),dp))/(n_meas*N)
    else 
      Mag = 0._dp
    endif

    ! calculate susceptibility
    if (method==1) then
      chi = sum(real(N_SW_2,dp)/N**2)/n_meas
      err_chi = std_err(real(N_SW_2,dp)/N**2)
    elseif (method==2) then
      !chi = 1._dp/L**2*sum(real(m,dp)**2)/(n_meas*L**2)
      chi = N_SW_mean/N !sum(real(N_SW_2,dp)/N)/n_meas !N_SW_mean/N !- Mag**2 
      err_chi = std_err(real(N_SW,dp)/N)
    endif

    err_Mag = std_err(real(m,dp)/N)
  end subroutine

  pure subroutine calc_spec_heat(BE,L,Cv,err_Cv)
    ! calculates specific heat, per particle
    real(dp), intent(in)     :: BE(:)
    integer(lng), intent(in) :: L
    real(dp), intent(out)    :: Cv, err_Cv
    
    real(dp)     :: mu_BE, mu_BE_2, mu_BE_4
    integer(lng) :: N

    N = L**2
    mu_BE = sum(BE)/n_meas
    mu_BE_2 = sum((BE-mu_BE)**2)/n_meas
    mu_BE_4 = sum(block_avg(BE-mu_BE)**4)/n_blocks**2

    Cv = mu_BE_2/N
    err_Cv = sqrt(mu_BE_4)/N
  end subroutine

  pure subroutine calc_corr_function(g,r,c_ss_fit,c_ss,nu)
    ! calculate correlation function 
    real(dp), intent(in)  :: g(:,:), r(:)
    real(dp), intent(out) :: c_ss(:), c_ss_fit(:), nu
    
    real(dp) :: offset, err_nu

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

    err_slope = s/(sqrt(ss_xx)*6._dp) ! error in slope calc
  end subroutine
end module
