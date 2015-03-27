module process_data 
  use constants
  implicit none
  private 
  public :: sim_proc_output

contains
  pure subroutine sim_proc_output(L,N_SWC,m,start_time,end_time,g,r,BE,&
      c_ss_fit,c_ss,nu,Mag,Cv,runtime,Chi)
    ! calculates various physical quantities from simulation
    integer, intent(in)   :: L, m(:), start_time, end_time, N_SWC(:)
    real(dp), intent(in)  :: g(:,:), BE(:), r(:)
    integer, intent(out)  :: runtime
    real(dp), intent(out) :: c_ss_fit(:), c_ss(:), Mag, Cv, nu, Chi

    real(dp)  :: N_SWC_mean, err_nu, offset
    integer   :: N

    ! initialize variables
    N = L**2
    N_SWC_mean = sum(real(N_SWC,dp))/n_meas
    Mag = 0._dp
    
    if (N_SWC_mean > N/2) then
      Mag = sum(real(abs(m),dp))/(n_meas*N)
    endif

    ! calculate susceptibility
    chi = N_SWC_mean/N - Mag**2
    !chi = 1._dp/L**2*sum(real(m,dp)**2)/(n_meas*L**2)
    !chi = sum(N_SWC**2/L**2)/n_meas 

    ! calculate specific heat, per particle
    Cv = sum(BE**2)/n_meas - sum(BE/n_meas)**2
    Cv = Cv/N

    ! calculate correlation function 
    c_ss = sum(g,1)/n_meas 
    call lin_fit(nu,err_nu,offset,-log(c_ss),log(r))
    c_ss_fit = exp(-offset)*r**(-nu)
    
    ! calculate runtime
    runtime = (end_time - start_time)/1000
  end subroutine

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
