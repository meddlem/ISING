module initialize
  use constants
  implicit none
  private
  public :: init_random_seed, init_lattice, init_LT

contains
  subroutine init_lattice(S, L)
    integer, intent(out)  :: S(:,:)
    real(dp), allocatable :: u(:,:)
    integer, intent(in)   :: L
    ! assign initial spins at random, corresponds to T=∞ 

    allocate(u(L,L))
    call random_number(u)
    
    S = -1
    where (u > 0.5_dp) S = 1
    deallocate(u)
  end subroutine 

  pure subroutine init_LT(L,L_s,BJ,T_s)
    real(dp), intent(out) :: BJ(:)
    integer, intent(out)  :: L(:), L_s, T_s

    integer :: i

    L_s = size(L)
    T_s = size(BJ)

    forall(i=1:L_s) L(i) = 2**i
    forall(i=1:T_s) BJ(i) = 0.44_dp + 0.0015_dp*(i-10)
  end subroutine
  
  ! initialize random seed, taken from ICCP github
  subroutine init_random_seed()
    integer, allocatable :: seed(:)
    integer :: i, m, un, istat, dtime(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = m)
    allocate(seed(m))
    open(newunit=un, file="/dev/urandom", access="stream",&
      form="unformatted", action="read", status="old", &
      iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      call system_clock(count)
      if (count /= 0) then
        t = transfer(count, t)
      else
        call date_and_time(values=dtime)
        tms = (dtime(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
          + dtime(2) * 31_8 * 24 * 60 * 60 * 1000 &
          + dtime(3) * 24 * 60 * 60 * 60 * 1000 &
          + dtime(5) * 60 * 60 * 1000 &
          + dtime(6) * 60 * 1000 + dtime(7) * 1000 &
          + dtime(8)
        t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (m >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (m > 3) then
          seed(4:) = s + 37 * (/ (i, i = 0, m - 4) /)
        end if
      else
        seed = s + 37 * (/ (i, i = 0, m - 1 ) /)
      end if
    end if
    call random_seed(put=seed)
  end subroutine
end module 
