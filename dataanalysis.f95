module dataanalysis
    use constants
    implicit none
    private
    public :: calc_specificheat
contains

    subroutine calc_specificheat(BE, startindex, Cv)
        real(dp), intent(in)    :: BE(:)
        integer, intent(in)     :: startindex
        integer                 :: stopindex
        real(dp), intent(out)   :: Cv
        stopindex = size(BE)
        !Specificheat : Cv = Kb*B**2(<BE**2> - (<BE>)**2), Boltzmann constant set to one.
        Cv = (dot_product(BE(startindex:stopindex),BE(startindex:stopindex))+sum(BE(startindex:stopindex)))/(stopindex-startindex)
    end subroutine

    subroutine calc_magneticsusceptability


        !susceptability : X = Kb*B(<m**2> - (<m>)**2), Boltzmann constant set to one.

    end subroutine

end module