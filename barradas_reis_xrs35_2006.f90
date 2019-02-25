!-----------------------------------------------------------------------
subroutine pileup(n, spectrum, f, Tw, Tp, time, pile) 
! Calculate pile-up according to Barradas & Reis XRS 35 (2006).
!  Tw     shapping time (s)
!  Tp     time required for the pulse to reach the maximum value
!  time   total pulse time t (s)
!  Nt     total pulse rate
!  Tpulse = f*Tw, f depends on the amplifier (ref f=2)
!-----------------------------------------------------------------------
    implicit none
    
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: n
    real(dp), intent(in) :: f, Tw, Tp, time
    real(dp), dimension(n), intent(in) :: spectrum
    real(dp), dimension(n), intent(out) :: pile
    real(dp) :: Nt, Tpulse, ixj, t1, t2, LG
    integer :: i, j, k, ipj

    Tpulse = f * Tw
    Nt = sum(spectrum) / time
    pile = 0
    do i = 1, n
        do j = i, n
            ipj = i+j
            ixj = i*j
            do k = j, ipj
                t1 = sqrt((ipj-k)*ipj/ixj)
                t2 = sqrt((ipj-k+1.0)*ipj/ixj)
                LG = spectrum(i) * spectrum(j) * (t2 - t1)
                pile(i) = pile(i) - LG
                pile(j) = pile(j) - LG
                if (k < n) then
                    pile(k) = pile(k) + LG
                end if
            end do
        end do
    end do
   pile = pile * exp(-Nt * Tpulse) * Tp / time

end subroutine pileup
