!-----------------------------------------------------------------------
subroutine pileup(n, spectrum, f, Tw, Tp, time, L, G)
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
    real(dp), dimension(2*n), intent(out) :: L, G
    real(dp) :: Nt, Tpulse, ixj, LG
    integer :: i, j, k, ipj

    Tpulse = f * Tw
    Nt = sum(spectrum) / time
    L = 0
    G = 0

    do i = 1, n
        if (spectrum(i) .le. 0) cycle
        do j = i, n
            if (spectrum(j) .le. 0) cycle
            ipj = i+j
            ixj = i*j
            do k = j, ipj
                LG = spectrum(i) * spectrum(j) * (sqrt((ipj-k+1.0)*ipj/ixj) - sqrt((ipj-k)*ipj/ixj))
                L(i) = L(i) - LG
                L(j) = L(j) - LG
                G(k) = G(k) + LG
            end do
        end do
    end do
   L = L * exp(-Nt * Tpulse) * Tp / time
   G = G * exp(-Nt * Tpulse) * Tp / time

end subroutine pileup
