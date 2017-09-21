!***************************************************************************
! synch_simulate_misalignment.f90
! -------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2015, LI-COR Biosciences
!
! This file is part of EddyPro (TM).
!
! EddyPro (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyPro (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!*******************************************************************************
!
! \brief       Simulate synchronization issues using w and sonic temperature
!              Synchronization issues are:
!              (i)  Jitter
!              (ii) Relative clocks drift
! \author      Gerardo Fratini
! \notes
!
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine sync_simulate_misalignment(Set, nrow, ncol, suffixOutString, PeriodRecords)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: PeriodRecords
    character(*), intent(in) :: suffixOutString
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variable
    integer :: i
    integer, parameter :: nj = 6
    integer, parameter :: nd = 6
    integer :: jits(nj)
    data jits(1:nj) /0, 1, 5, 10, 50, 100/
    integer :: drifts(nd)
    data drifts(1:nd) /0, 10, 30, 60, 90, 180/
    real(kind = dbl) :: SyncStats(5, nj+nd)
    real(kind = dbl) :: SetBak(nrow, ncol)


    SetBak = Set
    if (RPsetup%Sync%default_simulation) then
        print*, 'Performing default synching simulation'
        !> Simulate default jits
        do i = 1, nj
            call simulate_jitter(Set, nrow, ncol, jits(i))
            call BasicStats(Set, nrow, ncol, 1, .false.)
            call get_sync_stats(SyncStats(:, i), size(SyncStats, 2))
            Set = SetBak
        end do
        !> Simulate default drifts
        do i = 1, nd
            call simulate_drift(Set, nrow, ncol, drifts(i))
            call BasicStats(Set, nrow, ncol, 1, .false.)
            call get_sync_stats(SyncStats(:, nj+i), size(SyncStats, 2))
            Set = SetBak
        end do
        call write_sync_stats(usync, SyncStats, size(SyncStats, 1), &
            size(SyncStats, 2), suffixOutString, PeriodRecords)
        return
    end if

    !> ===== SIMULATE JITTER ===================================================
    if (RPsetup%Sync%simulate_jitter) then
        !> Apply jitter to sonic temperature data
        call simulate_jitter(Set, nrow, ncol, RPsetup%Sync%jitter)
    end if

    !> ===== SIMULATE RELATIVE CLOCKS DRIFT ====================================
    if (RPsetup%Sync%simulate_drift) then
        !> Apply clock drift to sonic temperature
        call simulate_drift(Set, nrow, ncol, RPsetup%Sync%drift)
    end if

    !> ===== CALCULATE AND OUTPUT BASIC STATS ==================================
    call BasicStats(Set, nrow, ncol, 1, .false.)
    call WriteOutStats(usync, Stats, suffixOutString, PeriodRecords)

end subroutine sync_simulate_misalignment


!*******************************************************************************
!*******************************************************************************

subroutine simulate_jitter(Set, nrow, ncol, jitter)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: jitter
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    real(kind = dbl) :: jit
    real(kind = dbl) :: off
    double precision, external :: interpolate


    !> Express jitter as a fraction (from 0 to 1) of the sampling interval
    jit = jitter * 1d-3 * Metadata%ac_freq

    do i = 2, nrow-1
        if (Set(i, ts) /= error &
            .and. Set(i-1, ts) /= error &
            .and. Set(i+1, ts) /= error) then

            !> Generate random number in the interval [-jit, jit)
            off = 2d0 * jit * rand(0) - jit

            !> Replace measured Ts with jittered value
            if (off > 0) then
                Set(i, ts) = interpolate(Set(i, ts), Set(i+1, ts), off)
            else
                Set(i, ts) = interpolate(Set(i-1, ts), Set(i, ts), 1+off)
            end if
        end if
    end do
end subroutine simulate_jitter


!*******************************************************************************
!*******************************************************************************

subroutine simulate_drift(Set, nrow, ncol, drift)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: drift
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    integer :: ii
    integer :: next
    real(kind = dbl) :: off
    real(kind = dbl) :: pts(nrow)
    real(kind = dbl) :: dts(nrow)
    real (kind = dbl) :: tts(nrow)
    double precision, external :: interpolate


    !> Build timestamp of precise and drifted time series, in microseconds
    do i = 1, nrow
        pts(i) = 1d6 / nint(Metadata%ac_freq) * (i-1)
        dts(i) = (1d6 + drift) / nint(Metadata%ac_freq) * (i-1)
    end do

    !> Replace measured Ts with drifted value
    next = 2
    do i = 2, nrow
        do ii = next, nrow
            if (dts(i) >= pts(ii-1) .and. dts(i) < pts(ii)) then
                next = ii
                off = (dts(i) - pts(ii-1)) &
                    / 1E6 * nint(Metadata%ac_freq)
                tts(i) = interpolate(Set(ii-1, ts), Set(ii, ts), off)
                exit
            end if
        end do
    end do
    Set(:, ts) = tts
end subroutine simulate_drift


!*******************************************************************************
!*******************************************************************************

subroutine get_sync_stats(SyncStat, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer :: ncol
    real(kind = dbl), intent(out) :: SyncStat(ncol)

    SyncStat(1) = Stats%Mean(w)
    SyncStat(2) = Stats%Mean(ts)
    SyncStat(3) = Stats%StDev(w)**2
    SyncStat(4) = Stats%StDev(ts)**2
    SyncStat(5) = Stats%Cov(w, ts)


end subroutine get_sync_stats


!*******************************************************************************
!*******************************************************************************

subroutine write_sync_stats(unt, SyncStats, nrow, ncol, string, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: N
    real (kind = dbl), intent(in) :: SyncStats(nrow, ncol)
    character(*), intent(in) :: string
    !> local variable
    integer :: i
    integer :: j
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum = ''


    call clearstr(dataline)

    !> add file info
    call AddDatum(dataline, trim(adjustl(string)), separator)
    call WriteDatumInt(N, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> Add stats for all simulations
    do j = 1, ncol
        do i = 1, nrow
            call WriteDatumFloat(SyncStats(i, j), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
    end do

    write(unt, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine write_sync_stats


!*******************************************************************************
!*******************************************************************************

double precision function interpolate(x1, x2, frac) result(xx)
    use m_numeric_kinds
    implicit none
    real(kind = dbl), intent(in) :: x1
    real(kind = dbl), intent(in) :: x2
    real(kind = dbl), intent(in) :: frac

    xx = x1 + (x2 - x1) * frac
end function interpolate


