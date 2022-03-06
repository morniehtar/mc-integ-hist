module config
    implicit none

! Threads count
    integer, parameter :: thr = 4
! Number of random points in each integral
    integer, parameter :: slw = 1d3 ! 1d9 takes around 8 seconds
! Bins count
    integer, parameter :: bins = 1d3
! Integrals count
    integer, parameter :: ilim = 1d5

! Beginning of the segment
    real(8), parameter :: fl = 0d0
! End of the segment
    real(8), parameter :: cl = 1d0

! Real comparison limiter
    real(8), parameter :: eps = 1d-9

contains

    pure real(8) function func(arg)
        real(8), intent(in) :: arg
        func = dcos(arg)
    end function

    pure real(8) function sol()
        sol = dsin(cl) - dsin(fl)
    end function

end module

program main
    use config
    implicit none
    abstract interface
        real(8) function fct(x)
        real(8), intent(in) :: x
        end function
    end interface

    procedure(fct), pointer :: ptr
    ptr => func
    call random_seed()

    call drawHist(ptr, fl, cl, slw, bins, ilim)

contains

    real(8) function integ(fptr, bt, tp, dots)

        integer, intent(in) :: dots
        real(8), intent(in) :: bt, tp
        procedure(fct), pointer, intent(in) :: fptr

        real(8) :: rnd
        integer :: j

        integ = 0
        !$omp parallel default(none) shared(dots, fptr, bt, tp) private(j, rnd) num_threads(thr)
            !$omp do schedule(static) reduction(+: integ)
            do j = 1, dots
                call random_number(rnd)
                rnd = bt + (tp-bt)*rnd
                integ = integ + fptr(rnd)
            end do
            !$omp end do
        !$omp end parallel

        integ = (tp-bt)*integ/dots

    end function integ

    subroutine drawHist(fptr, bt, tp, dots, nbins, icount)

        integer, intent(in) :: dots, nbins, icount
        real(8), intent(in) :: bt, tp
        procedure(fct), pointer, intent(in) :: fptr

        real(8), dimension(nbins) :: hits
        real(8), dimension(0:nbins+1) :: bnd
        real(8) :: tmp
        integer :: i, j

        hits = 0
        do j = 0, nbins
            bnd(j) = fl + (cl - fl) * j / nbins
        end do

        do i = 1, icount
            tmp = integ(fptr, bt, tp, dots)
            do j = 0, nbins
                if ((bnd(j) - tmp < eps).and.(tmp - bnd(j+1) < eps)) hits(j) = hits(j) + 1
            end do
        end do

        hits = hits / icount

        i=0
        open(unit = 1, file = 'graphics/prob.dat')
        do j = 0, nbins
            tmp = (bnd(j)+bnd(j+1))/2
            write(1, *) tmp, hits(j)
            i = i + hits(j)
        end do
        close(1)

    end subroutine drawHist

end program main

