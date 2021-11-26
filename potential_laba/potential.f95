program Potential_laba
    implicit none
    call graph
end program Potential_laba

SUBROUTINE potential(r, o, fi)
    real r, o, I1, I2, I3, fi
    I1 = simpson(0.0, r, r, o)
    I2 = 0.603154*f_2(0.322548, r, o)+0.357419*f_2(1.745761, r, o)+0.0388879*f_2(4.536620, r, o)+0.000539295*f_2(9.395071, r, o)
    I3 = 0.603154*f_3(0.322548, r, o)+0.357419*f_3(1.745761, r, o)+0.0388879*f_3(4.536620, r, o)+0.000539295*f_3(9.395071, r, o)
    !print *, "I1: ", I1
    !print *, "I2: ", I2
    !print *, "I3: ", I3
    fi = (I1 + I2 + I3)*3.0 / 8.0 / 3.1415926535
    return
end SUBROUTINE

real function f_1(x, r, o)
    real x, r, o
    f_1 = 4.0/3.0*3.1415926535*(x*x)*exp(-x)/r - 4.0/15.0*3.1415926535*(x*x*x*x)*exp(-x)/(r*r*r)*(1 - 3*cos(o)*cos(o))
    return
end function

real function f_2(x, r, o)
    real x, r, o
    f_2 = 4.0/3.0*3.1415926535*(x + r)*exp(-r)
    return
end function

real function f_3(x, r, o)
    real x, r, o
    f_3 = -4.0/15.0*3.1415926535*exp(-r)*r*r/(x + r)*(1 - 3*cos(o)*cos(o))
    return
end function

real function simpson(a, b, r, o)
    real xi, xii, h, d, a, b, r, o
    integer i, n
    n = 1001
    h = ((b - a) / (n - 1))
    simpson = f_1(a, r, o) - f_1(b, r, o)
    do i=1,n-1,2
        xi = a + i * h
        xii = a + (i + 1) * h
        simpson = simpson + (4.0 * f_1(xi, r, o)) + (2.0 * f_1(xii, r, o))
    end do
    simpson = ((simpson * h) / 3)
    return
end function

SUBROUTINE graph
    real x0, xn, x, y, xi, xj, yi, c
    integer i, j, k, n
    integer, parameter :: m = 1000
    real mass_1(1:m), mass_2(1:m), mass_3(1:m), mass_4(1:m), mass_5(1:m), mass_6(1:m)
    x0 = 0.01
    xn = 100.0
    open(2, file = "data.txt")
    do k=1,m
        x = x0 + (k - 1)*(xn - x0) / (m - 1)
        call potential (x, 0.0, mass_2(k))
        call potential (x, 3.1415926535 / 6, mass_3(k))
        call potential (x, 3.1415926535 / 4, mass_4(k))
        call potential (x, 3.1415926535 / 3, mass_5(k))
        call potential (x, 3.1415926535 / 2, mass_6(k))
        mass_1(k) = x
        write(2,*) mass_1(k), mass_2(k), mass_3(k), mass_4(k), mass_5(k)!, mass_6(k)
    end do
    close(2)
end SUBROUTINE graph
