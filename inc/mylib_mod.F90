! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
module mylib
  ! -------------------------------------------------------
  ! GENERAL FUNCTIONS AND SUBROUTINES
  ! -------------------------------------------------------
  use general
  implicit none
  ! =======================================================
contains
  ! -------------------------------------------------------
  ! computations
  ! -------------------------------------------------------
  subroutine generateRandomScal(r, n, seed)
    ! generates a random 2D array
    implicit none
    integer, intent(in):: n(1:2), seed
    real(rk), intent(out):: r(0:n(1)+1,0:n(2)+1)
    integer, dimension(33):: iseed ! 33: size of seed for macbookpro

    !  to not have same array on all nodes
    if(seed /= 0) then
       iseed(:) = seed ! regenerate the same array for each run
       call random_seed(put=iseed)
    end if
    call random_number(r)

    return
  end subroutine generateRandomScal
  ! -------------------------------------------------------
  subroutine generateRandomScal1D(r, n, seed)
    ! not in parallel
    implicit none
    integer, intent(in):: n, seed
    real(rk), intent(out):: r(1:n)
    integer, dimension(33):: iseed ! 33: size of seed for macbookpro
    
    if(seed /= 0) then
       iseed(:) = seed
       call random_seed(put=iseed)
    end if
    call random_number(r)
    
    return
  end subroutine generateRandomScal1D
  ! =======================================================
  ! MERGE SORT ALGORITHM
  ! =======================================================
  subroutine merge(a, na, b, nb, c, nc)
    ! merges vectors a and b into c
    implicit none
    integer, intent(in):: na, nb, nc ! requires nc=na+nb
    real(rk), dimension(na), intent(in):: a
    real(rk), dimension(nb), intent(in):: b
    real(rk), dimension(nc), intent(inout):: c
    integer:: i, j, k

    i = 1
    j = 1
    k = 1

    do while(i <= na .and. j <= nb)
       if(a(i) <= b(j)) then
          c(k) = a(i)
          i = i+1
       else
          c(k) = b(j)
          j = j+1
       end if
       k = k+1
    end do

    do while(i <= na)
       c(k) = a(i)
       i = i+1
       k = k+1
    end do

    return
  end subroutine merge
  ! ---------------------------------------------
  recursive subroutine mergeSort(a, n, t)
    ! orders real vector a(1:n); needs a vector t(1:n/2)
    implicit none
    integer, intent(in):: n
    real(rk), dimension(n), intent(inout):: a
    real(rk), dimension((n+1)/2), intent(out):: t
    integer:: na, nb
    real(rk):: val

    if(n < 2) return

    if(n == 2) then
       if(a(1) > a(2)) then
          val = a(1)
          a(1) = a(2)
          a(2) = val
       end if
       return
    end if

    na = (n+1) / 2
    nb = n - na

    call mergeSort(a, na, t)
    call mergeSort(a(na+1), nb, t)

    if(a(na) > a(na+1)) then
       t(1:na) = a(1:na)
       call merge(t, na, a(na+1), nb, a, n)
    end if

    return
  end subroutine mergeSort
  ! ---------------------------------------------
  subroutine sort(vec, n)
    ! sort the real vector vec(1:n)
    implicit none
    integer, intent(in):: n
    real(rk), dimension(n), intent(inout):: vec
    real(rk), dimension((n+1)/2):: buff

    call mergeSort(vec, n, buff)
    return
  end subroutine sort
  ! ---------------------------------------------
  ! =======================================================
  ! MISCELLANEOUS
  ! =======================================================
  subroutine makeName(fileName, stem, number, nameSuff, nameLevel)
    ! create fileName as "stem_{number}" or optionally,
    ! "stem_{nameSuff}{number}" and optionally also, with _nameLevel
    implicit none
    character(len=100), intent(in):: stem
    integer, intent(in):: number
    character(len=100), optional:: nameSuff
    integer, optional:: nameLevel
    character(len=100), intent(out):: fileName
    character(len=4):: cnumber
    character(len=2):: clevel

    call IntToChar4(cnumber, number)

    if(present(nameLevel)) then
       call IntToChar2(clevel, nameLevel)
       if(present(nameSuff)) then
          filename = &
               adjustl(trim(stem))//"_"//adjustl(trim(nameSuff))//cnumber//"_"//clevel
       else
          filename = &
               adjustl(trim(stem))//"_"//cnumber//"_"//clevel
       end if
    else
       if(present(nameSuff)) then
          filename = &
               adjustl(trim(stem))//"_"//adjustl(trim(nameSuff))//cnumber
       else
          filename = &
               adjustl(trim(stem))//"_"//cnumber
       end if
    end if

    return
  end subroutine makeName
  ! -------------------------------------------------------
  subroutine IntToChar2(cnumber, number)
    ! put number into a character of length 2
    implicit none
    integer, intent(in):: number
    character(len=2), intent(out):: cnumber

    write(cnumber, '(i2)') number

    if(number < 10) then
       cnumber = "0"//adjustl(trim(cnumber))
    else if(number < 100) then
       cnumber = adjustl(trim(cnumber))
    else
       write(*,*) &
            "ERROR in number for a file name: cannot handle numbers greater than 99"
       call endRun()
    end if

    return
  end subroutine IntToChar2
  ! =======================================================
  subroutine IntToChar4(cnumber, number)
    ! put number into a character of length 4
    implicit none
    integer, intent(in):: number
    character(len=4), intent(out):: cnumber

    write(cnumber, '(i4)') number

    if(number < 10) then
       cnumber = "000"//adjustl(trim(cnumber))
    else if(number < 100) then
       cnumber = "00"//adjustl(trim(cnumber))
    else if(number < 1000) then
       cnumber = "0"//adjustl(trim(cnumber))
    else if(number < 10000) then
       cnumber = adjustl(trim(cnumber))
    else
       write(*,*) "ERROR in number for a file name: cannot handle numbers greater than 9,999"
       call endRun()
    end if

    return
  end subroutine IntToChar4
  ! -------------------------------------------------------
  subroutine numberOfLines(nlines, fileName)
    ! gives the number of lines in file fileName
    ! doesn't check number of columns in each line
    implicit none
    character(len=100), intent(in):: fileName
    integer, intent(out):: nlines
    integer:: iostatus
    character(len=100):: buff

    nlines = 0
    open(unit=25, file=fileName, status='old')
    do
       read(25, *, iostat=iostatus) buff
       if(iostatus > 0) then
          write(*,'(2a)') "Problem while reading ", adjustl(trim(fileName))
          nlines = 0
          exit
       else if(iostatus < 0) then
          ! reached the end of the file
          exit
       else
          nlines = nlines + 1
       end if
    end do
    close(25)

    return
  end subroutine numberOfLines
  !-------------------------------------------------------
  subroutine existFile(filename, existing)
    ! inquire if filename exists
    implicit none
    character(len=100), intent(in):: filename
    logical, intent(out):: existing

    inquire(file=filename, exist=existing)

    return
  end subroutine existFile
  ! -------------------------------------------------------
  subroutine order(x1, x2, switch)
    ! reverse x1, x2 if x1 > x2
    implicit none
    real(rk), intent(inout):: x1, x2
    logical, intent(out):: switch
    real(rk):: buff

    switch = .false.
    if(x1 > x2) then
       switch = .true.
       buff = x1
       x1 = x2
       x2 = buff
    end if
    
    return
  end subroutine order
  ! -------------------------------------------------------
  ! =======================================================
  ! Functions
  ! =======================================================
  function myfloor(a)
    ! for rounding problems
    integer:: myfloor
    real(rk), intent(in):: a
    real(rk):: locsmall = 1.0E-12
    
    myfloor = floor(a + locsmall)
    
  end function myfloor
  ! =======================================================
  ! Routines to solve linear systems
  ! =======================================================
  subroutine det3x3(det, a11, a12, a13, a21, a22, a23, a31, a32, a33)
    implicit none
    real(rk), intent(in):: a11, a12, a13, a21, a22, a23, a31, a32, a33
    real(rk), intent(out):: det

    det = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a32 * a21 - &
         a11 * a23 * a32 - a22 * a13 * a31 - a33 * a12 * a21

    return
  end subroutine det3x3
  ! =============================================
  subroutine det4x4(det, a11, a12, a13, a14, a21, a22, a23, a24, &
       a31, a32, a33, a34, a41, a42, a43, a44)
    implicit none
    real(rk), intent(in):: a11, a12, a13, a14
    real(rk), intent(in):: a21, a22, a23, a24
    real(rk), intent(in):: a31, a32, a33, a34
    real(rk), intent(in):: a41, a42, a43, a44
    real(rk), intent(out):: det
    real(rk):: det11, det12, det13, det14

    call det3x3(det11, a22, a23, a24, a32, a33, a34, a42, a43, a44)
    call det3x3(det12, a21, a23, a24, a31, a33, a34, a41, a43, a44)
    call det3x3(det13, a21, a22, a24, a31, a32, a34, a41, a42, a44)
    call det3x3(det14, a21, a22, a23, a31, a32, a33, a41, a42, a43)

    det = a11 * det11 - a12 * det12 + a13 * det13 - a14 * det14

    return
  end subroutine det4x4
  ! =============================================
  subroutine solve2x2(sol, a11, a12, a21, a22, y)
    implicit none
    real(rk), dimension(1:2), intent(out):: sol
    real(rk), dimension(1:2), intent(in):: y
    real(rk), intent(in):: a11, a12, a21, a22
    real(rk):: det
    
    det = a11 * a22 - a12 * a21
    det = 1.0 / det

    sol(1) = det * (a22 * y(1) - a12 * y(2))
    sol(2) = det * (-a21 * y(1) + a11 * y(2))

    return
  end subroutine solve2x2
  ! =============================================
  subroutine solve3x3(sol, a11, a12, a13, a21, a22, a23, a31, a32, a33, y)
    ! -----------------------------
    ! Cramer's Method to solve
    ! [a11 a12 a13] [sol1]   [y1]
    ! [a21 a22 a23] [sol2] = [y2]
    ! [a31 a32 a33] [sol3]   [y3]
    ! -----------------------------
    implicit none
    real(rk), dimension(1:3), intent(out):: sol
    real(rk), dimension(1:3), intent(in):: y
    real(rk), intent(in):: a11, a12, a13, a21, a22, a23, a31, a32, a33
    real(rk):: det, det1, det2, det3

    call det3x3(det, a11, a12, a13, a21, a22, a23, a31, a32, a33)

    call det3x3(det1, &
         y(1), a12, a13, &
         y(2), a22, a23, &
         y(3), a32, a33)
    call det3x3(det2, &
         a11, y(1), a13, &
         a21, y(2), a23, &
         a31, y(3), a33)
    call det3x3(det3, &
         a11, a12, y(1), &
         a21, a22, y(2), &
         a31, a32, y(3))

    sol(1) = det1 / det
    sol(2) = det2 / det
    sol(3) = det3 / det

    return
  end subroutine solve3x3
  ! =============================================
  subroutine solve4x4(sol, a11, a12, a13, a14, a21, a22, a23, a24, &
       a31, a32, a33, a34, a41, a42, a43, a44, y)
    ! -----------------------------
    ! Cramer's Method to solve
    ! [a11 a12 a13 a14] [sol1]   [y1]
    ! [a21 a22 a23 a24] [sol2] = [y2]
    ! [a31 a32 a33 a34] [sol3]   [y3]
    ! [a41 a42 a43 a44] [sol4]   [y4]
    ! -----------------------------
    implicit none
    real(rk), dimension(1:4), intent(out):: sol
    real(rk), dimension(1:4), intent(in):: y
    real(rk), intent(in):: a11, a12, a13, a14
    real(rk), intent(in):: a21, a22, a23, a24
    real(rk), intent(in):: a31, a32, a33, a34
    real(rk), intent(in):: a41, a42, a43, a44
    real(rk):: det, det1, det2, det3, det4

    call det4x4(det, a11, a12, a13, a14, a21, a22, a23, a24, &
         a31, a32, a33, a34, a41, a42, a43, a44)

    call det4x4(det1, &
         y(1), a12, a13, a14, &
         y(2), a22, a23, a24, &
         y(3), a32, a33, a34, &
         y(4), a42, a43, a44)
    call det4x4(det2, &
         a11, y(1), a13, a14, &
         a21, y(2), a23, a24, &
         a31, y(3), a33, a34, &
         a41, y(4), a43, a44)
    call det4x4(det3, &
         a11, a12, y(1), a14, &
         a21, a22, y(2), a24, &
         a31, a32, y(3), a34, &
         a41, a42, y(4), a44)
    call det4x4(det4, &
         a11, a12, a13, y(1), &
         a21, a22, a23, y(2), &
         a31, a32, a33, y(3), &
         a41, a42, a43, y(4))

    sol(1) = det1 / det
    sol(2) = det2 / det
    sol(3) = det3 / det
    sol(4) = det4 / det

    return
  end subroutine solve4x4
  ! =============================================
  subroutine tridiag(sol, aa, bb, cc, y, nval)
    ! LU decomposition to solve M SOL = Y
    ! (M: tridiagonal matrix n*n with tridiag = aa, bb, cc)
    implicit none
    integer, intent(in):: nval
    real(rk), dimension(1:nval), intent(in):: aa, bb, cc, y
    real(rk), dimension(1:nval), intent(out):: sol
    real(rk), dimension(1:nval):: beta, u, v
    integer:: i

    beta(1) = bb(1)
    u(1) = cc(1) / beta(1)
    v(1) = y(1) / beta(1)
    do i=2, nval
       beta(i) = bb(i) - aa(i) * u(i-1)
       u(i) = cc(i) / beta(i)
       v(i) = (y(i) - aa(i) * v(i-1)) / beta(i)
    end do

    sol(nval) = v(nval)
    do i=nval-1, 1, -1
       sol(i) = v(i) - u(i) * sol(i+1)
    end do

    return
  end subroutine tridiag
  ! =============================================
  subroutine nearlyTridiag(sol, aa, bb, cc, y, nval)
    ! solve tridiagonal system with non-zero values in the corners
    ! LU decomposition
    implicit none
    integer, intent(in):: nval
    real(rk), dimension(1:nval), intent(in):: aa, bb, cc, y
    real(rk), dimension(1:nval), intent(out):: sol
    real(rk), dimension(1:nval):: beta, u, v, nu, ell
    real(rk):: som
    integer:: i

    !  diagonal of L and upper diag of U - - - - -
    ell(1) = bb(1)
    u(1) = cc(1) / ell(1)
    do i=2, nval-2
       ell(i) = bb(i) - aa(i) * u(i-1)
       u(i) = cc(i) / ell(i)
    end do
    ell(nval-1) = bb(nval-1) - aa(nval-1) * u(nval-2)
    
    ! last column of U - - - - - - - - - - - - -
    nu(1) = aa(1) / ell(1)
    do i=2, nval-2
       nu(i) = -aa(i) * nu(i-1) / ell(i)
    end do ! i
    nu(nval-1) = (cc(nval-1) - aa(nval-1) * nu(nval-2)) / ell(nval-1)

    ! last row of L - - - - - - - - - - - - - - 
    beta(1) = cc(nval)
    do i=2, nval-2
       beta(i) = -beta(i-1) * u(i-1)
    end do ! i
    beta(nval-1) = aa(nval) - beta(nval-2) * u(nval-2)

    som = 0.0
    do i=1, nval-1
       som = som + beta(i) * nu(i)
    end do
    beta(nval) = bb(nval) - som

    ! solving L v = y - - - - - - - - - - - - - -
    v(1) = y(1) / ell(1)
    som = beta(1) * v(1)
    do i=2, nval-1
       v(i) = (y(i) - aa(i) * v(i-1)) / ell(i)
       som = som + beta(i) * v(i)
    end do ! i
    v(nval) = (y(nval) - som) / beta(nval)

    ! solving U sol = v
    sol(nval) = v(nval)
    sol(nval-1) = v(nval-1) - nu(nval-1) * sol(nval)
    do i=nval-2, 1, -1
       sol(i) = v(i) - nu(i) * sol(nval) - u(i) * sol(i+1)
    end do ! i

    return
  end subroutine nearlyTridiag
  ! =============================================
  subroutine pentadiag(sol, ee, cc, dd, aa, bb, y, nval)
    ! algorithm PTRANS-II de Askar and Karawia, 2015, Math.
    ! Problems in Engineering, doi:10.1155/2015/232456
    ! (with correction of their wrong phi(4) to phi(3) in sig(2))
    ! 
    ! order of letters is changed compared to notation elsewhere here,
    ! to match the paper by Askar and Karawia.
    ! (ee, cc, dd, aa, bb)
    implicit none
    integer, intent(in):: nval
    real(rk), dimension(1:nval), intent(in):: aa, bb, cc, dd, ee, y
    real(rk), dimension(1:nval), intent(out):: sol
    real(rk), dimension(1:nval):: psi, sig, phi, ome, rho
    integer:: i

    ! step 3 of Askar and Karawia - - - - - -
    psi(nval) = dd(nval)
    sig(nval) = cc(nval) / psi(nval)
    phi(nval) = ee(nval) / psi(nval)
    ome(nval) = y(nval) / psi(nval)

    ! step 4 - - - - - -
    rho(nval-1) = aa(nval-1)
    psi(nval-1) = dd(nval-1) - sig(nval) * rho(nval-1)
    sig(nval-1) = (cc(nval-1) - phi(nval) * rho(nval-1)) / psi(nval-1)
    phi(nval-1) = ee(nval-1) / psi(nval-1)
    ome(nval-1) = (y(nval-1) - ome(nval) * rho(nval-1)) / psi(nval-1)

    ! step 5 - - - - - -
    do i=nval-2, 3, -1
       rho(i) = aa(i) - sig(i+2) * bb(i)
       psi(i) = dd(i) - phi(i+2) * bb(i) - sig(i+1) * rho(i)
       sig(i) = (cc(i) - phi(i+1) * rho(i)) / psi(i)
       phi(i) = ee(i) / psi(i)
       ome(i) = (y(i) - ome(i+2) * bb(i) - ome(i+1) * rho(i)) / psi(i)
    end do

    rho(2) = aa(2) - sig(4) * bb(2)
    psi(2) = dd(2) - phi(4) * bb(2) - sig(3) * rho(2)
    sig(2) = (cc(2) - phi(3) * rho(2)) / psi(2)

    rho(1) = aa(1) - sig(3) * bb(1)
    psi(1) = dd(1) - phi(3) * bb(1) - sig(2) * rho(1)
    
    ome(2) = (y(2) - ome(4) * bb(2) - ome(3) * rho(2)) / psi(2)
    ome(1) = (y(1) - ome(3) * bb(1) - ome(2) * rho(1)) / psi(1)

    ! step 6 - - - - -
    sol(1) = ome(1)
    sol(2) = ome(2) - sig(2) * sol(1)
    do i=3, nval
       sol(i) = ome(i) - sig(i) * sol(i-1) - phi(i) * sol(i-2)
    end do
    
    return
  end subroutine pentadiag
  ! =============================================
  subroutine nearlyPentadiag(sol, ee, aa, dd, cc, ff, y, nval)
    ! algorithm NPENTA of Neossi Nguetchue and Abelamn, 2008,
    ! Appl. Math and Computation, Vol.203, doi: 10.1016/j.amc.2008.05.012
    !
    ! letters changed to ee, aa, dd, cc, ff to match their writing.
    ! Changes to their notation: ee(1)=p_1, aa(1)=q_1, ee(2)=r_2,
    ! ff(n-1)=alpha_{n-1}, cc(n)=beta_n, ff(n)=gamma_n
    implicit none
    integer, intent(in):: nval
    real(rk), dimension(1:nval), intent(in):: ee, aa, dd, cc, ff, y
    real(rk), dimension(1:nval), intent(out):: sol
    real(rk), dimension(1:nval):: beta, gamm, h, k, l, v, w, z
    integer:: i
    real(rk):: som

    ! ----------------------------------
    ! step 1 - - - - - - - - - - - - - -
    gamm(1) = dd(1)
    beta(2) = aa(2) / gamm(1)
    h(1) = cc(1)
    k(1) = ff(nval-1) / gamm(1) ! ff(n-1)=alpha_{n-1} in their notation
    w(1) = aa(1)  ! q_1 in their notation
    v(1) = ee(1)  ! p_1 in their notation
    l(1) = cc(nval) / gamm(1)  ! cc(n)=beta_n in their notation

    gamm(2) = dd(2) - beta(2) * h(1)
    k(2) = -k(1) * h(1) / gamm(2)
    w(2) = ee(2) - beta(2) * w(1)  ! ee(2)=r_2
    v(2) = -beta(2) * v(1)
    l(2) = (ff(nval) - l(1) * h(1)) / gamm(2) ! ff(n)=gamma_n
    h(2) = cc(2) - beta(2) * ff(1)

    ! ----------------------------------
    ! step 2 - - - - - - - - - - - - - -
    do i=3, nval-3
       beta(i) = (aa(i) - ee(i) / gamm(i-2) * h(i-2)) / gamm(i-1)
       h(i) = cc(i) - beta(i) * ff(i-1)
       gamm(i) = dd(i) - ee(i) / gamm(i-2) * ff(i-2) - beta(i) * h(i-1)
    end do

    ! ----------------------------------
    ! step 3 - - - - - - - - - - - - - -
    beta(nval-2) = (aa(nval-2) - ee(nval-2) / gamm(nval-4) * h(nval-4)) / gamm(nval-3)
    gamm(nval-2) = dd(nval-2) - ee(nval-2) / gamm(nval-4) * ff(nval-4) - beta(nval-2) * h(nval-3)

    ! ----------------------------------
    ! step 4 - - - - - - - - - - - - - -
    do i=3, nval-4
       k(i) = -(k(i-2) * ff(i-2) + k(i-1) * h(i-1)) / gamm(i)
       v(i) = -ee(i) / gamm(i-2) * v(i-2) - beta(i) * v(i-1)
    end do ! i

    ! ----------------------------------
    ! step 5 - - - - - - - - - - - - - -
    k(nval-3) = (ee(nval-1) - k(nval-5) * ff(nval-5) - k(nval-4) * h(nval-4)) / gamm(nval-3)
    k(nval-2) = (aa(nval-1) - k(nval-4) * ff(nval-4) - k(nval-3) * h(nval-3)) / gamm(nval-2)
    v(nval-3) = ff(nval-3) - ee(nval-3) / gamm(nval-5) * v(nval-5) - beta(nval-3) * v(nval-4)
    v(nval-2) = cc(nval-2) - ee(nval-2) / gamm(nval-4) * v(nval-4) - beta(nval-2) * v(nval-3)
    som = 0.0
    do i=1, nval-2
       som = som + k(i) * v(i)
    end do
    gamm(nval-1) = dd(nval-1) - som

    ! ----------------------------------
    ! step 6 - - - - - - - - - - - - - -
    do i=3, nval-3
       w(i) = -ee(i) / gamm(i-2) * w(i-2) - beta(i) * w(i-1)
       l(i) = -(l(i-2) * ff(i-2) + l(i-1) * h(i-1)) / gamm(i)
    end do

    ! ----------------------------------
    ! step 7 - - - - - - - - - - - - - -
    w(nval-2) = ff(nval-2) - ee(nval-2) / gamm(nval-4) * w(nval-4) - beta(nval-2)  * w(nval-3)
    som = 0.0
    do i=1, nval-2
       som = som + k(i) * w(i)
    end do
    w(nval-1) = cc(nval-1) - som

    l(nval-2) = (ee(nval) - l(nval-4) * ff(nval-4) - l(nval-3) * h(nval-3)) / gamm(nval-2)
    som = 0.0
    do i=1, nval-2
       som = som + l(i) * v(i)
    end do
    l(nval-1) = (aa(nval) - som) / gamm(nval-1)

    som = 0.0
    do i=1, nval-1
       som = som + l(i) * w(i)
    end do
    gamm(nval) = dd(nval) - som

    ! ----------------------------------
    ! step 8 - - - - - - - - - - - - - -
    z(1) = y(1)
    z(2) = y(2) - beta(2) * z(1)

    ! step 9 - - - - - - - - - - - - - -
    do i=3, nval-2
       z(i) = y(i) - beta(i) * z(i-1) - ee(i) / gamm(i-2) * z(i-2)
    end do

    ! step 10 - - - - - - - - - - - - - -
    som = 0.0
    do i=1, nval-2
       som = som + k(i) * z(i)
    end do
    z(nval-1) = y(nval-1) - som

    som = 0.0
    do i=1, nval-1
       som = som + l(i) * z(i)
    end do
    z(nval) = y(nval) - som

    ! step 11 - - - - - - - - - - - - - -
    sol(nval) = z(nval) / gamm(nval)
    sol(nval-1) = (z(nval-1) - w(nval-1) * sol(nval)) / gamm(nval-1)
    sol(nval-2) = (z(nval-2) - v(nval-2) * sol(nval-1) - w(nval-2) * sol(nval)) / gamm(nval-2)
    sol(nval-3) = (z(nval-3) - h(nval-3) * sol(nval-2) - v(nval-3) * sol(nval-1) - w(nval-3) * sol(nval)) / gamm(nval-3)
    do i=nval-4, 1, -1
       sol(i) = (z(i) - h(i) * sol(i+1) - ff(i) * sol(i+2) - v(i) * sol(nval-1) - w(i) * sol(nval)) / gamm(i)
    end do

    return
  end subroutine nearlyPentadiag
  ! =============================================
  !-------------------------------------------------------
end module mylib
! =======================================================
