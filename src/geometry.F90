! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
logical function powerTwoOk(i)
  ! checks that i is a power of two
  implicit none
  integer, intent(in):: i
  real a
  a = log(real(i)) / log(2.0)
  if(2.0**a /= i) then
     powerTwoOk = .false.
  else
     powerTwoOk = .true.
  endif
end function powerTwoOk
! =======================================================
subroutine setupGeometry()
  ! computes grid spacing and number of cells on local grid
  use general
  implicit none
  logical:: okX, okZ, powerTwoOk
  integer:: n360

  curvedGeometry = .true.
  if(coordinates == "cylindrical") then
     ed = 1
     letterGeom = "C"
  else if(coordinates == "spherical") then
     ed = 2
     letterGeom = "S"
  else if(coordinates == "Cartesian" .or. coordinates == "cartesian") then
     curvedGeometry = .false.
  else
     write(*,'(a)') "ERROR: coordinates not recognized"
     call endRun()
  end if

  if(curvedGeometry) then ! ---------------
     ! checking fcurv and angledeg are okay
     if(fcurv <= 0.0 .or. fcurv >= 1.0)  then
        write(*,'(a)') "ERROR: fcurv has to be so that 0 < fcurv < 1"
call endRun()
     end if

     if(anglemax <= 0.0 .or. anglemax > 360.0) then
        write(*,'(a)') "WARNING: anglemax is not in [0;360] => modified"
        n360 = int(abs(anglemax) / 360.0)
        if(anglemax >= 0.0) then
           anglemax = anglemax - n360 * 360.0
        else
           anglemax = anglemax + n360 * 360.0
        end if
        write(*,'(a,f8.3)')"         new value = ", anglemax
     end if
  end if ! curvedGeometry ----------------
  
  ! checking nx >= nz
  if(nx < nz) then
     write(*,'(a)') "ERROR: for MG, we must have nphi >= nz"
     call endRun()
  end if

  okX = powerTwoOk(nx)
  okZ = powerTwoOk(nz)

  if(.not.okX .or. .not.okZ) then
     write(*,'(a)') "ERROR: nx and nz must be power of two"
     call endRun()
  end if

  n0(1) = nx
  n0(2) = nz

  if(curvedGeometry) then
     dx = anglemax / 180.0 * pi / n0(1)

     rbot = fcurv / (1.0 - fcurv)
     rtop = 1.0 + rbot
  else
     dx = lengthX / n0(1)
  end if
  invdx = 1.0 / dx
  
  dz = 1.0 / n0(2)
  invdz = 1.0 / dz
  
  dd(1) = dx
  dd(2) = dz
  invdd(1) = invdx
  invdd(2) = invdz

  ! variables to simplify writing of coefficients
  ! for biharmonic equation in Cartesian geometry
  r2_cart = dd(2)*dd(2) * invdd(1)*invdd(1)
  d2_cart = dd(1)*dd(1) * dd(2)*dd(2)
  invr2_cart = 1.0 / r2_cart

  ! for maps
  j25 = int(0.25 * n0(2))
  j50 = int(0.50 * n0(2))
  j75 = int(0.75 * n0(2))
  ! for NCells
  j375 = int(0.375 * n0(2))

  return
end subroutine setupGeometry
! =======================================================
subroutine getPos(position, i, idir, level, id)
  ! gives the position in idir-direction of cell located at i on the local grid.
  ! level is for the level of the multigrid considered (0 if finest grid)
  use general
  use grids
  implicit none
  integer, intent(in):: i, level, idir, id
  real(rk), intent(out):: position
  real(rk):: d(1:2)

  if(idir == 1) then
     position = grid(level)%pos(i, i0, idir)
  else if(idir == 2) then
     position = grid(level)%pos(i0, i, idir)
  end if
  
  if( id == 3 &
       .or. (id == 2 .and. idir == 1) &
       .or. (id == 1 .and. idir == 2) ) then
     ! add half d at cells' centers or middle of edges:
     call getD(d, level)  ! dx or dz
     position = position + half * d(idir)
  end if

  return
end subroutine getPos
! ==================================================
subroutine getPosX(posx, i, level, id)
  ! gives the position in horizontal direction of cell located at (i,:)
  ! on the local grid.
  use general
  use grids
  implicit none
  integer, intent(in):: i, level, id
  real(rk), intent(out):: posx

  call getPos(posx, i, 1, level, id)

  return
end subroutine getPosX
! ===================================================
subroutine getPosZ(posz, j, level, id)
  ! gives the position in vertical direction of cell located at (:,j)
  ! on the local grid.
  use general
  use grids
  implicit none
  integer, intent(in):: j, level, id
  real(rk), intent(out):: posz

  call getPos(posz, j, 2, level, id)

  return
end subroutine getPosZ
! ===================================================
subroutine getCoordX(i, posx, level, id)
  ! gives the coordinate i corresponding to the position posH
  use general
  use grids
  real(rk), intent(in):: posx
  integer, intent(in):: level, id
  integer, intent(out):: i
  real(rk):: invd(1:2)

  call getInvD(invd, level) ! invd(1) = 1/dx
  
  if( id == 1 .or. id == 4 ) then
     i = floor(posx * invd(1)) + i0
  else ! id=2 or 3
     i = floor(posx * invd(1) - half) + i0
  end if
  
  return
end subroutine getCoordX
! =======================================================
subroutine getCoordZ(j, posZ, level, id)
  ! gives the coordinate iz corresponding to the position posz
  use general
  use grids
  implicit none
  real(rk), intent(in):: posZ
  integer, intent(in):: level, id
  integer, intent(out):: j
  real(rk):: invd(1:2)

  call getInvD(invd, level) ! invd(2) = 1/dr
  
  if(id == 2 .or. id == 4) then
     j = floor(posZ * invd(2)) + i0  
  else ! id=1 or 3
     j = floor(posZ * invd(2) - half) + i0
  end if
  
  return  
end subroutine getCoordZ
! =======================================================
subroutine getCoords(i, j, posX, posZ, level, id)
  ! gives the horizontal and vertical coord. for the position posH-posV
  ! id is for the left bottom corner of the considered cell.
  use general
  use grids
  implicit none
  real(rk), intent(in):: posX, posZ
  integer, intent(in):: level, id
  integer, intent(out):: i, j
  
  call getCoordZ(j, posZ, level, id)
  call getCoordX(i, posX, level, id)
  
  return
end subroutine getCoords
! =======================================================
subroutine geometricalCoefficients()
  ! compute all the geometrical coefficients needed
  ! for boundary conditions and for the solver
  ! in stream function
  use general
  use grids
  implicit none
  integer:: level, iez, j
  real(rk):: a

  if(.not. curvedGeometry) then
     do level=0, nLevel
        grid(level)%topdirichlet(:) = one
        grid(level)%topfreeslip = one
        grid(level)%botdirichlet(:) = one
        grid(level)%botfreeslip = one
     end do ! level
     return  ! 
  end if

  ! ---------------------------------------------
  ! COEFFICIENTS FOR BOUNDARY CONDITIONS FOR SF
  ! FOR CURVED GEOMETRIES
  ! ---------------------------------------------
  do level=0, nLevel
     iez = iend(level, 2)
     ! TOP - - - - -
     a = grid(level)%alpha(iez+1)

     ! freeslip.......
     grid(level)%topfreeslip = &
          ((1.0 - 2.0*a) / (1.0 + 2.0*a))**(ed-1) * &
          ((1.0+a) /(1.0-a))**ed
     ! dirichlet........
     grid(level)%topdirichlet(1) = ((1.0+a) / (1.0+2.0*a))**(ed-1)
     grid(level)%topdirichlet(2) = ((1.0-2.0*a) / (1.0-a))**(ed-1)
     
     ! BOTTOM - - - - -
     a = grid(level)%alpha(i0)
     ! freeslip.............
     grid(level)%botfreeslip = &
          ((1.0 + 2.0*a) / (1.0 - 2.0*a))**(ed-1) * &
          ((1.0-a) /(1.0+a))**ed
     ! dirichlet...........
     grid(level)%botdirichlet(1) = ((1.0-a) / (1.0-2.0*a))**(ed-1)
     grid(level)%botdirichlet(2) = ((1.0+2.0*a) / (1.0+a))**(ed-1)

  end do ! level
  
  ! ---------------------------------------------
  ! COEFFICIENTS FOR THE SOLVER
  ! ---------------------------------------------
  
  do level=0, nLevel

     iez = iend(level, 2)
     do j=i0-1, iez+1
        a = grid(level)%alpha(j)
        
        grid(level)%gtop(j) = (1.0 + a)**(ed+1) + (1.0 + a)**(1-ed)        
        grid(level)%gbot(j) = (1.0 - a)**(ed+1) + (1.0 - a)**(1-ed)
        grid(level)%gmid(j) = (1.0 + a)**(-ed) + (1.0 - a)**(-ed)
        grid(level)%g1top(j) = (1.0 + 2.0 * a)**(ed+2) / (1.0 + a)**(2*ed)
        grid(level)%g1bot(j) = (1.0 - 2.0 * a)**(ed+2) / (1.0 - a)**(2*ed)

        grid(level)%gcentertop(j) = (1.0 + a)**(ed+1) / (1.0 + 2.0 * a) + &
             (1.0 + a)**(1-ed) * (1.0 + 2.0 * a)**(ed-2)
        grid(level)%gcenterbot(j) = (1.0 - a)**(ed+1) / (1.0 - 2.0 * a) + &
             (1.0 - a)**(1-ed) * (1.0 - 2.0 * a)**(ed-2)

        grid(level)%gcornertop(j) = (1.0 + 2.0 * a)**(ed-1) / (1.0 + a)**ed
        grid(level)%gcornerbot(j) = (1.0 - 2.0 * a)**(ed-1) / (1.0 - a)**ed

        grid(level)%g4(j) = (1.0 + 2.0 * a)**(2*ed+1) / (1.0 + a)**ed * &
             ( (1.0 + a)**(-ed) + (1.0 + 3.0 * a)**(-ed) ) 
        grid(level)%g5(j) = (1.0 - 2.0 * a)**(2*ed+1) / (1.0 - a)**ed * &
             ( (1.0 - a)**(-ed) + (1.0 - 3.0 * a)**(-ed) )

        grid(level)%g12(j) = (1.0 + 2.0 * a)**(ed+2) * (1.0 + 4.0 * a)**(ed-1) * &
             (1.0 + a)**(-ed) * (1.0 + 3.0 * a)**(-ed)
        grid(level)%g13(j) = (1.0 - 2.0 * a)**(ed+2) * (1.0 - 4.0 * a)**(ed-1) * &
             (1.0 - a)**(-ed) * (1.0 - 3.0 * a)**(-ed)
     end do ! j 
  end do ! level
  
  return
end subroutine geometricalCoefficients
! =======================================================
subroutine poissonCoefficients()
  ! Compute coefficient for simple Poisson solver (Lap(f) = RHS)
  use general
  use grids
  implicit none
  integer:: level, iez, j
  real(rk), pointer:: pa(:) => null()
  real(rk), pointer:: pr(:) => null()
  real(rk), pointer:: pp(:) => null()
  real(rk), pointer:: pptop(:) => null()
  real(rk), pointer:: ppbot(:) => null()
  real(rk), pointer:: ppmid(:) => null()
  real(rk):: invd2(1:2), invr2

  call getInvD2(invd2, level)

  call associatePtrLine(pp, poisson, level)
  call associatePtrLine(pptop, poissontop, level)
  call associatePtrLine(ppbot, poissonbot, level)
  call associatePtrLine(ppmid, poissonmid, level)

  if(curvedGeometry) then !______________________
     do level = 0, nLevel

        call associatePtrLine(pa, alpha, level)
        call associatePtrLine(pr, r, level)

        iez = iend(level, 2)

        do j=i0-1, iez+1
           invr2 = 1.0 / pr(j)**2

           pp(j) = 2.0 * invd2(1) * invr2 + &
                invd2(2) * ( (1.0 + pa(j))**ed + (1.0 - pa(j))**ed )

           pptop(j) = invd2(2) * (1.0 + pa(j))**ed
           ppbot(j) = invd2(2) * (1.0 - pa(j))**ed
           ppmid(j) = invd2(1) * invr2
        end do ! j

        nullify(pa, pr)

     end do ! level
  else ! Cartesian geometry _________________________
     ! only needed to have a uniform writing with a
     ! curved or cartesian geometry
     do level=0, nLevel

        iez = iend(level, 2)

        do j=i0-1, iez+1
           pp(j) = 2.0 * (invd2(1) + invd2(2)) ! 2 * (1/dx^2 + 1/dz^2)

           pptop(j) = invd2(2)
           ppbot(j) = invd2(2)
           ppmid(j) = invd2(1)
        end do ! j
        
     end do ! level
  end if ! Curved geometry or not ___________________

  nullify(pp, pptop, ppbot, ppmid)
  
  return
end subroutine poissonCoefficients
! =======================================================
