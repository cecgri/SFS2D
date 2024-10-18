! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS_2D
! =======================================================
subroutine ADI()
  ! Alternating Direction Implicit method
  !  for the advection and diffusion of the temperature.
  use general
  use grids
  implicit none
  integer:: id

  do id=1,2
     call fluxLimiters(temp, id, TVDscheme_temp) ! computes flim_m,.._p if needed
     call coeffAdv(vel, id) ! computes cadv_m, .._c and .._p
  end do

  do id=1,2
     call advdiff_1d(temp, id)       
     call applyBcs(temp, 0, 3)
  end do

  return
end subroutine ADI
! =======================================================
subroutine advdiff_1d(what, idir)
  ! advection AND diffusion in the direction idir of what->pt
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: idir
  real(rk), pointer:: pt(:,:) => null()
  real(rk), dimension(i0-1:iend0(1)+1, i0-1:iend0(2)+1):: aa, bb, cc, ff
  real(rk), dimension(i0-1:iend0(1)+1, i0-1:iend0(2)+1):: sol
  real(rk):: coeffCC, opp, firstVal, lastVal
  integer:: i, j, nsize

  call associatePtrScal(pt, what, 0, 3)
  
  ! uniform heating
  heat(:,:) = internalHeating 
  
  if(idir == 1) then  ! phi - - -
     do j=i0, iend0(2)
        do i=i0, iend0(1)
           aa(i,j) = -invdh(j) + cadv_m(i,j,1)
           bb(i,j) = 4.0 * dh(j) / timestep + 2.0 * invdh(j) + cadv_c(i,j,1)
           cc(i,j) = -invdh(j) + cadv_p(i,j,1)

           coeffCC = 4.0 * dh(j) / timestep - 2.0 * invdh(j) - cadv_c(i,j,1) - &
                2.0 * dh(j) * invdv(j) * ( &
                invdvp(j) + invdz + cadv_c(i,j,2) )
           opp = pt(i,j-1) * (invdz - cadv_m(i,j,2)) + &
                pt(i,j+1) * (invdvp(j) - cadv_p(i,j,2))
           
           ff(i,j) = 2.0 * heat(i,j) * dh(j) + &
                pt(i-1,j) * (invdh(j) - cadv_m(i,j,1)) + &
                pt(i+1,j) * (invdh(j) - cadv_p(i,j,1)) + &
                pt(i,j) * coeffCC + &
                2.0 * dh(j) * invdv(j) * opp
        end do ! i
     end do ! j
  else ! idir == 2  ! r - - -
     do j=i0, iend0(2)
        do i=i0, iend0(1)
           aa(i,j) = -invdz + cadv_m(i,j,2)
           bb(i,j) = 4.0 * dv(j) / timestep + invdvp(j) + invdz + cadv_c(i,j,2)
           cc(i,j) = -invdvp(j) + cadv_p(i,j,2)

           coeffCC = 4.0 * dv(j) / timestep - invdvp(j) - invdz - cadv_c(i,j,2) - &
                2.0 * dv(j) * invdh(j) * ( 2.0 * invdh(j) + cadv_c(i,j,1) )
           opp = pt(i-1,j) * (invdh(j) - cadv_m(i,j,1)) + &
                pt(i+1,j) * (invdh(j) - cadv_p(i,j,1))
           
           ff(i,j) = 2.0 * heat(i,j) * dv(j) + &
                pt(i,j-1) * (invdz - cadv_m(i,j,2)) + &
                pt(i,j+1) * (invdvp(j) - cadv_p(i,j,2)) + &
                pt(i,j) * coeffCC + &
                2.0 * dv(j) * invdh(j) * opp
        end do ! i
     end do ! j
  end if ! idir=1,2
    
  if(idir == 1) then
     nsize = iend0(1) - i0 + 1
     if(leftT == dirichlet) then
        firstVal = leftValueT
     else if(leftT == neumann) then
        firstVal = 0.0
     end if
     if(rightT == dirichlet) then
        lastVal = rightValueT
     else if(rightT == neumann) then
        lastVal = 0.0
     end if
     ! firstVal and lastVal are not used if periodic
     do j=i0, iend0(2)
        call solve_tridiag( sol(i0:iend0(1), j), &
             aa(i0:iend0(1), j), bb(i0:iend0(1), j), cc(i0:iend0(1), j), &
             ff(i0:iend0(1), j), nsize, firstVal, lastVal, idir )
        !first and last val are not used for idir=1
     end do ! j
  else if(idir == 2) then
     nsize = iend0(2) - i0 + 1
     lastVal = topValueT ! doesnt cover special topT boundary conditions
     if(botT == dirichlet) then
        firstVal = botValueT
     else if(botT == neumann) then
        firstVal = botValueT * dz
     end if
     do i=i0, iend0(1)
        call solve_tridiag( sol(i, i0:iend0(2)), &
             aa(i, i0:iend0(2)), bb(i, i0:iend0(2)), cc(i, i0:iend0(2)), &
             ff(i, i0:iend0(2)), nsize, firstVal, lastVal, idir )
     end do ! i
  end if

  do j=i0, iend0(2)
     do i=i0, iend0(1)
        pt(i, j) = sol(i, j)
     end do !i
  end do !j
  nullify(pt)
  
  return
end subroutine advdiff_1d
! =======================================================
! flux limiter and coefficients for advection
! =======================================================
subroutine fluxLimiters(what, idir, scheme)
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: idir
  character(len=20):: scheme
  real(rk), pointer:: pt(:,:) => null()
  real(rk), dimension(i0-1:iend0(1)+1, i0-1:iend0(2)+1):: dif, rp, rm
  integer:: kx, kz, i, j

  call associatePtrScal(pt, what, 0, 3)
  
  if(scheme == "upwind") then
     flim_m(:,:,:) = 0.0
     flim_p(:,:,:) = 0.0
  else if(scheme == "centered") then
     flim_m(:,:,:) = 1.0
     flim_p(:,:,:) = 1.0
  else
     kx = kr(1,idir)
     kz = kr(2,idir)

     ! grad of pt - - - - - - - - -
     do j=i0, iend0(2)+1
        do i=i0, iend0(1)+1
           dif(i,j) = pt(i,j) - pt(i-kx,j-kz)
        end do
     end do
     ! bcs for grad - - - - - - - -
     ! zero flux at the bottom
     dif(:,i0-1) = dif(:,i0)
     ! vertical walls
     if(leftT == periodic) then ! then rightV periodic too
        dif(i0-1,:) = dif(iend0(1),:)
        dif(iend0(1)+1,:) = dif(i0,:)
     else
        ! MAYBE TO RECHECK AND MODIFY....
        dif(i0-1,:) = dif(i0,:)
        dif(iend0(1)+1,:) = dif(iend0(1),:)
     end if

     ! ratio of grads - - - - - - - -
     do j=i0, iend0(2)+1
        do i=i0, iend0(1)+1
           if(dif(i,j) /= 0.0) then
              rp(i,j) = dif(i-kx,j-kz) / dif(i,j)
              if(i+kx < iend0(1)+1 .and. j+kz < iend0(2)+1) &
                   rm(i,j) = dif(i+kx,j+kz) / dif(i,j)
           else
              rp(i,j) = 0.0
              rm(i,j) = 0.0
           end if
        end do ! i 
     end do ! j
     ! topbcs
     rm(:,iend0(2)+1) = rm(:,iend0(2))
     ! vertical wall
     if(leftT == periodic) then
        rm(iend0(1)+1,:) = rm(i0,:)
     else
        ! MAYBE TO MODIFY?
        rm(iend0(1)+1,:) = rm(iend0(1),:)
     end if

     ! compute flux limiters
     if(scheme == "minmod") then
        do j=i0, iend0(2)+1
           do i=i0, iend0(1)+1
              flim_p(i,j,idir) = max(0.0, min(1.0, rp(i,j)))
              flim_m(i,j,idir) = max(0.0, min(1.0, rm(i,j)))
           end do
        end do
     else if(scheme == "superbee") then
        do j=i0, iend0(2)+1
           do i=i0, iend0(1)+1
              flim_p(i,j,idir) = max( 0.0, min(1.0, 2.0*rp(i,j)), min(2.0, rp(i,j)) )
              flim_m(i,j,idir) = max( 0.0, min(1.0, 2.0*rm(i,j)), min(2.0, rm(i,j)) )
           end do
        end do
     else if(scheme == "koren") then
        do j=i0, iend0(2)+1
           do i=i0, iend0(1)+1
              flim_p(i,j,idir) = max( 0.0, min(2.0*rp(i,j), 0.5 + 0.5*rp(i,j), 2.0) )
              flim_m(i,j,idir) = max( 0.0, min(2.0*rm(i,j), 0.5 + 0.5*rm(i,j), 2.0) )
           end do
        end do
     end if ! koren / superbee / minmod

  end if ! scheme upwind / centered / others

  nullify(pt)
  
  return
end subroutine fluxLimiters
! =======================================================
subroutine coeffAdv(what, idir)
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: idir
  integer:: kx, kz, i, j
  real(rk):: vc, vp, vcp, vcm, vpp, vpm
  real(rk), pointer:: pa(:) => null()
  real(rk), pointer:: pv(:,:,:) => null()
  real(rk):: coeffgeom
  
  kx = kr(1,idir)
  kz = kr(2,idir)

  call associatePtrVec(pv, what, 0)
  call associatePtrLine(pa, alpha, 0)
  coeffgeom = 1.0
     
  do j=i0, iend0(2)
     if(idir == 2 .and. curvedGeometry) then
        coeffgeom = (1.0 + 2.0 * pa(j))**ed
     end if
     do i=i0, iend0(1)
        vc = pv(i,j,idir)
        vp = pv(i+kx,j+kz,idir)
        vcp = 0.5 * (vc + abs(vc))
        vcm = 0.5 * (vc - abs(vc))
        vpp = 0.5 * (vp + abs(vp))
        vpm = 0.5 * (vp - abs(vp))

        cadv_m(i,j,idir) = &
             - vcp * (1.0 - 0.5*flim_p(i,j,idir)) &
             - vcm * 0.5*flim_m(i,j,idir)
        cadv_c(i,j,idir) = &
             coeffgeom * vpp * (1.0 - 0.5*flim_p(i+kx,j+kz,idir)) + &
             coeffgeom * vpm * 0.5*flim_m(i+kx,j+kz,idir) - &
             vcp * 0.5*flim_p(i,j,idir) - &
             vcm * (1.0 - 0.5*flim_m(i,j,idir))
        cadv_p(i,j,idir) = &
             coeffgeom * vpp * 0.5*flim_p(i+kx,j+kz,idir) + &
             coeffgeom * vpm * (1.0 - 0.5*flim_m(i+kx,j+kz,idir))
     end do ! i
  end do ! j

  nullify(pv, pa)
  
  return
end subroutine coeffAdv
! =======================================================
subroutine solve_tridiag(sol, a, b, c, rhs, n, firstVal, lastVal, idir)
  use general
  use mylib
  implicit none
  integer, intent(in):: n, idir
  real(rk), intent(in):: firstVal, lastVal
  real(rk), dimension(1:n):: sol, a, b, c, rhs
  logical:: truetridiag

  truetridiag = .true.

  if(idir == 1) then
     if(leftT == periodic) then  ! rightT also periodic
        truetridiag = .false.
        if(a(1)==0 .and. c(n)==0) &
             truetridiag = .true.
     else
        if(leftT == dirichlet) then
           b(1) = b(1) - a(1)
           rhs(1) = rhs(1) - 2.0 * a(1) * firstVal
        else if(leftT == neumann) then
           b(1) = a(1) + b(1)
           rhs(1) = rhs(1) + a(1) * firstVal
        end if

        if(rightT == dirichlet) then
           b(n) = b(n) - c(n)
           rhs(n) = rhs(n) - 2.0 * c(n) * lastVal
        else if(rightT == neumann) then
           b(n) = b(n) + c(n)
           rhs(n) = rhs(n) - c(n) * lastVal
        end if
     end if     
  else
     if(botT == dirichlet) then
        b(1) = b(1) - a(1)
        rhs(1) = rhs(1) - 2.0 * a(1) * firstVal
     else if(botT == neumann) then
        b(1) = a(1) + b(1)
        rhs(1) = rhs(1) + a(1) * firstVal
     end if
     if(topT == dirichlet) then
        b(n) = b(n) - c(n)
        rhs(n) = rhs(n) - 2.0 * c(n) * lastVal
     else if(topT == neumann) then
        b(n) = b(n) + c(n)
        rhs(n) = rhs(n) - c(n) * lastVal
     end if
  end if

  if(truetridiag) then
     call tridiag(sol, a, b, c, rhs, n)
  else
     call nearlyTridiag(sol, a, b, c, rhs, n)
  end if

  return
end subroutine solve_tridiag
! =======================================================
