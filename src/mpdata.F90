! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine mpdata(what)
  ! Multidimensional Positive Definite Advection
  ! Transport Algorithm, by Smolarkiewicz (1984)
  ! - - - - - - - 
  ! advects field that can only be in position "3" (center of cell)
  ! (used for temperature)
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, parameter:: iord=3
  integer:: io

  call diffusion(what) ! check if order changes something...
  call copyVec(work1, vel, 0) ! copy v in work1

  ! MPDATA classique - order=iord
  do io=1, iord

     call upwind(what, work1)

     if(io < iord) then
        call antiveloc(work2, what, work1)
        call copyVec(work1, work2, 0)
     end if

  end do ! io=1, iord

  return
end subroutine mpdata
! =======================================================
subroutine diffusion(what)
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  real(rk), pointer:: ptr(:,:) => null()
  real(rk), dimension(i0-2:iend0(1)+2, i0-2:iend0(2)+2):: diff
  integer:: i,j
  real(rk), pointer:: pr(:) => null()
  real(rk), pointer:: pa(:) => null()
  real(rk):: invd2(1:2), coeffgeom(1:2)
  
  call associatePtrScal(ptr, what, 0, 3) !level 0, id=3 for cells' centers
  call associatePtrLine(pr, r, 0)
  call associatePtrLine(pa, alpha, 0)
  call getInvD2(invd2, 0) ! 1/dx^2 and 1/dz^2
 
  diff = 0.0
  ! uniform conductivity
  if(curvedGeometry) then ! ______________________________
     do j=i0, iend0(2)
        coeffgeom(1) = invd2(1) / (pr(j) * (1.0 + pa(j)))**2
        coeffgeom(2) = invd2(2) / (1.0 + pa(j))**ed
        do i=i0, iend0(1)
           diff(i,j) = &
                coeffgeom(1) * ( ptr(i+1,j) + ptr(i-1,j) - 2.0 * ptr(i,j) ) + &
                coeffgeom(2) * ( (1.0 + 2.0 * pa(j))**ed * ptr(i,j+1) + ptr(i,j-1) - &
                ((1.0 + 2.0 * pa(j))**ed + 1.0) * ptr(i,j) )
        end do
     end do
  else ! Cartesian ______________________________________
     do j=i0, iend0(2)
        do i=i0, iend0(1)
           diff(i,j) = &
                invd2(1) * ( ptr(i+1,j) + ptr(i-1,j) - 2.0 * ptr(i,j) ) + &
                invd2(2) * ( ptr(i,j+1) + ptr(i,j-1) - 2.0 * ptr(i,j) )
        end do
     end do
  end if ! Curved geometry or not _______________________
  
  heat(:,:) = internalHeating

  do j=i0, iend0(2)
     do i=i0, iend0(1)
        ptr(i,j) = ptr(i,j) + timestep * (diff(i,j) + heat(i, j))
     end do
  end do
  
  call applyBcs(what, 0, 3)

  nullify(ptr, pa, pr)
  return
end subroutine diffusion
! =======================================================
subroutine upwind(what, whatvel)
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what, whatvel
  real(rk), dimension(1:n0(1)+2, 1:n0(2)+2):: adv
  real(rk), pointer:: ptr(:,:) => null()
  real(rk), pointer:: pv(:,:,:) => null()
  real(rk), pointer:: pa(:) => null()
  real(rk), pointer:: pr(:) => null()
  integer::i, j
  real(rk):: invd(1:2), coeffgeom(1:2)

  call associatePtrScal(ptr, what, 0, 3)
  call associatePtrVec(pv, whatvel, 0)
  call getInvd(invd, 0) ! 1/dx and 1/dz
  
  adv(:,:) = 0.0

  if(curvedGeometry) then ! _______________________________
     
     call associatePtrLine(pr, r, 0)
     call associatePtrLine(pa, alpha, 0)

     do j=i0, iend0(2)
        coeffgeom(1) = invd(1) / (pr(j) * (1.0 + pa(j)))
        coeffgeom(2) = invd(2) / (1.0 + pa(j))**ed
        do i=i0, iend0(1)
           adv(i,j) = coeffgeom(1) * ( &
                flux1d(pv(i+1,j,1), ptr(i,j), ptr(i+1,j)) - &
                flux1d(pv(i,j,1), ptr(i-1,j), ptr(i,j)) ) + &
                coeffgeom(2) * ( &
                (1.0 + 2.0*pa(j))**ed * flux1d(pv(i,j+1,2), ptr(i,j), ptr(i,j+1)) - &
                flux1d(pv(i,j,2), ptr(i,j-1), ptr(i,j)) )
        end do ! i
     end do ! j

     nullify(pr, pa)
     
  else ! Cartesian ______________________________________
     
     do j=i0, iend0(2)
        do i=i0, iend0(1)
           adv(i,j) = invd(1) * ( &
                flux1d(pv(i+1,j,1), ptr(i,j), ptr(i+1,j)) - &
                flux1d(pv(i,j,1), ptr(i-1,j), ptr(i,j)) ) + &
                invd(2) * ( &
                flux1d(pv(i,j+1,2), ptr(i,j), ptr(i,j+1)) - &
                flux1d(pv(i,j,2), ptr(i,j-1), ptr(i,j)) )
        end do ! i
     end do ! j
     
  end if ! Curved geometry or not _______________________
     
  do j=i0, iend0(2)
     do i=i0, iend0(1)
        ptr(i,j) = ptr(i,j) - timestep * adv(i,j)
     end do
  end do

  call applyBcs(what, 0, 3)

  nullify(ptr, pv)

  return

contains
  function flux1d(veloc, a1, a2)
    ! a1: a(i), a2: a(i+1)
    use general
    implicit none
    real(rk):: flux1d
    real(rk), intent(in):: veloc, a1, a2
    flux1d = 0.5 * &
         ( (veloc + abs(veloc))*a1 + (veloc - abs(veloc))*a2 )
  end function flux1d

end subroutine upwind
! =======================================================
subroutine antiveloc(vmod, what, vin)
  ! vmod: computed/modified from vin and what
  ! WARNING: method quickly implemented for cylindrical
  !   and spherical annulus version. To check!!!
  !   Recommanded to use ADI and not mpdata
  ! - - - - - - - - -
  use general
  use grids
  implicit none
  character(len=10), intent(in):: vmod, vin, what
  real(rk), pointer:: pv(:,:,:) => null()
  real(rk), pointer:: pvmod(:,:,:) => null()
  real(rk), pointer:: ptr(:,:) => null()
  integer:: i,j
  real(rk):: den, dhh

  call zeroingVec(vmod, 0)
  call associatePtrVec(pv, vin, 0)
  call associatePtrVec(pvmod, vmod, 0)
  call associatePtrScal(ptr, what, 0, 3)

  ! horizontal direction
  do j=i0, iend0(2)
     do i=i0, iend0(1)
        den = ptr(i,j) + ptr(i-1,j)
        if(den /= 0.0) then
           pvmod(i,j,1) = &
                ( abs(pv(i,j,1)) * dh(j) - timestep * pv(i,j,1)**2 ) * &
                (ptr(i,j) - ptr(i-1,j)) / den
        end if
        
        ! cross-term
        den = ptr(i,j+1) + ptr(i,j-1) + ptr(i-1,j+1) + ptr(i-1,j-1)
        if(den /= 0.0) then
           pvmod(i,j,1) = pvmod(i,j,1) &
                - 0.5 * timestep * pv(i,j,1) * &
                0.25 * dz * &
                (pv(i,j,2) + pv(i,j+1,2) + pv(i-1,j,2) + pv(i-1,j+1,2)) * &
                (ptr(i,j+1) - ptr(i,j-1) + ptr(i-1,j+1) - ptr(i-1,j-1)) / den
        end if
     end do ! i
  end do ! j

  ! vertical direction
  do j=i0, iend0(2)
     if(curvedGeometry) then
        dhh = dx * r0(j)
     else
        dhh = dx
     end if

     do i=i0, iend0(1)
        den = ptr(i,j) + ptr(i,j-1)
        if(den /= 0.0) then
           pvmod(i,j,2) = &
                ( abs(pv(i,j,2)) * dz - timestep * pv(i,j,2)**2 ) * &
                (ptr(i,j) - ptr(i,j-1)) / den
        end if

        ! cross-term
        den = ptr(i+1,j) + ptr(i-1,j) + ptr(i+1,j-1) + ptr(i-1,j-1)
        if(den /= 0.0) then
           pvmod(i,j,2) = pvmod(i,j,2) &
                - 0.5 * timestep * pv(i,j,2) * &
                0.25 * dhh * &
                (pv(i,j,1) + pv(i+1,j,1) + pv(i,j-1,1) + pv(i+1,j-1,1)) * &
                (ptr(i+1,j) - ptr(i-1,j) + ptr(i+1,j-1) - ptr(i-1,j-1)) / den
        end if
     end do ! i
  end do ! j
   
  nullify(pv, pvmod, ptr)

  return
end subroutine antiveloc
! =======================================================
