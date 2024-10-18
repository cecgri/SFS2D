! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine computeHeatFlux(qmeanTop, qmeanBot)
  use general
  use grids
  implicit none
  real(rk), intent(out):: qmeanTop, qmeanBot
  real(rk), pointer:: pt(:,:) => null()
  real(rk), pointer:: pq(:,:) => null()
  integer:: i, iex, iez
  character(len=10):: flux=work3 ! temporary

 ! take temperature on x-faces for simple trapeze method
  call centerToFace(temp, 0, 1) 
  call associatePtrScal(pt, temp, 0, 1)

  call associatePtrScal(pq, flux, 0, 4)

  iex = iend0(1)
  iez = iend0(2)

  ! top flux - - - -
  do i=i0, iex+1
     pq(i, iez+1) = invdz * (pt(i, iez) - pt(i, iez+1))
  end do ! i
  
  call getMeanPlan(qmeanTop, flux, 0, 4, 2, iez+1)

  ! bottom flux - - - - - - 
  do i=i0, iex+1
     pq(i, i0) = invdz * (pt(i, i0-1) - pt(i, i0))
  end do ! i
  
  call getMeanPlan(qmeanBot, flux, 0, 4, 2, i0)
  
  nullify(pt, pq)

  return
end subroutine computeHeatFlux
! =======================================================
subroutine computeVsurf(vsurfRms, vsurfMax)
  use general
  use grids
  implicit none
  real(rk), intent(out):: vsurfRms, vsurfMax
  real(rk), pointer:: pv(:,:) => null()
  real(rk):: vsurf
  integer:: i, iez

  call associatePtrScal(pv, vel, 0, 1)
  vsurfMax = 0.0
  vsurfRms = 0.0

  iez = iend0(2)
  do i=i0, iend0(1)
     vsurf = 0.5 * (pv(i, iez) + pv(i, iez+1))
     vsurfRms = vsurfRms + vsurf*vsurf
     vsurf = abs(vsurf)
     vsurfMax = max(vsurf, vsurfMax)
  end do

  vsurfRms = sqrt(vsurfRms / dble(n0(1)))

  nullify(pv)

  return
end subroutine computeVsurf
! =======================================================
subroutine computeVmid(vrms, vmax, percentiles)
  use general
  use grids
  implicit none
  real(rk), intent(out):: vrms, vmax
  real(rk), dimension(1:5), intent(out):: percentiles
  real(rk), pointer:: pv(:,:,:) => null()
  real(rk), pointer:: pvmid(:,:) => null()
  character(len=10), parameter:: vmid=work1
  integer:: i, j

  call associatePtrVec(pv, vel, 0)
  call associatePtrScal(pvmid, vmid, 0, 3)
  do j=i0, iend0(2)+1
     do i=i0, iend0(1)+1
        pvmid(i, j) = &
             ( pv(i, j, 1) + pv(i+1, j, 1) )**2 + &
             ( pv(i, j, 2) + pv(i, j+1, 2) )**2
        pvmid(i, j) = sqrt(fourth * pvmid(i, j))
     end do
  end do

  nullify(pv, pvmid)

  call getRms(vrms, vmid, 0)
  call getMax(vmax, vmid, 0, 3, borders=.true.)
  call getPercentiles(vmid, percentiles, 3)

  return
end subroutine computeVmid
! =======================================================
subroutine computeNCells(ncells)
  ! ---------------------------------------------
  ! Computes the number of cells, using the
  ! stream function and how many times it changes
  ! sign at mid depth.
  ! Uses j50 for the nodes at mid-depth
  ! ---------------------------------------------
  use general
  use grids
  implicit none
  real(rk), intent(out):: ncells ! could be integer, but easier for writings
  real(rk), pointer:: psf(:, :) => null()
  integer:: i
  real(rk):: currentsign
  
  call associatePtrScal(psf, sf, 0, 4)  ! at cells' corners
  ncells = zero
  currentsign = sign(one, psf(i0, j50))
  i = i0+1
  do while(i <= iend0(1) + 1)
     if(sign(one, psf(i, j50)) /= currentsign) then
        currentsign = -currentsign
        ncells = ncells + one
     end if
     i = i+1
  end do

  nullify(psf)
  return
end subroutine computeNCells
! =======================================================
subroutine computeVsurfEdot(vsurfrms, vsurfmax, percentiles, edotmax)
  use general
  use grids
  implicit none
  real(rk), intent(out):: vsurfrms, vsurfmax, edotmax
  real(rk), dimension(1:5):: percentiles
  real(rk), pointer:: pvx(:,:) => null()
  real(rk), dimension(n0(1)):: vec
  real(rk):: coeffE
  integer:: i, iez, n

  call faceToCorner(vel, 0, 1, plan=iend0(2)+1)
  ! vel(:,:,4) is now v(:,:,1) at the surface
  vsurfRms = 0.0
  vsurfMax = 0.0
  call getMeanPlan(vsurfRms, vel, 0, 4, 2, iend0(2)+1, rms=.true.)
  call getMaxPlan(vsurfMax, vel, 0, 4, 2, iend0(2)+1, absval=.true.)

  call associatePtrScal(pvx, vel, 0, 1)
  iez = iend0(2)
  if(curvedGeometry) then
     coeffE = invdx / rtop
  else ! Cartesian
     coeffE = invdx
  end if
  do i=i0, iend0(1)
     ! xval : epsilon_xx
     xval(1, i) =  coeffE * ( &
             pvx(i+1, iez) + pvx(i+1, iez+1) - &
             pvx(i, iez) - pvx(i, iez+1) )
     xval(1, i) = abs(xval(1, i))
  end do ! i

  nullify(pvx)

  edotMax = 0.0
  do n=1, n0(1)
     vec(n) = xval(1, n+i0-1)
     edotMax = max(edotMax, vec(n))
  end do
  call sort(vec, n0(1))
  n = 1
  i = 1
  do
     if(dble(n) / dble(n0(1)) > percValues(i)) then
        percentiles(i) = vec(n)
        i = i+1
     end if
     if(i > 5) exit
     n = n+1
  end do

  return
end subroutine computeVsurfEdot
! =======================================================
subroutine computeVcells(vtopmax, vbotmax, vupmax, vdownmax, vtoprms, vbotrms)
  ! March 2022: adding the computation of the maximum of horizontal
  ! velocity in the top half and bottom half of the layer, and
  ! of the max and min of the vertical velocity (up and down plumes);
  ! March 2024: added vtoprms and vbotrms
  use general
  use grids
  implicit none
  real(rk), intent(out):: vtopmax, vbotmax, vupmax, vdownmax, vtoprms, vbotrms
  real(rk), pointer:: pvx(:,:) => null()
  real(rk), pointer:: pvxtop(:,:) => null()
  real(rk), pointer:: pvxbot(:,:) => null()
  character(len=10), parameter:: vtop=work1, vbot=work2  ! tempo
  real(rk):: posZ, zmid
  integer:: i, j
  
  call associatePtrScal(pvx, vel, 0, 1)
  call associatePtrScal(pvxtop, vtop, 0, 1)  ! work1
  call associatePtrScal(pvxbot, vbot, 0, 1)  ! work2
  pvxtop = 0.0
  pvxbot = 0.0

  zmid = half  ! middle of the model if Cartesian
  if(curvedGeometry) then
     zmid = half * (rtop + rbot)
  end if

  do j=i0, iend0(2)+1
     call getPosZ(posZ, j, 0, 1)
     if(posZ > zmid) then
        do i=i0, iend0(1)+1
           pvxtop(i, j) = pvx(i, j)
        end do
     else
        do i=i0, iend0(1)+1
           pvxbot(i, j) = pvx(i, j)
        end do
     end if
  end do

  call getMax(vtopmax, vtop, 0, 1, absVal=.true.)
  call getMax(vbotmax, vbot, 0, 1, absVal=.true.)
  call getMax(vupmax, vel, 0, 2)
  call getMin(vdownmax, vel, 0, 2)

  call getMeanPlan(vtoprms, vel, 0, 1, 2, iend0(2), rms=.true.)
  call getMeanPlan(vbotrms, vel, 0, 1, 2, i0, rms=.true.)
  
  nullify(pvx, pvxtop, pvxbot)
  
  return
end subroutine computeVcells
! =======================================================
subroutine computeExtra(etamin, etamax, work, dissipation)
  ! ------------------------------------------------------
  ! computes some extra diagnosis that
  ! are included in the benchmark by Tosi et al., 2015, G3
  ! -----------------------------------------------------
  use general
  use grids
  implicit none
  real(rk), intent(out):: etamin, etamax, work, dissipation
  real(rk), pointer:: pvz(:, :) => null()
  real(rk), pointer:: pt(:, :) => null()
  real(rk), pointer:: peta(:, :) => null()
  real(rk), pointer:: psist(:, :) => null() !Second Invariant of the Strainrate Tensor
  real(rk), pointer:: pwork(:, :) => null()
  real(rk), pointer:: pdiss(:, :) => null()
  integer:: i, j
  character(len=10):: wrk=work4, diss=work5  ! temporary

  call getMin(etamin, etaEff, 0, 4, borders=.true.)
  call getMax(etamax, etaEff, 0, 4, borders=.true.)

!!$  ! all computed at cells' corners:
  call associatePtrScal(pt, temp, 0, 4)
  call faceToCorner(vel, 0, 2)  ! vz put at cells' corners
  call associatePtrScal(pvz, vel, 0, 4)
  call associatePtrScal(peta, etaEff, 0, 4)
  ! second invariant of the strainrate tensor: 
  call associatePtrScal(psist, sist, 0, 4) ! must have been computed before

  call associatePtrscal(pwork, wrk, 0, 4)
  call associatePtrscal(pdiss, diss, 0, 4)
  ! work = vz * temp
  ! dissipation = tauij:eij = 2 * eta * eij:eij
  !             = 2 * eta * (2 * sist) = 4 * eta * sist
  do j=i0, iend0(2)+1
     do i=i0, iend0(1)+1
        pwork(i, j) = pt(i, j) * pvz(i, j)
        pdiss(i, j) = four * peta(i, j) * psist(i, j)
     end do
  end do
  
  call getMean(work, wrk, 0)
  call getMean(dissipation, diss, 0)

  nullify(pvz, pt, peta, psist, pwork, pdiss)
  return
end subroutine computeExtra
! =======================================================
