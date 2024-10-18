! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine newtonianVisco(what, id)
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: id
  real(rk), pointer:: peta(:,:) => null()
  real(rk), pointer:: ptemp(:,:) => null()
  integer:: i, j, ibord

  if(isoviscous) return 

  call associatePtrScal(peta, what, 0, id)
  call associatePtrScal(ptemp, temp, 0, id)

  ibord = 0
  if(id == 4) ibord = 1
  
  ! thermal viscosity
  if(viscoLaw == "Arrhenius") then
     !note: coeffViscoArr_z(j) = Ea + Va*(1-z)
     do j=i0, iend(0,2) + ibord
        do i=i0, iend(0,1) + ibord
           peta(i, j) = coeffVisco * &
                exp(coeffViscoArr_z(j, id) / (ptemp(i, j) + surfTAdim + ToffAdim))
        end do !i
     end do !j
  else if(viscoLaw == "FK") then
     !note: coeffViscoFK_z(j) = exp(-gammaVisco * z)
     do j=i0, iend(0,2) + ibord
        do i=i0, iend(0,1) + ibord
           peta(i,j) = coeffViscoFK_z(j, id) * coeffVisco * &
                exp(-thetaVisco * ptemp(i, j))
        end do ! i
     end do ! j
  end if ! viscolaw

  if(cutOffMax) then
     if(.not.smoothCutOff) then
        do j=i0, iend(0,2) + ibord
           do i=i0, iend(0,1) + ibord
              peta(i,j) = min(peta(i,j), viscoCutOffMax)
           end do !i
        end do !j
     else
        do j=i0, iend(0,2) + ibord
           do i=i0, iend(0,1) + ibord
              peta(i,j) = 1.0 / (1.0 / peta(i,j) + 1.0 / viscoCutOffMax)
           end do !i
        end do !j
     end if
  end if
  if(cutOffMin) then
     do j=i0, iend(0,2) + ibord
        do i=i0, iend(0,1) + ibord
           peta(i,j) = max(peta(i,j), viscoCutOffMin)
        end do !i
     end do !j
  end if

  call applyBcs(etaT, 0, id)

  nullify(peta, ptemp)
  
  return
end subroutine newtonianVisco
! =======================================================
subroutine computeYielding(id)
  ! computes etayield(:,:, id)
  use general
  use grids
  implicit none
  integer, intent(in):: id
  real(rk), pointer:: psist(:,:) => null()
  real(rk), pointer:: pexx(:,:) => null()
  real(rk):: edot, yieldVal
  integer:: ibord, i, j
  
  call associatePtrScal(psist, sist, 0,  id)

  ibord = 0
  if(id == 4) ibord = 1

  etayield(:,:,id) = infinity

  if(.not. anderson) then
     do  j=i0, iend(0, 2) + ibord
        do i=i0, iend(0, 1) + ibord
           edot = 2.0 * sqrt(psist(i, j))
           if(edot > 0.0) then
              etayield(i, j, id) = yieldT_z(j, id) / edot + etaStar
           end if
        end do !i
     end do !j
  else ! anderson model
     ! anderson: harder to reach failure in compression
     ! determined by 'compressionFactor' (normally > 1, but possible to
     ! test compression easier than tension)
     ! NOTE: exx has been updated with sf_to_sist() done before

     !!!! WARNING: TO CHECK AND MODIFY FOR CURVED GEOMETRY ?? !!!!!!!!!!!!!
     call associatePtrScal(pexx, exx, 0, id)
     do j=i0, iend(0, 2) + ibord
        do i=i0, iend(0, 1) + ibord              
           edot = 2.0 * sqrt(psist(i, j))
           if(edot > 0.0) then
              yieldVal = &
                   (sign(half, pexx(i, j)) + sign(half, abs(pexx(i, j)))) * &
                   yieldT_z(j, id) - &
                   (sign(half, pexx(i, j)) - sign(half, abs(pexx(i, j)))) * &
                   yieldC_z(j, id)
              
              etayield(i, j, id) = yieldVal / edot + etaStar
           end if
        end do !i
     end do !j
     nullify(pexx)
  end if  ! anderson or not

  nullify(psist)
  return
end subroutine computeYielding
! =======================================================
subroutine effectiveVisco(whatEff, whatLin, id)
  ! effective viscosity with yielding.
  ! whatLin (etaT) and sist (with good id) must be already computed.
  use general
  use grids
  implicit none
  character(len=10), intent(in):: whatEff, whatLin
  integer, intent(in):: id
  real(rk), pointer:: petaT(:,:) => null()
  real(rk), pointer:: petaEff(:,:) => null()
  integer:: i, j, ibord

  call computeYielding(id)
  ! updates etayield(:,:,id)
  
  call associatePtrScal(petaT, whatLin, 0, id)
  call associatePtrScal(petaEff, whatEff, 0, id)
  ibord = 0
  if(id == 4) ibord = 1

  if(methodEtaEff == "harmonic") then
     do  j=i0, iend(0, 2) + ibord
        do i=i0, iend(0, 1) + ibord
           petaEff(i,j) = etaEffHarmonicCoeff / &
                (1.0 / etayield(i, j, id) + 1.0 / petaT(i,j))
        end do ! i
     end do !j
  else if(methodEtaEff == "minharmonic") then
     do  j=i0, iend(0, 2) + ibord
        do i=i0, iend(0, 1) + ibord
           petaEff(i,j) = min(petaT(i, j), &
                1.0 / (1.0 / etayield(i, j, id) + 1.0 / petaT(i, j)))
        end do ! i
     end do !j
  else if(methodEtaEff == "min") then
     do  j=i0, iend(0, 2) + ibord
        do i=i0, iend(0, 1) + ibord
           petaEff(i,j) = min(petaT(i,j), etayield(i, j, id))
        end do ! i
     end do !j
  end if ! methodEtaEff

  nullify(petaEff, petaT)
  return
end subroutine effectiveVisco
! =======================================================
subroutine yieldZ()
  ! initiate yield as a function of depth: yield_z
  use grids
  implicit none
  integer:: j, id
  real(rk):: posZ, postop

  postop = 1.0
  if(curvedGeometry) postop = rtop

  open(unit=22, file="yield_z.out", status='unknown')
  do j=i0-1, iend(0,2)+1
     do id=3, 4
        call getPosZ(posZ, j, 0, id)
        yieldT_z(j, id) = min(yieldFix, &
             yieldSurf + slopeYield * (postop - posZ))

        if(anderson) then
           yieldC_z(j, id) = min(yieldFix, &
                yieldSurf + compressionFactor * slopeYield * (postop - posZ))
           if(id == 4) then
              write(22, '(f6.3, 2e12.4)') posZ, yieldT_z(j, id), yieldC_z(j, id)
           end if ! id=4                      
        else
           yieldC_z(j, id) = yieldT_z(j, id)
           if(id == 4) then
              write(22, '(f6.3, e12.4)') posZ, yieldT_z(j, id)
           end if ! id=4           
        end if

     end do ! id=3, 4
  end do ! j
  close(22)
  
  return
end subroutine yieldZ
! =======================================================
subroutine coeffViscoFKZ()
  ! computes the effect of pressure (z) on the viscosity
  ! when a Frank-Kamenetskii law is used
  ! (effect of z is independent on effect of T then, no need
  ! to be re-computed).
  use grids
  implicit none
  integer:: j, id
  real(rk):: posZ, postop

  postop = 1.0
  if(curvedGeometry) postop = rtop

  do j=i0-1, iend(0, 2)+1
     do id=3, 4
        call getPosZ(posZ, j, 0, id) 
        coeffViscoFK_z(j, id) = exp(-gammaVisco * (postop - rbot))
     end do
  end do
  
  return
end subroutine coeffViscoFKZ
! =======================================================
subroutine coeffViscoArrZ()
  ! computes the prefactor for the effect of pressure (z) on the viscosity
  ! when an Arrhenius law is used
  use grids
  implicit none
  integer:: j, id
  real(rk):: posZ, postop

  postop = 1.0
  if(curvedGeometry) postop = rtop

  do j=i0-1, iend(0, 2)+1
     do id=3, 4
        call getPosZ(posZ, j, 0, id) 
        coeffViscoArr_z(j, id) = Ea + Va * (postop - posZ)
     end do
  end do
  
  return
end subroutine coeffViscoArrZ
! =======================================================
subroutine initVisco()
  ! set if isoviscous or not,
  ! computes the prefactor if reference is bottom or top,
  ! and computes etaTop or etaBot (just for display)
  use general
  use grids
  use io
  implicit none

  if(.not. isoviscous) then
     if(viscoLaw == "FK") then
        if(thetaVisco <= 0.0 .and. gammaVisco <= 0.0) then
           isoviscous = .true.
        end if
     else if(viscoLaw == "Arrhenius") then
        if(Ea <= 0.0 .and. Va <= 0.0) then
           isoviscous = .true.
        end if
     end if
     cutOffMax = .false.
     if(viscoCutOffMax > 0.0) then
        cutOffMax = .true.
     end if
     cutOffMin = .false.
     if(viscoCutOffMin > 0.0) then
        cutOffMin = .true.
     end if     
  end if
  
  if(isoviscous) then

     etaTop = 1.0
     etaBot = 1.0
     call setValueVec(eta, 0, one)
     call setValueVec(etaT, 0, one)
     call setValueVec(etaEff, 0, one) 
     call updateBiHarmCoeff(eta) !-> coeff computed at all levels, only once
     
  else

     ! different viscosity laws - - - - - - 
     if(viscoLaw == "Arrhenius") then

        if(reference == "bottom") then
           coeffVisco = exp(-(Ea + Va) / (botTref + surfTAdim + ToffAdim))
           etaTop = exp(Ea / (topTref + surfTAdim + ToffAdim) - &
                (Ea + Va) / (botTref + surfTAdim + ToffAdim))
           etaBot = 1.0
        else if(reference == "top") then
           coeffVisco = exp(-Ea / (topTref + surfTAdim + ToffAdim))
           etaTop = 1.0
           etaBot = exp((Ea + Va) / (botTref + surfTAdim + ToffAdim) - &
               Ea / (topTref + surfTAdim + ToffAdim))
        else ! reference == "T_ref"
           coeffVisco = exp(-Ea / (TrefAdim + surfTAdim + ToffAdim))
           etaTop = exp(Ea / (topTref + surfTAdim + ToffAdim) - &
                Ea / (TrefAdim + surfTAdim + ToffAdim))
           etaBot = exp((Ea + Va) / (botTref + surfTAdim + ToffAdim) - &
                Ea / (TrefAdim + surfTAdim + ToffAdim))
        end if

        call coeffViscoArrZ()

     else if(viscoLaw == "FK") then

        if(reference == "bottom") then
           coeffVisco = exp(thetaVisco * botTref)
           etaTop = exp(thetaVisco * (botTref - topTref) - gammaVisco)
           etaBot = 1.0
        else if (reference == "top") then
           coeffVisco = exp(thetaVisco * topTref + gammaVisco)
           etaTop = 1.0
           etaBot = exp(-thetaVisco * (botTref - topTref) + gammaVisco)
        else ! reference == "T_ref"
           coeffVisco = exp(thetaVisco * TrefAdim + gammaVisco)
           etaTop = exp(thetaVisco * (TrefAdim - topTref))
           etaBot = exp(thetaVisco * (TrefAdim - botTref) + gammaVisco)
        end if

        call coeffViscoFKZ()  ! computes the effect of pressure for FK

     end if ! viscolaw

     if(cutOffMax) then
        if(.not. smoothCutOff) then
           etaTop = min(viscoCutOffMax, etaTop)
        else
           etaTop = 1.0 / (1.0 / viscoCutOffMax + 1.0 / etaTop)
        end if
     end if
     
     if(cutOffMin) then
        etaBot = max(viscoCutOffMin, etaBot)
     end if
     
     call newtonianVisco(etaT, 3)  ! computed at cells' center (id=3)
     
  end if

  if(isoviscous) yielding = .false.
  
  ! yielding - - - - -
  if(yielding) then
     if(yieldFix < 0.0 .and. yieldSurf < 0.0 &
          .and. slopeYield < 0.0) then
        yielding = .false.
     end if
  end if

  if(yielding) then
     ! effectiveVisco: etayield=min(yieldFix, yieldSurf + slopeYield * (1-y))
     ! if param put to neg => infinity
     if(yieldFix < 0.0) then
        yieldFix = infinity
     end if
     if(slopeYield < 0.0) then
        slopeYield = infinity
     end if
     if(yieldSurf < 0.0) then
        yieldSurf = infinity
     end if
     if(etaStar < 0.0) then
        etaStar = 0.0
     end if
     if(compressionFactor > 0.0) then
        anderson = .true.
     end if
     
     write(*, '(a)') "---------------"
     write(*, '(a)') "Yielding with: "
     write(*,'(a,e10.2e3,a,2e10.3,a,e10.2)') &
          "yieldFix = ", yieldFix, &
          "     yieldSurf, slopeYield = ", yieldSurf, slopeYield, &
          "     etaStar = ", etaStar
     if(anderson) then
        write(*, '(a, f5.2)') "    compressionFactor = ", compressionFactor
     end if

     call yieldZ()

     ! for underrelaxations in the Picard iterations for yielding
     if(alphaVisco < 0.0 .or. alphaVisco > 1.0) then
        write(*, '(a)') "Error in the value of alphaVisco: &
             & must be between 0.0 and 1.0"
        write(*, '(a)') "-> set to alphaVisco = 0.5"
        alphaVisco = 0.5
     end if

     if(.not.etaRead) then
        call modifyEffVisco(one) ! computes eta
     end if

  end if ! yielding

  return
end subroutine initVisco
! =======================================================
subroutine updateBiHarmCoeff(whatEta)
  ! computes eta on cells' corners and coeff
  ! needed for the coeff in the solver of form Lap(Lap(sf)) + RHS = 0
  use general
  use grids
  implicit none
  character(len=10), intent(in):: whatEta
  integer:: level
  
  do level = 0, nLevel
     ! bcs should be done before, and whatEta(4) computed
     ! directly too, from effectiveVisco(id=4)
     call computeBiHarmCoeff(whatEta, level) ! -> computes cbh
     
     if(level < nLevel) then
        call restriction(whatEta, level+1, 3)
        call applyBcs(whatEta, level+1, 3)
        call centerToCorner(whatEta, level+1) ! from id=3 to id=4        
        call applyBcs(whatEta, level+1, 4)
     end if ! level

  end do ! level

  ! at coarsest level, for the matrix for exact solution:
  call computeCoarsestCoeff(cbh)
  
  return
end subroutine updateBiHarmCoeff
! =======================================================
subroutine modifyEffVisco(exponent)
  use general
  use grids
  implicit none
  real(rk), intent(in):: exponent
  character(len=10), parameter:: etanew=work4 
  
  call sf_to_sist(0, 3)

  if(exponent < 1.0) then
     call effectiveVisco(etaEff, etaT, 3)
     call coeffProduct(eta, etaEff, eta, 0, 3, exponent)
  else
     call effectiveVisco(eta, etaT, 3)
  end if
  
  call applyBcs(eta, 0, 3)
  call centerToCorner(eta, 0, geometric=viscoGeom)
  call applyBcs(eta, 0, 4)

  return
end subroutine modifyEffVisco
! =============================================
subroutine computeBiHarmCoeff(whatEta, level)
  ! whatEta can be etaT or etaEff
  use general
  use grids
  implicit none
  character(len=10), intent(in):: whatEta
  integer, intent(in):: level
  real(rk), pointer:: pc(:, :, :) => null()
  real(rk), pointer:: peta(:, :, :) => null()
  integer:: i, j
  real(rk):: invd(1:2), invr2, invr4, invdz2, invdx2, invdz4, invdx4
  real(rk), pointer:: pr(:) => null()
  real(rk), pointer:: pgtop(:) => null()
  real(rk), pointer:: pgbot(:) => null()
  real(rk), pointer:: pgmid(:) => null()
  real(rk), pointer:: pg1top(:) => null()
  real(rk), pointer:: pg1bot(:) => null()
  real(rk), pointer:: pgcentertop(:) => null()
  real(rk), pointer:: pgcenterbot(:) => null()
  real(rk), pointer:: pgcornertop(:) => null()
  real(rk), pointer:: pgcornerbot(:) => null()
  real(rk), pointer:: pg4(:) => null()
  real(rk), pointer:: pg5(:) => null()
  real(rk), pointer:: pg12(:) => null()
  real(rk), pointer:: pg13(:) => null()
  
  !-------- stencil ---------
  !            12
  !            |
  !       7----4----6
  !       |    |    |
  !  11---3----1----2---10
  !       |    |    |
  !       9----5----8
  !            |
  !            13
  !--------------------------  
  
  call associatePtrVec(peta, whatEta, level)
  call associatePtrVec(pc, cbh, level)

  call getInvD(invd, level)

  invdx2 = invd(1)**2  ! actually = 1/dphi^2 for curved Geometry (angle, not length)
  invdx4 = invd(1)**4
  invdz2 = invd(2)**2
  invdz4 = invd(2)**4

  ! ____________________________________________________________________
  if(curvedGeometry) then ! Spherical annulus or cylindrical geometry
  ! ____________________________________________________________________     
     call associatePtrLine(pr, r, level)

     call associatePtrLine(pgtop, gtop, level)
     call associatePtrLine(pgbot, gbot, level)
     call associatePtrLine(pgmid, gmid, level)
     call associatePtrLine(pg1top, g1top, level)
     call associatePtrLine(pg1bot, g1bot, level)
     call associatePtrLine(pgcentertop, gcentertop, level)
     call associatePtrLine(pgcenterbot, gcenterbot, level)
     call associatePtrLine(pgcornertop, gcornertop, level)
     call associatePtrLine(pgcornerbot, gcornerbot, level)
     call associatePtrLine(pg4, g4, level)
     call associatePtrLine(pg5, g5, level)
     call associatePtrLine(pg12, g12, level)
     call associatePtrLine(pg13, g13, level)

     if(isoviscous) then
        ! the other writing works too, but simpler here and computed only once
        do j=i0, iend(level, 2)

           invr2 = 1.0 / pr(j)**2
           invr4 = 1.0 / pr(j)**4

           do i=i0, iend(level, 1)

              pc(1, i, j) = 4.0 * invr2 * invdx2 * invdz2 * &
                   ( pgtop(j) + pgbot(j) - pgmid(j) ) + &
                   invdz4 * ( pgmid(j)**2 + pg1top(j) + pg1bot(j) ) + &
                   6.0 * invr4 * invdx4

              pc(2, i, j) = - 2.0 * invr2 * invdx2 * invdz2 * &
                   ( pgtop(j) + pgbot(j) - pgmid(j) ) &
                   - 4.0 * invr4 * invdx4 
              pc(3, i, j) = pc(2, i, j)

              pc(4, i, j) = -4.0 * invr2 * invdx2 * invdz2 * &
                   ( pgcentertop(j) - pgcornertop(j) ) - &
                   invdz4 * ( pgcornertop(j) * pgmid(j) + pg4(j) )
              pc(5, i, j) = -4.0 * invr2 * invdx2 * invdz2 * &
                   ( pgcenterbot(j) - pgcornerbot(j) ) - &
                   invdz4 * ( pgcornerbot(j) * pgmid(j) + pg5(j) )

              pc(6, i, j) = 2.0 * invr2 * invdx2 * invdz2 * &
                   ( pgcentertop(j) - pgcornertop(j) )
              pc(7, i, j) = pc(6, i, j)

              pc(8, i, j) = 2.0 * invr2 * invdx2 * invdz2 * &
                   ( pgcenterbot(j) - pgcornerbot(j) )
              pc(9, i, j) = pc(8, i, j)

              pc(10, i, j) = invr4 * invdx4
              pc(11, i, j) = pc(10, i, j)

              pc(12, i, j) = invdz4 * pg12(j)
              pc(13, i, j) = invdz4 * pg13(j)

           end do  ! i 
        end do ! j

     else ! not isoviscous = = = = = = = = = = = = = = =

        do j=i0, iend(level, 2)

           invr2 = 1.0 / pr(j)**2
           invr4 = 1.0 / pr(j)**4

           do i=i0, iend(level, 1)

              pc(1, i, j) = 2.0 * invr2 * invdx2 * invdz2 * ( &
                   (peta(i, j, 3) + peta(i-1, j, 3)) * pgtop(j) + &
                   (peta(i, j-1, 3) + peta(i-1, j-1, 3)) * pgbot(j) - &
                   2.0 * peta(i, j, 4) * pgmid(j) &
                   ) + invdz4 * ( &
                   peta(i, j, 4) * pgmid(j)**2 &
                   + peta(i, j+1, 4) * pg1top(j) + peta(i, j-1, 4) * pg1bot(j) ) &
                   + invr4 * invdx4 * &
                   ( 4.0 * peta(i, j, 4) + peta(i+1, j, 4) + peta(i-1, j, 4) )

              pc(2, i, j) = - invr2 * invdx2 * invdz2 * ( &
                   2.0 * peta(i, j, 3) * pgtop(j) + &
                   2.0 * peta(i, j-1, 3) * pgbot(j) &
                   - (peta(i, j, 4) + peta(i+1, j, 4)) * pgmid(j) ) &
                   - 2.0 * invr4 * invdx4 * (peta(i, j, 4) + peta(i+1, j, 4))
              pc(3, i, j) = - invr2 * invdx2 * invdz2 * ( &
                   2.0 * peta(i-1, j, 3) * pgtop(j) + &
                   2.0 * peta(i-1, j-1, 3) * pgbot(j) &
                   - (peta(i, j, 4) + peta(i-1, j, 4)) * pgmid(j) ) &
                   - 2.0 * invr4 * invdx4 * (peta(i, j, 4) + peta(i-1, j, 4))

              pc(4, i, j) = - 2.0 * invr2 * invdx2 * invdz2 * ( &
                   (peta(i, j, 3) + peta(i-1, j, 3)) * pgcentertop(j) - &
                   (peta(i, j, 4) + peta(i, j+1, 4)) * pgcornertop(j) ) &
                   - invdz4 * ( &
                   peta(i, j, 4) * pgcornertop(j) * pgmid(j) + &
                   peta(i, j+1, 4) * pg4(j) )
              pc(5, i, j) = - 2.0 * invr2 * invdx2 * invdz2 * ( &
                   (peta(i, j-1, 3) + peta(i-1, j-1, 3)) * pgcenterbot(j) - &
                   (peta(i, j, 4) + peta(i, j-1, 4)) * pgcornerbot(j) ) &
                   - invdz4 * ( &
                   peta(i, j, 4) * pgcornerbot(j) * pgmid(j) + &
                   peta(i, j-1, 4) * pg5(j) )

              pc(6, i, j) = invr2 * invdx2 * invdz2 * ( &
                   2.0 * peta(i, j, 3) * pgcentertop(j) - &
                   (peta(i+1, j, 4) + peta(i, j+1, 4)) * pgcornertop(j) )
              pc(7, i, j) = invr2 * invdx2 * invdz2 * ( &
                   2.0 * peta(i-1, j, 3) * pgcentertop(j) - &
                   (peta(i-1, j, 4) + peta(i, j+1, 4)) * pgcornertop(j) )

              pc(8, i, j) = invr2 * invdx2 * invdz2 * ( &
                   2.0 * peta(i, j-1, 3) * pgcenterbot(j) - &
                   (peta(i+1, j, 4) + peta(i, j-1, 4)) * pgcornerbot(j) )
              pc(9, i, j) = invr2 * invdx2 * invdz2 * ( &
                   2.0 * peta(i-1, j-1, 3) * pgcenterbot(j) - &
                   (peta(i-1, j, 4) + peta(i, j-1, 4)) * pgcornerbot(j) )

              pc(10, i, j) = invr4 * invdx4 * peta(i+1, j, 4)
              pc(11, i, j) = invr4 * invdx4 * peta(i-1, j, 4)

              pc(12, i, j) = invdz4 * pg12(j) * peta(i, j+1, 4)
              pc(13, i, j) = invdz4 * pg13(j) * peta(i, j-1, 4)
           end do ! i
        end do ! j

        nullify(pr)
        nullify(pgtop, pgbot, pgmid, pg1top, pg1bot)
        nullify(pgcentertop, pgcenterbot, pgcornertop, pgcornerbot)
        nullify(pg4, pg5, pg12, pg13)

     end if ! isoviscous or not
  
     !______________________________________________________________________
  else  ! Cartesian geometry 
     !______________________________________________________________________

     if(isoviscous) then

        do j=i0, iend(level,2)
           do i=i0, iend(level, 1)

              pc(1, i, j) = 6.0 * (invdx4 + invdz4) + 8.0 * invdx2 * invdz2

              pc(2, i, j) = 4.0 * invdx2 * (invdz2 - invdx2) - 8.0 * invdx2 * invdz2
              pc(3, i, j) = pc(2, i, j)
              pc(4, i, j) = 4.0 * invdz2 * (invdx2 - invdz2) - 8.0 * invdx2 * invdz2
              pc(5, i, j) = pc(4, i, j)

              pc(6, i, j) = 2.0 * invdx2 * invdz2
              pc(7, i, j) = pc(6, i, j)
              pc(8, i, j) = pc(6, i, j)
              pc(9, i, j) = pc(6, i, j)

              pc(10, i, j) = invdx4
              pc(11, i, j) = invdx4
              pc(12, i, j) = invdz4
              pc(13, i, j) = invdz4

           end do ! i
        end do ! j

     else ! not isoviscous = = = = = = = = = = = = = = =

        do j=i0, iend(level,2)
           do i=i0, iend(level, 1)

              pc(1, i, j) = &
                   4.0 * (invdx4 + invdz4 - 2.0 * invdx2 * invdz2) * peta(i, j, 4) + &
                   invdx4 * (peta(i+1, j, 4) + peta(i-1, j, 4)) + &
                   invdz4 * (peta(i, j+1, 4) + peta(i, j-1, 4)) + &
                   4.0 * invdx2 * invdz2 * &
                   (peta(i, j, 3) + peta(i-1, j, 3) + peta(i, j-1, 3) + peta(i-1, j-1, 3))

              pc(2, i, j) = &
                   2.0 * invdx2 * (invdz2 - invdx2) * (peta(i, j, 4) + peta(i+1, j, 4)) - &
                   4.0 * invdx2 * invdz2 * (peta(i, j, 3) + peta(i, j-1, 3))
              pc(3, i, j) = &
                   2.0 * invdx2 * (invdz2 - invdx2) * (peta(i-1, j, 4) + peta(i, j, 4)) - &
                   4.0 * invdx2 * invdz2 * (peta(i-1, j, 3) + peta(i-1, j-1, 3))

              pc(4, i, j) = &
                   2.0 * invdz2 * (invdx2 - invdz2) * (peta(i, j, 4) + peta(i, j+1, 4)) - &
                   4.0 * invdx2 * invdz2 * (peta(i, j, 3) + peta(i-1, j, 3))
              pc(5, i, j) = &
                   2.0 * invdz2 * (invdx2 - invdz2) * (peta(i, j-1, 4) + peta(i, j, 4)) - &
                   4.0 * invdx2 * invdz2 * (peta(i, j-1, 3) + peta(i-1, j-1, 3))
              
              pc(6, i, j) = invdx2 * invdz2 * &
                   (4.0 * peta(i, j, 3) - peta(i, j+1, 4) - peta(i+1, j, 4))
              pc(7, i, j) = invdx2 * invdz2 * &
                   (4.0 * peta(i-1, j, 3) - peta(i, j+1, 4) - peta(i-1, j, 4))
              pc(8, i, j) = invdx2 * invdz2 * &
                   (4.0 * peta(i, j-1, 3) - peta(i, j-1, 4) - peta(i+1, j, 4))
              pc(9, i, j) = invdx2 * invdz2 * &
                   (4.0 * peta(i-1, j-1, 3) - peta(i, j-1, 4) - peta(i-1, j, 4))

              pc(10, i, j) = invdx4 * peta(i+1, j, 4)
              pc(11, i, j) = invdx4 * peta(i-1, j, 4)
              pc(12, i, j) = invdz4 * peta(i, j+1, 4)
              pc(13, i, j) = invdz4 * peta(i, j-1, 4)
              
           end do ! i
        end do ! j

     end if ! isoviscous or not
     !______________________________________________________________
  end if ! curvedGeometry or not
  !_________________________________________________________________

  nullify(pc, peta)

  return
end subroutine computeBiHarmCoeff
! =======================================================
subroutine computeCoarsestCoeff(whatCoeff)
  use general
  use grids
  implicit none
  character(len=10), intent(in):: whatCoeff  ! type of coeff (cbh, cbhy or cbhPrec)
  real(rk), pointer:: pc(:, :, :) => null()
  integer:: i, nval, iez

  iez = iend(nLevel, 2)
  call associatePtrVec(pc, whatCoeff, nLevel)

  do i=i0, iend(nLevel, 1)
     coarsestGamma(i) = pc(1, i, iez)
     if(topV == 'dirichlet') then
        coarsestGamma(i) = coarsestGamma(i) + pc(12, i, iez)
     else
        coarsestGamma(i) = coarsestGamma(i) - pc(12, i, iez)
     end if
     if(botV == 'dirichlet') then
        coarsestGamma(i) = coarsestGamma(i) + pc(13, i, iez)
     else
        coarsestGamma(i) = coarsestGamma(i) - pc(13, i, iez)
     end if     
  end do

  if(coarsestMode == 1) then
     ! at coarsest level: nphi=nr=2 and side not periodic -> only sf(i0+1) to compute
     mcc(1, 1) = kleft * pc(11, i0+1, iez) + coarsestGamma(i0+1) &
          + kright * pc(10, i0+1, iez)
  else if(coarsestMode == 2) then
     ! nphi=nr=2 and side='periodic'
     mcc(1, 1) = pc(10, i0, iez) + coarsestGamma(i0) + pc(11, i0, iez)
     mcc(1, 2) = pc(2, i0, iez) + pc(3, i0, iez)
     mcc(2, 1) = pc(2, i0+1, iez) + pc(3, i0+1, iez)
     mcc(2, 2) = pc(10, i0+1, iez) + coarsestGamma(i0+1) + pc(11, i0+1, iez)
  else if(coarsestMode == 3) then
     !nr=2, nphi=4 and side not periodic
     mcc(1, 1) = kleft * pc(11, i0+1, iez) + coarsestGamma(i0+1)
     mcc(1, 2) = pc(2, i0+1, iez)
     mcc(1, 3) = pc(10, i0+1, iez)
     !
     mcc(2, 1) = pc(3, i0+2, iez)
     mcc(2, 2) = coarsestGamma(i0+2)
     mcc(2, 3) = pc(2, i0+2, iez)
     !
     mcc(3, 1) = pc(11, i0+3, iez)
     mcc(3, 2) = pc(3, i0+3, iez)
     mcc(3, 3) = coarsestGamma(i0+3) + kright * pc(10, i0+3, iez)
  else if(coarsestMode == 4) then
     !nr=2, nphi=4 and side='periodic'
     mcc(1, 1) = coarsestGamma(i0)
     mcc(1, 2) = pc(2, i0, iez)
     mcc(1, 3) = pc(10, i0, iez) + pc(11, i0, iez)
     mcc(1, 4) = pc(3, i0, iez)
     !
     mcc(2, 1) = pc(3, i0+1, iez)
     mcc(2, 2) = coarsestGamma(i0+1)
     mcc(2, 3) = pc(2, i0+1, iez)
     mcc(2, 4) = pc(10, i0+1, iez) + pc(11, i0+1, iez)
     !
     mcc(3, 1) = pc(10, i0+2, iez) + pc(11, i0+2, iez)
     mcc(3, 2) = pc(3, i0+2, iez)
     mcc(3, 3) = coarsestGamma(i0+2)
     mcc(3, 4) = pc(2, i0+2, iez)
     !
     mcc(4, 1) = pc(2, i0+3, iez)
     mcc(4, 2) = pc(10, i0+3, iez) + pc(11, i0+3, iez)
     mcc(4, 3) = pc(3, i0+3, iez)
     mcc(4, 4) = coarsestGamma(i0+3)
  else
     ! mcc are no more directly the elements of the matrix
     ! but the 5 elements of a line in the (nearly)pentadiagonal matrix
     ! [c d e ..... a b]
     ! [b c d e ..... a]
     ! [a b c d e ....0]
     ! [0 a b c d e...0]
     ! :
     ! coarsest(:,1)=a, coarsest(:,2)=b, coarsest(:,3)=c, coarsest(:,4)=d, coarsest(:,5)=e 

     nval = coarsestMode ! for clarity
     mcc(:, :) = 0.0
     
     if(mod(nval, 2) == 0.0) then
        ! side = 'periodic' -> nearly pentadiagonal
        do i=1, nval
           mcc(i, 1) = pc(11, i0+i-1, iez)
           mcc(i, 2) = pc(3, i0+i-1, iez)
           mcc(i, 3) = coarsestGamma(i0+i-1)
           mcc(i, 4) = pc(2, i0+i-1, iez)
           mcc(i, 5) = pc(10, i0+i-1, iez)
        end do
     else
        ! side not periodic -> perfectly pentadiagonal
        mcc(1, 3) = kleft * pc(11, i0+1, iez) + coarsestGamma(i0+1)
        mcc(1, 4) = pc(2, i0+1, iez)
        mcc(1, 5) = pc(10, i0+1, iez)
        !
        mcc(2, 2) = pc(3, i0+2, iez)
        mcc(2, 3) = coarsestGamma(i0+2)
        mcc(2, 4) = pc(2, i0+2, iez)
        mcc(2, 5) = pc(10, i0+2, iez)
        !
        do i=3, nval - 2
           mcc(i, 1) = pc(11, i0+i, iez)
           mcc(i, 2) = pc(3, i0+i, iez)
           mcc(i, 3) = coarsestGamma(i0+i)
           mcc(i, 4) = pc(2, i0+i, iez)
           mcc(i, 5) = pc(10, i0+i, iez)
        end do
        !
        mcc(nval-1, 1) = pc(11, i0+nval-1, iez)
        mcc(nval-1, 2) = pc(3, i0+nval-1, iez)
        mcc(nval-1, 3) = coarsestGamma(i0+nval-1)
        mcc(nval-1, 4) = pc(2, i0+nval-1, iez)
        ! 
        mcc(nval, 1) = pc(11, i0+nval, iez)
        mcc(nval, 2) = pc(3, i0+nval, iez)
        mcc(nval, 3) = coarsestGamma(i0+nval) + kright * pc(10, i0+nval, iez)
     end if
  endif

  nullify(pc)

  return
end subroutine computeCoarsestCoeff
! =======================================================
subroutine changeEta(rmsChange, maxChange, whatEta, whatEtaPrec, id)
  ! relative change of the viscosity
  use general
  use grids
  implicit none
  character(len=10), intent(in):: whatEta, whatEtaPrec
  integer, intent(in):: id
  real(rk), intent(out):: rmsChange, maxChange
  real(rk), pointer:: peta(:,:) => null()
  real(rk), pointer:: petaPrec(:,:) => null()
  real(rk), pointer:: pdiff(:,:) => null()
  real(rk):: l0, l1
  integer:: i, j
 
  call associatePtrScal(peta, whatEta, 0, 3)
  call associatePtrScal(petaPrec, whatEtaPrec, 0, 3)
  call associatePtrScal(pdiff, work1, 0, 3) ! work1: tempo

  do j=i0, iend(0, 2)
     do i=i0, iend(0, 1)
        l0 = log(petaPrec(i, j))
        l1 = log(peta(i, j))
        pdiff(i, j) = 2.0 * abs(l1 - l0) / (l1 + l0)
     end do
  end do
  nullify(peta, petaPrec, pdiff)
  
  call getRmsSimple(rmsChange, work1, 0, id)
  call getMax(maxChange, work1, 0, id, absVal=.true.)
    
  return
end subroutine changeEta
! =======================================================
