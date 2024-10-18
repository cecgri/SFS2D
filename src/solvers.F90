! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine solveVelocities()
  use general
  use grids
  use io
  use relax
  implicit none
  logical:: solver_converge, converge, almostConverge
  integer, parameter:: nloopmax=100
  integer:: iloop
  real(rk):: rmsRes, rmsSf, etamin
  real(rk):: rmsCorrSf, rmsChange, tolerance
  character(len=20), parameter:: methodYield="newton"

  if(continueRun .and. counter == 0) then
     ! sf is already read. Just computes vel:
     call applyBcs(sf, 0, 4)
     call sf_to_vel(sf, vel, 0)
     return
  end if
  
  call cpu_time(t_start)

  ! DEFINE RHS =============================
  call zeroingScal(rhs, 0, 4)
  
  ! Temp on z-faces
  if(.not.uniformTemp) then
     call centerToFace(temp, 0, 2)
     call applyBcs(temp, 0, 2)
     call op_dx(rhs, 4, temp, 2, 1, 0, mult=rayleigh)
  end if
  
  ! COMPUTE COEFFICIENTS FOR BIHARMONIC EQUATION ===============
  if(.not.isoviscous) then
     call newtonianVisco(etaT, 3)  ! computed at cells' centers (id=3)
     
     if(yielding) then
        call modifyEffVisco(alphaVisco)  ! does id=3 and id=4 ! one or alphaVisco?
        call updateBiHarmCoeff(eta)       
     else
        call centerToCorner(etaT, 0, geometric=viscoGeom)
        call applyBcs(etaT, 0, 4)
        call updateBiHarmCoeff(etaT)
     end if 
  end if
     
  if(yielding) then
     call copyScal(sfPrec, 4, sf, 4, 0)
     call solver(relax_biharmonic, sf, rhs, res, 4, errorRes, &
          solver_converge)

     iloop = 0

     if(methodYield == "picard") then 
        !===================!
        ! PICARD ITERATIONS !
        !===================!
        do iloop = 0, nloopmax
           call modifyEffVisco(alphaVisco)

           ! residue: computed with etaEff if alphaVisco < 1 
           if(alphaVisco < 1.0) then
              call updateBiHarmCoeff(etaEff)
           else
              call updateBiHarmCoeff(eta)
           end if
           call residue_biharmonic(rmsRes, sf, rhs, res, 4, do_relax=.true.)

           call getRmsChange(rmsChange, sf, sfPrec, 0, 4)
           call getRmsSimple(rmsSf, sf, 0, 4)
           !! TEMPO
           if(rmsSf > 1.0E12) then
              print*, 'rmsSf = ', rmsSf, ' => too high'
              call endRun()
           end if

           converge = (rmsRes < errorRes * rmsSf .and. &
                rmsChange < errorChange * rmsSf)

 !          almostConverge = (rmsRes < 5.0 * errorRes * rmsSf .and. iloop > 10)
           almostConverge = .false.
           
           if(verbose >= 3 .or. (verbose >= 2 .and. converge)) then
              call getMin(etamin, eta, 0, 4)
              write(*,'(a, 2i4, 4e14.4, 2l4)') &
                   "       yield picard: ", iloop, nrelax, &
                   rmsSf, rmsRes / rmsSf, rmsChange / rmsSf,  &
                   etamin, converge, almostConverge
           end if
           if(converge) exit

           ! for solver: using eta and not etaEff if alphaVisco < 1
           if(alphaVisco < 1.0) then
              call updateBiHarmCoeff(eta)
           end if
           
           call copyScal(sfPrec, 4, sf, 4, 0)
           if(.not. almostConverge) then
              call solver(relax_biharmonic, sf, rhs, res, 4, &
                   errorRes, solver_converge)
           else ! almostConverge
              call solver(relax_biharmonic, sf, rhs, res, 4, &
                   errorRes, solver_converge, levelMax=3)              
           end if

        end do ! iloop

        if(.not.converge) then
           write(*, '(a, 3e13.3, i7)') "WARNING: yielding not converging", &
                rmsSf, rmsRes / rmsSf, rmsChange, itime
        end if

     else if(methodYield == "newton") then
        !===============!
        ! NEWTON METHOD !
        !===============!
        call modifyEffVisco(one)  ! do one modification of eta first
        call updateBiHarmCoeff(eta) ! computes the new coeff. cbh
        call solver(relax_biharmonic, sf, rhs, res, 4, &
             errorRes, solver_converge)

        converge = .false.
        
        call zeroingScal(corrSf, 0, 4) ! check for position
        call getRmsSimple(rmsSf, sf, 0, 4)
        
        do iloop=0, nloopmax
           call modifyEffVisco(one)
           call updateBiHarmCoeff(eta)
           call residue_biharmonic(rmsRes, sf, rhs, res, 4, do_relax=.true.)
           
           call copyScal(corrRhs, 4, res, 4, 0)
           call zeroingScal(corrSf, 0, 4) ! check for position
           call solver(relax_biharmonic, corrSf, corrSf, res, 4, &
                errorRes, solver_converge, givenRms=rmsSf)
           
           call getRmsSimple(rmsSf, sf, 0, 4)
           call getRmsSimple(rmsCorrSf, corrSf, 0, 4)
           
           converge = (rmsRes < errorRes * rmsSf) .and. &
                (rmsCorrSf < errorChange * rmsSf)
              
           if(verbose >= 3 .or. (verbose >= 2 .and. converge)) then
              call getMin(etamin, eta, 0, 4)
              write(*,'(a, 2i4, 5e14.4, l4)') &
                   "       yield newton: ", iloop, nrelax, &
                   rmsSf, rmsCorrSf, rmsRes / rmsSf, rmsCorrSf / rmsSf, &
                   etamin, converge
           end if
           if(converge) exit
           
           call addScal(sf, corrSf, 0, 4, mult=alphaCorrSf)
           call applyBcs(sf, 0, 4)
        end do ! iloop=0, nloopmax

        if(.not.converge) then
           write(*, '(a, 3e13.3, i7)') "WARNING: yielding not converging", &
                rmsSf, rmsRes / rmsSf, rmsCorrSf / rmsSf, itime
        end if

     else
        ! =============================================
        ! mixed method: Picard until certain error and
        ! then Newton approximate
        ! =============================================
        ! -> does not seem to work
        tolerance = 10.0
        call modifyEffVisco(alphaVisco) ! computes eta and etaEff
        call updateBiHarmCoeff(eta)

        ! PICARD = = = = = = = =
        do iloop = 0, nloopmax
           call residue_biharmonic(rmsRes, sf, rhs, res, 4, do_relax=.true.)
           !call getRmsChange(rmsChange, sf, sfPrec, 0, 4)
           ! TEST
           call getRmsChange(rmsChange, sf, sfPrec, 0, 4, relative=.true.)
           call getRmsSimple(rmsSf, sf, 0, 4)

           converge = (rmsRes < tolerance * errorRes * rmsSf .and. &
                rmsChange < tolerance * errorChange * rmsSf)
           
           if(verbose >= 3 .or. (verbose >= 2 .and. converge)) then
              call getMin(etamin, eta, 0, 4)
              write(*,'(a, 2i4, 4e14.4, l4)') &
                   "       yield picard: ", iloop, nrelax, &
                   rmsSf, rmsRes / rmsSf, rmsChange,  &
                   etamin, converge
           end if
           if(converge) exit
           
           call copyScal(sfPrec, 4, sf, 4, 0)
           call solver(relax_biharmonic, sf, rhs, res, 4, &
                errorRes, solver_converge)
           
           call modifyEffVisco(alphaVisco) ! computes eta and etaEff           
           call updateBiHarmCoeff(eta)
        end do ! iloop
        
        ! NEWTON = = = = = = = =
        converge = .false.
        call zeroingScal(corrSf, 0, 4)

        do iloop=0, nloopmax
           call modifyEffVisco(one)
           call updateBiHarmCoeff(eta)
           call residue_biharmonic(rmsRes, sf, rhs, res, 4, do_relax=.true.)
           
           call copyScal(corrRhs, 4, res, 4, 0)
           call zeroingScal(corrSf, 0, 4) ! check for position
           call solver(relax_biharmonic, corrSf, corrRhs, res, 4, &
                errorRes, solver_converge, givenRms=rmsSf)
           
           call getRmsSimple(rmsSf, sf, 0, 4)
           call getRmsSimple(rmsCorrSf, corrSf, 0, 4)
           
           converge = (rmsRes < errorRes * rmsSf) .and. &
                (rmsCorrSf < errorChange * rmsSf)
              
           if(verbose >= 3 .or. (verbose >= 2 .and. converge)) then
              call getMin(etamin, eta, 0, 4)
              write(*,'(a, 2i4, 5e14.4, l4)') &
                   "       yield newton: ", iloop, nrelax, &
                   rmsSf, rmsCorrSf, rmsRes / rmsSf, rmsCorrSf / rmsSf, &
                   etamin, converge
           end if
           if(converge) exit
           
           call addScal(sf, corrSf, 0, 4, mult=alphaCorrSf)
           call applyBcs(sf, 0, 4)
        end do ! iloop=0, nloopmax

        if(.not.converge) then
           write(*, '(a, 3e13.3, i7)') "WARNING: yielding not converging", &
                rmsSf, rmsRes / rmsSf, rmsCorrSf / rmsSf, itime
        end if
        
     end if  ! method Picard or approximate Newton
        
  else ! not yielding
     
     call solver(relax_biharmonic, sf, rhs, res, 4, errorRes, &
          solver_converge)
     
  end if
  
  call applyBcs(sf, 0, 4)
  call sf_to_vel(sf, vel, 0)  ! includes bcs
  
  call cpu_time(t_end)
  cpuTimeSolver = cpuTimeSolver + t_end - t_start
  
  return
end subroutine solveVelocities
! =======================================================
