! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine startRun()
  use general
  implicit none

  running = .true.
  write(*,*) " ======================================================="
  write(*,*) "|                    starting run                       | "
  write(*,*) " ======================================================="

  call initCpuTime() 

end subroutine startRun
! =======================================================
subroutine endRun()
  use grids
  use io
  implicit none

  call freeGrids()
  call clearLogs()

  write(*,*) " ======================================================="
  write(*,*) "|                      ending run                       | "
  write(*,*) " ======================================================="

  stop
end subroutine endRun
! =======================================================
subroutine initCpuTime()
  use general
  use io
  implicit none

  if(continueRun) then
     cpuTimeAll = logDiag%data(2)
     cpuTimeSolver = logDiag%data(3)
     cpuTimeTemperature = logDiag%data(4)
     cpuTimeOutput = logDiag%data(5)
     nLoopTot = int(logDiag%data(6))
  else
     cpuTimeAll = 0.0
     cpuTimeSolver = 0.0
     cpuTimeTemperature = 0.0
     cpuTimeOutput = 0.0
     nLoopTot = 0
  end if ! continueRun

  return
end subroutine initCpuTime
! =======================================================
subroutine nexttimestep()
  use general
  implicit none

  call computeTimestep()
  
  call heatConservation()

  time = time + timeStep
  itime = itime + 1 
  counter = counter + 1

  if( (endType == "steps" .and. (counter > endSteps)) .or. &
      (endType == "time" .and. (time > endTime + 10.0*timeStep)) ) then
     running = .false.
  end if

  return
end subroutine nexttimestep
! ==================================================
subroutine computeTimestep()
  use general
  use grids
  implicit none
  real(rk):: dt
  
  call courant(dt)
  if(.not. uniformTemp) then
     timestep = min(dt, minTimestep)
  else
     timestep = dt
  end if

  return
end subroutine computeTimestep
! =======================================================
subroutine computeMinTimestep(mindt)
  use general
  implicit none
  real(rk), intent(out):: mindt
  real(rk):: dmin, dxx

  dxx = dx ! Cartesian geometry
  if(curvedGeometry) dxx = dx * rbot ! dx is an angle

  dmin = min(dz, dxx)
  mindt = fourth * dmin * dmin

  return
end subroutine computeMinTimestep
! =======================================================
subroutine courant(dt)
  use general
  use grids
  implicit none
  real(rk), intent(out):: dt
  real(rk):: maxV(1:2), dmin, dxx
  integer:: id

  do id=1, 2
     call getMax(maxV(id), vel, 0, id, absVal=.true.)
  end do

  dxx = dx ! Cartesian geometry
  if(curvedGeometry) dxx = dx * rbot ! dx is an angle

  dmin = min(dz, dxx)
  dt = half * dmin / maxval(maxV) * timeFactor
  
  return
end subroutine courant
! =======================================================
subroutine heatConservation()
  use general
  use grids
  implicit none

  if(uniformTemp) return ! special case (Rayleigh-Taylor)
  
  call cpu_time(t_start)

  ! choice of method for advection-diffusion --------
  if(heat_method == "mpdata") then
     
     call mpdata(temp)
     
  else if(heat_method == "ADI") then
     
     call ADI() ! advection/diffusion of aglo using vglo
     call correct_negativeTemp() ! -> does not seem needed
     call applyBcs(temp, 0, 3)
     
  end if

  call cpu_time(t_end)
  cpuTimeTemperature = cpuTimeTemperature + t_end - t_start

  return
end subroutine heatConservation
! =======================================================
subroutine correct_negativeTemp()
  use grids
  implicit none
  real(rk), pointer:: pt(:,:) => null()
  integer:: i, j

  call associatePtrScal(pt, temp, 0, 3)
  do j=i0, iend(0, 2)
     do i=i0, iend(0, 1)
        if(pt(i,j) < topValueT) then
           print*, 'neg: ', i, j, pt(i,j)
           pt(i,j) = topValueT
        end if
     end do ! i
  end do ! j

  nullify(pt)

  return
end subroutine correct_negativeTemp
! =======================================================
