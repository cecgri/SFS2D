!=======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
!=======================================================
program SFS2D

  use general
  use grids
  use io
  implicit none

  ! initial readings --------
  call defaultValues()
  call readParameters()
  call check()

  ! init geometry and grid
  call setupGeometry()
  call initGrids()
  call geometricalCoefficients()

 ! prepare output files, check existing ones 
 ! and init temperature, velocity, etc.
  call initValues()
  call initLogs()
  call initMaps()

  call startRun()
  
  do while(running)
     
     call cpu_time(t_allStart)

     call solveVelocities()
     
     call writings()
     
     call nexttimestep() 
     
     call cpu_time(t_allEnd)
     cpuTimeAll = cpuTimeAll + t_allEnd - t_allStart

  end do ! running
  
  ! finalize ----------------
  call endRun()

end program SFS2D
!=======================================================
