! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine check()
  ! check different read parameters
  implicit none

  ! heat_method - - - - -
  call checkHeatMethod()

  ! viscosity - - -
  call checkViscosity()

  ! bcs - - - - - - 
  call checkBoundaryConditions()

  ! writings - - - - -
  call checkWritings()
  
  return
end subroutine check
! =======================================================
subroutine checkHeatMethod()
  use general
  implicit none

  if(uniformTemp) return
  
  if(.not.(heat_method == "mpdata" .or. &
       heat_method == "ADI")) then
     write(*,*) "ERROR: heat_method is not recognized"
     call endRun()
  end if

  if(heat_method == "ADI") then
     if(.not. (TVDscheme_temp == "upwind" .or. &
          TVDscheme_temp == "centered" .or. &
          TVDscheme_temp == "minmod" .or. &
          TVDscheme_temp == "superbee" .or. &
          TVDscheme_temp == "koren") ) then
        write(*,*) "ERROR: TVDscheme_temp is not recognized"
        call endRun()
     end if
  end if
  if(heat_method == "mpdata" .and. timeFactor > 1.0) then
     write(*,*) "WARNING: &
          & timeFactor must be equal or less than one for mpdata -> changed"
     timeFactor = 1.0
  end if
  
  return
end subroutine checkHeatMethod
! =======================================================
subroutine checkViscosity()
  use general
  implicit none

  if(viscoLaw == "isoviscous" .or. viscoLaw == "") then
     isoviscous = .true.
     return
  end if

  if(uniformTemp) return
  
  if(.not.(viscoLaw == "FK" .or. viscoLaw == "Arrhenius" &
       .or. viscoLaw == "FKC")) then
     write(*,*) "ERROR: viscoLaw is not recognized"
     call endRun()
  end if

  if(.not.(methodEtaEff == "harmonic" .or. methodEtaEff == "min" .or. &
       methodEtaEff == "minharmonic")) then
     write(*,*) "ERROR: methodEtaEff is not recognized"
     call endRun()
  end if

  if(.not.(reference == "bottom" .or. reference == "top" .or. &
       reference == "T_ref")) then
     write(*,*) "ERROR: reference is not recognized"
     call endRun()
  end if
  
  return
end subroutine checkViscosity
! =======================================================
subroutine checkBoundaryConditions()
  ! checks that boundary conditions given in param.in are okay
  use general
  implicit none

  !! CHECK PERIODICITY ==================================
  if((leftT == "periodic" .and. .not. rightT == "periodic") .or. &
       (.not. leftT == "periodic" .and. rightT == "periodic")) then
     write(*,*) &
          "ERROR: if leftT or rightT is periodic, both have to be"
     call endRun()
  end if

  if((leftV == "periodic" .and. .not. rightV == "periodic") .or. &
       (.not. leftV == "periodic" .and. rightV == "periodic")) then
     write(*,*) &
          "ERROR: if leftV or rightV is periodic, both have to be"
     call endRun()
  end if

  !! TEMPERATURE ========================================
  if(.not. (topT=="dirichlet" .or. topT=="neumann")) then
     write(*,*) "ERROR in boundary conditions: topT not recognized"
     call endRun()
  end if

  if(.not. (botT=="dirichlet" .or. botT=="neumann")) then
     write(*,*) "ERROR in boundary conditions: botT not recognized"
     call endRun()
  end if

  if(.not. (leftT=="dirichlet" .or. leftT=="neumann" .or. leftT=="periodic")) then
     write(*,*) "ERROR in boundary conditions: leftT not recognized"
     call endRun()
  end if
  if(leftT=="neumann" .and. leftValueT/=0.0 .and. heat_method=="ADI") then
     write(*,*) "ERROR in boundary conditions: leftT=neumann only &
          & implemented with leftValueT=0 (zero gradient) for ADI"
     call endRun()
  end if

  if(.not. (rightT=="dirichlet" .or. rightT=="neumann" .or. &
       rightT=="periodic")) then
     write(*,*) "ERROR in boundary conditions: leftT not recognized"
     call endRun()
  end if
  if(rightT=="neumann" .and. rightValueT/=0.0 .and. heat_method=="ADI") then
     write(*,*) "ERROR in boundary conditions: rightT=neumann only & 
          & implemented with rightValueT=0 (zero gradient) for ADI"
     call endRun()
  end if

  !! VELOCITY ============================================
  if(.not. (topV=="freeslip" .or. topV=="dirichlet" .or. topV=="file")) then
     write(*,*) "ERROR in boundary conditions: topV not recognized"
     call endRun()
  end if

  if(.not. (botV=="freeslip" .or. botV=="dirichlet") ) then
     write(*,*) "ERROR in boundary conditions: botV not recognized"
     call endRun()
  end if

  if(.not. (leftV=="dirichlet" .or. leftV=="neumann" .or. &
       leftV=="periodic")) then
     write(*,*) "ERROR in boundary conditions: leftV not recognized"
     call endRun()
  end if

  if(.not. (rightV=="dirichlet" .or. rightV=="neumann" .or. &
       rightV=="periodic")) then
     write(*,*) "ERROR in boundary conditions: rightV not recognized"
     call endRun()
  end if

  if(leftV == "neumann" .and. leftValueV /= 0.0) then
     write(*,*) &
          "WARNING: leftV=neumann is actually freeslip -> leftValueV is set to zero"
     leftValueV = 0.0
  end if
  
  if(rightV == "neumann" .and. rightValueV /= 0.0) then
     write(*,*) &
          "WARNING: rightV=neumann is actually freeslip -> rightValueV is set to zero"
     rightValueV = 0.0
  end if

  return
end subroutine checkBoundaryConditions
! =======================================================
subroutine checkWritings()
  use general
  implicit none

  if(.not. (writeType == "time" .or. writeType == "steps")) then
     write(*,*) &
          "ERROR in writing settings: writeTimeType not recognized"
     call endRun()
  end if

  if(.not. (endType == "time" .or. endType == "steps")) then
     write(*,*) &
          "ERROR in writing settings: endType not recognized"
     call endRun()     
  end if
  
  ! check that logs are not written less often than full files
  if(freqLogs > freqFiles) then
     write(*,*) "WARNING: freqFiles < freqLogs & 
          & -> setting both to freqLogs"
     freqFiles = freqLogs
  end if
  
  ! check that maps are not written more often than logs
  if(freqLogs > freqMaps) then
     write(*,*) "WARNING: freqMaps < freqLogs & 
          & -> setting both to freqLogs"
     freqMaps = freqLogs     
  end if
  
  ! check that numbers are > 1 if writeType or endType = 'steps'
  if(writeType == 'steps') then
     if(freqLogs < 1) then
        write(*,*) "ERROR: freqLogs is too small for writeType == 'steps'"
        call endRun()
     end if
  end if
  if(endType == 'steps') then
     if(endValue < 1.0) then
        write(*,*) "ERROR: endValue is too small for endType == 'steps'"
        call endRun()
     end if
  end if
  
  return
end subroutine checkWritings
! =======================================================
subroutine initValues()
  use general
  use grids
  use io
!  use relax
  implicit none

  !time:
  if(endType == "time") then
     endTime = endValue
  else ! endType == "steps"
     endSteps = int(endValue)
  end if
  if(writeType == "time") then
     timeWriteLogs = freqLogs
     timeWriteFiles = freqFiles
     timeWriteMaps = freqMaps
  else ! writeType == "steps"
     nWriteLogs = int(freqLogs)
     nWriteFiles = int(freqFiles)
     nWriteMaps = int(freqMaps)
  end if
  
  counter = 0
  itime = 0
  time = 0.0
  filenumber = 0
  nextTimeWriteLogs = 0.0
  nextTimeWriteFiles = 0.0
  nextTimeWriteMaps = 0.0
  
  call computeMinTimestep(minTimestep)

  call checkExistingFiles(startingFromRead)
  !! FILES ARE ALSO READ HERE !!!!!!!!!!!!!!!
  !! CONTINUERUN and ETAREAD ARE SET !!!!!!!!

  if(.not.startingFromRead) then
     initNumber = 0
     write(*,'(a)') &
          "Run is starting with a linear profile of temperature"
     call conductiveTemp()
  end if

  ! ------------------------------------------------
  ! special top boundary conditions for the velocity
  ! ------------------------------------------------
  if(topV == "file") then
     ! read imposed top velocity vh=uphi if topV == "file"
     call read1DFile(fileVtop, vxtop, n0(1), 1) ! (.., .., idir, idplan)
  end if

  ! init kright and kleft for boundary conditions in the multigrid -------
  if(leftV == dirichlet) then
     kleft = 1.0
  else if(leftV == neumann) then
     kleft = -1.0
  end if
  
  if(rightV == dirichlet) then
     kright = 1.0
  else if(rightV == neumann) then
     kright = -1.0
  end if
  
  ! init the viscosity ------------------------------
  call initVisco() ! defines if isoviscous, yielding, or not +
  ! computes etaT and eta if needed
  
  write(*,'(a)') "========================================" 
  if(isoviscous) then
     write(*, '(a)') "isoviscous"
  else
     write(*,'(a,e16.8)') 'Eta top: ', etaTop
     write(*,'(a,e16.8)') 'Eta bot: ', etaBot
     if(reference == "bottom") then
        write(*,'(a,e16.8)') 'Rayleigh top: ', rayleigh * etaBot / etaTop
        write(*,'(a,e16.8)') 'Rayleigh bot (REFERENCE): ', rayleigh
     else if(reference == "top") then
        write(*,'(a,e16.8)') 'Rayleigh top (REFERENCE): ', rayleigh
        write(*,'(a,e16.8)') 'Rayleigh bot: ', rayleigh * etaTop / etaBot
     else  ! reference == "T_ref"
        write(*,'(a,e16.8)') 'Rayleigh top: ', rayleigh / etaTop
        write(*,'(a,e16.8)') 'Rayleigh T_ref (REFERENCE): ', rayleigh
        write(*,'(a,e16.8)') 'Rayleigh bot: ', rayleigh / etaBot
     end if
  end if

  call applyBcs(temp, 0, 3)
  
  ! init thermal conductivity -----------------------
  call setValueVec(k, 0, one)  ! k=1 everywhere except tracers

  ! COMPUTE parameters -------------------------------
  relax_on_error = .false. ! init, DO NOT CHANGE!!!
  relax_corrSf = .false. ! IDEM ! DO NOT CHANGE
  nrelax = nrelaxMin
  counterNRelax = 0
  nextIncrease = 10
  
  ! write parameters 
  call writeParameters()

  return
end subroutine initValues
! =======================================================
subroutine checkExistingFiles(readFile)
  ! check if the files that should be read are present (beginMethod=="read"),
  ! or see if there are existing files to continue the run from
  ! (beginMethod="continue")
  use general
  use grids
  use io
  implicit none
  character(len=100):: fileName
  logical, intent(out):: readFile
  logical:: existOutput, existInput, existLog
  integer:: number

  continueRun = .false.
  readFile = .false.
  etaRead = .false.
  
  existOutput = .false.
  existInput = .false.
  existLog = .false.
  
  if(beginMethod(1:8) == "continue") then  ! - - - - - - - - - - - - - - - - - 

     ! searching for file with largest number
     number = 9999

     do while(.not. existOutput .and. number >= 0)
        call makeName(fileName, outputStem, number)
        call existFile(fileName, existOutput)
        if(existOutput) then
           readFile = .true.
           write(*,'(2a)') "beginMethod=continue: Found file ", &
                trim(adjustl(fileName))
        end if
        number = number - 1
     end do ! while .not. existOutput
     ! read temperature if one file was found
     if(readFile) then
        continueRun = .true.
        call readInitFile(fileName, temp, 3, sf, 4)

        ! read also the viscosity if exists
        suffixe = "eta"
        call makeName(filename, outputStem, number+1, suffixe) ! +1 because of l.311
        call existFile(filename, existOutput)
        if(existOutput) then
           call readInitFile(filename, eta, 3, sist, 3, &
                setTime=.false., mandatory=.false.)
           etaRead = .true.
        end if
        
        ! initializing time and itime is done in readInitFile()
        initNumber = number + 1
        filenumber = initNumber + 1
     else
        write(*, '(a)') "ERROR: beginMethod=continue & 
             & but could not find any file."
        call endRun()
     end if

  else

     ! check if log files already exist here
     fileName = trim(adjustl(outputStem))//"_temp.out"
     call existFile(fileName, existLog)
     if(.not.existLog) then
        fileName = trim(adjustl(outputStem))//"_veloc.out"
        call existFile(fileName, existLog)
        if(.not.existLog) then
           fileName = trim(adjustl(outputStem))//"_surf.out"
           call existFile(fileName, existLog)
        end if
     end if
     ! check if other files exist:
     if(.not.existLog) then
        number = 0
        do while(.not.existOutput .and. number < 10000)
           ! checking that no output files already exists here
           call makeName(fileName, outputStem, number)
           call existFile(fileName, existOutput)
           number = number + 1
        end do
     end if
     if(existLog .or. existOutput) then
        write(*,'(a)') &
             "ERROR: output files already exist here & 
             & -> change outputStem before launching run again."
        call endRun()
     end if

     if(beginMethod(1:4) == "read") then ! - - - - - - - - - - - - -

        call existFile(fileToRead, existInput)
        if(existInput) then
           write(*,'(2a)') "beginMethod=Read: Found file ", &
                trim(adjustl(fileToRead))
           call readInitFile(fileToRead, temp, 3, sf, 4)

           readFile = .true.
           initNumber = 0

           ! read also the viscosity if exists
           call existFile(fileEtaToRead, existOutput)
           if(existOutput) then
              call readInitFile(fileEtaToRead, eta, 3, sist, 3, &
                   setTime=.false., mandatory=.false.)
              etaRead = .true.
           end if
           
        else ! input file does not exist
           write(*,'(3a)') "ERROR when trying to start from ", &
                trim(adjustl(fileToRead)), " -> not found"
           call endRun()
        end if

     end if ! beginMethod == read

  end if ! beginMethod

  return
end subroutine checkExistingFiles
! =======================================================
subroutine lithoPressure()
  use grids
  implicit none
  real(rk):: posZ
  integer:: j, id
  
  do j=i0-1, iend(0,2)+1
     do id=3, 4
        call getPosZ(posZ, j, 0, id)
        plitho(j, id) = rayleigh / alphaDt * (rtop - posZ)
     end do
  end do
  
  return
end subroutine lithoPressure
! =======================================================
subroutine conductiveTemp()
  ! makes a simple conductive profile for the temperature
  use general
  use grids
  implicit none
  integer j
  real(rk) posZ, val

  do j=i0-2, iend0(2)+2
     
     call getPosZ(posZ, j, 0, 3) ! r-position of the node j

     if(coordinates == "cylindrical") then
        val = topValueT + (botValueT - topValueT) / log(rtop/rbot) * &
             log(rtop / posZ)
     else ! spherical
        val = botValueT + posZ * (topValueT - botValueT) * &
             rbot / posZ * (rtop - posZ) / (rtop-rbot)    
     end if

     call setValuePlan(temp, 0, 3, 2, j, val)
  end do

  call addNoise()

  return
end subroutine conductiveTemp
! =======================================================
subroutine addNoise()
  ! add random noise to the temperature field, 
  ! with a max amplitude equal to "amplitudeNoise"
  use general
  use grids
  use mylib
  implicit none
  real(rk), pointer:: ptemp(:,:) => null()
  real(rk), allocatable, dimension(:,:):: randnoise
  real(rk), parameter:: amplitudeNoise = 0.05
  integer, parameter:: seed = 7164
  integer:: i,j

  call associatePtrScal(ptemp, temp, 0, 3)

  allocate(randnoise(i0-1:iend0(1)+1, i0-1:iend0(2)+1))
  ! make random field
  call generateRandomScal(randnoise, n0, seed)

  do j=i0, iend0(2)
     do i=i0, iend0(1)
        
        ptemp(i,j) = &
             ptemp(i,j) + &
             amplitudeNoise * (2.0 * randnoise(i,j) - 1.0)

        if(ptemp(i,j) > botValueT .and. botT == "dirichlet") &
             ptemp(i,j) = botValueT

        if(ptemp(i,j) < topValueT) &
             ptemp(i,j) = topValueT

     end do
  end do

  deallocate(randnoise)
  nullify(ptemp)

  return
end subroutine addNoise
! =============================================
