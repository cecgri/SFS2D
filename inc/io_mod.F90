! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
module io

  use general
  use grids
  use mylib
  implicit none

  logical:: continueRun
  logical:: etaRead
  logical:: writeLogsNow, writeFilesNow, writeMapsNow
  logical:: writeSurfNow
  character(len=100):: suffixe
  character(len=100), parameter:: &
       fmtHeaderCurv='(a8,i8,a9,e16.6e3,a18,2i5,a10,f10.6,a12,f8.2,a8,a1)'
  character(len=50), parameter:: &
       fmtHeaderCart='(a8,i8,a9,e16.6e3,a18,2i5,a6,f10.6)'
  character(len=50), parameter:: fmtOneCol='(e16.6e3)'
  character(len=50), parameter:: fmtTwoCols='(2e16.6e3)'
  
  type logOutput 
     character(len=100):: filename
     character(len=50):: fmt
     character(len=200):: header
     integer:: ndata, unit, ncount
     real(rk), allocatable:: data(:)
  end type logOutput

  type(logOutput):: logTemp, logVeloc, logSurf, logConc, logVcell
  type(logOutput):: logDiag, logExtra

  type mapOutput
     character(len=100):: filename
     integer:: unit
  end type mapOutput

  character(len=50):: mapFmt
  character(len=200):: mapHeader
  type(mapOutput):: mapVsurf, mapQtop, mapCont, mapTemp25, mapTemp50, mapTemp75

  
contains
  ! =======================================================
  subroutine initLogs()
    implicit none
    integer:: ncount, nunit

    ncount = 1
    nunit = 60
    
    ! temperature - - - - -
    logTemp%filename = trim(adjustl(outputStem))//"_temp.out"
    logTemp%ncount = ncount  ! identifier
    ncount = ncount + 1
    logTemp%header = &
         "#Step(1) Time(2) Tmean(3) Qtop(4) Qbot(5) NCells(6)"
    logTemp%ndata = 4
    logTemp%fmt = '(i10, 5e16.8)'
    logTemp%unit = nunit
    allocate(logTemp%data(1:logTemp%ndata))

    ! velocity - - - - - -
    logVeloc%filename = trim(adjustl(outputStem))//"_veloc.out"
    logVeloc%ncount = ncount
    ncount = ncount + 1
    logVeloc%header = &
         "#Step(1) Time(2) Vrms(3) Vmax(4) & 
         & V10(5) V25(6) V50(7) V75(8) V90(9)"
    logVeloc%ndata = 7
    logVeloc%fmt = '(i10, 8e16.8)'
    nunit = nunit + 1
    logVeloc%unit = nunit
    allocate(logVeloc%data(1:logVeloc%ndata))
    
    ! surface velocity and strain rate - - - - -
    logSurf%filename = trim(adjustl(outputStem))//"_surf.out"
    logSurf%ncount = ncount
    ncount = ncount + 1
    logSurf%header = &
         "#Step(1) Time(2) Vrms(3) Vmax(4) &
         & edot10(5) edot25(6) edot50(7) edot75(8) &
         & edot90(9) edotMax(10) mobility(11)"
    logSurf%ndata = 9
    logSurf%fmt = '(i10, 10e16.8)'
    nunit = nunit + 1
    logSurf%unit = nunit
    allocate(logSurf%data(1:logSurf%ndata))

    ! velocities around the cell - - - - - - -
    logVcell%filename = trim(adjustl(outputStem))//"_vcell.out"
    logVcell%ncount = ncount
    ncount = ncount + 1
    logVcell%header = &
         "#Step(1) Time(2) & 
         & Utopmax(3) Ubotmax(4) Vup(5) Vdown(6) &
         & Utoprms(7) Ubotrms(8)"
    logVcell%ndata  = 6
    logVcell%fmt = '(i10, 7e16.8)'
    nunit = nunit + 1
    logVcell%unit = nunit
    allocate(logVcell%data(1:logVcell%ndata))
    
    ! diagnosis on computation time - - - - - -
    logDiag%filename = trim(adjustl(outputStem))//"_diag.out"
    logDiag%ncount = ncount
    ncount = ncount + 1
    logDiag%header = &
         "#Step(1) Time(2) Timestep(3) & 
         & cpuTime:_all(4) _solver(5) _heatcons(6) _writing(7) _nloop(8) _nrelax(9)"
    logDiag%ndata = 7
    logDiag%fmt = '(i10, 8e16.8)'
    nunit = nunit + 1
    logDiag%unit = nunit
    allocate(logDiag%data(1:logDiag%ndata))

    ! Extra diagnosis for benchmark by Tosi et al., 2015, G3,  - - - -
    logExtra%filename = trim(adjustl(outputStem))//"_extraDiag.out"
    logExtra%ncount = ncount
    ncount = ncount + 1
    logExtra%header = &
         "#Step(1) Time(2) etaMin(3) etaMax(4) Work(5) Dissipation(6)"
    logExtra%ndata = 4
    logExtra%fmt = '(i10, 5e16.8)'
    nunit = nunit + 1
    logExtra%unit = nunit
    allocate(logExtra%data(1:logExtra%ndata))
    
    call initFile(logTemp)
    call initFile(logVeloc)
    call initFile(logSurf)
    call initFile(logVcell)
    call initFile(logDiag)
    call initFile(logExtra)

    return
  end subroutine initLogs
  ! =======================================================
  subroutine clearLogs()
    implicit none

    if(allocated(logTemp%data)) deallocate(logTemp%data)
    if(allocated(logVeloc%data)) deallocate(logVeloc%data)
    if(allocated(logSurf%data)) deallocate(logSurf%data)
    if(allocated(logConc%data)) deallocate(logConc%data)
    if(allocated(logVcell%data)) deallocate(logVcell%data)
    if(allocated(logDiag%data)) deallocate(logDiag%data)
    if(allocated(logExtra%data)) deallocate(logExtra%data)

  end subroutine clearLogs
  ! =======================================================
  subroutine readParameters()
    ! read the file paramFile that contains input parameters
    implicit none

    ! read file with input parameters
    open(unit=10, file=paramFile, status='old')
    read(10, physics)
    read(10, geometry)
    read(10, inout)
    read(10, boundaries)
    read(10, timer)
    read(10, compute)
    read(10, advdiff)

    write(*,*) '   '
    write(*,*) '--------- physics ---------'
    write(*, physics)
    write(*,*) '--------- geometry --------'
    write(*, geometry)
    write(*,*) '---------- inout ----------'
    write(*, inout)
    write(*,*) '-------- boundaries -------'
    write(*, boundaries)
    write(*,*) '---------- timer ----------'
    write(*, timer)
    write(*,*) '--------- compute ---------'
    write(*, compute)
    write(*,*) '--------- advdiff ---------'
    write(*, advdiff)
    write(*,*) '---------------------------'
    
    write(*,*) '      '       
    
    close(10)  

    return
  end subroutine readParameters
  ! =======================================================
  subroutine readInitFile(filename, what1, id1, what2, id2, setTime, mandatory)
    ! read file=filename, that has two columns (T and phi for instance)
    !  or only T
    implicit none
    character(len=10), intent(in):: what1, what2
    integer, intent(in):: id1, id2
    character(len=100), intent(in):: filename
    logical, optional:: setTime, mandatory
    integer:: nlines, nexpected, i,j
    real(rk), dimension(0:n0(1)+1, 0:n0(2)+1):: a1, a2
    real(rk), pointer:: ptr(:,:) => null()
    character(len=100):: c_buff, r_geom
    real(rk):: r_time, r_lx, r_fcurv, r_angle
    character(len=1):: r_letter
    integer:: r_n0(1:2), r_itime, ibord, iostatus
    integer, parameter:: nmax_cols = 2
    real(rk), dimension(1:nmax_cols):: test_array
    integer:: ncols
    logical:: okayFile, readBorders, twoCols
    logical:: doSetTime, doMandatory

    doSetTime = .true.
    doMandatory = .true.
    if(present(setTime)) then
       doSetTime = setTime
    end if
    if(present(mandatory)) then
       doMandatory = mandatory

    end if
    
    okayFile = .true.
    readBorders = .true.
    ibord = 1
    twoCols = .true.
    
    !check the number of lines in filename
    call numberOfLines(nlines, filename)
    nexpected = (n0(1)+2) * (n0(2)+2) + 1 ! +1 for header line

    if(nlines /= nexpected) then
       ! check if it could be a file without ghost points:
       if(nlines == n0(1)*n0(2) + 1) then
          nexpected = nlines
          readBorders = .false.
          ibord = 0
          write(*,'(3a)') "Warning: file ", trim(adjustl(filename)), &
               " does not include borders."
       end if
    end if

    if(nlines == 0 .or. nlines /= nexpected) then 
       ! reading in numberOfLines went wrong.
       write(*,'(a)') "ERROR in read ", trim(adjustl(filename)),  " -> number of lines does not match param.in"
       if(doMandatory) then
          call endRun()
       end if
    endif
    
    ! reading = = = = = = 
    open(10, file=filename, status='old')
    ! first line (header)
    
    if(curvedGeometry) then
       ! spherical annulus or cylindrical geometry ------------------------
       read(10, fmtHeaderCurv) c_buff, r_itime, c_buff, &
            r_time, c_buff, r_n0(1), r_n0(2), c_buff, r_fcurv, &
            c_buff, r_angle, c_buff, r_letter
       
       if(r_letter /= letterGeom) then
          r_geom = "spherical"
          if(r_letter == "C") r_geom = "cylindrical"
          write(*, '(4a)') "WARNING: geometry is set to ", &
               trim(adjustl(coordinates)), &
               " while read file was done with ", &
               trim(adjustl(r_geom))
       end if
    else ! CARTESIAN ------------------------------
       read(10, fmtHeaderCart) c_buff, r_itime, c_buff, &
            r_time, c_buff, r_n0(1), r_n0(2), c_buff, r_lx
    end if ! curvedGeometry or not ----------------

    if(r_n0(1) /= n0(1) .or. r_n0(2) /= n0(2)) then
       okayFile = .false.
       write(*,'(4a,2i6)') "ERROR: ", &
            trim(filename), " has grid:" , r_n0(1), r_n0(2)
       write(*,'(a,2i6,a)') "Grid ", n0(1), n0(2), " was expected."
    end if

    ! check the number of columns
    if(okayFile) then
       ! read second line to check if one or two columns
       read(10, '(a)', iostat=iostatus) c_buff
       do ncols=1, nmax_cols
          read(c_buff, *, iostat=iostatus) test_array(1:ncols)
          if(iostatus /= 0) exit 
       end do
       ncols = ncols - 1
       if(ncols == 1) twoCols = .false.
    end if

    if(okayFile) then
       ! restart at the beginning of the file
       rewind(10)

       ! first line (header)
       if(curvedGeometry) then ! -------------
          read(10, fmtHeaderCurv) c_buff, r_itime, c_buff, &
               r_time, c_buff, r_n0(1), r_n0(2), c_buff, r_fcurv, &
               c_buff, r_angle, c_buff, r_letter
       else ! Cartesian ---------------------
          read(10, fmtHeaderCart) c_buff, r_itime, c_buff, &
               r_time, c_buff, r_n0(1), r_n0(2), c_buff, r_lx
       end if ! curvedGeometry ---------------

       if(twoCols) then
          do j=1-ibord, n0(2)+ibord
             do i=1-ibord, n0(1)+ibord
                read(10, fmtTwoCols) a1(i,j), a2(i, j)
             end do ! i
          end do ! j
       else
          do j=1-ibord, n0(2)+ibord
             do i=1-ibord, n0(1)+ibord
                read(10, fmtOneCol) a1(i,j)
             end do ! i
          end do ! j
       end if

    end if ! okayFile

    close(10)
    write(*,'(a)')  "-------- readInitFile ---------"  
    write(*,'(3a)') " File ", trim(filename), " read"
    
    if(.not. okayFile .and. doMandatory) &
         call endRun() ! file didn't have the right size
    
    if(continueRun .and. doSetTime) then
       itime = r_itime
       time = r_time
       write(*,'(a,e13.4)') "...starting time is: ",  time
       nextTimeWriteLogs = time + timeWriteLogs
       nextTimeWriteFiles = time + timeWriteFiles
       nextTimeWriteMaps = time + timeWriteMaps
    end if

    ! FIRST COLUMN (temperature for normal output file)
    call associatePtrScal(ptr, what1, 0, id1)
    ptr(i0-1:iend0(1)+1,i0-1:iend0(2)+1) = a1(0:n0(1)+1,0:n0(2)+1)
    nullify(ptr)
    
    if(twoCols) then
       ! SECOND COLUMN 
       call associatePtrScal(ptr, what2, 0, id2)
       ptr(i0-1:iend0(1)+1, i0-1:iend0(2)+1) =a2(0:n0(1)+1,0:n0(2)+1)
       nullify(ptr)
       write(*, '(3a)') "-> ", what1, what2
    else
       write(*, '(2a)') "-> ", what1
    end if  ! twoCols
    
    if(.not.readBorders) then
       call applyBcs(what1, 0, id1)
       call applyBcs(what2, 0, id2)
    end if

    return
  end subroutine readInitFile
  ! =======================================================
  subroutine initFile(log)
    implicit none
    type(logOutput), intent(inout):: log
    logical:: existLog, cutFile
    
    call existFile(log%filename, existLog)

    ! Used in the beginning, existLog can be T or F
    if(.not.existLog) then
       open(unit=log%unit, file=log%filename, status='new')
       write(log%unit, '(a)') trim(adjustl(log%header))
       close(log%unit)
    else
       if(continueRun) then
          ! Continuing a run: log file needs to end at "itime"
          call prepareFile(log, cutFile)
          if(log%ncount == 1) then
             write(*, '(a)') "========================================"
             write(*,'(a)') "Warning: log files already exists -> append"
             if(cutFile) then
                write(*,'(a)') &
                     "         files are going past time in read file &
                     & -> cut"
             end if
          end if
       else
          ! one logFile is present while there was no output file to read: 
          ! must be error -> abort
          write(*,'(3a)') "ERROR: ", trim(log%filename), &
               " already exists here -> CANCELLED"
          call endRun()
       end if
    end if ! .not.existLog

    return
  end subroutine initFile
  ! =======================================================
  subroutine prepareFile(log, cutFile)
    ! check the line in a log file where "time"
    ! is reached. Re-write the log file until this 
    ! line only.
    implicit none
    type(logOutput), intent(inout):: log
    logical, intent(out):: cutFile
    character(len=200):: read_header
    integer:: iostatus, r_itime, nline, n, k
    real(rk):: r_time
    real(rk), dimension(1:log%ndata):: r_buff
    integer, allocatable:: r_idata(:)
    real(rk), allocatable:: r_data(:,:)

    iostatus = 0
    nline = 0
    cutFile = .false.

    ! find lines where "time" is reached in the logFile
    open(unit=log%unit, file=log%filename, status='old')
    read(log%unit, '(a)') read_header
    nline = 1

    if(read_header /= log%header) then
       write(*,'(3a)') "ERROR: ", trim(adjustl(log%filename)), &
            " does not have the correct header"
       write(*,'(3a)') "--> move the file ", &
            trim(log%filename), " before continuing"
       call endRun()
    end if
       
    do 
       read(log%unit, log%fmt, iostat=iostatus) &
            r_itime, r_time, (r_buff(k), k=1, log%ndata)
       if(iostatus==0 .and. r_time > time) then
          cutFile = .true.
          exit
       end if
       if(iostatus /= 0) then
          if(iostatus > 0) then
             write(*,'(2a)') "ERROR: problem while reading ", &
                  trim(adjustl(log%filename))
          end if
          ! if iostatus<0: EOF
          exit
       end if
       nline = nline + 1
    end do
    
    close(log%unit)

    if(iostatus > 0) then ! there was a problem while reading
       call endRun()
    end if

    nline = nline - 1
    allocate(r_idata(1:nline))
    allocate(r_data(1:nline, 0:log%ndata)) !0 to add the time
       
    ! read data in logFile
    open(unit=log%unit, file=log%filename, status='old')
    read(log%unit, '(a)') read_header
       
    do n=1, nline 
       read(log%unit, log%fmt) &
            r_idata(n), (r_data(n,k), k=0, log%ndata)
    end do ! n=1,nline
    close(log%unit)
    
    ! give last values to log%data(:) (useful for log=diag)
    log%data(1:log%ndata) = r_data(nline, 1:log%ndata)
    
    if(cutFile) then
       ! rewrite logFile with data only until "time"
       open(unit=log%unit, file=log%filename, status='replace')
       write(log%unit, '(a)') trim(adjustl(log%header))
       do n=1, nline
          write(log%unit, log%fmt) r_idata(n), (r_data(n,k), k=0, log%ndata)
       end do
       close(log%unit)
    end if ! cutFile
    
    deallocate(r_data, r_idata)
  
    return
  end subroutine prepareFile
  ! =======================================================
  subroutine writeScal(what, id, number, suff, level)
    ! writes "what(:,:,id)" in a file with filename "stem_{suffixe}{number}"
    ! e.g. iso_t0003 if what=temp, suffixe=t and number=3.
    ! Writes only (i0:n0(1), i0:n0(2)) if id == 3
    ! If id = 4 : writes (i0:n0(1)+1, i0:n0(2)+1)
    ! If id=1 or id=2: faceToCorner and then writing id=4 (corners)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: id, number
    character(len=100), optional:: suff
    integer, optional:: level
    character(len=100):: fileName
    real(rk), pointer:: ptr(:,:) => null()
    real(rk), allocatable, dimension(:,:):: a0
    integer:: i,j, lev, nl(1:2)

    lev = 0
    if(present(level)) then
       lev = level
    end if

    if(present(suff)) then
       if(present(level)) then
          call makeName(filename, outputStem, number, nameSuff=suff, &
               nameLevel=lev)
       else
          call makeName(filename, outputStem, number, nameSuff=suff)
       end if
    else
       if(present(level)) then
          call makeName(filename, outputStem, number, nameLevel=lev)
       else
          call makeName(filename, outputStem, number)
       end if
    end if
    
    call associatePtrScal(ptr, what, lev, id)
    nl = shape(ptr) - 4
    allocate(a0(0:nl(1)+1, 0:nl(2)+1))

    do j=0,nl(2)+1
       do i=0,nl(1)+1
          a0(i,j) = ptr(i-1+i0,j-1+i0)
       end do
    end do
    
    open(unit=10, file=fileName, status='unknown')

    if(curvedGeometry) then ! -----------------------------------
       write(10, fmtHeaderCurv) "#itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", nl(1), nl(2), &
            "  fcurv: ", fcurv, "  anglemax: ", anglemax, &
            "  geom: ", letterGeom
    else ! Cartesian --------------------------------------------
       write(10, fmtHeaderCart) "itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", nl(1), nl(2), &
            "  lx: ", lengthX
    end if ! curvedGeometry -------------------------------------

    do j=0, nl(2)+1
       do i=0, nl(1)+1
          write(10, '(e16.8)') a0(i, j)
       end do
    end do
    close(10)
    if(verbose > 0) &
         write(*,'(2a)') trim(adjustl(fileName)), " written"
    
    deallocate(a0)
    nullify(ptr)

    return
  end subroutine writeScal
  ! =======================================================
  subroutine writeHorizontalWalls(number)
    ! write velocity and heatflux on top and bottom surface
    implicit none
    integer, intent(in):: number
    character(len=100):: fileName
    real(rk), dimension(i0:iend0(1)):: vtop, vbot, nutop, nubot
    real(rk), pointer:: pv(:,:) => null()
    real(rk), pointer:: pt(:,:) => null()
    integer:: i, iez

    suffixe = "horiz"
    call makeName(filename, outputStem, number, nameSuff=suffixe)

    call associatePtrScal(pv, vel, 0, 1) ! vh
    call associatePtrScal(pt, temp, 0, 3)

    iez = iend0(2)
    do i=i0, iend0(1)
       vtop(i) = 0.5 * (pv(i, iez) + pv(i, iez+1))
       vbot(i) = 0.5 * (pv(i, i0) + pv(i, i0-1))
       nutop(i) = (pt(i, iez) - pt(i, iez+1)) * invdz
       nubot(i) = (pt(i, i0-1) - pt(i, i0)) * invdz
    end do

    nullify(pv, pt)

    open(unit=15, file=filename, status='new')
    do i=i0, iend0(1)
       write(15, '(4e16.8)') vtop(i), vbot(i), nutop(i), nubot(i)
    end do ! i
    close(15)
    if(verbose > 0) &
         write(*,'(2a)') trim(adjustl(fileName)), " written"
    
    return
  end subroutine writeHorizontalWalls
  ! =======================================================
  subroutine writeOutputFile(number)
    ! writes Temp, SF in stem_{number}
    implicit none
    integer, intent(in):: number
    real(rk), pointer:: pt(:,:) => null()
    real(rk), pointer:: psf(:,:) => null()
    integer:: i, j
    character(len=100):: filename
    
    call makename(filename, outputStem, number)

    ! temperature at grid center (id=3)
    call associatePtrScal(pt, temp, 0, 3)
    ! stream function at grid corner (id=4)
    call associatePtrScal(psf, sf, 0, 4)
    
    open(unit=10, file=filename, status='new')
    if(curvedGeometry) then ! -----------------------------------------
       write(10,fmtHeaderCurv) "#itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  fcurv: ", fcurv, "  anglemax: ", anglemax, &
            "  geom: ", letterGeom
    else ! Cartesian ---------------------------------------------------
       write(10,fmtHeaderCart) "#itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  lx: ", lengthX
    end if ! curvedGeometry or not --------------------------------------

    do j=i0-1, iend0(2)+1
       do i=i0-1, iend0(1)+1
          write(10, fmtTwoCols) &
               pt(i,j), psf(i, j)
       end do ! i
    end do ! j
    close(10)

    nullify(pt, psf)

    if(verbose > 0) then
       write(*,'(a)') " --------------------------------------------------------"
       write(*,'(2a, e12.4)') trim(fileName), " written;  time = ", time
    end if
 
    if(yielding .or. writeEta) then
       call writeEtaFile(number)
    end if
    
    if(verbose > 0) then
       write(*,'(a)') " --------------------------------------------------------"
    end if
    
    return
  end subroutine writeOutputFile
  ! =======================================================
  subroutine writeEtaFile(number)
    ! writes etaEff in one column and the second invariant
    ! of the strainrate tensor in the second column.
    ! WARNING: phi_to_sist() must have been done before
    implicit none
    integer, intent(in):: number
    real(rk), pointer:: peta(:,:) => null()
    real(rk), pointer:: psist(:,:) => null()
    integer:: i, j, id
    character(len=100):: filename

    suffixe = "eta"
    call makename(filename, outputStem, number, suffixe)

    id = 3
    if(yielding) then
       ! eta 
       call associatePtrScal(peta, eta, 0, id)
    else
       ! etaT
       call associatePtrScal(peta, etaT, 0, id)
    end if
    call associatePtrScal(psist, sist, 0, id)

    open(unit=23, file=filename, status='new')
    if(curvedGeometry) then ! -----------------------------------------
       write(23,fmtHeaderCurv) "#itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  fcurv: ", fcurv, "  anglemax: ", anglemax, &
            "  geom: ", letterGeom
    else ! Cartesian ---------------------------------------------------
       write(23,fmtHeaderCart) "#itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  lx: ", lengthX
    end if ! curvedGeometry or not --------------------------------------

    do j=i0-1, iend0(2)+1
       do i=i0-1, iend0(1)+1
          write(23, fmtTwoCols) &
               peta(i,j), psist(i,j)
       end do
    end do

    nullify(peta, psist)

    close(23)
    if(verbose > 0) &
         write(*,'(2a, e12.4)') trim(fileName), " written; time = ", time

    return
  end subroutine writeEtaFile
  ! =======================================================
  subroutine writeParameters()
    ! write parameters at the start of a run (itime)
    implicit none
    character(len=100):: filename

    suffixe = "param"
    call makeName(filename, outputStem, initNumber, suffixe)

    open(unit=10, file=filename, status='unknown')
    write(10,*) '--------- physics ---------'
    write(10, physics)
    write(10,*) '--------- geometry --------'
    write(10, geometry)
    write(10,*) '---------- inout ----------'
    write(10, inout)
    write(10,*) '-------- boundaries -------'
    write(10, boundaries)
    write(10,*) '---------- timer ----------'
    write(10, timer)
    write(10,*) '--------- compute ---------'
    write(10, compute)
    write(10,*) '--------- advdiff ---------'
    write(10, advdiff)
    write(10,*) '---------------------------'

    return
  end subroutine writeParameters
  ! =======================================================
  subroutine setWriteNow()
    implicit none

    writeLogsNow = .false.
    writeFilesNow = .false.
    writeMapsNow = .false.
    
    if(continueRun .and. counter == 0) return

    if(writeType == "time") then
       if(time >= nextTimeWriteLogs) then
          writeLogsNow = .true.
          nextTimeWriteLogs = nextTimeWriteLogs + timeWriteLogs
       end if
       if(time >= nextTimeWriteFiles) then
          writeFilesNow = .true.
          nextTimeWriteFiles = nextTimeWriteFiles + timeWriteFiles
       end if
       if(time >= nextTimeWriteMaps) then
          writeMapsNow = .true.
          nextTimeWriteMaps = nextTimeWriteMaps + timeWriteMaps
       end if
    else ! if(writeType == "steps") then
       if(modulo(itime, nWriteLogs) == 0) writeLogsNow = .true.
       if(modulo(itime, nWriteFiles) == 0) writeFilesNow = .true.
       if(modulo(itime, nWriteMaps) == 0) writeMapsNow = .true.
    end if

    return
  end subroutine setWriteNow
  ! =======================================================
  subroutine writings()
    implicit none

    call setWriteNow()
    
    if(writeFilesNow .or. writeLogsNow) then
       call applyBcs(temp, 0, 3)
       call centerToCorner(temp, 0) ! put values at cells corners
    end if

    if(writeLogsNow) call writeLogs()

    if(writeMapsNow) call writeMaps()
    
    if(writeFilesNow) then
       call writeFiles()
       if(writeProfiles) then
          call writeSurfProf()
       end if
       filenumber = filenumber + 1
    end if

    return
  end subroutine writings
  ! =======================================================
  subroutine writeLogs()
    implicit none
    real(rk):: meanT, &
         qmeanTop, qmeanBot, ncells, &
         vrms, vmax, edotmax, vsurfrms, vsurfmax, &
         etaMin, etaMax, work, dissipation, &
         vtopmax, vbotmax, vupmax, vdownmax, vtoprms, vbotrms
    real(rk):: percentiles(1:5)
    
    call cpu_time(t_start)

    ! compute data for log = logTemp - - - - - 
    call getMean(meanT, temp, 0) ! done for id=3
    call computeHeatFlux(qmeanTop, qmeanBot)
    call computeNCells(ncells)
    logTemp%data(1) = meanT
    logTemp%data(2) = qmeanTop
    logTemp%data(3) = qmeanBot
    logTemp%data(4) = ncells
    
    ! compute data for veloc file (log=logVeloc) - - - - - 
    call computeVmid(vrms, vmax, percentiles)
    logVeloc%data(1) = vrms
    logVeloc%data(2) = vmax
    logVeloc%data(3:7) = percentiles(1:5)

    ! vsurf and edot at the surface - - - - - - - - 
    call computeVsurfEdot(vsurfrms, vsurfmax, percentiles, edotmax)
    !warning: the outputs here are known only on nodezero
    logSurf%data(1) = vsurfrms
    logSurf%data(2) = vsurfmax
    logSurf%data(3:7) = percentiles(1:5)
    logSurf%data(8) = edotmax
    logSurf%data(9) = vsurfrms / vrms
    
    ! special velocities
    !(horizontal max in top half, bottom half, max verticals) - - - -
    call computeVcells(vtopmax, vbotmax, vupmax, vdownmax, &
         vtoprms, vbotrms)
    logVcell%data(1) = vtopmax
    logVcell%data(2) = vbotmax
    logVcell%data(3) = vupmax
    logVcell%data(4) = -vdownmax
    logVcell%data(5) = vtoprms
    logVcell%data(6) = vbotrms

    ! diagnosis file - - - - - - - - - - - - -
    logDiag%data(1) = timestep
    logDiag%data(2) = cpuTimeAll
    logDiag%data(3) = cpuTimeSolver
    logDiag%data(4) = cpuTimeTemperature
    logDiag%data(5) = cpuTimeOutput
    logDiag%data(6) = dble(nLoopTot)
    logDiag%data(7) = dble(nrelax)

    ! extra diagnosis (for benchmark with Tosi et al., 2015) - - - - 
    if (writeExtraDiagnosis) then
       call computeExtra(etaMin, etaMax, work, dissipation)
       logExtra%data(1) = etaMin
       logExtra%data(2) = etaMax
       logExtra%data(3) = work
       logExtra%data(4) = dissipation
    end if
    
    ! writing - - - - - - - -
    call writeLogFile(logTemp)
    call writeLogFile(logVeloc)
    call writeLogFile(logSurf)
    call writeLogFile(logDiag)
    call writeLogFile(logVcell)
    call writeLogFile(logExtra)

    if(verbose > 0) then
       write(*,'(i8, e14.6, 3f10.4, 2e13.4, f5.1, i4)') itime, time, &
            meanT, qmeanTop, qmeanBot, vrms, vsurfRms, ncells, nrelax
    end if

    call cpu_time(t_end)
    cpuTimeOutput = cpuTimeOutput + t_end - t_start

    return
  end subroutine writeLogs
  ! =======================================================
  subroutine writeFiles()
    implicit none

    call cpu_time(t_start)

    call writeOutputFile(filenumber)
    
    call cpu_time(t_end)
    cpuTimeOutput = cpuTimeOutput + t_end - t_start

    return
  end subroutine writeFiles
  ! =======================================================
  subroutine writeLogFile(log)
    implicit none
    type(logOutput), intent(inout):: log
    integer:: k

    open(unit=log%unit, file=log%filename, status='old', position='append')
    write(log%unit, log%fmt) itime, time, &
         (log%data(k), k=1, log%ndata)
    close(log%unit)

  end subroutine writeLogFile
  ! =======================================================
  subroutine writeSurfProf()
    ! writes vsurf and qtop at the surface,
    ! and average vertical profile of temperature
    implicit none
    integer:: i, j, iex, iez
    character(len=100):: filename
    real(rk), pointer:: pvx(:,:) => null()
    real(rk), pointer:: pt(:,:) => null()
    
    call cpu_time(t_start)

    call associatePtrScal(pvx, vel, 0, 1)
    call associatePtrScal(pt, temp, 0, 3)

    iex = iend0(1)
    iez = iend0(2)

    ! surf: vsurf and qtop = = = = = = = = 
    suffixe = "surf"
    call makeName(filename, outputStem, filenumber, suffixe)
    
    do i=i0-1, iex+1 
       ! surface velocity
       xval(1, i) = 0.5 * (pvx(i, iez) + pvx(i, iez+1))
       ! heat flux
       xval(2, i) = (pt(i, iez) - pt(i, iez+1)) * invdz
    end do

    nullify(pt, pvx)

    open(unit=20, file=filename, status='unknown')
    if(curvedGeometry) then ! -----------------------------------------
       write(20, fmtHeaderCurv) "#itime =", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  fcurv: ", fcurv, "  anglemax: ", anglemax, &
            "  geom: ", letterGeom
    else ! Cartesian --------------------------------------------------
       write(20, fmtHeaderCart) "itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  lx: ", lengthX
    end if ! curvedGeometry -------------------------------------------

    do i=i0-1, iend0(1)+1
       write(20, '(2e16.8)') xval(1, i), xval(2, i)
    end do
    close(20)
    if(verbose > 0) &
         write(*, '(3a)') "      ", trim(filename), " written"

    ! prof: <T>, <vh> et <vv>
    suffixe = "prof"
    call makeName(filename, outputStem, filenumber, suffixe)

    do j=i0, iez+1
       ! temperature <T> (centerToCorner must be done before)
       call getMeanPlan(zval(1, j), temp, 0, 4, 2, j)
       ! horizontal velocity <vh>
       call faceToCorner(vel, 0, 1)
       call getMeanPlan(zval(2, j), vel, 0, 4, 2, j, rms=.true.)
       ! vertical velocity <vv>
       call getMeanPlan(zval(3, j), vel, 0, 2, 2, j, rms=.true.)
    end do

    open(unit=30, file=filename, status='unknown')
    if(curvedGeometry) then ! -----------------------------------------
       write(30, fmtHeaderCurv) "#itime =", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  fcurv: ", fcurv, "  anglemax: ", anglemax, &
            "  geom: ", letterGeom
    else ! Cartesian --------------------------------------------------
       write(30, fmtHeaderCart) "itime = ", itime, "  time = ", &
            time, "  grid(nx,nz): ", n0(1), n0(2), &
            "  lx: ", lengthX
    end if ! curvedGeometry -------------------------------------------

    do j=i0, iend0(2)+1 ! written on faces (i0-1 not useful)
       write(30, '(3e16.8)') zval(1, j), &
            zval(2, j), zval(3, j)
    end do
    if(verbose > 0) &
         write(*, '(3a)') "      ", trim(filename), " written"
    close(30)

    call cpu_time(t_end)
    cpuTimeOutput = cpuTimeOutput + t_end - t_start
    
    return
  end subroutine writeSurfProf
  ! =============================================
  subroutine read1DFile(filename, a, nl, idir)
    implicit none
    character(len=100), intent(in):: filename
    integer, intent(in):: nl, idir
    real(rk), dimension(i0-1:i0+nl), intent(inout):: a
    integer:: i, iostatus
    logical:: exist

    if(.not.(idir == 1 .or. idir == 2)) then
       write(*,*) "ERROR: read1DFile with idir/= 1 or 2"
       call endRun()
    end if
    
    call existFile(filename, exist)
    if(.not. exist) then
       write(*, '(3a)') "ERROR in read1DFile: ", &
            trim(adjustl(filename)), " cannot be found"
       call endRun()
    end if

    open(unit=90, file=filename, status='old')
    do i=i0-1, i0+nl
       read(90, fmtOneCol, iostat=iostatus) a(i)
    end do
    if(iostatus > 0) then
       write(*,'(2a)') "ERROR: problem while reading ", &
            trim(adjustl(filename))
    end if
    close(90)
    write(*,'(3a)') " File ", trim(filename), " read"
    
    return
  end subroutine read1DFile
   ! =============================================
  subroutine initMaps()
    implicit none
    integer:: nunit
    character(len=10):: cnh, cnhp

    write(cnh, '(i5)') n0(1)
    write(cnhp, '(i5)') n0(1) + 1 ! +1 for time
    
    mapFmt = '(i7,*(e14.4e3))'
    mapHeader = "#Step(1) Time(2)  grid(nx): "//trim(adjustl(cnh))

    nunit = 70
    ! vsurf - - - - - -
    mapVsurf%filename = trim(adjustl(outputStem))//"_mapvsurf.out"
    mapVsurf%unit = nunit
    nunit = nunit + 1
    call initOneMap(mapVsurf)
    
    ! heat flux - - - - - -
    mapQtop%filename = trim(adjustl(outputStem))//"_mapqtop.out"
    mapQtop%unit = nunit
    nunit = nunit + 1
    call initOneMap(mapQtop)

    ! temperature at height 0.25, 0.50 and 0.75
    mapTemp25%filename = trim(adjustl(outputStem))//"_maptemp25.out"
    mapTemp25%unit = nunit
    nunit = nunit + 1
    call initOneMap(mapTemp25)
    
    mapTemp50%filename = trim(adjustl(outputStem))//"_maptemp50.out"
    mapTemp50%unit = nunit
    nunit = nunit + 1
    call initOneMap(mapTemp50)

    mapTemp75%filename = trim(adjustl(outputStem))//"_maptemp75.out"
    mapTemp75%unit = nunit
    nunit = nunit + 1
    call initOneMap(mapTemp75)
    
    return
  end subroutine initMaps
  ! =============================================
  subroutine initOneMap(map)
    implicit none
    type(mapOutput), intent(inout):: map
    logical:: exist
    
    call existFile(map%filename, exist)
    if(.not.exist) then
       ! create new file
       open(unit=map%unit, file=map%filename, status='new')
       write(map%unit, '(a)') trim(adjustl(mapHeader))
       close(map%unit)
    else  !
       if(continueRun) then
          ! Continuing a run: log file needs to end at "itime"
          call prepareMapFile(map)
       else
          ! one logFile is present while there was no output file to read: 
          ! must be error -> abort
          write(*,'(3a)') "ERROR: ", trim(map%filename), &
               " already exists here -> CANCELLED"
          call endRun()
       end if
    end if
    
    return
  end subroutine initOneMap
  ! =============================================
  subroutine prepareMapFile(map)
    ! check the line in a map file where "time"
    ! is reached. Re-write the map file until this 
    ! line only.
    implicit none
    type(mapOutput), intent(inout):: map
    character(len=200):: read_header
    logical:: cutFile
    integer:: iostatus, r_itime, nline, n, k
    real(rk):: r_time, r_buff
    integer, allocatable:: r_idata(:)
    real(rk), allocatable:: r_data(:,:)

    iostatus = 0
    nline = 0
    cutFile = .false.
    
    open(unit=map%unit, file=map%filename, status='old')
    read(map%unit, '(a)') read_header
    nline = 1

    if(read_header /= mapHeader) then
       write(*,'(3a)') "ERROR: ", trim(adjustl(map%filename)), &
            " does not have the correct header"
       write(*,'(3a)') "--> move the file ", &
            trim(map%filename), " before continuing"
       call endRun()
    end if

    do 
       read(map%unit, mapFmt, iostat=iostatus) &
            r_itime, r_time, (r_buff, k=1, n0(1))
       if(iostatus==0 .and. r_time > time) then
          cutFile = .true.
          exit
       end if
       if(iostatus /= 0) then
          if(iostatus > 0) then
             write(*,'(2a)') "ERROR: problem while reading ", &
                  trim(adjustl(map%filename))
          end if
          ! if iostatus<0: EOF
          exit
       end if
       nline = nline + 1
    end do

    close(map%unit)

    if(iostatus > 0) then ! there was a problem while reading
       call endRun()
    end if

    if(cutFile) then
       nline = nline - 1

       allocate(r_idata(1:nline))
       allocate(r_data(1:nline, 0:n0(1))) !0 to add time

       ! read data in mapFile
       open(unit=map%unit, file=map%filename, status='old')
       read(map%unit, '(a)') read_header

       do n=1, nline 
          read(map%unit, mapFmt) &
               r_idata(n), (r_data(n,k), k=0, n0(1))
       end do ! n=1,nline
       close(map%unit)

       ! rewrite logFile with data only until "time"
       open(unit=map%unit, file=map%filename, status='replace')
       write(map%unit, '(a)') mapHeader
       do n=1, nline
          write(map%unit, mapFmt) r_idata(n), (r_data(n,k), k=0, n0(1))
       end do
       close(map%unit)
       
       deallocate(r_data, r_idata)
    end if ! cutFile
    
    return
  end subroutine prepareMapFile
  ! =============================================
  subroutine writeMaps()
    implicit none
    integer:: i, iez
    real(rk), pointer:: pt(:,:) => null()
    real(rk), pointer:: pv(:,:) => null()

    call associatePtrScal(pt, temp, 0, 3) ! for temp and qtop
    
    ! temperature, at height 0.25, 0.50, and 0.75
    open(mapTemp25%unit, file=mapTemp25%filename, status='old', position='append')
    write(mapTemp25%unit, mapFmt) itime, time, &
         (pt(i, j25), i=i0, iend(0, 1))
    close(mapTemp25%unit)

    open(mapTemp50%unit, file=mapTemp50%filename, status='old', position='append')
    write(mapTemp50%unit, mapFmt) itime, time, &
         (pt(i, j50), i=i0, iend(0, 1))
    close(mapTemp50%unit)

    open(mapTemp75%unit, file=mapTemp75%filename, status='old', position='append')
    write(mapTemp75%unit, mapFmt) itime, time, &
         (pt(i, j75), i=i0, iend(0, 1))
    close(mapTemp75%unit)
    
    ! qtop - - - - - - - - - - - - - - - - - - - - - - -
    iez = iend(0, 2)
    do i=i0, iend(0, 1)
       qtop(i) = (pt(i, iez) - pt(i, iez+1)) * invdz
    end do
    open(unit=mapQtop%unit, file=mapQtop%filename, status='old', position='append')
    write(mapQtop%unit, mapFmt) itime, time, &
         (qtop(i), i=i0, iend(0, 1))
    close(mapQtop%unit)

    nullify(pt)
    
    ! vhsurf - - - - - - - - - - - - - - - - - - -
    call faceToCenter(vel, 0, 1)
    call associatePtrScal(pv, vel, 0, 3)
    open(unit=mapVsurf%unit, file=mapVsurf%filename, status='old', &
         position='append')
    write(mapVsurf%unit, mapFmt) itime, time, &
         (pv(i, iend(0, 2)), i=i0, iend(0, 1))
    close(mapVsurf%unit)
    nullify(pv)

    return
  end subroutine writeMaps
  ! =============================================
end module io
! =========================================================
