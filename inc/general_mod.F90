! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
module general

  implicit none

  integer, public, parameter:: rk=kind(1.0d0)
  integer, parameter:: i0=3 !starting index in grids
  integer, allocatable, dimension(:,:):: iend

  ! geometry ----------------------------------------
  character(len=100):: coordinates
  logical:: curvedGeometry
  integer, public:: nx, nz
  integer, dimension(1:2):: n0, iend0 ! at level 0
  integer, public:: ed !! ed=1 if cylindrical, ed=2 if spherical annulus
  real(rk):: anglemax, lengthX
  real(rk):: d2_cart, r2_cart, invr2_cart  ! d2 = dx^2*dz^2, r2 = dz^2 / dx^2, invr2=1/r2
  character(len=1):: letterGeom
  real(rk):: dx, dz, invdx, invdz 
  real(rk), dimension(1:2):: dd, invdd ! dd(1)=dx, dd(2)=dz
  ! for curved geometries:
  real(rk):: rbot, rtop
  real(rk), allocatable, dimension(:):: r0 ! at level 0
  real(rk), allocatable, dimension(:):: phi0 ! at level 0
  real(rk), allocatable, dimension(:):: alpha0 ! at level 0
  real(rk), allocatable, dimension(:):: extent0 ! at level 0
  real(rk), allocatable, dimension(:):: dv, dh, dvp ! level 0, for ADI
  real(rk), allocatable, dimension(:):: invdv, invdh, invdvp

  ! timer -------------------------------------------
  character(len=100):: endType, writeType
  real(rk):: endValue
  real(rk):: freqLogs, freqFiles, freqMaps
  integer:: nWriteLogs, nWriteFiles, nWriteMaps, endSteps
  real(rk):: timeWriteLogs, timeWriteFiles, timeWriteMaps
  real(rk):: endTime
  
  ! I/O --------------------------------------------
  character(len=100), parameter:: paramFile='param.in'
  character(len=100):: outputStem, fileToRead, fileEtaToRead
  integer:: initNumber
  character(len=100):: beginMethod
  logical:: writeEta, writeMisc, writeExtraDiagnosis, writeProfiles
  
  ! compute/solver ----------------------------------
  logical:: extraRelaxInMG, gauss_seidel
  real(rk):: alphaRelax, errorRes, errorChange
  integer:: nrelaxMin, nrelaxMax, nrelax, counterNRelax, nextIncrease
  character(len=20):: methodEtaEff
  real(rk):: etaEffHarmonicCoeff, alphaVisco, alphaCorrSf
  logical:: retrySolver
  logical:: Fcycles
  integer:: coarsestMode
  logical, parameter:: viscoGeom = .true.

  ! advection/diffusion -----------------
  real(rk):: timeFactor
  character(len=20):: heat_method, TVDscheme_temp
  character(len=20):: TVDscheme_conc

  ! boundary conditions ---------------------------
  real(rk):: topValueT, botValueT, topValueV, botValueV
  real(rk):: leftValueT, rightValueT, leftValueV, rightValueV
  character(len=20):: topT, botT, topV, botV, leftT, rightT, leftV, rightV
  character(len=100):: fileVtop
  real(rk):: kleft, kright
  character(len=20), parameter:: & !define char arrays for easy use of bcs
       top="top", bot="bot", &
       left="left", right="right", &
       dirichlet="dirichlet", neumann="neumann", &
       freeslip="freeslip", &
       secondDerivative="secondDerivative", periodic="periodic"
  character(len=20), parameter:: & ! options
       surface="surface", horizontal="horizontal"
  integer:: nplates
  type plate
     integer:: nleft, nright, npoints
     real(rk):: v, sf
  end type plate
  type(plate), allocatable, dimension(:):: plates

  ! physics --------------------------------------
  real(rk), parameter:: gasConstant = 8.314
  real(rk), parameter:: topTref = 0.0, botTref = 1.0
  real(rk):: rayleigh, internalHeating
  real(rk):: fcurv 
  character(len=20):: viscoLaw, reference
  logical:: isoviscous, yielding, anderson, uniformTemp
  logical:: cutOffMax, cutOffMin, smoothCutOff
  real(rk):: Ea, Va, &
       surfTAdim, thetaVisco, gammaVisco, &
       yieldFix, yieldSurf, slopeYield, etaStar, &
       alphaDt, viscoCutOffMax, viscoCutOffMin, &
       coeffVisco, compressionFactor, TrefAdim, ToffAdim

  ! miscellaneous --------------------------------------------
  integer verbose ! used to define how much output is written
  ! during a run (from 0 to 5... )
  real(rk):: pi
   
  real(rk), parameter:: zero=0.0
  real(rk), parameter:: one=1.0
  real(rk), parameter:: two=2.0
  real(rk), parameter:: four=4.0
  real(rk), parameter:: eight=8.0
  real(rk), parameter:: sixteen=16.0
  real(rk), parameter:: minusone=-1.0
  real(rk), parameter:: half=0.500
  real(rk), parameter:: fourth=0.250
  real(rk), parameter:: eighth=0.125
  real(rk), parameter:: infinity=huge(zero)
  real(rk), parameter:: small=1.0E-30
  real(rk), parameter:: limitStagnantLid=0.05
  
  integer:: j25, j50, j75, j375  ! heights
  
  ! kronecker symbol, its "opposite" and opposit direction
  integer:: kr(1:2, 1:2), nkr(1:2, 1:2), opposit(1:4)
  ! for gauss-seidel biharmonic:
  integer:: gs_idx(0:7, 1:2)
  ! computation
  logical:: relax_on_error
  logical:: relax_corrSf

   ! time
  real(rk):: time, timestep, minTimestep
  real(rk):: nextTimeWriteLogs, nextTimeWriteFiles, nextTimeWriteTracers, &
       nextTimeWriteMaps
  integer:: itime, counter, nLoopTot
  integer:: filenumber
  logical:: running

  ! io
  logical:: startingFromRead
  ! inout diagnosis:
  real(rk), dimension(1:5), parameter:: percValues = (/ 0.10, 0.25, 0.50, 0.75, 0.90 /)
  
  ! physics
  real(rk):: etaTop, etaBot

  ! computation time
  real(rk):: cpuTimeAll, cpuTimeSolver, &
       cpuTimeTemperature, cpuTimeOutput
  real(rk):: t_start, t_end, t_allStart, t_allEnd

  !========================================
  ! NAMELISTS
  !========================================       
  namelist /physics/ rayleigh, internalHeating, fcurv, &
       viscoLaw, Ea, Va, thetaVisco, surfTAdim, ToffAdim, gammaVisco, &
       yielding, yieldFix, yieldSurf, slopeYield, etaStar, &
       alphaDt, viscoCutOffMin, viscoCutOffMax, smoothCutOff, &
       reference, compressionFactor, TrefAdim, uniformTemp

  namelist /geometry/ coordinates, anglemax, &
       lengthX, nx, nz

  namelist /inout/ verbose, outputStem, &
       beginMethod, fileToRead, fileEtaToRead, &
       writeEta, writeMisc, writeExtraDiagnosis, writeProfiles
  
  namelist /boundaries/ topT, botT, topV, botV, &
       leftT, rightT, leftV, rightV, &
       topValueT, botValueT, topValueV, botValueV, &
       leftValueT, rightValueT, leftValueV, rightValueV, &
       fileVtop
  
  namelist /timer/ endType, endValue, &
       writeType, freqLogs, freqFiles, freqMaps

  namelist /compute/ errorRes, errorChange, &
       nrelaxMin, nrelaxMax, alphaRelax, extraRelaxInMG, &
       gauss_seidel, methodEtaEff, &
       etaEffHarmonicCoeff, &
       Fcycles, alphaVisco, alphaCorrSf

  namelist /advdiff/ heat_method, TVDscheme_temp, &
       timeFactor, TVDscheme_conc

  ! -------------------------------------------------------

  !=======================================================
contains
  !-------------------------------------------------------
  subroutine defaultValues()
    implicit none
    integer:: i

    ! geometry -----------------
    coordinates = 'cartesian'  ! cartesian, spherical or cylindrical
    curvedGeometry = .false.
    nx = 512
    nz = 64
    lengthX = 4.0
    anglemax = 360.0 ! if curved, full annulus as default

    ! timer -------------------
    writeType = "steps"
    nWriteLogs = 10
    nWriteFiles = 1000
    timeWriteLogs = 0.001
    timeWriteFiles = 0.01
    
    endType = "steps"
    endSteps = 1000
    endTime = 1.0
    
    ! inout --------------------
    verbose = 1
    outputStem = "iso"
    beginMethod = "continue"
    fileToRead = "temperature_init.dat"
    fileEtaToRead = "eta_init.dat"
    writeEta = .false.
    writeMisc = .false.
    writeExtraDiagnosis = .false.
    writeProfiles = .false.
    
    ! compute and advdiff ------------------
    errorRes = 1.0E-4
    errorChange = 1.0E-2
    nrelaxMin = 2
    nrelaxMax = 10
    alphaRelax = 1.50
    heat_method = "ADI"
    TVDscheme_temp = "minmod"
    TVDscheme_conc = "minmod"
    timeFactor = 2.0
    extraRelaxInMG = .false.
    gauss_seidel = .true.
    methodEtaEff = "min"
    etaEffHarmonicCoeff = 1.0
    retrySolver = .false.
    alphaVisco = 0.5
    alphaCorrSf = 0.5
    Fcycles = .true.
    coarsestMode = 2
    
    ! boundaries ---------------
    topT = "dirichlet"
    botT = "dirichlet"
    topValueT = 0.0
    botValueT = 1.0

    topV = "neumann"
    botV = "neumann"
    topValueV = 0.0
    botValueV = 0.0
    fileVtop = "vtop.dat"
    nplates = 1
    
    leftT = "periodic"
    rightT = "periodic"
    leftValueT = 0.0
    rightValueT = 0.0
    
    leftV = "periodic"
    rightV = "periodic"
    leftValueV = 0.0
    rightValueV = 0.0
    kleft = 0.0
    kright = 0.0
    
    ! physics -----------------
    rayleigh = 1.0E7
    internalHeating = 0.0
    fcurv = 0.5448  ! default: Earth (D=2900km and R=6371km)
    isoviscous = .false.
    anderson = .false.
    uniformTemp = .false.
    viscoLaw = "Arrhenius"
    reference = "bottom"
    Ea = 0.0
    thetaVisco = 50.0
    surfTAdim = 0.2037     ! 275K / 1350K
    ToffAdim = 0.0 
    TrefAdim = 1.0  ! not used unless reference == "T_ref"
    Va = 0.0
    gammaVisco = 0.0
    viscoCutOffMax = -1.0     ! not used if negative value
    viscoCutOffMin = -1.0     ! not used if negative value
    smoothCutOff = .false. 
    yielding = .false.
    yieldFix = -1.0  ! infinity if Neg
    slopeYield = -1.0 ! infinity if Neg
    compressionFactor = -1.0 ! not used if Neg
    yieldSurf = -1.0 ! infinity if Neg
    etaStar = 0.0 ! (effective viscosity at high stresses (like in Stein et al., 2014)
    alphaDt = 0.02
    ! ---------------------------------
    ! miscellaneous -------------------
    ! stuff not read in param.in
    pi = 4.0 * atan(1.0)
    
    ! kroneker delta and its opposite
    kr(:,:) = 0
    nkr(:,:) = 1
    do i=1, 2
       kr(i,i) = 1
       nkr(i,i) = 0
    end do

    ! opposite directions
    opposit(1) = 2
    opposit(2) = 1
    opposit(3) = 4
    opposit(4) = 3

    ! gauss_seidel positions in grid for 13-point stencil:
    gs_idx(0, 1:2) = (/ 0, 0 /)
    gs_idx(1, 1:2) = (/ 1, 0 /)
    gs_idx(2, 1:2) = (/ 0, 1 /)
    gs_idx(3, 1:2) = (/ 1, 2 /)
    gs_idx(4, 1:2) = (/ 2, 1 /)
    gs_idx(5, 1:2) = (/ 1, 1 /)
    gs_idx(6, 1:2) = (/ 1, 3 /)
    gs_idx(7, 1:2) = (/ 0, 2 /)
    
    return
  end subroutine defaultValues
  ! -------------------------------------------------------
end module general
! =======================================================
