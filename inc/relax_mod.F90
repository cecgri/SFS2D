! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
module relax

  use general
  use grids
  use mylib
  implicit none

  abstract interface
     subroutine relax_interface(what, whatrhs, whatres, level, id, resOnly)
       character(len=10), intent(in):: what, whatrhs, whatres
       integer, intent(in):: level, id
       logical, optional:: resOnly
     end subroutine relax_interface
  end interface

  ! =======================================================
contains
  ! =======================================================
  subroutine relax_poisson(what, whatrhs, whatres, level, id, resOnly)
    implicit none
    character(len=10), intent(in):: what, whatrhs, whatres
    integer, intent(in):: level, id
    logical, optional:: resOnly
    logical:: update, redblack
    integer:: i,j, irb, irbmax, istep, ix0, iz0
    real(rk), pointer:: pf(:,:) => null()
    real(rk), pointer:: prhs(:,:) => null()
    real(rk), pointer:: pres(:,:) => null()
    real(rk), pointer:: pp(:) => null() ! poisson
    real(rk), pointer:: pptop(:) => null() ! poissontop
    real(rk), pointer:: ppbot(:) => null() ! poissonbot
    real(rk), pointer:: ppmid(:) => null() ! poissonmid

    update = .true.
    redblack = .true.

    if(present(resOnly)) then
       if(resOnly) then
          update = .false.
          redblack = .false.
       end if
    end if

    istep = 1
    irbmax = 0
    if(redblack) then
       istep = 2
       irbmax = 1
    end if

    call associatePtrScal(pf, what, level, id)
    call associatePtrScal(pres, whatres, level, id)
    call associatePtrScal(prhs, whatrhs, level, id)

    call associatePtrLine(pp, poisson, level)
    call associatePtrLine(pptop, poissontop, level)
    call associatePtrLine(ppbot, poissonbot, level)
    call associatePtrLine(ppmid, poissonmid, level)

    iz0 = i0
    if((id == 4 .or. id == 2)) iz0 = i0+1

    do irb = 0, irbmax

       do j=iz0, iend(level,2)

          ix0 = i0
          if(redblack) &
               ix0 = i0 + mod(j+irb+1, 2)

          do i=ix0, iend(level,1), istep
             pres(i,j) = prhs(i,j) + &
                  ppmid(j) * (pf(i+1,j) + pf(i-1,j)) + &
                  pptop(j) * pf(i,j+1) + ppbot(j) * pf(i,j-1) &
                  - pp(j) * pf(i,j)

             if(redblack .and. update) then
                pf(i,j) = pf(i,j) + &
                     pres(i,j) / pp(j) * alphaRelax
             end if
          end do !i

       end do !j

       if(update) then
          call applyBcs(what, level, id)
       end if

    end do ! irb

    if(update .and. .not.redblack) then
       do j=iz0, iend(level,2)
          do i=i0, iend(level,1)
             pf(i,j) = pf(i,j) + &
                  pres(i,j) / pp(j) * alphaRelax
          end do
       end do
       call applyBcs(what, level, id)
    end if

    nullify(pf, prhs, pres, pp, pptop, ppbot, ppmid)

    return
  end subroutine relax_poisson
  ! =======================================================
  subroutine relax_biharmonic(what, whatrhs, whatres, level, id, resOnly)
    implicit none
    character(len=10), intent(in):: what, whatrhs, whatres
    integer, intent(in):: level, id
    logical, optional:: resOnly
    logical:: update, gauss
    real(rk), pointer:: pf(:,:) => null()
    real(rk), pointer:: prhs(:,:) => null()
    real(rk), pointer:: pres(:,:) => null()
    real(rk), pointer:: pc(:, :, :) => null()
    real(rk), dimension(1:13):: p
    integer:: i, j, ic, ix0, iz0, ixg0, ngs, nngs, igs(0:7)
    real(rk):: sum
    !------ stencil ---------
    !            12
    !            |
    !       7----4----6
    !       |    |    |
    !  11---3----1----2---10
    !       |    |    |
    !       9----5----8
    !            |
    !            13
    !-----------------------  

    update = .true.
    gauss = gauss_seidel
    if(present(resOnly)) then
       if(resOnly) then
          update = .false.
          gauss = .false.
       end if
    end if

    if(id /= 4) then
       write(*,'(a)') "ERROR: relax_biharmonic only for id=4"
       call endRun()
    end if

    ! TEST (should not be needed)
    call zeroingScal(whatres, level, id) 

    call associatePtrScal(pf, what, level, id)
    call associatePtrScal(pres, whatres, level, id)
    call associatePtrScal(prhs, whatrhs, level, id)
    if(.not.relax_corrSf) then ! so far relax_corrSf = .false. always
       call associatePtrVec(pc, cbh, level) ! (13, :, :)
    else
       call associatePtrVec(pc, cbhy, level)
    end if

    iz0 = i0 + 1
    ix0 = i0
    if(leftV /= periodic) then
       ix0 = i0 + 1
    end if

    if(gauss .and. update) then

       igs = (/ 0, 7, 5, 6, 2, 4, 1, 3 /)
       do nngs = 0, 7
          ngs = igs(nngs)
          do j = iz0 + mod(gs_idx(ngs, 2), 2), iend(level, 2), 2
             ixg0 = ix0 + mod((gs_idx(ngs, 1) + mod((j-i0 - gs_idx(ngs, 2)), 4)), 4)
             do i = ixg0, iend(level, 1), 4
                p(1) = pf(i, j)
                p(2) = pf(i+1, j)
                p(3) = pf(i-1, j)
                p(4) = pf(i, j+1)
                p(5) = pf(i, j-1)
                p(6) = pf(i+1, j+1)
                p(7) = pf(i-1, j+1)
                p(8) = pf(i+1, j-1)
                p(9) = pf(i-1, j-1)
                p(10) = pf(i+2, j)
                p(11) = pf(i-2, j)
                p(12) = pf(i, j+2)
                p(13) = pf(i, j-2)
                sum = 0.0
                do ic = 1, 13
                   sum = sum + pc(ic, i, j) * p(ic)
                end do ! ic=1, 13

                pres(i, j) = prhs(i, j) + sum

                pf(i, j) = pf(i, j) - pres(i, j) / pc(1, i, j) * alphaRelax

             end do ! i
          end do ! j

          call applyBcs(what, level, id)

       end do ! nngs

    else !! Jacobi relaxations or computes residue only

       do j=iz0, iend(level, 2)
          do i=ix0, iend(level, 1)
             p(1) = pf(i, j)
             p(2) = pf(i+1, j)
             p(3) = pf(i-1, j)
             p(4) = pf(i, j+1)
             p(5) = pf(i, j-1)
             p(6) = pf(i+1, j+1)
             p(7) = pf(i-1, j+1)
             p(8) = pf(i+1, j-1)
             p(9) = pf(i-1, j-1)
             p(10) = pf(i+2, j)
             p(11) = pf(i-2, j)
             p(12) = pf(i, j+2)
             p(13) = pf(i, j-2)
             sum = 0.0
             do ic=1, 13
                sum = sum + pc(ic, i, j) * p(ic)
             end do

             pres(i, j) = prhs(i, j) + sum
          end do ! i
       end do ! j

       if(update) then
          do j=iz0, iend(level, 2)
             do i=ix0, iend(level, 1)
                pf(i, j) = pf(i, j) - pres(i, j) / pc(1, i, j) * alphaRelax
             end do ! i
          end do ! j

          call applyBcs(what, level, id)
       end if ! update

    end if ! Gauss-Seidel or Jacobi

    nullify(pf, prhs, pres, pc)

    return
  end subroutine relax_biharmonic
  ! =======================================================
  subroutine solver(relax, what, whatrhs, whatres, id, limitRes, &
       converge, multigrid, givenRms, levelMax)
    use io

    implicit none
    procedure(relax_interface):: relax
    character(len=10), intent(in):: what, whatrhs, whatres
    integer, intent(in):: id
    real(rk), intent(in):: limitRes
    logical, intent(out):: converge
    logical, optional:: multigrid
    real(rk), optional:: givenRms
    integer, optional:: levelMax
    integer:: iloop, nrel
    integer, parameter:: nloopmax=100
    real(rk):: rmsRes, rms, resPrec(0:1), rmsChange
    character(len=10), parameter:: whatFirst=work2, whatprec=work3
    logical:: firstTry, doMG, fixRms, diverge

    call copyScal(whatfirst, id, what, id, 0)
    call copyScal(whatprec, id, what, id, 0)
    firstTry = (nrelax /= nrelaxMax)
    converge = .false.

    doMG = .true.
    if(present(multigrid)) doMG = multigrid

    fixRms = .false.
    if(present(givenRms)) then
       fixRms = .true.
       rms = givenRms
    end if

    resPrec(:) = 0.0
    nextIncrease = 10

1   iloop = 0
    do while(iloop <= nloopmax)       
2      do nrel=1, nrelax
          call relax(what, whatrhs, whatres, 0, id)
       end do
       call relax(what, whatrhs, whatres, 0, id, resOnly=.true.)
       call applyBcs(whatres, 0, id)

       if(iloop > 0) then ! do at least one loop before checking
          resPrec(mod(iloop, 2)) = rmsRes

          if(.not.fixRms) then
             call getRmsSimple(rms, what, 0, id)
          end if
          call residue_biharmonic(rmsRes, what, whatrhs, whatres, id)
          call getRmsChange(rmsChange, what, whatprec, 0, 4)

          diverge = (rmsRes > 0.5 * (resPrec(0) + resPrec(1)))
          if(what == sf) then
             converge = (rmsRes < limitRes * rms .and. &
                  rmsChange < errorChange * rms)
          else  ! corrSf
             converge = (rmsRes < 0.1 * limitRes * rms)
          end if

          if(verbose >= 5 .or. (verbose >= 4 .and. converge)) then
             write(*,'(3a, 2i5, 3e16.6, l4)') &
                  "          relax ", trim(adjustl(what)), " : ", iloop, &
                  nrelax, rms, rmsRes / rms, rmsChange / rms, converge
          end if

          ! adapt nrelax if iloop very high or low
          if(converge .and. iloop < 4 .and. counterNRelax > 20) then
             nrelax = nrelax - 2
             nextIncrease = max(10, nextIncrease - 10)
             if(nrelax < nrelaxMin) then
                nrelax = nrelaxMin
             end if
             counterNRelax = 0
          else if(iloop > nextIncrease) then
             nrelax = min(nrelaxMax, nrelax + 2)
             counterNRelax = 0
             nextIncrease = nextIncrease + 10
             if(diverge) then
                write(*, '(a, i4, 3e14.4)') &
                     "warning: diverge -> reset  ", nrelax, &
                     rmsRes, resPrec(0), resPrec(1)
                call copyScal(what, id, whatFirst, id, 0)
                go to 2
             end if
          end if
          if(nrelax > nrelaxMin) counterNRelax = counterNRelax + 1

          if(converge) exit
       end if ! iloop > 0

       ! --- MULTIGRID ---
       if(doMG) then
          if(Fcycles) then
             call Fcycle(relax, what, whatres, id)
          else
             if(present(levelMax)) then
                call Vcycle(relax, what, whatres, id, levelMax=levelMax)
             else
                call Vcycle(relax, what, whatres, id)
             end if
          end if
       end if
       ! -----------------
       call copyScal(whatprec, id, what, id, 0)

       iloop = iloop + 1

    end do ! iloop

    call applyBcs(what, 0, id) ! useless?

    nLoopTot = nLoopTot + iloop

    if(.not.converge) then
       call getRmsSimple(rms, what, 0, id)
       write(*,'(i6, i4, 2a, i3, 3e16.6)') itime, nrelax, &
            "  WARNING: solver did not converge   ", &
            trim(what), id, rms, rmsRes / rms, rmsChange / rms

       if(firstTry) then
          ! increase nrelax to nrelaxMax
          firstTry = .false.
          nrelax = nrelaxMax
          call copyScal(what, id, whatFirst, id, 0)
          write(*, '(a, i5)') "RETRY with nrelax = ", nrelax
          go to 1
       end if

       if((nrelax == nrelaxMax .and. rmsRes > 1.0E4 * limitRes) .or. &
            isnan(rmsRes)) then
          write(*,'(a)') &
               "rmsRes too large or NaN -> end run."
          call endRun()
       end if

    end if ! not converge

    return
  end subroutine solver
  ! =======================================================
  subroutine Vcycle(relax, what, whatres, id, levelMax)
    use io
    implicit none
    procedure(relax_interface):: relax
    character(len=10), intent(in):: what, whatres
    integer, intent(in):: id
    integer, optional:: levelMax
    integer:: nlevelMax
    integer:: level, nrel, nrelaxMG
    character(len=10), parameter:: prolong=work1
    real(rk), parameter:: restrictionFactor = 1.0 !! =16 in Cartesian - TEST !!!

    nLevelMax = nLevel
    if(present(levelMax)) then
       nLevelMax = min(nLevel, levelMax)
    end if

    relax_on_error = .true.  ! useless in present version (?)

    call zeroingScal(error, 0, id)

    call copyScal(errorRhs, id, whatres, id, 0)
    do nrel=1, nrelax
       call relax(error, errorRhs, whatres, 0, id)
    end do
    call relax(error, errorRhs, whatres, 0, id, resOnly=.true.)
    call applyBcs(whatres, 0, id)

    ! descending - - - - - - - - - - - - - -
    do level=1, nLevelMax

       call zeroingScal(error, level, id)
       call restriction(whatres, level, id) !from level-1 to level
       call copyScal(errorRhs, id, whatres, id, level, mult=restrictionFactor)

       if(level < nLevelMax) then
          nrelaxMG = nrelax
          if(extraRelaxInMG) nrelaxMG = nrelax * level

          do nrel=1, nrelaxMG
             call relax(error, errorRhs, whatres, level, id)
          end do
          call relax(error, errorRhs, whatres, level, id, resOnly=.true.)
          call applyBcs(whatres, level, id)

       end if

    end do !level

    if(nLevelMax == nLevel) then
       ! exact solution at coarsest level (nz=2)
       call exactCoarsest(error, errorRhs, id)
    else
       ! do extra relaxations at the coarsest level
       do nrel=1, nrelax * 10
          call relax(error, errorRhs, whatres, nLevelMax, id)
       end do
    end if

    ! climbing back up - - - - - - - - - - -
    do level=nLevelMax, 1, -1

       call prolongation(prolong, error, level, id) ! from level to level-1
       call addScal(error, prolong, level-1, id)
       call applyBcs(error, level-1, id)

       do nrel=1, nrelax
          call relax(error, errorRhs, whatres, level-1, id)
       end do

    end do !level

    ! add error to solution at level 0
    call addScal(what, error, 0, id)
    call applyBcs(what, 0, id)

    relax_on_error = .false. ! useless in present version

    return
  end subroutine Vcycle
  ! =======================================================
  subroutine Fcycle(relax, what, whatres, id)
    implicit none
    procedure(relax_interface):: relax
    character(len=10), intent(in):: what, whatres
    integer, intent(in):: id
    integer:: level, nrel
    character(len=10), parameter:: prolong=work1
    integer:: nFine
    real(rk), parameter:: restrictionFactor = 16.0

    relax_on_error = .true.

    call zeroingScal(error, 0, id)

    call copyScal(errorRhs, id, whatres, id, 0)
    do nrel=1, nrelax
       call relax(error, errorRhs, whatres, 0, id)
    end do
    call relax(error, errorRhs, whatres, 0, id, resOnly=.true.)
    call applyBcs(whatres, 0, id)

    ! descending ==========================================
    do level=1, nLevel

       call zeroingScal(error, level, id)
       call restriction(whatres, level, id) !from level-1 to level
       call copyScal(errorRhs, id, whatres, id, level, mult=restrictionFactor)

       if(level < nLevel) then
          do nrel=1, nrelax
             call relax(error, errorRhs, whatres, level, id)
          end do
          call relax(error, errorRhs, whatres, level, id, resOnly=.true.)
          call applyBcs(whatres, level, id)
       end if

    end do !level

    ! exact solution at coarsest level (nz=2) =========
    call exactCoarsest(error, errorRhs, id) 

    ! climbing back up in several stages ======================
    do nFine = nLevel-1, 0, -1
       ! climbing back until nFine - - - - - - - - 
       do level=nLevel, nFine+1, -1

          call prolongation(prolong, error, level, id) ! from level to level-1
          call addScal(error, prolong, level-1, id)
          call applyBcs(error, level-1, id)

          do nrel=1, nrelax
             call relax(error, errorRhs, whatres, level-1, id)
          end do

       end do ! level=nLevel, nFine+1, -1

       ! descending again until nLevel - - - - - -
       do level=nFine, nLevel-1
          call relax(error, errorRhs, whatres, level, id, resOnly=.true.)
          call applyBcs(whatres, level, id)

          call restriction(whatres, level+1, id) ! from level to level+1
          call copyScal(errorRhs, id, whatres, id, level+1, mult=restrictionFactor)

          call zeroingScal(error, level+1, id) ! better without?
          if(level+1 < nLevel) then 
             do nrel=1, nrelax
                call relax(error, errorRhs, whatres, level+1, id)
             end do
          else
             ! exact solution at coarsest level (nz=2) =========
             call exactCoarsest(error, errorRhs, id) 
          end if

       end do ! level=nFine, nLevel

    end do ! nFine = nLevel-1, 1, -1

    ! Final climbing back, until level = 0 =========
    do level=nLevel, 1, -1
       call prolongation(prolong, error, level, id)
       call addScal(error, prolong, level-1, id)
       call applyBcs(error, level-1, id)

       do nrel=1, nrelax
          call relax(error, errorRhs, whatres, level-1, id)
       end do

    end do ! level

    ! add error to solution at level 0
    call addScal(what, error, 0, id)
    call applyBcs(what, 0, id)

    relax_on_error = .false.

    return
  end subroutine Fcycle
  ! =======================================================
  subroutine residue_biharmonic(rmsRes, what, whatrhs, whatres, id, do_relax)
    implicit none
    character(len=10), intent(in):: what, whatrhs, whatres
    integer, intent(in):: id
    real(rk), intent(out):: rmsRes
    logical, optional:: do_relax
    character(len=10), parameter:: res_on_C1 = work4 ! TEMPO
    real(rk), pointer:: pc(:,:,:) => null()
    real(rk), pointer:: pres(:,:) => null()
    real(rk), pointer:: presdiv(:,:) => null()
    integer:: i, j
    logical:: divideC1

    if(present(do_relax)) then
       if(do_relax) &
            call relax_biharmonic(what, whatrhs, whatres, 0, id, resOnly=.true.)
    end if

    divideC1 = .true.  ! TEST

    if(.not.divideC1) then
       call getRmsSimple(rmsRes, whatres, 0, id)
    else
       if(.not.relax_corrSf) then  ! so far relax_corrSf = .false. always
          call associatePtrVec(pc, cbh, 0)
       else
          call associatePtrVec(pc, cbhY, 0)
       end if
       call associatePtrScal(pres, res, 0, id)
       call associatePtrScal(presdiv, res_on_C1, 0, id)
       do j=i0, iend(0, 2) + 1
          do i=i0, iend(0, 1) + 1
             if(pc(1, i, j) > 1.0) then
                presdiv(i, j) = pres(i, j) / pc(1, i, j)
             else
                presdiv(i, j) = pres(i, j)
             end if
          end do
       end do
       nullify(pc, pres, presdiv)
       call getRmsSimple(rmsRes, res_on_C1, 0, id)
    end if

    return
  end subroutine residue_biharmonic
  ! =============================================
  subroutine exactCoarsest(what, whatrhs, id)
    ! boundary conditions are included
    ! WARNING: parallel not taken into account at all yet
    implicit none
    character(len=10), intent(in):: what, whatrhs
    integer, intent(in):: id
    real(rk), pointer:: prhs(:,:) => null()
    real(rk), pointer:: pf(:,:) => null()
    integer:: nval, iez, iex

    iex = iend(nLevel, 1)
    iez = iend(nLevel, 2) ! should be = 4 for i0=3
    nval = coarsestMode ! for clarity

    call associatePtrScal(prhs, whatrhs, nLevel, id)
    call associatePtrScal(pf, what, nLevel, id)

    if(coarsestMode == 1) then
       pf(i0+1, iez) = -prhs(i0+1, iez) / mcc(1, 1)
    else if(coarsestMode == 2) then
       ycc(1:2) = -prhs(i0:i0+1, iez)
       call solve2x2(xcc, &
            mcc(1,1), mcc(1,2), &
            mcc(2,1), mcc(2,2), ycc)
       pf(i0:i0+1, iez) = xcc(1:2)
    else if(coarsestMode == 3) then
       ycc(1:3) = -prhs(i0+1:i0+3, iez)
       call solve3x3(xcc, &
            mcc(1,1), mcc(1,2), mcc(1,3), &
            mcc(2,1), mcc(2,2), mcc(2,3), &
            mcc(3,1), mcc(3,2), mcc(3,3), ycc)
       pf(i0+1:i0+3, iez) = xcc(1:3)
    else if(coarsestMode == 4) then
       ycc(1:4) = -prhs(i0:i0+3, iez)
       call solve4x4(xcc, &
            mcc(1,1), mcc(1,2), mcc(1,3), mcc(1,4), &
            mcc(2,1), mcc(2,2), mcc(2,3), mcc(2,4), &
            mcc(3,1), mcc(3,2), mcc(3,3), mcc(3,4), &
            mcc(4,1), mcc(4,2), mcc(4,3), mcc(4,4), ycc)
       pf(i0:i0+3, iez) = xcc(1:4)
    else ! pentadiagonal system
       if(mod(coarsestMode, 2) == 0) then
          ! leftV=rightV='periodic' -> nearly-pentadiagonal system
          ycc(1:nval) = -prhs(i0:i0+nval-1, iez)
          call nearlyPentadiag(xcc, &
               mcc(:,1), mcc(:,2), mcc(:,3), mcc(:,4), mcc(:,5), ycc, nval)
          pf(i0:i0+nval-1, iez) = xcc(1:nval)
       else
          ! --> pentadiagonal system
          ycc(1:nval) = -prhs(i0+1:i0+nval, iez)
          call pentadiag(xcc, &
               mcc(:,1), mcc(:,2), mcc(:,3), mcc(:,4), mcc(:,5), ycc, nval)
          pf(i0+1:i0+nval, iez) = xcc(1:nval)
       end if
    end if

    ! Horizontal bcs, top and bottom ------------------
    pf(i0:iex, i0) = 0.0
    pf(i0:iex, iez+1) = 0.0
    if(topV == dirichlet .or. topV == "file") then  ! (? for "file")
       pf(i0:iex, iez+2) = pf(i0:iex, iez) * grid(nLevel)%topdirichlet(1) * &
            grid(nLevel)%topdirichlet(2)
    else ! neumann 
       pf(i0:iex, iez+2) = -pf(i0:iex, iez) * grid(nLevel)%topfreeslip
    end if
    if(botV == dirichlet) then
       pf(i0:iex, i0-1) = pf(i0:iex, i0+1) * grid(nLevel)%botdirichlet(1) * &
            grid(nLevel)%botdirichlet(2)
    else ! neumann
       pf(i0:iex, i0-1) = -pf(i0:iex, i0+1) * grid(nLevel)%botfreeslip
    end if

    ! Vertical bcs, left and right ---------------
    if(leftV == periodic) then ! rightV == periodic also then
       pf(i0-2, :) = pf(iex-1, :)
       pf(i0-1, :) = pf(iex, :)
       pf(iex+1, :) = pf(i0, :)
       pf(iex+2, :) = pf(i0+1, :)
    else  ! not periodic
       pf(i0, :) = 0.0
       if(leftV == dirichlet) then
          pf(i0-1, :) = pf(i0+1, :)
       else if(leftV == neumann) then
          pf(i0-1, :) = -pf(i0+1, :)
       else 
          pf(i0-1, :) = -pf(i0+1, :)
       end if

       pf(iex+1, :) = 0.0
       if(rightV == dirichlet) then
          pf(iex+2, :) = pf(iex, :)
       else if(rightV == neumann) then
          pf(iex+2, :) = -pf(iex, :)
       else
          pf(iex+2, :) = -pf(iex, :)
       end if
    end if

    nullify(pf, prhs)

    return
  end subroutine exactCoarsest
  ! =============================================
end module relax
