! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
module grids

  use general
  use mylib
  implicit none

  integer:: nLevel  ! number of levels of multigrid

  ! cell: id is 1: x-face, 2: z-face, 3: cell center, 4: cell corner (vertex)
  type celltype
     integer:: n(1:2), ntot
     real(rk):: dx, dz
     real(rk), allocatable:: pos(:,:,:) ! positions of nodes
     real(rk), allocatable:: r(:), alpha(:), extent(:), phi(:)
     ! GEOMETRICAL PARAMETERS --------
     ! for solvers (for curved geometries)
     real(rk), allocatable:: gtop(:), gbot(:), gmid(:), g1top(:), g1bot(:)
     real(rk), allocatable:: gcentertop(:), gcenterbot(:)
     real(rk), allocatable:: gcornertop(:), gcornerbot(:)
     real(rk), allocatable:: g4(:), g5(:), g12(:), g13(:)
     real(rk), allocatable:: poisson(:), poissontop(:)
     real(rk), allocatable:: poissonbot(:), poissonmid(:)
     ! for boundary conditions
     real(rk):: topfreeslip, botfreeslip, topdirichlet(1:2), botdirichlet(1:2)
     ! PHYSICAL PARAMETERS ----------
     real(rk), allocatable:: temp(:,:,:)
     real(rk), allocatable:: etaT(:,:,:)
     real(rk), allocatable:: etaEff(:,:,:)
     real(rk), allocatable:: etaPrec(:,:,:)
     real(rk), allocatable:: eta(:,:,:)  ! used in the solver
     real(rk), allocatable:: sf(:,:,:) ! streamlines
     real(rk), allocatable:: corrSf(:,:,:)
     real(rk), allocatable:: sfPrec(:,:,:)
     real(rk), allocatable:: k(:,:,:)
     real(rk), allocatable:: cpSf(:,:,:)
     real(rk), allocatable:: vel(:,:,:) !1: vh, 2: vr, 3:p
     real(rk), allocatable:: rhs(:,:,:)
     real(rk), allocatable:: errorRhs(:,:,:)
     real(rk), allocatable:: error(:,:,:)
     real(rk), allocatable:: corrRhs(:,:,:)
     real(rk), allocatable:: res(:,:,:)
     real(rk), allocatable:: omega(:,:,:) ! vorticity
     real(rk), allocatable:: exz(:,:,:) !4 only needed
     real(rk), allocatable:: exx(:,:,:) !3,4 only needed     
     real(rk), allocatable:: ezz(:,:,:) !3,4 only needed
     real(rk), allocatable:: sist(:,:,:) !Second Invariant Strainrate Tensor
     real(rk), allocatable:: work1(:,:,:) ! temporary
     real(rk), allocatable:: work2(:,:,:) 
     real(rk), allocatable:: work3(:,:,:)  
     real(rk), allocatable:: work4(:,:,:)  
     real(rk), allocatable:: work5(:,:,:)  
     real(rk), allocatable:: coeff(:,:,:)
     real(rk), allocatable:: cbh(:,:,:) ! coeff for biharmonic (13 points stencil)
     real(rk), allocatable:: cbhY(:,:,:) ! coeff for yielding (newton method??)
     real(rk), allocatable:: cbhPrec(:,:,:)
  end type celltype
  
  type(celltype), allocatable, target:: grid(:)

  ! define character chains, to simplify the finding of objects in grids
  character(len=10), parameter:: pos="pos", &
       vel="vel", rhs="rhs", error="error", errorRhs="errorRhs", res="res", &
       temp="temp", etaT="etaT", etaEff="etaEff", eta="eta", &
       sf="sf", omega="omega", &
       exx="exx", ezz="ezz", exz="exz", &
       work1="work1", work2="work2", &
       work3="work3", work4="work4", &
       work5="work5", &
       sist="sist", coeff="coeff", &
       etaPrec="etaPrec", &
       corrSf="corrSf", corrRhs="corrRhs", k="k", &
       cpSf="cpSf", cbh="cbh", cbhY="cbhY", cbhPrec="cbhPrec", &
       sfPrec="sfPrec", &
       extent='extent', alpha='alpha', r='r', phi='phi', &
       gtop='gtop', gbot='gbot', &
       gmid='gmid', gcentertop='gcentertop', gcenterbot='gcenterbot', &
       g1top='g1top', g1bot='g1bot', gcornertop='gcornertop', &
       gcornerbot='gcornerbot', g4='g4', g5='g5', g12='g12', g13='g13', &
       poisson="poisson", poissontop="poissontop", poissonbot="poissonbot", &
       poissonmid="poissonmid"

  ! coeff for ADI advection of temperature:
  real(rk), dimension(:,:,:), allocatable:: flim_m, flim_p, &
       cadv_m, cadv_c, cadv_p
  
  ! for multigrid coarsest level:
  real(rk), dimension(:,:), allocatable:: mcc ! Matrix of Coarsest Coefficients
  real(rk), dimension(:), allocatable:: coarsestGamma
  real(rk), dimension(:), allocatable:: xcc, ycc ! Solution and RHS on Coarsest Grid

  ! yielding:
  real(rk), dimension(:,:,:), allocatable:: etayield
  ! lithostatic pressure:
  real(rk), dimension(:,:), allocatable:: plitho
  ! yield as a function of depth:
  real(rk), dimension(:,:), allocatable:: yieldT_z, yieldC_z
  ! coeffs for the effect of pressure on the viscosity:
  real(rk), dimension(:,:), allocatable:: coeffViscoFK_z, coeffViscoArr_z
  ! any horizontal or vertical vector at level 0
  ! (first index is to identify several profiles,
  !  second is the coordinate) 
  real(rk), dimension(:,:), allocatable:: xval, zval
  ! 2D scalar
  real(rk), dimension(:,:), allocatable:: scalar
  ! 2d vector
  real(rk), dimension(:,:,:), allocatable:: vector
  ! top velocity when topV = 'file' or 'ridge'
  real(rk), dimension(:), allocatable:: vxtop, sftop
  ! bot and side velocity when topV = 'ridge'
  real(rk), dimension(:), allocatable:: vztop, vzbot
  real(rk), dimension(:), allocatable:: vxbot, sfbot
  real(rk), dimension(:), allocatable:: vxleft, sfleft
  real(rk), dimension(:), allocatable:: vxright, sfright
  real(rk), dimension(:), allocatable:: vzleft, vzright
  ! others
  real(rk), dimension(:,:,:), allocatable:: ktherm
  real(rk), dimension(:,:), allocatable:: heat, map
  
  ! for maps (heat flux)
  real(rk), dimension(:), allocatable:: qtop
! =======================================================
contains
  ! -------------------------------------------------------
  subroutine initGrids()
    ! allocate all the fields in grid() at all possible levels of multigrid (MG)
    implicit none
    integer:: ntempo(1:2)
    integer:: nl, idir
    integer:: i, j
    real(rk):: dxtempo, dztempo
    
    ! general grid
    iend0(:) = n0(:) + i0 - 1  ! i0=3 (points 1 and 2 are ghostpoints)

    ! for integrals (mean, rms... )
    allocate(extent0(1:iend0(2)+2))
    ! for ADI:
    allocate(dh(1:iend0(2)+2))
    allocate(dv(1:iend0(2)+2))
    allocate(dvp(1:iend0(2)+2))
    allocate(invdh(1:iend0(2)+2))
    allocate(invdv(1:iend0(2)+2))
    allocate(invdvp(1:iend0(2)+2))

    if(curvedGeometry) then ! _________________________________
       ! for geometry:
       allocate(r0(1:iend0(2)+2))
       allocate(alpha0(1:iend0(2)+2))
       allocate(phi0(1:iend0(1)+2))

       do j=i0-2, iend0(2) + 2
          r0(j) = (j-i0) * dz + rbot
          alpha0(j) = half * dz / r0(j)
          ! coeff for ADI and mpdata:
          dh(j) = r0(j) * (1.0 + alpha0(j)) * dx
          invdh(j) = 1.0 / dh(j)
          dv(j) = (1.0 + alpha0(j))**ed * dz
          invdv(j) = 1.0 / dv(j)
          dvp(j) = dz / (1.0 + 2.0 * alpha0(j))**ed
          invdvp(j) = 1.0 / dvp(j)
       end do ! j

       do i=i0-2, iend0(1)+2
          phi0(i) = (i-i0) * dx ! not used?
       end do

       if(ed == 1) then  ! (cylindrical)
          do j=i0-2, iend0(2) + 2
             extent0(j) = 2.0 * r0(j) * (1.0 + alpha0(j)) * dz / &
                  (rtop**2 * (1.0 - fcurv**2))
          end do
       else !! ed == 2 (spherical annulus)
          do j=i0-2, iend0(2) + 2
             extent0(j) = r0(j)**2 * &
                  (3.0 + 6.0 * alpha0(j) + 4.0 * alpha0(j)**2) * dz / &
                  (rtop**3 * (1.0 - fcurv**3))
          end do
       end if
    else !! CARTESIAN _________________________________________
       ! still define extent0, dh, dv etc. for simple unified writing elsewhere
       extent0(:) = dz
       dh(:) = dx
       dv(:) = dz
       invdh(:) = invdx
       invdv(:) = invdz
       dvp(:) = dz
       invdvp(:) = invdz
    end if ! curvedGeometry or not _______________________________

    ! setup multigrid ------------------------------
    nLevel = 0
    ntempo = n0
    do 
       ntempo = ntempo / 2
       nLevel = nLevel + 1
       if(ntempo(2) < 2) then
          exit
       end if
    end do
    nLevel = nLevel - 1

    ! allocate different levels of MG
    allocate(grid(0:nLevel))

    ! put back ntempo to finest grid values
    ntempo(:) = n0(:)
    ! size of cells at level 0
    dxtempo = dx
    dztempo = dz
    
    write(*,'(a)') "========================================"

    do nl = 0, nLevel
       
       grid(nl)%n(:) = ntempo(:)
       grid(nl)%dx = dxtempo
       grid(nl)%dz = dztempo

       allocate(grid(nl)%extent(1:ntempo(2)+4))
       
       if(curvedGeometry) then  ! ---------------

          allocate(grid(nl)%r(1:ntempo(2)+4))
          allocate(grid(nl)%alpha(1:ntempo(2)+4))
          allocate(grid(nl)%phi(1:ntempo(1)+4))

          ! compute positions and geometrical parameters of each node
          do j=i0-2, i0 + grid(nl)%n(2) + 1
             grid(nl)%r(j) = ( j - i0 ) * dztempo + rbot
             grid(nl)%alpha(j) = half * dztempo / grid(nl)%r(j)
          end do ! j

          ! position of each point:
          do i=i0-2, i0 + grid(nl)%n(1) + 1
             grid(nl)%phi(i) = ( i-i0 ) * dxtempo ! useful?
          end do ! i         
          
          ! relative surface area of each cell
          if (ed == 1) then ! (cylindrical)
             do j=i0-2, i0 + grid(nl)%n(2) + 1
                grid(nl)%extent(j) = 2.0 * grid(nl)%r(j) * &
                     (1.0 + grid(nl)%alpha(j)) * dztempo / &
                     (rtop**2 * (1.0 - fcurv**2))
             end do ! j
          else ! ed==2 (spherical annulus)
             do j=i0-2, i0 + grid(nl)%n(2) + 1
                grid(nl)%extent(j) = grid(nl)%r(j)**2 * &
                     (3.0 + 6.0 * grid(nl)%alpha(j) + 4.0 * grid(nl)%alpha(j)**2) * &
                     dztempo / (rtop**3 * (1.0 - fcurv**3))
             end do ! j
          end if

       else ! CARTESIAN 
          grid(nl)%extent(:) = dztempo
       end if ! curvedGeometry ------------------
       
       write(*,'(a, i2, a, 2i5)') "MG level:", nl, &
            "   nx, nz = ", (grid(nl)%n(idir), idir=1, 2)
       write(*,'(a, 2f10.5)')     "                dx, dz = ", dxtempo, dztempo

       ! all arrays are allocated from 1 to n+4 
       ! 1 - 2 - [ 3 - 4 - ...- n+1 - n+2 ] - n+3 - n+4
       ! Two ghost points on each side (1,2 and n+3,n+4)
       ! Real grid goes from 3 to n+2.
       ! Notations: i0 = 3 and iend = n+2
       allocate(grid(nl)%pos(1:ntempo(1)+4, 1:ntempo(2)+4, 1:2))
       allocate(grid(nl)%temp(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%etaT(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%etaEff(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%etaPrec(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%eta(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%corrSf(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%sfPrec(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%corrRhs(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%k(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%cpSf(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%sf(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%omega(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%vel(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%rhs(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%errorRhs(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%error(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%res(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%exx(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%ezz(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%exz(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%work1(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%work2(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%work3(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%work4(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%work5(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%coeff(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       allocate(grid(nl)%sist(1:ntempo(1)+4, 1:ntempo(2)+4, 1:4))
       !! SPECIAL FORMATS:
       allocate(grid(nl)%gtop(1:ntempo(2)+4))
       allocate(grid(nl)%gbot(1:ntempo(2)+4))
       allocate(grid(nl)%gmid(1:ntempo(2)+4))
       allocate(grid(nl)%g1top(1:ntempo(2)+4))
       allocate(grid(nl)%g1bot(1:ntempo(2)+4))
       allocate(grid(nl)%gcentertop(1:ntempo(2)+4))
       allocate(grid(nl)%gcenterbot(1:ntempo(2)+4))
       allocate(grid(nl)%gcornertop(1:ntempo(2)+4))
       allocate(grid(nl)%gcornerbot(1:ntempo(2)+4))
       allocate(grid(nl)%g4(1:ntempo(2)+4))
       allocate(grid(nl)%g5(1:ntempo(2)+4))
       allocate(grid(nl)%g12(1:ntempo(2)+4))
       allocate(grid(nl)%g13(1:ntempo(2)+4))
       allocate(grid(nl)%poisson(1:ntempo(2)+4))
       allocate(grid(nl)%poissontop(1:ntempo(2)+4))
       allocate(grid(nl)%poissonbot(1:ntempo(2)+4))
       allocate(grid(nl)%poissonmid(1:ntempo(2)+4))
       allocate(grid(nl)%cbh(1:13, 1:ntempo(1)+4, 1:ntempo(2)+4))
       allocate(grid(nl)%cbhY(1:13, 1:ntempo(1)+4, 1:ntempo(2)+4))
       allocate(grid(nl)%cbhPrec(1:13, 1:ntempo(1)+4, 1:ntempo(2)+4))
       call zeroingAll(nl) ! put all to zero at level=nl

       if(curvedGeometry) then
          do j=i0-2, i0 + grid(nl)%n(2) + 1
             do i=i0-2, i0 + grid(nl)%n(1) + 1
                grid(nl)%pos(i,j,1) = grid(nl)%phi(i)
                grid(nl)%pos(i,j,2) = grid(nl)%r(j)
             end do ! i
          end do ! j
       else
          do j=i0-2, i0 + grid(nl)%n(2) + 1
             do i=i0-2, i0 + grid(nl)%n(1) + 1
                grid(nl)%pos(i,j,1) = ( i-i0 ) * grid(nl)%dx
                grid(nl)%pos(i,j,2) = ( j-i0 ) * grid(nl)%dz
             end do ! i
          end do ! j
       end if ! curvedGeometry or not  ---------------

       ! get to the coarser level: 
       ! twice less nodes, and size twice larger.
       ntempo(:) = ntempo(:) / 2
       dxtempo = dxtempo * 2.0
       dztempo = dztempo * 2.0

    end do ! nl=0, nLevel
    
    allocate(iend(0:nLevel, 1:2))
    do nl=0, nLevel
       iend(nl,:) = grid(nl)%n(:) + i0 - 1
       grid(nl)%ntot = grid(nl)%n(1) * grid(nl)%n(2)
    end do

    ! for coarsest level of multigrid:
    if(leftV == "periodic") then  ! rightV is also periodic then
       coarsestMode = grid(nLevel)%n(1)
    else ! dirichlet or neumann
       coarsestMode = grid(nLevel)%n(1) - 1
    end if
    
    write(*,'(a,i3)') "coarsestMode = ", coarsestMode

    if(coarsestMode < 5) then
       allocate(mcc(1:coarsestMode, 1:coarsestMode)) ! square matrix
    else
       allocate(mcc(1:coarsestMode, 1:5)) ! square pentadiagonal matrix
    end if
    allocate(coarsestGamma(i0:iend(nLevel, 1))) ! diagonal term (taking into account topV and botV)
    allocate(xcc(1:coarsestMode)) ! solution on coarsest grid
    allocate(ycc(1:coarsestMode)) ! rhs on coarsest grid

    ! others:
    allocate(etayield(i0-1:iend0(1)+1, i0-1:iend0(2)+1, 1:4))
    allocate(plitho(i0-1:iend0(2)+1, 3:4))  ! alloc for id=3 and 4... 
    allocate(yieldT_z(i0-1:iend0(2)+1, 3:4)) ! ...(for different depths at cells' centers and corners)
    allocate(yieldC_z(i0-1:iend0(2)+1, 3:4)) ! ...(for different depths at cells' centers and corners)
    allocate(coeffViscoFK_z(i0-1:iend0(2)+1, 3:4))
    allocate(coeffViscoArr_z(i0-1:iend0(2)+1, 3:4))

    !! For maps: surface heat flux and horizontal veloc.
    allocate(qtop(i0:iend0(1)))
            
    allocate(xval(1:4, i0-1:iend0(1)+1)) ! 1:4 -> to cover different data
    allocate(zval(1:3, i0-1:iend0(2)+1))

    allocate(scalar(i0-1:iend0(1)+1, i0-1:iend0(2)+1))
    allocate(vector(i0-1:iend0(1)+1, i0-1:iend0(2)+1, 1:2))

    allocate(vxtop(i0-1:iend0(1)+1))
    allocate(vztop(i0-1:iend0(1)+1))
    allocate(sftop(i0-2:iend0(1)+2))
    
    allocate(vxbot(i0-1:iend0(1)+1))
    allocate(vzbot(i0-1:iend0(1)+1))
    allocate(sfbot(i0-2:iend0(1)+2))

    allocate(vxleft(i0-1:iend0(2)+1))
    allocate(vzleft(i0-1:iend0(2)+1))
    allocate(sfleft(i0-2:iend0(2)+2))

    allocate(vxright(i0-1:iend0(2)+1))
    allocate(vzright(i0-1:iend0(2)+1))
    allocate(sfright(i0-2:iend0(2)+2))
    
    if(heat_method == "ADI") then
       allocate(flim_m(1:n0(1)+4, 1:n0(2)+4, 1:2))
       allocate(flim_p(1:n0(1)+4, 1:n0(2)+4, 1:2))
       allocate(cadv_m(1:n0(1)+4, 1:n0(2)+4, 1:2))
       allocate(cadv_c(1:n0(1)+4, 1:n0(2)+4, 1:2))
       allocate(cadv_p(1:n0(1)+4, 1:n0(2)+4, 1:2))
    end if

    ! for heterogeneous internal heating and thermal conductivity:
    allocate(heat(i0-1:iend0(1)+1, i0-1:iend0(1)+2))
    allocate(ktherm(i0-1:iend0(1)+1, i0-1:iend0(1)+2, 1:2))
    
    write(*,'(a)') "========================================"

    return
  end subroutine initGrids
  ! -------------------------------------------------------
  subroutine freeGrids()
    ! clean and deallocate grid(:)%cell() and grid()
    implicit none
    integer:: n

    if(allocated(iend)) deallocate(iend)
    if(allocated(extent0)) deallocate(extent0)
    
    if(curvedGeometry) then
       if(allocated(r0)) deallocate(r0)
       if(allocated(phi0)) deallocate(phi0)
       if(allocated(alpha0)) deallocate(alpha0)
       if(allocated(dh)) deallocate(dh)
       if(allocated(dv)) deallocate(dv)
       if(allocated(dvp)) deallocate(dvp)
       if(allocated(invdh)) deallocate(invdh)
       if(allocated(invdv)) deallocate(invdv)
       if(allocated(invdvp)) deallocate(invdvp)
    end if ! curvedGeom

    if(allocated(grid)) then

       if(curvedGeometry) then !---------------      
          do n=0, nLevel
             if(allocated(grid(n)%r)) deallocate(grid(n)%r)
             if(allocated(grid(n)%phi)) deallocate(grid(n)%phi)
             if(allocated(grid(n)%alpha)) deallocate(grid(n)%alpha)
          end do ! n=0, nLevel
       end if ! curvedGeometry -----------------
       
       do n=0, nLevel
          if(allocated(grid(n)%extent)) deallocate(grid(n)%extent)
          if(allocated(grid(n)%pos)) deallocate(grid(n)%pos)       
          if(allocated(grid(n)%temp)) deallocate(grid(n)%temp)
          if(allocated(grid(n)%etaT)) deallocate(grid(n)%etaT)
          if(allocated(grid(n)%etaEff)) deallocate(grid(n)%etaEff)
          if(allocated(grid(n)%eta)) deallocate(grid(n)%eta)
          if(allocated(grid(n)%corrSf)) deallocate(grid(n)%corrSf)
          if(allocated(grid(n)%sfPrec)) deallocate(grid(n)%sfPrec)
          if(allocated(grid(n)%corrRhs)) deallocate(grid(n)%corrRhs)
          if(allocated(grid(n)%k)) deallocate(grid(n)%k)
          if(allocated(grid(n)%cpSf)) deallocate(grid(n)%cpSf)
          if(allocated(grid(n)%sf)) deallocate(grid(n)%sf)
          if(allocated(grid(n)%omega)) deallocate(grid(n)%omega)
          if(allocated(grid(n)%rhs)) deallocate(grid(n)%rhs)
          if(allocated(grid(n)%error)) deallocate(grid(n)%error)
          if(allocated(grid(n)%errorRhs)) deallocate(grid(n)%errorRhs)
          if(allocated(grid(n)%res)) deallocate(grid(n)%res)
          if(allocated(grid(n)%vel)) deallocate(grid(n)%vel)
          if(allocated(grid(n)%exx)) deallocate(grid(n)%exx)
          if(allocated(grid(n)%ezz)) deallocate(grid(n)%ezz)
          if(allocated(grid(n)%exz)) deallocate(grid(n)%exz)
          if(allocated(grid(n)%work1)) deallocate(grid(n)%work1)
          if(allocated(grid(n)%work2)) deallocate(grid(n)%work2)
          if(allocated(grid(n)%work3)) deallocate(grid(n)%work3)
          if(allocated(grid(n)%work4)) deallocate(grid(n)%work4)
          if(allocated(grid(n)%work5)) deallocate(grid(n)%work5)
          if(allocated(grid(n)%coeff)) deallocate(grid(n)%coeff)
          if(allocated(grid(n)%sist)) deallocate(grid(n)%sist)
          if(allocated(grid(n)%gtop)) deallocate(grid(n)%gtop)
          if(allocated(grid(n)%gbot)) deallocate(grid(n)%gbot)
          if(allocated(grid(n)%gmid)) deallocate(grid(n)%gmid)
          if(allocated(grid(n)%g1top)) deallocate(grid(n)%g1top)
          if(allocated(grid(n)%g1bot)) deallocate(grid(n)%g1bot)
          if(allocated(grid(n)%gcentertop)) deallocate(grid(n)%gcentertop)
          if(allocated(grid(n)%gcenterbot)) deallocate(grid(n)%gcenterbot)
          if(allocated(grid(n)%gcornertop)) deallocate(grid(n)%gcornertop)
          if(allocated(grid(n)%gcornerbot)) deallocate(grid(n)%gcornerbot)
          if(allocated(grid(n)%g4)) deallocate(grid(n)%g4)
          if(allocated(grid(n)%g5)) deallocate(grid(n)%g5)
          if(allocated(grid(n)%g12)) deallocate(grid(n)%g12)
          if(allocated(grid(n)%g13)) deallocate(grid(n)%g13)
          if(allocated(grid(n)%poisson)) deallocate(grid(n)%poisson)
          if(allocated(grid(n)%poissontop)) deallocate(grid(n)%poissontop)
          if(allocated(grid(n)%poissonbot)) deallocate(grid(n)%poissonbot)
          if(allocated(grid(n)%poissonmid)) deallocate(grid(n)%poissonmid)
          if(allocated(grid(n)%cbh)) deallocate(grid(n)%cbh)
          if(allocated(grid(n)%cbhY)) deallocate(grid(n)%cbhY)
          if(allocated(grid(n)%cbhPrec)) deallocate(grid(n)%cbhPrec)
       end do ! n=0, nLevel

       deallocate(grid)

    end if ! allocated(grid)

    if(allocated(flim_m)) deallocate(flim_m)
    if(allocated(flim_p)) deallocate(flim_p)
    if(allocated(cadv_m)) deallocate(cadv_m)
    if(allocated(cadv_c)) deallocate(cadv_c)
    if(allocated(cadv_p)) deallocate(cadv_p)
    if(allocated(etayield)) deallocate(etayield)
    if(allocated(plitho)) deallocate(plitho)
    if(allocated(yieldT_z)) deallocate(yieldT_z)
    if(allocated(yieldC_z)) deallocate(yieldC_z)
    if(allocated(coeffViscoFK_z)) deallocate(coeffViscoFK_z)
    if(allocated(coeffViscoArr_z)) deallocate(coeffViscoArr_z)
    if(allocated(xval)) deallocate(xval)
    if(allocated(zval)) deallocate(zval)
    if(allocated(scalar)) deallocate(scalar)
    if(allocated(vector)) deallocate(vector)
    if(allocated(qtop)) deallocate(qtop)
    if(allocated(vxtop)) deallocate(vxtop)
    if(allocated(vztop)) deallocate(vztop)
    if(allocated(sftop)) deallocate(sftop)
    if(allocated(vxbot)) deallocate(vxbot)
    if(allocated(vzbot)) deallocate(vzbot)
    if(allocated(sfbot)) deallocate(sfbot)
    if(allocated(vxleft)) deallocate(vxleft)
    if(allocated(vzleft)) deallocate(vzleft)
    if(allocated(sfleft)) deallocate(sfleft)
    if(allocated(vxright)) deallocate(vxright)
    if(allocated(vzright)) deallocate(vzright)
    if(allocated(sfright)) deallocate(sfright)
    
    if(allocated(mcc)) deallocate(mcc)
    if(allocated(xcc)) deallocate(xcc)
    if(allocated(ycc)) deallocate(ycc)
    if(allocated(coarsestGamma)) deallocate(coarsestGamma)

    if(allocated(ktherm)) deallocate(ktherm)
    if(allocated(heat)) deallocate(heat)

    return
  end subroutine freeGrids
  ! -------------------------------------------------------
  subroutine getN(nl, level)
    ! return the grid size nl(1), nl(2) at level=level
    implicit none
    integer, intent(in):: level
    integer, dimension(1:2), intent(out):: nl
    integer:: idir

    do idir = 1, 2
       nl(idir) = grid(level)%n(idir)
    end do

    return
  end subroutine getN
  !-------------------------------------------------------
  subroutine getN_1D(nl, level, idir)
    ! return the grid size nl(idir) at level=level
    implicit none
    integer, intent(in):: level, idir
    integer, intent(out):: nl

    nl = grid(level)%n(idir)

    return
  end subroutine getN_1D
  !-------------------------------------------------------
  subroutine getD(d, level)
    ! returns grid spacing dx (idir=1) and dz (idir=2)
    ! at level=level
    implicit none
    integer, intent(in):: level
    real(rk), intent(out):: d(1:2)

    d(1) = grid(level)%dx
    d(2) = grid(level)%dz

    return
  end subroutine getD
  !-------------------------------------------------------
  subroutine getInvD(invd, level)
    ! return indv(1) = 1/dx and indv(2)=1/dz 
    ! at level=level
    integer, intent(in)::level
    real(rk), dimension(1:2), intent(out):: invd

    invd(1) = 1.0 / grid(level)%dx
    invd(2) = 1.0 / grid(level)%dz

    return
  end subroutine getInvD
  !-------------------------------------------------------
  subroutine getInvD2(invd2, level)
    ! return indv2(1) = 1/dx^2 etc.
    ! at level=level
    integer, intent(in)::level
    real(rk), dimension(1:2), intent(out):: invd2

    invd2(1) = 1.0 / grid(level)%dx**2
    invd2(2) = 1.0 / grid(level)%dz**2

    return
  end subroutine getInvD2
  !-------------------------------------------------------
  subroutine zeroingAll(level)
    ! put all objects in grids at level=level to zero
    implicit none
    integer, intent(in):: level

    ! all following variables are defined on vectors
    call zeroingVec(temp, level)
    call zeroingVec(etaT, level)
    call zeroingVec(etaEff, level)
    call zeroingVec(etaPrec, level)
    call zeroingVec(eta, level)
    call zeroingVec(corrSf, level)
    call zeroingVec(sfPrec, level)
    call zeroingVec(corrRhs, level)
    call zeroingVec(k, level)
    call zeroingVec(cpSf, level)
    call zeroingVec(sf, level)
    call zeroingVec(omega, level)
    call zeroingVec(vel, level)
    call zeroingVec(rhs, level)
    call zeroingVec(res, level)
    call zeroingVec(error, level)
    call zeroingVec(errorRhs, level)
    call zeroingVec(exx, level)
    call zeroingVec(ezz, level)
    call zeroingVec(exz, level)
    
    call zeroingVec(work1, level)
    call zeroingVec(work2, level)
    call zeroingVec(work3, level)
    call zeroingVec(work4, level)
    call zeroingVec(work5, level)
    call zeroingVec(coeff, level)
    call zeroingVec(sist, level)

    ! special formats:
    call zeroingVec(cbh, level)
    call zeroingVec(cbhY, level)
    call zeroingVec(cbhPrec, level)
    
    return
  end subroutine zeroingAll
  !-------------------------------------------------------
    subroutine zeroingVec(what, level)
    ! put all objects of type "what" in grids at level=level to zero
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level

    call setValueVec(what, level, zero)

    return
  end subroutine zeroingVec
  !-------------------------------------------------------
  subroutine zeroingScal(what, level, id)
    ! put grid(level)%what(:,:,id) at zero
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id

    call setValueScal(what, level, id, zero)

    return
  end subroutine zeroingScal
  !-------------------------------------------------------
  subroutine setValueScal(what, level, id, val)
    ! gives the same value "val" to grid(level)%what(:,:,id)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), intent(in):: val
    real(rk), pointer:: ptr(:,:) => null()

    call associatePtrScal(ptr, what, level, id)
    ptr(:,:) = val

    nullify(ptr)

    return
  end subroutine setValueScal
  !-------------------------------------------------------
  subroutine setValueVec(what, level, val)
    ! gives the same value "val" to grid(level)%what(:,:,:)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    real(rk), intent(in):: val
    real(rk), pointer:: ptr(:,:,:) => null()

    call associatePtrVec(ptr, what, level)
    ptr(:,:,:) = val

    nullify(ptr)

    return
  end subroutine setValueVec
  !-------------------------------------------------------
  function getValue(what, level, id, i,j)
    ! gets the value at grid(level)%what(i,j,id)
    ! warning: not a function that should be included inside a big
    ! do loop -> all the if.. else if... are costly
    implicit none
    integer, intent(in):: level, id, i, j
    character(len=10), intent(in):: what
    real(rk):: getValue

    getvalue = infinity

    if(what == pos) then
       getValue = grid(level)%pos(i,j,id)
    else if(what == temp) then
       getValue = grid(level)%temp(i,j,id)
    else if(what == etaT) then
       getValue = grid(level)%etaT(i,j,id)
    else if(what == etaEff) then
       getValue = grid(level)%etaEff(i,j,id)
    else if(what == etaPrec) then
       getValue = grid(level)%etaPrec(i,j,id)       
    else if(what == eta) then
       getValue = grid(level)%eta(i,j,id)
    else if(what == corrSf) then
       getValue = grid(level)%corrSf(i,j,id)
    else if(what == sfPrec) then
       getValue = grid(level)%sfPrec(i,j,id)
    else if(what == corrRhs) then
       getValue = grid(level)%corrRhs(i,j,id)
    else if(what == k) then
       getValue = grid(level)%k(i,j,id)
    else if(what == cpSf) then
       getValue = grid(level)%cpSf(i,j,id)
    else if(what == vel) then
       getValue = grid(level)%vel(i,j,id)
    else if(what == exx) then
       getValue = grid(level)%exx(i,j,id)
    else if(what == ezz) then
       getValue = grid(level)%ezz(i,j,id)
    else if(what == exz) then
       getValue = grid(level)%exz(i,j,id)
    else if(what == rhs) then
       getValue = grid(level)%rhs(i,j,id)
    else if(what == res) then
       getValue = grid(level)%res(i,j,id)
    else if(what == error) then
       getValue = grid(level)%error(i,j,id)
    else if(what == errorRhs) then
       getValue = grid(level)%errorRhs(i,j,id)
    else if(what == work1) then
       getValue = grid(level)%work1(i,j,id)
    else if(what == work2) then
       getValue = grid(level)%work2(i,j,id)
    else if(what == work3) then
       getValue = grid(level)%work3(i,j,id)
    else if(what == work4) then
       getValue = grid(level)%work4(i,j,id)
    else if(what == work5) then
       getValue = grid(level)%work5(i,j,id)
    else if(what == coeff) then
       getValue = grid(level)%coeff(i,j,id)
    else if(what == sf) then
       getValue = grid(level)%sf(i,j,id)
    else if(what == omega) then
       getValue = grid(level)%omega(i,j,id)
    else if(what == sist) then
       getValue = grid(level)%sist(i,j,id)
    else if(what == r) then
       getValue = grid(level)%r(j)  ! i and id not needed
    else if(what == phi) then
       getValue = grid(level)%phi(i)  ! j and id not needed
    end if

    return
  end function getValue
  !-------------------------------------------------------
  subroutine setValue(val, what, level, id, i,j)
    ! puts val in grid(level)%what(i, j ,id),
    ! same warning as above.
    implicit none
    integer, intent(in):: level, id, i, j
    character(len=10), intent(in):: what
    real(rk), intent(in):: val

    if(what == pos) then
       grid(level)%pos(i,j,id) = val
    else if(what == temp) then
       grid(level)%temp(i,j,id) = val
    else if(what == etaT) then
       grid(level)%etaT(i,j,id) = val
    else if(what == etaEff) then
       grid(level)%etaEff(i,j,id) = val
    else if(what == etaPrec) then
       grid(level)%etaPrec(i,j,id) = val
    else if(what == eta) then
       grid(level)%eta(i,j,id) = val
    else if(what == corrSf) then
       grid(level)%corrSf(i,j,id) = val
    else if(what == sfPrec) then
       grid(level)%sfPrec(i,j,id) = val
    else if(what == corrRhs) then
       grid(level)%corrRhs(i,j,id) = val
    else if(what == k) then
       grid(level)%k(i,j,id) = val
    else if(what == sf) then
       grid(level)%sf(i,j,id) = val
    else if(what == vel) then
       grid(level)%vel(i,j,id) = val
    else if(what == exx) then
       grid(level)%exx(i,j,id) = val
    else if(what == ezz) then
       grid(level)%ezz(i,j,id) = val
    else if(what == exz) then
       grid(level)%exz(i,j,id) = val
    else if(what == rhs) then
       grid(level)%rhs(i,j,id) = val
    else if(what == res) then
       grid(level)%res(i,j,id) = val
    else if(what == error) then
       grid(level)%error(i,j,id) = val
    else if(what == errorRhs) then
       grid(level)%errorRhs(i,j,id) = val
    else if(what == work1) then
       grid(level)%work1(i,j,id) = val
    else if(what == work2) then
       grid(level)%work2(i,j,id) = val
    else if(what == work3) then
       grid(level)%work3(i,j,id) = val
    else if(what == work4) then
       grid(level)%work4(i,j,id) = val
    else if(what == work5) then
       grid(level)%work5(i,j,id) = val
    else if(what == coeff) then
       grid(level)%coeff(i,j,id) = val
    else if(what == sf) then
       grid(level)%sf(i,j,id) = val
    else if(what == omega) then
       grid(level)%omega(i,j,id) = val
    else if(what == sist) then
       grid(level)%sist(i,j,id) = val
    end if

    return
  end subroutine setValue
  !-------------------------------------------------------
  subroutine associatePtrScal(ptr, what, level, id)
    ! associate ptr to grid(level)%what(:,:, id)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), pointer:: ptr(:,:)

    if(what == temp) then
       ptr => grid(level)%temp(:,:,id)
    else if(what == etaT) then
       ptr => grid(level)%etaT(:,:,id)
    else if(what == etaEff) then
       ptr => grid(level)%etaEff(:,:,id)
    else if(what == etaPrec) then
       ptr => grid(level)%etaPrec(:,:,id)
    else if(what == eta) then
       ptr => grid(level)%eta(:,:,id)
    else if(what == corrSf) then
       ptr => grid(level)%corrSf(:,:,id)
    else if(what == sfPrec) then
       ptr => grid(level)%sfPrec(:,:,id)
    else if(what == corrRhs) then
       ptr => grid(level)%corrRhs(:,:,id)
    else if(what == k) then
       ptr => grid(level)%k(:,:,id)
    else if(what == cpSf) then
       ptr => grid(level)%cpSf(:,:,id)
    else if(what == pos) then
       ptr => grid(level)%pos(:,:,id)       
    else if(what == vel) then
       ptr => grid(level)%vel(:,:,id)       
    else if(what == exx) then
       ptr => grid(level)%exx(:,:,id)       
    else if(what == ezz) then
       ptr => grid(level)%ezz(:,:,id)       
    else if(what == exz) then
       ptr => grid(level)%exz(:,:,id)       
    else if(what == rhs) then
       ptr => grid(level)%rhs(:,:,id)
    else if(what == res) then
       ptr => grid(level)%res(:,:,id)
    else if(what == error) then
       ptr => grid(level)%error(:,:,id)
    else if(what == errorRhs) then
       ptr => grid(level)%errorRhs(:,:,id)
    else if(what == work1) then
       ptr => grid(level)%work1(:,:,id)
    else if(what == work2) then
       ptr => grid(level)%work2(:,:,id)
    else if(what == work3) then
       ptr => grid(level)%work3(:,:,id)
    else if(what == work4) then
       ptr => grid(level)%work4(:,:,id)
    else if(what == work5) then
       ptr => grid(level)%work5(:,:,id)
    else if(what == coeff) then
       ptr => grid(level)%coeff(:,:,id)
    else if(what == sf) then
       ptr => grid(level)%sf(:,:,id)
    else if(what == omega) then
       ptr => grid(level)%omega(:,:,id)
    else if(what == sist) then
       ptr => grid(level)%sist(:,:,id)
    end if

    return
  end subroutine associatePtrScal
  !-------------------------------------------------------
  subroutine associatePtrVec(ptr, what, level)
    ! associate ptr to grid(level)%what(:,:,:)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    real(rk), pointer:: ptr(:,:,:)

    if(what == temp) then
       ptr => grid(level)%temp
    else if(what == etaT) then
       ptr => grid(level)%etaT
    else if(what == etaEff) then
       ptr => grid(level)%etaEff
    else if(what == etaPrec) then
       ptr => grid(level)%etaPrec
    else if(what == eta) then
       ptr => grid(level)%eta
    else if(what == corrSf) then
       ptr => grid(level)%corrSf
    else if(what == sfPrec) then
       ptr => grid(level)%sfPrec
    else if(what == corrRhs) then
       ptr => grid(level)%corrRhs
    else if(what == k) then
       ptr => grid(level)%k
    else if(what == cpSf) then
       ptr => grid(level)%cpSf
    else if(what == pos) then
       ptr => grid(level)%pos
    else if(what == vel) then
       ptr => grid(level)%vel
    else if(what == exx) then
       ptr => grid(level)%exx
    else if(what == ezz) then
       ptr => grid(level)%ezz
    else if(what == exz) then
       ptr => grid(level)%exz
    else if(what == rhs) then
       ptr => grid(level)%rhs
    else if(what == res) then
       ptr => grid(level)%res
    else if(what == errorRhs) then
       ptr => grid(level)%errorRhs
    else if(what == error) then
       ptr => grid(level)%error
    else if(what == work1) then
       ptr => grid(level)%work1
    else if(what == work2) then
       ptr => grid(level)%work2
    else if(what == work3) then
       ptr => grid(level)%work3
    else if(what == work4) then
       ptr => grid(level)%work4
    else if(what == work5) then
       ptr => grid(level)%work5
    else if(what == coeff) then
       ptr => grid(level)%coeff
    else if(what == sf) then
       ptr => grid(level)%sf
    else if(what == omega) then
       ptr => grid(level)%omega
    else if(what == sist) then
       ptr => grid(level)%sist
    else if(what == cbh) then
       ptr => grid(level)%cbh
    else if(what == cbhY) then
       ptr => grid(level)%cbhY
    else if(what == cbhPrec) then
       ptr => grid(level)%cbhPrec
    end if

    return
  end subroutine associatePtrVec
  !-------------------------------------------------------
  subroutine givePlan(what, level, id, idir, iget, iset, mult)
    ! give 1D (plan(:)) to another plan, in the direction idir
    ! if idir==1: grid(level)%what(iset,:, id) = grid(level)%what(iget,:, id)
    ! if idir==2: grid(level)%what(:,iset, id) = grid(level)%what(:,iget, id)
    ! if mult is present: set = mult * get
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id, idir, iget, iset
    real(rk), optional:: mult
    real(rk), pointer:: ptr(:,:) => null()

    if(iget == iset) return !nothing to do here

    call associatePtrScal(ptr, what, level, id)
    if(idir == 1) then
       if(present(mult)) then
          ptr(iset,:) = mult * ptr(iget,:)
       else
          ptr(iset,:) = ptr(iget,:)
       end if
    else if(idir == 2) then
       if(present(mult)) then
          ptr(:,iset) = mult * ptr(:,iget)
       else
          ptr(:,iset) = ptr(:,iget)
       end if
    end if

    nullify(ptr)

    return
  end subroutine givePlan
  !-------------------------------------------------------
  subroutine setValuePlan(what, level, id, idir, iset, val)
    ! if idir==1: grid(level)%what(iset,:, id) = val
    ! if idir==2: grid(level)%what(:,iset, id) = val
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id, idir, iset
    real(rk), intent(in):: val
    real(rk), pointer:: ptr(:,:) => null()

    call associatePtrScal(ptr, what, level, id)
    if(idir == 1) then
       ptr(iset,:) = val
    else if(idir == 2) then
       ptr(:,iset) = val
    end if

    nullify(ptr)

    return
  end subroutine setValuePlan
  !-------------------------------------------------------
  subroutine addValuePlan(what, level, id, idir, iset, add)
    ! if idir==1: grid(level)%what(iset,:, id) += add
    ! if idir==2: grid(level)%what(:,iset, id) += add
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id, idir, iset
    real(rk), intent(in):: add
    real(rk), pointer:: ptr(:,:) => null()

    call associatePtrScal(ptr, what, level, id)
    if(idir == 1) then
       ptr(iset,:) = ptr(iset,:) + add
    else if(idir == 2) then
       ptr(:,iset) = ptr(:,iset) + add
    end if

    nullify(ptr)

    return
  end subroutine addValuePlan
  !-------------------------------------------------------
  subroutine setMiddleValuePlan(what, level, id, idir, iget, iset, val)
    ! e.g. if idir==1: makes grid(level)%what(iset,:, id) = 
    ! 2.0*val - grid(level)%what(iget,:,:, id)
    ! Useful for setting for instance 
    ! v(i0-1,:,:) = 2*val - v(i0,:,:) (boundary conditions)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id, idir, iget, iset
    real(rk), intent(in):: val
    real(rk), pointer:: ptr(:,:) => null()

    call associatePtrScal(ptr, what, level, id)
    if(idir == 1) then
       ptr(iset,:) = 2.0 * val - ptr(iget,:)
    else if(idir == 2) then
       ptr(:,iset) = 2.0 * val - ptr(:,iget)
    end if

    nullify(ptr)

    return
  end subroutine setMiddleValuePlan
  !-------------------------------------------------------
  subroutine multiplyVec(what, level, mult)
    ! multiply grid(level)%what(:,:,:) by mult
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    real(rk), intent(in):: mult
    real(rk), pointer:: ptr(:,:,:) => null()

    call associatePtrVec(ptr, what, level)
    ptr(:,:,:) = ptr(:,:,:) * mult

    nullify(ptr)

    return
  end subroutine multiplyVec
  !-------------------------------------------------------
  subroutine multiplyScal(what, level, id, mult)
    ! multiply grid(level)%what(:,:,id) by mult
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), intent(in):: mult
    real(rk), pointer:: ptr(:,:) => null()

    call associatePtrScal(ptr, what, level, id)
    ptr(:,:) = ptr(:,:) * mult

    nullify(ptr)

    return
  end subroutine multiplyScal
  !-------------------------------------------------------
  subroutine addVec(added, source, level, mult)
    ! adds grid(level)%source(:,:,:) to grid(level)%added(:,:,:)
    ! i.e. added = added + source
    ! if mult is present: added = added + mult * source
    implicit none
    character(len=10), intent(in):: added, source
    integer, intent(in):: level
    real(rk), optional:: mult
    real(rk), pointer:: ptrAd(:,:,:) => null()
    real(rk), pointer:: ptrSrc(:,:,:) => null()

    call associatePtrVec(ptrSrc, source, level)
    call associatePtrVec(ptrAd, added, level)

    if(present(mult)) then
       ptrAd(:,:,:) = ptrAd(:,:,:) + mult * ptrSrc(:,:,:)
    else
       ptrAd(:,:,:) = ptrAd(:,:,:) + ptrSrc(:,:,:)
    end if

    nullify(ptrAd, ptrSrc)

    return
  end subroutine addVec
  !-------------------------------------------------------
  subroutine addScal(added, source, level, id, mult)
    ! adds grid(level)%source(:,:,id) to grid(level)%added(:,:,id)
    ! i.e. added = added + source
    ! if mult is present: added = added + mult * source
    implicit none
    character(len=10), intent(in):: added, source
    integer, intent(in):: level, id
    real(rk), optional:: mult
    real(rk), pointer:: ptrAd(:,:) => null()
    real(rk), pointer:: ptrSrc(:,:) => null()

    call associatePtrScal(ptrSrc, source, level, id)
    call associatePtrScal(ptrAd, added, level, id)

    if(present(mult)) then
       ptrAd(:,:) = ptrAd(:,:) + mult * ptrSrc(:,:)
    else
       ptrAd(:,:) = ptrAd(:,:) + ptrSrc(:,:)
    end if

    nullify(ptrAd, ptrSrc)

    return
  end subroutine addScal
  !-------------------------------------------------------
  subroutine addFactorScal(what, whatAdd, level, id, factor)
    ! what = (1-factor)*what + factor*whatAdd
    character(len=10), intent(in):: what, whatAdd
    integer, intent(in):: level, id
    real(rk), intent(in):: factor
    real(rk), pointer:: ptrAdd(:,:) => null()
    real(rk), pointer:: ptr(:,:) => null()
    integer:: i, j
    real(rk):: beta
    
    call associatePtrScal(ptr, what, level, id)
    call associatePtrScal(ptrAdd, whatAdd, level, id)

    beta = 1.0 - factor
    do j=i0, iend(level, 2)
       do i=i0, iend(level, 1)
          ptr(i, j) = beta * ptr(i, j) + factor * ptrAdd(i, j)
       end do
    end do
    nullify(ptr, ptrAdd)
    
    return
  end subroutine addFactorScal
  !-------------------------------------------------------
  subroutine addValueScal(what, level, id, val)
    ! adds a scalar value "val" to all grid(level)%what(:,:,id)
    implicit none
    character(len=10), intent(in)::what
    integer, intent(in):: level, id
    real(rk), intent(in):: val
    real(rk), pointer:: ptr(:,:) => null()

    call associatePtrScal(ptr, what, level, id)
    ptr(:,:) = ptr(:,:) + val

    nullify(ptr)

    return
  end subroutine addValueScal
  !-------------------------------------------------------
  subroutine copyVec(copy, source, level, mult)
    ! copy grid(level)%source(:,:,:) in
    ! grid(level)%copy(:,:,:). 
    ! if mult is present : copy *= mult
    implicit none
    character(len=10), intent(in):: source, copy
    integer, intent(in):: level
    real(rk), optional:: mult
    real(rk), pointer:: ptrSource(:,:,:) => null()
    real(rk), pointer:: ptrCopy(:,:,:) => null()

    call associatePtrVec(ptrCopy, copy, level)
    call associatePtrVec(ptrSource, source, level)

    if(present(mult)) then
       ptrCopy(:,:,:) = mult * ptrSource(:,:,:)
    else
       ptrCopy(:,:,:) = ptrSource(:,:,:)
    end if

    nullify(ptrCopy, ptrSource)

    return
  end subroutine copyVec
  !-------------------------------------------------------
  subroutine copyScal(copy, icopy, source, isource, level, mult)
    ! copies grid(level)%source(:,:,isource) in 
    ! grid(level)%copy(:,:,icopy).
    ! "copy" and "source" must be defined as a "what".
    ! if mult is present : copy *= mult
    implicit none
    character(len=10), intent(in):: copy, source
    integer, intent(in):: icopy, isource, level
    real(rk), optional:: mult
    real(rk), pointer:: ptrSource(:,:) => null()
    real(rk), pointer:: ptrCopy(:,:) => null()

    call associatePtrScal(ptrSource, source, level, isource)
    call associatePtrScal(ptrCopy, copy, level, icopy)

    if(present(mult)) then
       ptrCopy(:,:) = mult * ptrSource(:,:)
    else
       ptrCopy(:,:) = ptrSource(:,:)
    end if

    nullify(ptrCopy, ptrSource)

    return
  end subroutine copyScal
  !-------------------------------------------------------
  subroutine copyPlan(copy, icopy, source, isource, level, idir, iplan, mult)
    ! copies grid(level)%source(:,iplan,isource) in
    ! grid(level)%copy(:,iplan,icopy) if idir=2 (other dir if idir=1)
    ! if mult is present : copy *= mult
    implicit none
    character(len=10), intent(in):: copy, source
    integer, intent(in):: icopy, isource, level, idir, iplan
    real(rk), optional:: mult
    real(rk), pointer:: ptrSource(:,:) => null()
    real(rk), pointer:: ptrCopy(:,:) => null()

    call associatePtrScal(ptrCopy, copy, level, icopy)
    call associatePtrScal(ptrSource, source, level, isource)

    if(idir == 1) then
       if(present(mult)) then
          ptrCopy(iplan,:) = mult * ptrSource(iplan,:)
       else
          ptrCopy(iplan,:) = ptrSource(iplan,:)
       end if
    else if(idir == 2) then
       if(present(mult)) then
          ptrCopy(:,iplan) = mult * ptrSource(:,iplan)
       else
          ptrCopy(:,iplan) = ptrSource(:,iplan)
       end if
    end if

    nullify(ptrCopy, ptrSource)

    return
  end subroutine copyPlan
  !-------------------------------------------------------
  subroutine copyId(what, level, idsource, idcopy)
    ! copies grid(level)%what(:,:,idsource) to
    ! grid(level)%what(:,:,idcopy)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, idsource, idcopy
    real(rk), pointer:: ptrSource(:,:) => null()
    real(rk), pointer:: ptrCopy(:,:) => null()

    call associatePtrScal(ptrSource, what, level, idsource)
    call associatePtrScal(ptrCopy, what, level, idcopy)

    ptrCopy(:,:) = ptrSource(:,:)
    
    nullify(ptrSource, ptrCopy)
    
    return
  end subroutine copyId
  !-------------------------------------------------------
  subroutine getRms(rms, what, level)
    ! compute the rms with considering normal cells, with the
    ! value at the centers id=3 (the interpolation
    ! must have been done before).

    ! For curvedGeometry:
    !    The "extent" (volume when ed=2, area when ed=1) of each 
    !    cell is taken into account.
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    real(rk), intent(out):: rms
    real(rk), pointer:: ptr(:,:) => null()
    real(rk), pointer:: pe(:) => null()
    integer, dimension(1:2):: nl
    integer:: i, j
    
    call associatePtrScal(ptr, what, level, 3)
    call associatePtrLine(pe, extent, level)
    call getN(nl, level)
    
    rms = 0.0
    do j=i0, iend(level, 2)
       do i=i0, iend(level, 1)
          rms = rms + pe(j) * ptr(i,j) * ptr(i,j) 
       end do
    end do
    rms = sqrt(rms / grid(level)%n(1))

    nullify(ptr, pe)
    return
  end subroutine getRms
  !-------------------------------------------------------
  subroutine getMean(meanVal, what, level, absVal)
    ! compute the mean with considering normal cells, with the
    ! value at the centers id=3 (the interpolation
    ! must have been done before).
    ! For curvedGeometry: 
    !    The "extent" (area if ed=1 and volume if ed=2) of each 
    !    cell is taken into account.
    character(len=10), intent(in):: what
    integer, intent(in):: level
    real(rk), intent(out):: meanVal
    logical, optional:: absVal
    logical:: doAbs
    integer, dimension(1:2):: nl
    real(rk), pointer:: ptr(:,:) => null()
    real(rk), pointer:: pe(:) => null()
    integer:: i, j
    
    doAbs = .false.
    if(present(absVal)) doAbs = absVal

    call associatePtrScal(ptr, what, level, 3)
    call associatePtrLine(pe, extent, level)
    call getN(nl, level)
    
    meanVal = 0.0
    if(doAbs) then
       do j=i0, iend(level, 2)
          do i=i0, iend(level, 1)
             meanVal = meanVal + pe(j) * abs(ptr(i,j))
          end do ! i
       end do ! j
    else 
       do j=i0, iend(level, 2)
          do i=i0, iend(level, 1)
             meanVal = meanVal + pe(j) * ptr(i,j)
          end do ! i
       end do ! j
    end if ! doAbs

    meanVal = meanVal / dble(grid(level)%n(1))
    
    nullify(ptr, pe)
    
    return 
  end subroutine getMean
  !-------------------------------------------------------
  subroutine getRmsSimple(rms, what, level, id)
    ! simply add together values, whatever their position
    ! (less accurate than getRms(), which considers the different
    ! sizes of cells)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), intent(out):: rms
    real(rk), pointer:: ptr(:,:) => null()
    integer:: i, j

    call associatePtrScal(ptr, what, level, id)
    rms = 0.0
   
    do j=i0, iend(level,2)
       do i=i0, iend(level,1)
          rms = rms + ptr(i,j)*ptr(i,j)
       end do
    end do

    rms = sqrt(rms / dble(grid(level)%ntot))

    nullify(ptr)

    return
  end subroutine getRmsSimple
  !-------------------------------------------------------
  subroutine getMeanSimple(meanVal, what, level, id, absVal)
    ! simply add together values, whatever their position
    ! (less accurate than getMean(), which considers the different
    ! sizes of cells)
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), intent(out):: meanVal
    logical, optional:: absVal
    logical:: doAbs
    real(rk), pointer:: ptr(:,:) => null()
    integer:: i, j

    doAbs = .false.
    if(present(absVal)) doAbs = absVal

    call associatePtrScal(ptr, what, level, id)
    meanVal = 0.0
    
    if(doAbs) then
       do j=i0, iend(level, 2)
          do i=i0, iend(level, 1)
             meanVal = meanVal + abs(ptr(i, j))
          end do ! i
       end do ! j
    else 
       do j=i0, iend(level, 2)
          do i=i0, iend(level, 1)
             meanVal = meanVal + ptr(i, j)
          end do ! i
       end do ! j
    end if ! doAbs

    meanVal = meanVal / dble(grid(level)%ntot)

    nullify(ptr)

    return
  end subroutine getMeanSimple
  !-------------------------------------------------------
  subroutine getMax(maxValue, what, level, id, absVal, borders)
    !returns max value of grid(level)%what(:,:,id)
    ! absVal: optional
    ! if .not.present: simple max
    ! if present: check abs(max)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), intent(out):: maxValue
    logical, optional:: absVal, borders
    logical:: doAbs
    integer:: ibord
    real(rk), pointer:: ptr(:,:) => null()
    integer:: i,j

    call associatePtrScal(ptr, what, level, id)

    doAbs = .false.
    if(present(absVal)) doAbs = absVal

    ibord = 0
    if(present(borders)) then
       if(borders) ibord = 1 ! includes borders (should be used mainly for id=4)
    end if

    maxValue = -1.0E-10
    if(doAbs) then
       do j=i0, iend(level,2) + ibord
          do i=i0, iend(level,1) + ibord
             maxValue = max( abs(ptr(i, j)), maxValue )
          end do
       end do
    else
       do j=i0, iend(level,2) + ibord
          do i=i0, iend(level,1) + ibord
             maxValue = max( ptr(i, j), maxValue )
          end do
       end do
    end if

    nullify(ptr)

    return
  end subroutine getMax
  !-------------------------------------------------------
  subroutine getMin(minValue, what, level, id, absVal, borders)
    !returns min value of grid(level)%what(:,:,id)
    ! absVal: optional
    ! if .not.present: simple min
    ! if present: check abs(min)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), intent(out):: minValue
    logical, optional:: absVal, borders
    logical:: doAbs
    integer:: ibord
    real(rk), pointer:: ptr(:,:) => null()
    integer:: i,j

    call associatePtrScal(ptr, what, level, id)

    doAbs = .false.
    if(present(absVal)) doAbs = absVal

    ibord = 0
    if(present(borders)) then
       if(borders) ibord = 1
    end if

    minValue = +1.0E30
    if(doAbs) then
       do j=i0, iend(level,2) + ibord
          do i=i0, iend(level,1) + ibord
             minValue = min( abs(ptr(i, j)), minValue )
          end do
       end do
    else
       do j=i0, iend(level,2) + ibord
          do i=i0, iend(level,1) + ibord
             minValue = min( ptr(i, j), minValue )
          end do
       end do
    end if

    nullify(ptr)

    return
  end subroutine getMin
  !-------------------------------------------------------
  subroutine getMeanPlan(meanVal, what, level, id, idir, iplan, rms)
    ! returns meanVal, the mean value of 
    ! grid(level)%what(iplan,:, id) if idir=1
    ! grid(level)%what(:,iplan, id) if idir=2
    ! rms: optional
    ! WARNING: for idir=1, computation is done at centers
    !    of cells (interpolation to id=3 must have been done before)
    implicit none
    real(rk), intent(out):: meanVal
    character(len=10), intent(in):: what
    integer, intent(in):: level, id, idir, iplan
    logical, optional:: rms
    logical:: do_rms
    real(rk), pointer:: ptr(:,:) => null()
    real(rk), pointer:: pe(:) => null()
    integer:: i, j, nl(1:2), iex, iez

    if(idir == 1 .and. .not. id == 3) then
       write(*,*) &
            "ERROR: getMeanPlan() not implemented for idir=1 and id=/=3"
       call endRun()
    end if
    
    call associatePtrScal(ptr, what, level, id)
    call associatePtrLine(pe, extent, level)
    
    call getN(nl, level)

    do_rms = .false.
    if(present(rms)) do_rms = rms

    meanVal = 0.0
    if(idir == 1) then
       ! --------------------------------------------------------------
       ! sum along a vertical plan (the horizontal coordinate is fixed)
       ! takes into account the different extent of cells 
       ! --------------------------------------------------------------
       iez = iend(level, 2)
       if(do_rms) then
          do j=i0, iez
             meanVal = meanVal + pe(j) * ptr(iplan, j) * ptr(iplan, j)
          end do !j
       else  ! .not.do_rms
          do j=i0, iez
             meanVal = meanVal + pe(j) * ptr(iplan, j)
          end do
       end if ! rms or not
    else ! idir == 2
       ! --------------------------------------------------------------
       ! sum along a horizontal plan (the vertical coordinate is fixed)
       ! --------------------------------------------------------------
       iex = iend(level, 1)
       if(do_rms) then
          do i=i0, iex
             meanVal = meanVal + ptr(i, iplan) * ptr(i, iplan)
          end do
       else ! not rms
          do i=i0, iex
             meanVal = meanVal + ptr(i, iplan)
          end do
       end if ! rms or not

       meanVal = meanVal / dble(nl(1))

    end if ! idir
    
    if (do_rms) meanVal = sqrt(meanVal)
    
    nullify(ptr, pe)

    return
  end subroutine getMeanPlan 
  ! -------------------------------------------------------
  subroutine getMaxPlan(maxVal, what, level, id, idir, iplan, absval)
    ! returns maxVal, the maximum value of 
    ! grid(level)%what(iplan,:, id) if idir=1
    ! grid(level)%what(:,iplan, id) if idir=2
    ! abs: optional
    implicit none
    real(rk), intent(out):: maxVal
    character(len=10), intent(in):: what
    integer, intent(in):: level, id, idir, iplan
    logical, optional:: absval
    logical:: do_abs
    real(rk), pointer:: ptr(:,:) => null()
    integer:: i, j, nl(1:2), iex, iez, ibord

    call associatePtrScal(ptr, what, level, id)

    call getN(nl, level)

    do_abs = .false.
    if(present(absval)) do_abs = absval
    
    maxVal = -infinity
    if(idir == 1) then
       ! vertical plan (the horizontal coordinate is fixed)
       iez = iend(level, 2)
       ibord = 0 ! for id=3 or 1 (centers of cells)
       if(id == 2 .or. id == 4) then ! borders of cells
          ibord = 1
       end if
       if(do_abs) then
          do j=i0, iez+ibord
             maxval = max(maxval, abs(ptr(iplan, j)))
          end do
       else
          do j=i0, iez+ibord
             maxval = max(maxval, ptr(iplan, j))
          end do
       end if
    else if(idir == 2) then
       ! horizontal plan (vertical coordinate is fixed)
       iex = iend(level, 1)
       ibord = 0  ! for id=2 or id=3 (centers of cells)
       if(id == 1 .or. id == 4) then ! borders
          ibord = 1
       end if
       if(do_abs) then
          do i=i0, iex+ibord
             maxval = max(maxval, abs(ptr(i, iplan)))
          end do ! i
       else
          do i=i0, iex+ibord
             maxval = max(maxval, ptr(i, iplan))
          end do ! i
       end if
    end if ! idir
    nullify(ptr)
    
    return
  end subroutine getMaxPlan
  ! -------------------------------------------------------
  subroutine getPercentiles(what, percentiles, id)
    ! returns P10, P25(=Q1), P50(=median), P75(=Q3) and P90
    ! of "what".
    ! if id=1 or id=3: does not include top and right walls
    ! if id=2 or id=4: walls are included.
    ! Only level zero is implemented.
    ! WARNING: does not take into account the different
    ! sizes of cells.
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: id
    real(rk), dimension(1:5), intent(out):: percentiles
    real(rk), pointer:: ptr(:,:) => null()
    real(rk), dimension(1:(nx+1)*(nz+1)):: vec
    integer:: ntot, ibord, i, j, n
    
    call associatePtrScal(ptr, what, 0, id)
    ibord = 0
    ntot = nx * nz
    if(id == 2 .or. id == 4) then
       ibord = 1
       ntot = (nx+1) * (nz+1)
    end if

    n = 0
    vec(:) = 0.0
    do j=i0, iend0(2) + ibord
       do i=i0, iend0(1) + ibord
          n = n + 1
          vec(n) = ptr(i, j)
       end do ! i
    end do ! j
    call sort(vec, ntot)  ! vec is ptr values in increasing order
    n = 1
    i = 1
    do
       if(dble(n) / dble(ntot) > percValues(i)) then
          percentiles(i) = vec(n)
          i = i+1
       end if
       if(i > 5) exit ! all 5 percentiles are computed
       n = n + 1
    end do

    nullify(ptr)
    
    return
  end subroutine getPercentiles
  !=======================================================
  ! subroutines useful for multigrid
  !=======================================================
  subroutine restriction(what, level, id, geometric)
    ! restriction of "what" from finer level "level-1" 
    ! to coarser level "level".
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    logical, optional:: geometric
    integer:: iendf(1:2)
    integer:: i, j, ic, jc, kx, kz, nkx, nkz
    real(rk), pointer:: fine(:,:) => null()
    real(rk), pointer:: coarse(:,:) => null()
    real(rk):: c1, c2, c3
    logical:: geom
    
    geom = .false.
    if(present(geometric)) then
       geom = geometric
    end if
    if(geom .and. id /= 3) then
       write(*,*) &
            'ERROR: geometric=.true. implemented only for id=3 in restriction()'
       call endRun()
    end if
    
    if(level > nLevel) then
       write(*,*) &
            'ERROR: level larger than nLevel in restriction()'
       call endRun()
    end if

    ! put coarser level to zero
    call zeroingScal(what, level, id)

    call associatePtrScal(fine, what, level-1, id)
    call associatePtrScal(coarse, what, level, id)

    iendf(:) = iend(level-1, :) !finer level

    if(id == 3) then ! cell center

       if(.not. geom) then
          do j=i0, iendf(2)-1, 2

             jc = (j+i0) / 2 

             do i=i0, iendf(1)-1, 2

                ic = (i+i0) / 2

                coarse(ic,jc) = fourth * ( &
                     fine(i,j) + fine(i+1,j) + &
                     fine(i,j+1) + fine(i+1,j+1) )

             end do
          end do
       else  ! geom
          do j=i0, iendf(2)-1, 2

             jc = (j+i0) / 2 

             do i=i0, iendf(1)-1, 2

                ic = (i+i0) / 2

                coarse(ic,jc) = ( &
                     fine(i,j) * fine(i+1,j) * &
                     fine(i,j+1) * fine(i+1,j+1) )**fourth

             end do
          end do
       end if ! .not.geom or geom
       
    else if(id == 4) then ! corner
       c1 = 0.25 ! center (same corner)
       c2 = 0.125 ! closest 
       c3 = 0.0625 ! farthest

       if(.not. geom) then
          do j=i0, iendf(2)-1, 2 
             !check boundaries (was done until iendf+1 before???)
             
             jc = (j+i0) / 2
             
             do i=i0, iendf(1)-1, 2
                
                ic = (i+i0) / 2
                
                coarse(ic,jc) = c1 * fine(i,j) + &
                     c2 * &
                     (fine(i+1,j)+fine(i-1,j)+fine(i,j+1)+fine(i,j-1)) + &
                     c3 * &
                     (fine(i-1,j-1)+fine(i+1,j-1)+fine(i-1,j+1)+fine(i+1,j+1))
                
             end do
          end do
       else ! geom
          do j=i0, iendf(2)-1, 2 
             !check boundaries (was done until iendf+1 before???)
             
             jc = (j+i0) / 2
             
             do i=i0, iendf(1)-1, 2
                
                ic = (i+i0) / 2
                
                coarse(ic,jc) = fine(i,j)**c1 * &
                     (fine(i+1,j)*fine(i-1,j)*fine(i,j+1)*fine(i,j-1))**c2 * &
                     (fine(i-1,j-1)*fine(i+1,j-1)*fine(i-1,j+1)*fine(i+1,j+1))**c3
                
             end do
          end do
       end if ! geom or not
       
    else if(id == 1 .or. id == 2) then ! faces

       c1 = 0.25 ! same line
       c2 = 0.125 ! other line

       ! for velocities 
       kx = kr(1, id) ! kx=1 if id=1, and kx=0 if id=2
       kz = kr(2, id) ! kz=1 if id=2, and kz=0 if id=1
       nkx = nkr(1, id) ! nkx=0 if id=1, and nkx=1 if id=2
       nkz = nkr(2, id) ! nkz=0 if id=2, and nkz=1 if id=1

       do j=i0, iendf(2)+kz, 2

          jc = (j+i0) / 2

          do i=i0, iendf(1)+kx, 2

             ic = (i+i0) / 2

             coarse(ic,jc) = &
                  c1 * (fine(i,j) + fine(i+nkx, j+nkz)) + &
                  c2 * (fine(i-kx,j-kz) + fine(i-kx+nkx,j-kz+nkz) + &
                  fine(i+kx,j+kz) + fine(i+kx+nkx,j+kz+nkz)) 

          end do
       end do

    end if ! id

    nullify(fine, coarse)

    return
  end subroutine restriction
  !-------------------------------------------------------
  subroutine prolongation(whatOut, whatIn, level, id)
    ! interpolation of "what" from coarser level "level"
    ! to finer level "level-1".
    ! result is put in whatOut.
    implicit none
    character(len=10), intent(in):: whatIn, whatOut
    integer, intent(in):: level, id
    integer:: iendc(1:2) !coarse grid
    integer:: ic, jc, i, j, kh, nkh, kv, nkv
    real(rk), pointer:: fine(:,:) => null()
    real(rk), pointer:: coarse(:,:) => null()
    real(rk):: c1, c2, c3, cpc, cpp

    if(level < 1) then
       write(*,*) &
            'ERROR: level==0 in interpolation()'
       call endRun()
    end if

    call associatePtrScal(fine, whatOut, level-1, id)
    call associatePtrScal(coarse, whatIn, level, id)

    iendc(:) = iend(level, :) ! coarser level

    if(id == 3) then ! cell center
       !note: prolongation at cell center
       ! include bcs

       c1 = 0.5625 ! closest   9/16
       c2 = 0.1875 ! mid       3/16
       c3 = 0.0625 ! farthest  1/16

       do jc=i0-1, iendc(2)

          j = 2*jc + 1 - i0

          do ic=i0-1, iendc(1)

             i = 2*ic + 1 - i0

             fine(i,j) = c1 * coarse(ic,jc) + &
                  c2 * ( coarse(ic+1,jc) + coarse(ic,jc+1) ) + &
                  c3 * coarse(ic+1,jc+1)

             fine(i+1,j) = c1 * coarse(ic+1,jc) + &
                  c2 * ( coarse(ic,jc) + coarse(ic+1,jc+1) ) + &
                  c3 * coarse(ic,jc+1)

             fine(i,j+1) = c1 * coarse(ic,jc+1) + &
                  c2 * ( coarse(ic,jc) + coarse(ic+1,jc+1) ) + &
                  c3 * coarse(ic+1,jc)

             fine(i+1,j+1) = c1 * coarse(ic+1,jc+1) + &
                  c2 * ( coarse(ic+1,jc) + coarse(ic,jc+1) ) + &
                  c3 * coarse(ic,jc)

          end do ! i

       end do ! j

    else if(id == 4) then ! cell corner

       do jc=i0-1, iendc(2) 

          j = 2*jc + 1 - i0

          do ic=i0-1, iendc(1)

             i = 2*ic + 1 - i0

             !center
             fine(i,j) = fourth * ( &
                  coarse(ic,jc) + coarse(ic+1,jc) + &
                  coarse(ic,jc+1) + coarse(ic+1,jc+1) &
                  )

             !z-wall
             fine(i+1,j) = half * ( &
                  coarse(ic+1,jc) + coarse(ic+1,jc+1) )

             !x-wall
             fine(i,j+1) = half * ( &
                  coarse(ic,jc+1) + coarse(ic+1,jc+1) )

             !right-up corner
             fine(i+1,j+1) = coarse(ic+1,jc+1)

          end do !i

       end do !j

    else if(id == 1 .or. id == 2) then ! faces

       kh = kr(1, id)
       kv = kr(2, id)
       nkh = nkr(1, id)
       nkv = nkr(2, id)

       do jc=i0-1, iendc(2)

          j = 2*jc + 1 - i0

          do ic=i0-1, iendc(1)

             i = 2*ic + 1 - i0

             cpc = coarse(ic,jc) + coarse(ic+kh,jc+kv)
             cpp = coarse(ic+nkh,jc+nkv) + &
                  coarse(ic+kh+nkh,jc+kv+nkv)

             fine(i,j) = 0.375 * cpc + 0.125 * cpp

             fine(i+nkh,j+nkv) = 0.125 * cpc + 0.375 * cpp             

             fine(i+kh,j+kv) = &
                  0.75 * coarse(ic+kh,jc+kv) + &
                  0.25 * coarse(ic+kh+nkh,jc+kv+nkv)

             fine(i+kh+nkh,j+kv+nkv) = &
                  0.25 * coarse(ic+kh,jc+kv) + &
                  0.75 * coarse(ic+kh+nkh,jc+kv+nkv)

          end do !i

       end do !j

    end if ! type

    nullify(fine, coarse)

    return
  end subroutine prolongation
  !-------------------------------------------------------
  ! =======================================================
  ! other grid interpolations
  ! =======================================================
  subroutine centerToFace(what, level, id)
    ! put interpolation of grid(level)%what(:,:,3) (center)
    ! on face grid(level)%what(:,:,id) 
    ! (id = 1 or 2, for x- and z-face)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), pointer:: pt(:,:,:) => null()
    integer:: i, j, kx, kz, nkx, nkz

    call associatePtrVec(pt, what, level)

    kx = kr(1, id)
    kz = kr(2, id)
    nkx = nkr(1, id)
    nkz = nkr(2, id)

    do j=i0-nkz, iend(level,2)+1
       do i=i0-nkx, iend(level,1)+1
          pt(i,j,id) = half * &
               ( pt(i, j, 3) + pt(i-kx, j-kz, 3) )
       end do
    end do
    nullify(pt)

    return
  end subroutine centerToFace
  !-------------------------------------------------------
  subroutine cornerToCenter(what, level, geometric)
    ! put interpolation of grid(level)%what(:,:,4) (corner)
    ! to center grid(level)%what(:,:,3)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    logical, optional:: geometric
    real(rk), pointer:: pt(:,:,:) => null()
    integer:: i, j
    logical:: geom

    geom = .false.
    if(present(geometric)) then
       geom = geometric
    end if
    
    call associatePtrVec(pt, what, level)

    if(.not.geom) then
       do j=i0-1, iend(level,2)+1
          do i=i0-1, iend(level,1)+1
             pt(i,j,3) = fourth * ( &
                  pt(i,j,4) + pt(i,j+1,4) + &
                  pt(i+1,j,4) + pt(i+1,j+1,4) )
          end do
       end do
    else
       do j=i0-1, iend(level,2)+1
          do i=i0-1, iend(level,1)+1
             pt(i,j,3) = ( &
                  pt(i,j,4) * pt(i,j+1,4) * &
                  pt(i+1,j,4) * pt(i+1,j+1,4) )**fourth
          end do
       end do
    end if
       
    nullify(pt)

    return
  end subroutine cornerToCenter
  ! -------------------------------------------------------
  subroutine centerToCorner(what, level, geometric)
    ! put interpolation of grid(level)%what(:,:,3) (center)
    ! to corner grid(level)%what(:,:,4).
    ! Includes top and right walls.
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    logical, optional:: geometric
    real(rk), pointer:: pt(:,:,:) => null()
    integer:: i, j
    logical:: geom

    geom = .false.
    if(present(geometric)) then
       geom = geometric
    end if

    call associatePtrVec(pt, what, level)

    if(.not.geom) then
       do j=i0, iend(level,2)+1
          do i=i0, iend(level,1)+1
             pt(i,j,4) = fourth * ( &
                  pt(i-1,j-1,3) + pt(i,j-1,3) + &
                  pt(i-1,j,3) + pt(i,j,3) )
          end do
       end do
    else
       do j=i0, iend(level,2)+1
          do i=i0, iend(level,1)+1
             pt(i,j,4) = ( &
                  pt(i-1,j-1,3) * pt(i,j-1,3) * &
                  pt(i-1,j,3) * pt(i,j,3) )**fourth
          end do
       end do
    end if
       
    nullify(pt)

    return
  end subroutine centerToCorner
  !--------------------------------------------------
  subroutine faceToCorner(what, level, id, plan)
    !put interpolation of grid(level)%what(:,:,id) 
    ! in corners (id=4)
    ! plan: optional-> if present, do interpol only on this plan
    ! (for vsurf mainly)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    integer, optional:: plan
    real(rk), pointer:: pt(:,:,:) => null()
    integer:: i,j, kx, kz, nkx, nkz

    kx = kr(1, id)
    kz = kr(2, id)
    nkx = nkr(1, id)
    nkz = nkr(2, id)

    call associatePtrVec(pt, what, level)

    if(.not.present(plan)) then
       do j=i0-kz, iend(level,2)+1
          do i=i0-kx, iend(level,1)+1
             pt(i,j,4) = half * &
                  ( pt(i-nkx,j-nkz,id) + pt(i,j,id) )
          end do
       end do
    else ! present(plan)
       if(id == 1) then
          do i=i0-kx, iend(level,1)+1
             pt(i, plan, 4) = half * &
                  ( pt(i, plan, id) + pt(i, plan-1, id) )
          end do
       else ! id=2
          do j=i0-kz, iend(level,2)+1
             pt(plan, j, 4) = half * &
                  ( pt(plan, j, id) + pt(plan-1, j, id) )
          end do
       end if
    end if ! present plan

    nullify(pt)

    return
  end subroutine faceToCorner
  !--------------------------------------------------
  subroutine faceToCenter(what, level, id)
    ! interpolation of grid(level)%what(:,:,id)
    ! at cells' centers (id=3)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level, id
    real(rk), pointer:: pv(:,:,:) => null()
    integer:: i, j, kx, kz, nkx, nkz

    call associatePtrVec(pv, what, level)

    kx = kr(1, id)
    kz = kr(2, id)
    nkx = nkr(1, id)
    nkz = nkr(2, id)

    do j=i0-1, iend(level, 2) + nkz
       do i=i0-1, iend(level, 1) + nkx
          pv(i, j, 3) = half * &
               ( pv(i+kx, j+kz, id) + pv(i, j, id) )
       end do
    end do
    nullify(pv)

    return
  end subroutine faceToCenter
  ! =======================================================
  ! OPERATORS 
  ! =======================================================
  subroutine op_dx(whatOut, idout, whatIn, idin, idir, level, mult)
    ! compute the derivative of whatIn(idin) in the direction idir
    ! WARNING: not optimized (two steps/loops for idir=1)
    ! not to be used in solver loops
    implicit none
    character(len=10), intent(in):: whatOut, whatIn
    integer, intent(in):: idout, idin, idir, level
    real(rk), optional:: mult
    real(rk), pointer:: pin(:,:) => null()
    real(rk), pointer:: pout(:,:) => null()
    real(rk), pointer:: pr(:) => null()
    real(rk), pointer:: pa(:) => null()
    real(rk), dimension(1:2):: invd
    real(rk):: coeff
    logical:: cas1, cas2
    integer:: i,j, kx, kz

    call getInvD(invd, level) ! invd(1)=1/dx and invd(2)=1/dz 

    call associatePtrScal(pin, whatIn, level, idin)
    call associatePtrScal(pout, whatOut, level, idout)
    
    kx = kr(idir,1)
    kz = kr(idir,2)

    if(present(mult)) invd(:) = invd(:) * mult

    cas1 = &
         (idin == 1 .and. idout == 3) .or. &
         (idin == 2 .and. idout == 3) .or. &
         (idin == 4 .and. idout == 1) .or. &
         (idin == 4 .and. idout == 2)

    cas2 = & 
         (idin == 1 .and. idout == 4) .or. &
         (idin == 2 .and. idout == 4) .or. &
         (idin == 3 .and. idout == 1) .or. &
         (idin == 3 .and. idout == 2)

    ! do only d/dx and d/dz here:
    if(cas1) then
       do j=i0-1, iend(level,2)
          do i=i0-1, iend(level,1)
             pout(i,j) = invd(idir) * &
                  (pin(i+kx,j+kz) - pin(i,j))
          end do
       end do       
    else if(cas2) then
       do j=i0, iend(level,2)+1
          do i=i0, iend(level,1)+1
             pout(i,j) = invd(idir) * &
                  (pin(i,j) - pin(i-kx,j-kz))
          end do
       end do
    else
       write(*,'(a)') "ERROR in op_dx: &
            & idin/idout not appropriate"
       call endRun()
    end if

    if(curvedGeometry) then ! -----------------------

       call associatePtrLine(pa, alpha, level)
       call associatePtrLine(pr, r, level)

       ! for idir=1, divide results by r (at the good position)
       if(idir == 1) then
          if(idout == 4 .or. idout == 2) then
             do j=i0-1, iend(level,2)+1
                coeff = 1.0 / pr(j)
                do i=i0-1, iend(level,1)+1
                   pout(i,j) = pout(i,j) * coeff
                end do ! i
             end do ! j
          else ! idout=1 or idout=3
             do j=i0-1, iend(level,2)+1
                coeff = 1.0 / (pr(j) * (1.0 + pa(j)))
                do i=i0-1, iend(level,1)+1
                   pout(i,j) = pout(i,j) * coeff
                end do ! i
             end do ! j
          end if
       end if ! idir == 1

       nullify(pin, pout, pa, pr)

    end if ! curvedGeometry --------------------------
 
    return
  end subroutine op_dx
  ! ===================================================
  subroutine sf_to_vel(whatsf, whatvel, level, where)
    ! Compute velocities vel from streamfunction sf
    !
    ! vx = -1 / r^(ed-1) * d(r^(ed-1) * sf) / dz
    ! and  vz = 1/r * d(sf)/dphi in curvedGeometry
    implicit none
    integer, intent(in):: level
    character(len=10), intent(in):: whatsf, whatvel
    character(len=20), optional:: where
    real(rk), pointer:: psf(:,:) => null()
    real(rk), pointer:: pv(:,:,:) => null()
    real(rk), pointer:: pr(:) => null()
    real(rk), pointer:: pa(:) => null()
    real(rk):: invd(1:2), coeffDen, coeffNum
    integer:: i,j, iez
    logical:: surf, horiz, all

    all = .true.
    surf = .false.
    horiz = .false.
    if(present(where)) then
       if(where == surface) then
          all = .false.
          surf = .true.
       else if(where == horizontal) then
          all = .false.
          horiz = .true.
       end if
    end if

    call getInvD(invd, level)

    call associatePtrScal(psf, whatsf, level, 4)
    call associatePtrVec(pv, whatvel, level)

    if(curvedGeometry) then ! ---------------------------------------

       call associatePtrLine(pa, alpha, level)
       call associatePtrLine(pr, r, level)

       if(all) then
          ! Vx ------------
          do j=i0-2, iend(level,2)+1
             coeffDen = - invd(2) / (1.0 + pa(j))**(ed-1)
             coeffNum = (1.0 + 2.0 * pa(j))**(ed-1)
             do i=i0-2, iend(level,1)+2
                pv(i,j,1) = coeffDen * (coeffNum * psf(i,j+1) - psf(i,j))
             end do
          end do
          ! Vz ------------
          do j=i0-2, iend(level,2)+2
             coeffDen = invd(1) / pr(j)
             do i=i0-2, iend(level,1)+1
                pv(i,j,2) = coeffDen * (psf(i+1,j) - psf(i,j))
             end do
          end do
       else if(surf) then
          iez = iend(level, 2)
          coeffDen = - invd(2) / (1.0 + pa(iez))**(ed-1)
          coeffNum = (1.0 + 2.0 * pa(iez))**(ed-1)
          do i=i0-2, iend(level, 1)+2
             pv(i, iez, 1) = coeffDen * &
                  (coeffNum * psf(i,iez+1) - psf(i,iez))
          end do
       else if(horiz) then
          do j=i0-1, iend(level,2)+2
             coeffDen = - invd(2) / (1.0 + pa(j))**(ed-1)
             coeffNum = (1.0 + 2.0 * pa(j))**(ed-1)
             do i=i0-2, iend(level,1)+2
                pv(i,j,1) = coeffDen * (coeffNum * psf(i,j+1) - psf(i,j))
             end do
          end do
       end if ! all or not

       nullify(pa, pr)

    else ! CARTESIAN ------------------------------------

       if(all) then
          ! Vx ------------
          do j=i0-2, iend(level,2)+1
             do i=i0-2, iend(level,1)+2
                pv(i,j,1) = -invd(2) * (psf(i,j+1) - psf(i,j))
             end do
          end do
          ! Vz ------------
          do j=i0-2, iend(level,2)+2
             do i=i0-2, iend(level,1)+1
                pv(i,j,2) = invd(1) * (psf(i+1,j) - psf(i,j))
             end do
          end do
       else if(surf) then
          iez = iend(level, 2)
          do i=i0-2, iend(level, 1)+2
             pv(i, iez, 1) = -invd(2) * &
                  (psf(i, iez+1) - psf(i, iez)) 
          end do
       else if(horiz) then
          do j=i0-1, iend(level,2)+2
             do i=i0-1, iend(level,1)+1
                pv(i,j,1) = -invd(2) * (psf(i,j+1) - psf(i,j))
             end do
          end do
       end if ! all or not

    end if ! curvedGeometry or not ----------------------

    nullify(psf, pv)

    return
  end subroutine sf_to_vel
  ! =======================================================
  subroutine sf_to_sist(level, id)
        ! computes the opposite of the 
    ! Second Invariant of the Strainrate Tensor (SIST) from sf.
    ! sf at cells' corners (id=4)
    ! and sist at cells' centers (id=3) or corners (id=4)
    ! SIST is real definition : = 1/2 * (Tr(E)^2 - Tr(EE))
    ! Also updates exz, exx, and ezz at cells' corners.
    ! Cylindrical: exx = - ezz, only ezz and exz needed for sist.
    ! Spherical annulus: exx must also be computed, and exp is different
    implicit none
    integer, intent(in):: level, id
    real(rk), pointer:: psf(:,:) => null()
    real(rk), pointer:: pezz(:,:) => null()
    real(rk), pointer:: pexz(:,:) => null()
    real(rk), pointer:: pexx(:,:) => null()
    real(rk), pointer:: psist(:,:) => null()
    real(rk), pointer:: pr(:) => null()
    real(rk), pointer:: pa(:) => null()
    real(rk), dimension(1:2):: invd, invd2
    real(rk):: coeff1, coeff2, c1, c4, c5
    integer:: i, j, ibord
    
    if(.not.(id == 3 .or. id == 4)) then
       write(*,'(a,i3,a)') "ERROR: id = ", id, &
            " not implemented in sf_to_sist()"
       call endRun()
    end if

    call getInvD(invd, level)
    call getInvD2(invd2, level)
    call associatePtrScal(psf, sf, level, 4)
    ! BCS on sf must be up-to-date !!!!!!!!!!

    if(curvedGeometry) then ! -------------------
       
       call associatePtrLine(pr, r, level)
       call associatePtrLine(pa, alpha, level)
    
       ! Compute Ezz at cells' centers (id=3)
       call associatePtrScal(pezz, ezz, level, 3)
       do j=i0-1, iend(level, 2)+1
          coeff1 = 2.0 / pr(j) * invd(1) * invd(2) ! 2/(r*dphi*dz)
          coeff2 = 1.0 / (1.0 + 2.0 * pa(j))
          do i=i0-1, iend(level, 1)+1
             pezz(i, j) = coeff1 * (coeff2 * (psf(i+1, j+1) - psf(i, j+1)) - &
                  psf(i+1, j) + psf(i, j))
          end do ! i
       end do ! j
    
       ! Exx, also at cells' centers (id=3)
       ! If cylindrical: Exx = - Ezz
       call associatePtrScal(pexx, exx, level, 3)
       if(coordinates == "cylindrical") then
          pexx = -pezz
       else ! spherical annulus
          do j=i0-1, iend(level, 2)+1
             ! coeff1 = -2 / (r*dphi*dz*(1+alpha))
             coeff1 = -2.0 / (pr(j) * (1.0 + pa(j))) * invd(1) * invd(2)
             do i=i0-1, iend(level, 1)+1
                pexx(i, j) = coeff1 * (psf(i+1,j+1) - psf(i,j+1) &
                     - psf(i+1,j) + psf(i,j))
             end do ! i
          end do ! j
       end if ! cylindrical or spherical annulus
    
       ! if id=4: interpolation on cells' corners of Ezz and Exx
       if(id == 4) then
          nullify(pexx, pezz)  ! were defined at id=3
          call centerToCorner(exx, level)
          call centerToCorner(ezz, level)
       end if

       ! Compute Exz at cells' corners (id=4)
       call associatePtrScal(pexz, exz, level, 4)
       coeff2 = - invd2(2)  ! - 1/dr^2
       do j=i0, iend(level, 2) + 1  ! required for cornerToCenter if id=3
          coeff1 = invd2(1) / pr(j)**2  ! 1/(r^2*dphi^2)
          ! ed=1 for cylindrical and ed=2 for spherical annulus
          c1 = - (1.0 / (1.0+pa(j))**ed + 1.0 / (1.0-pa(j))**ed) ! coeff for sf(i,j)
          c4 = (1.0 + 2.0*pa(j))**(ed-1) / ((1.0 + pa(j)))**ed ! coeff for sf(i,j+1)
          c5 = (1.0 - 2.0*pa(j))**(ed-1) / ((1.0 - pa(j)))**ed ! coeff for sf(i,j-1)
          do i=i0, iend(level, 1) + 1          
             pexz(i,j) = coeff1 * (psf(i+1,j) + psf(i-1,j) - 2.0 * psf(i,j)) + &
                  coeff2 * (c4*psf(i,j+1) + c5*psf(i,j-1) + c1*psf(i,j))
          end do ! i
       end do ! j
    
       ! if id=3, interpolation at cells' centers
       if(id == 3) then
          nullify(pexz) ! was defined at id=4
          call cornerToCenter(exz, level)
       end if

       ! compute -sist
       if(id == 3) then
          call associatePtrScal(pexz, exz, level, 3)
       else ! id=4
          call associatePtrScal(pexx, exx, level, 4)
          call associatePtrScal(pezz, ezz, level, 4)
       end if
    
       call associatePtrScal(psist, sist, level, id)
       ibord = 0
       if(id == 4) then
          ibord = 1
       end if

       do j=i0, iend(level, 2) + ibord
          do i=i0, iend(level, 1) + ibord
             psist(i, j) = pezz(i,j)*pezz(i,j) + pexz(i,j)*pexz(i,j)
          end do ! i
       end do ! j
       if(coordinates == "spherical") then
          do j=i0, iend(level, 2) + ibord
             do i=i0, iend(level, 1) + ibord
                psist(i, j) = psist(i,j) + pexx(i,j)*pexx(i,j) + &
                     pezz(i,j)*pexx(i,j)
             end do ! i
          end do ! j       
       end if
    
       nullify(pexx, pezz, pexz, pa, pr, psist)
       
    else ! Cartesian -----------------------------
       
       ! compute exz = 0.5 * (d^2(sf)/dx^2-d^2(sf)/dz^2)
       ! at cells' corners (id=4)
       call associatePtrScal(pexz, exz, level, 4)
       do j=i0, iend(level, 2) + 1  ! required for cornerToCenter if id=3
          do i=i0, iend(level, 1) + 1
             pexz(i, j) = &
                  invd2(1) * (psf(i+1, j) + psf(i-1, j) - 2.0 * psf(i, j)) - &
                  invd2(2) * (psf(i, j+1) + psf(i, j-1) - 2.0 * psf(i, j))
             pexz(i, j) = half * pexz(i, j)
          end do ! i
       end do ! j
       if(id == 3) then
          nullify(pexz)
          call cornerToCenter(exz, level)  ! from id=4 to id=3
       end if
       
       ! compute exx = - d^2(sf)/(dx.dz)
       ! at cells' centers (id=3)
       call associatePtrScal(pexx, exx, level, 3)
       ! bcs on sf must be up-to-date
       do j=i0-1, iend(level, 2)+1
          do i=i0-1, iend(level, 1)+1
             pexx(i, j) = - invd(1)*invd(2) * &
                  (psf(i+1, j+1) - psf(i, j+1) - psf(i+1, j) + psf(i, j))
          end do
       end do
       if(id == 4) then
          nullify(pexx)
          !call applyBcs(exx, level, 3)
          call centerToCorner(exx, level) ! from id=3 to id=4
          ! check if should not redo applyBcs
       end if
       
       ! compute -sist = (exx^2 + exz^2) 
       if (id == 3) then
          call associatePtrScal(pexz, exz, level, 3)
       else
          call associatePtrScal(pexx, exx, level, 4)
       end if
       
       call associatePtrScal(psist, sist, level, id)
       ibord = 0
       if (id == 4) then
          ibord = 1
       end if
       do j=i0, iend(level, 2) + ibord
          do i=i0, iend(level, 1) + ibord
             psist(i, j) = pexx(i,j)*pexx(i,j) + pexz(i,j)*pexz(i,j)
          end do
       end do

       nullify(pexz, pexx, psist)

    end if ! curvedGeometry or Cartesian -------------------

    nullify(psf)
    
    return
  end subroutine sf_to_sist
  ! ===================================================
    subroutine divergence(whatout, whatin, level)
    ! Divergence on any vector, written as the velocity
    ! (x-component with id=1 and z-component with id=2,
    ! result put at cells' center: id=3)
    implicit none
    character(len=10), intent(in):: whatin, whatout
    integer, intent(in):: level
    real(rk), pointer:: pin(:,:,:) => null()
    real(rk), pointer:: pout(:,:) => null()
    real(rk), pointer:: pa(:) => null()
    real(rk), pointer:: pr(:) => null()
    real(rk), dimension(1:2):: invd
    real(rk):: ch, cv, cvp
    integer:: i,j

    call getInvD(invd, level)
    call associatePtrVec(pin, whatin, level)
    call associatePtrScal(pout, whatout, level, 3)

    if(curvedGeometry) then ! --------------------------
       
       call associatePtrLine(pr, r, level)
       call associatePtrLine(pa, alpha, level)
    
       do j=i0, iend(level,2)
          ch = invd(1) / (pr(j) * (1.0 + pa(j))) ! 1/(r*(1+alpha)*dphi)
          cv = invd(2) / (1.0 + pa(j))**ed ! 1/((1+alpha)^d * dr)
          cvp = (1.0 + 2.0*pa(j))**ed  ! multiplies only pin(i,j+1)
          do i=i0, iend(level,1)
             pout(i,j) = ch * (pin(i+1,j,1) - pin(i,j,1)) + &
                  cv * (cvp * pin(i,j+1,2) - pin(i,j,2))
          end do
       end do
    
       nullify(pa, pr)

    else ! Cartesian ------------------------------------

       do j=i0, iend(level,2)
          do i=i0, iend(level,1)
             pout(i,j) = &
                  invd(1) * (pin(i+1,j,1) - pin(i,j,1)) + &
                  invd(2) * (pin(i,j+1,2) - pin(i,j,2))
          end do
       end do
       
    end if  ! curvedGeometry or Cartesian -----------------

    nullify(pin, pout)
    
    return
  end subroutine divergence
  ! ========================================================
    subroutine laplacian(whatout, whatin, level)
    ! laplacian of a scalar, computed only for id=4
    ! for both input and output
    implicit none
    character(len=10), intent(in):: whatout, whatin
    integer, intent(in):: level
    real(rk), pointer:: pin(:,:) => null()
    real(rk), pointer:: pout(:,:) => null()
    real(rk), pointer:: pa(:) => null()
    real(rk), pointer:: pr(:) => null()
    real(rk), dimension(1:2):: invd2
    real(rk):: ch, cv, cvp, cvm
    integer:: i,j

    call getInvD2(invd2, level)
    call associatePtrScal(pin, whatin, level, 4)
    call associatePtrScal(pout, whatout, level, 4)

    if(curvedGeometry) then ! ------------------------
       
       call associatePtrLine(pr, r, level)
       call associatePtrLine(pa, alpha, level)

       do j=i0, iend(level, 2)+1
          ch = invd2(1) / pr(j)**2  ! 1/(r^2 * dphi^2)
          cv = 2.0 * invd2(2) ! 2/dr^2
          cvp = (1.0 + ed * pa(j)) * cv
          cvm = (1.0 - ed * pa(j)) * cv
          do i=i0, iend(level, 1)+1
             pout(i, j) = ch * (pin(i+1, j) + pin(i-1, j) - 2.0 * pin(i, j)) + &
                  cvp * pin(i, j+1) + cvm * pin(i, j-1) - cv * pin(i, j)
          end do ! i
       end do !j
       
       nullify(pr, pa)
       
    else ! Cartesian -------------------------------

       do j=i0, iend(level, 2)+1
          do i=i0, iend(level, 1)+1
             pout(i, j) = invd2(1) * (pin(i+1, j) + pin(i-1, j) - 2.0 * pin(i, j)) + &
                  invd2(2) * (pin(i, j+1) + pin(i, j-1) - 2.0 * pin(i, j)) 
          end do ! i
       end do !j

    end if ! curvedGeometry or Cartesian -----------------
    
    nullify(pin, pout) 
    
    return
  end subroutine laplacian
  ! ========================================================
  subroutine coeffProduct(what, what1, what2, level, id, alpha)
    implicit none
    character(len=10), intent(in):: what, what1, what2
    integer, intent(in):: level, id
    real(rk), intent(in):: alpha
    real(rk), pointer:: p(:, :) => null()
    real(rk), pointer:: p1(:, :) => null()
    real(rk), pointer:: p2(:, :) => null()
    integer:: i,j
    real(rk):: beta

    call associatePtrScal(p, what, level, id)
    call associatePtrScal(p1, what1, level, id)
    call associatePtrScal(p2, what2, level, id)

    beta = 1.0 - alpha
    do j=i0, iend(level, 2) + 1
       do i=i0, iend(level, 1) + 1
          p(i, j) = p1(i, j)**alpha * p2(i, j)**beta
       end do
    end do
    
    nullify(p, p1, p2)
    
    return
  end subroutine coeffProduct
  ! ========================================================
  subroutine maxDifference(maxval, what1, what2, level, id)
    ! computes the maximum abs difference between what1 and what2
    implicit none
    character(len=10), intent(in):: what1, what2
    integer, intent(in):: level, id
    real(rk), intent(out):: maxval
    real(rk), pointer:: p1(:, :) => null()
    real(rk), pointer:: p2(:, :) => null()
    integer:: i, j

    call associatePtrScal(p1, what1, level, id)
    call associatePtrScal(p2, what2, level, id)
    
    maxval = 0.0
    do j=i0, iend(level, 2)
       do i=i0, iend(level, 1)
          maxval = max(maxval, abs(p1(i, j) - p2(i, j)))
       end do ! i
    end do ! j

    nullify(p1, p2)

    return
  end subroutine maxDifference
  ! =======================================================
  subroutine getRmsChange(rmsChange, what1, what0, level, id, relative)
    ! check the rms of the modification of "what".
    ! WARNING: does not take into account the different sizes of cells
    implicit none
    character(len=10), intent(in):: what0, what1
    integer, intent(in):: level, id
    real(rk), intent(out):: rmsChange
    logical, optional:: relative
    real(rk), pointer:: p0(:,:) => null()
    real(rk), pointer:: p1(:,:) => null()
    integer:: i, j
    logical:: doRelative
    real(rk):: den
    
    call associatePtrScal(p0, what0, level, id)
    call associatePtrScal(p1, what1, level, id)

    doRelative = .false.
    if(present(relative)) then
       doRelative = relative
    end if
    
    rmsChange = 0.0
    if(doRelative) then
       do j=i0, iend(level, 2)
          do i=i0, iend(level, 1)
             den = (p1(i,j) + p0(i,j))
             if(den /= 0.0) then
                rmsChange = rmsChange + &
                     (2.0 / den * (p1(i,j) - p0(i,j)))**2
             end if
          end do ! i
       end do ! j
    else
       do j=i0, iend(level, 2)
          do i=i0, iend(level, 1)
             rmsChange = rmsChange + (p1(i,j) - p0(i,j))**2
          end do ! i 
       end do ! j
    end if
    nullify(p0, p1)

    rmsChange = sqrt(rmsChange / grid(level)%ntot)

    return
  end subroutine getRmsChange
  !
  ! =======================================================
  ! Coefficients for geometric parameters
  ! ========================================================
  !
  subroutine associatePtrLine(ptr, what, level)
    implicit none
    character(len=10), intent(in):: what
    integer, intent(in):: level
    real(rk), pointer:: ptr(:)

    if(what == alpha) then
       ptr => grid(level)%alpha
    else if(what == r) then
       ptr => grid(level)%r
    else if(what == phi) then
       ptr => grid(level)%phi
    else if(what == extent) then
       ptr => grid(level)%extent
    else if(what == gtop) then
       ptr => grid(level)%gtop
    else if(what == gbot) then
       ptr => grid(level)%gbot
    else if(what == gmid) then
       ptr => grid(level)%gmid
    else if(what == g1top) then
       ptr => grid(level)%g1top
    else if(what == g1bot) then
       ptr => grid(level)%g1bot
    else if(what == gcentertop) then
       ptr => grid(level)%gcentertop
    else if(what == gcenterbot) then
       ptr => grid(level)%gcenterbot
    else if(what == gcornertop) then
       ptr => grid(level)%gcornertop
    else if(what == gcornerbot) then
       ptr => grid(level)%gcornerbot
    else if(what == g4) then
       ptr => grid(level)%g4
    else if(what == g5) then
       ptr => grid(level)%g5
    else if(what == g12) then
       ptr => grid(level)%g12
    else if(what == g13) then
       ptr => grid(level)%g13
    end if
    
    return
  end subroutine associatePtrLine
  ! ========================================================
end module grids
! =======================================================
