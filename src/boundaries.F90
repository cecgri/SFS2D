! =======================================================
! author:cecile.grigne@univ-brest.fr
! this file is part of SFS2D
! =======================================================
subroutine bcs(what, type, id, where, level, val)
  ! ----------------------------------------
  ! General boundary conditions.
  ! -type=dirichlet: imposes the fixed value "val".
  ! -type=neumann: imposes the fixed derivative "val".
  ! -id: what(id) is modified in grid(level)%what(:,:,id)
  ! -where: can be "top", "bot", "left" or "right"
  ! -level: level in the multigrid.
  ! ----------------------------------------
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  character(len=20), intent(in):: type, where
  real(rk), intent(in):: val
  integer, intent(in):: id, level
  integer:: iex, iez
  real(rk):: d(1:2), buff

  if(type /= dirichlet) then
     call getD(d, level) ! needed for neumann
  end if
  iex = iend(level,1)
  iez = iend(level,2)

  if(where == top) then

     if(type == dirichlet) then
        
        if(id == 2 .or. id == 4) then !z-face or corners
           call setValuePlan(what, level, id, 2, iez+1, val)
        else if(id == 1 .or. id == 3) then !x-face or center
           call setMiddleValuePlan(what, level, id, 2, iez, iez+1, val)
        end if ! id
        
     else if(type == neumann) then
        
        ! same thing for all positions in the cell
        call givePlan(what, level, id, 2, iez, iez+1)
        if(val /= zero) then
           ! e.g v(:,n+1) = v(:n) + dr*val
           buff = d(2) * val
           call addValuePlan(what, level, id, 2, iez+1, buff)
        end if
        
     else if(type == secondDerivative) then
        
        if(id == 2 .or. id == 4) then ! r-faces or corners
           ! what(iez+2) = 2*val - what(iez)
           call setMiddleValuePlan(what, level, id, 2, iez, iez+2, val)
        else
           ! what(iez+2) =  2*val - what(iez-1)
           call setMiddleValuePlan(what, level, id, 2, iez-1, iez+2, val)
        end if
        
     else           
       
        write(*,*) "ERROR: type not recognized in bcs()"
        write(*,*) "   top ", what, "   ", type
        
        call endRun()
     end if

  else if(where == bot) then

     if(type == dirichlet) then

        if(id == 2 .or. id == 4) then ! z-face or corner
           call setValuePlan(what, level, id, 2, i0, val)
           call setValuePlan(what, level, id, 2, i0-1, val)
        else if(id == 1 .or. id == 3) then
           call setMiddleValuePlan(what, level, id, 2, i0, i0-1, val)
        end if ! id

     else if(type == neumann) then

        ! same thing for all positions in the cell
        call givePlan(what, level, id, 2, i0, i0-1)
        if(val /= zero) then
           buff = - d(2) * val
           call addValuePlan(what, level, id, 2, i0-1, buff)
        end if

     else if(type == secondDerivative) then

        if(id == 2 .or. id == 4) then ! z-faces or corners
           ! what(i0-1) = 2*val - what(i0+1)
           call setMiddleValuePlan(what, level, id, 2, i0+1, i0-1, val)
        else
           ! what(i0-2) =  2*val - what(i0-1)
           call setMiddleValuePlan(what, level, id, 2, i0+1, i0-2, val)
        end if

     else 
        write(*,*) "ERROR: type not recognized in bcs()"
        write(*,*) "    bot ", what, "   ", type
        call endRun()
     end if

  else if(where == left) then
     ! NOTE: periodic boundary conditions must be taken care of 
     ! with periodicSides() before

     if(type == dirichlet) then

        if(id == 1 .or. id == 4) then !phi-face or corner
           call setValuePlan(what, level, id, 1, i0, val)
           call setValuePlan(what, level, id, 1, i0-1, val) !???? 
        else if(id == 2 .or. id == 3) then !r-face or center
           call setMiddleValuePlan(what, level, id, 1, i0, i0-1, val)
        end if ! id

     else if(type == neumann) then

        call givePlan(what, level, id, 1, i0, i0-1)
        if(val /= zero) then
           write(*,*) &
                "WARNING in bcs(): Neumann condition horizontally with non-zero val is not implemented."
        end if

     end if ! type

  else if(where == right) then
     
     if(type == dirichlet) then

        if(id == 1 .or. id == 4) then !x-face or corner
           call setValuePlan(what, level, id, 1, iex+1, val)
        else if(id == 2 .or. id == 3) then !z-face or center
           call setMiddleValuePlan(what, level, id, 1, iex, iex+1, val)
        end if ! id

     else if(type == neumann) then
        call givePlan(what, level, id, 1, iex, iex+1)
        if(val /= zero) then
           write(*,*) &
                "WARNING in bcs(): Neumann condition horizontally with non-zero val is not implemented."
        end if

     end if ! type

  end if !! where == top, bot, left or right

  return
end subroutine bcs
! ==================================================
subroutine periodicSides(what, level, id)
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level, id
  real(rk), pointer:: ptr(:,:) => null()
  integer:: iex

  call associatePtrScal(ptr, what, level, id)
  iex = iend(level, 1)
  ptr(i0-1, :) = ptr(iex, :)
  ptr(i0-2, :) = ptr(iex-1, :)
  ptr(iex+1, :) = ptr(i0, :)
  ptr(iex+2, :) = ptr(i0+1, :)
  nullify(ptr)

  return
end subroutine periodicSides
! ==================================================
subroutine bcsTemp(what, level, id)
  ! imposes boundary conditions for temperature
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level, id

  call bcs(what, topT, id, top, level, topValueT)
  call bcs(what, botT, id, bot, level, botValueT)

  if(leftT == periodic) then
     call periodicSides(what, level, id)
  else
     call bcs(what, leftT, id, left, level, leftValueT)
     call bcs(what, rightT, id, right, level, rightValueT)
  end if

  return
end subroutine bcsTemp
! ==================================================
subroutine bcsEta(what, level, id)
  ! imposes boundary conditions for viscosity.
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level, id

  if (id == 3 .or. id == 4) then
     call bcs(what, neumann, id, top, level, zero)
     call bcs(what, neumann, id, bot, level, zero)
     !! TEST
!!$     call bcs(what, dirichlet, id, top, level, etaTop)
!!$     call bcs(what, dirichlet, id, bot, level, etaBot)
     !!!! 
  else 
     write(*,*) "ERROR: bcsEta implemented only for id=3 or 4"
     call endRun()
  end if

  if(leftT == periodic) then
     call periodicSides(what, level, id)
  else
     call bcs(what, neumann, id, left, level, zero)
     call bcs(what, neumann, id, right, level, zero)
  end if

  return
end subroutine bcsEta
! =======================================================
subroutine bcsVel(what, level, id)
  ! imposes boundary conditions for vx and vz
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level, id
  real(rk), pointer:: pv(:,:) => null()
  integer:: i, iex, iez
  real(rk):: a, coeffA

  if(.not. (id == 1 .or. id == 2)) then
     write(*,*) "ERROR: bcsVel only implemented for id=1 or 2"
     call endRun()
  end if   ! id=1, 2
   
  call associatePtrScal(pv, what, level, id) ! 

  iex = iend(level, 1)
  iez = iend(level, 2)
  
  coeffA = 1.0  ! default value for Cartesian, computed otherwise
  ! ----------
  ! TOP
  ! ----------
  if(id == 1) then  ! horizontal velocity
     if(topV == freeslip) then
        if(curvedGeometry) then
           a = grid(level)%alpha(iez+1)
           coeffA = (1.0+a) / (1.0-a)
        end if ! curvedGeometry
        do i=i0-1, iex+1                 
           pv(i, iez+1) = pv(i, iez) * coeffA
        end do ! i
     else if(topV == dirichlet .and. level == 0) then
        do i=i0-1, iex+1
           pv(i, iez+1) = 2.0 * topValueV - pv(i, iez)
        end do ! i
     else if(topV == "file" .and. level == 0) then
        do i=i0-1, iex+1
           pv(i, iez+1) = 2.0 * vxtop(i) - pv(i, iez)
        end do  ! i
     else ! can only be level>0 and not freeslip
        do i=i0-1, iex+1
           pv(i, iez+1) = -pv(i, iez)
        end do
     end if
  else if (id == 2) then  ! vertical velocity
     ! = 0 always
     call bcs(what, dirichlet, id, top, level, zero)
  end if
  ! ----------
  ! BOTTOM
  ! ----------
  if(id == 1) then  ! horizontal velocity
     if(botV == freeslip) then
        if(curvedGeometry) then
           a = grid(level)%alpha(i0)
           coeffA = (1.0-a) / (1.0+a)
        end if ! curvedGeometry
        do i=i0-1, iex+1                 
           pv(i, i0-1) = pv(i, i0) * coeffA 
        end do ! i
     else if(botV == dirichlet .and. level == 0) then
        do i=i0-1, iex+1
           pv(i, i0-1) = 2.0 * botValueV - pv(i, i0)
        end do ! i
     else ! can only be level>0 and not freeslip
        do i=i0-1, iex+1
           pv(i, i0-1) = -pv(i, i0)
        end do
     end if
  else if (id == 2) then  ! vv, ur, vertical velocity
     ! = 0 always
     call bcs(what, dirichlet, id, bot, level, zero)
  end if

  ! ----------
  ! SIDES
  ! ----------
  if(leftV == periodic) then
     call periodicSides(what, level, id)
  end if

  ! ----------
  ! LEFT
  ! ----------
  if(leftV /= periodic) then
     if (id == 1) then ! horizontal velocity
        call bcs(what, dirichlet, id, left, level, zero)
     else if (id == 2) then
        if(level > 0.0 .or. leftV == neumann .or. leftV == dirichlet) then
           call bcs(what, leftV, id, left, level, leftValueV)
           !warning: works only if leftValueV=0 for neumann
        else if(level == 0) then
           call bcs(what, leftV, id, left, level, zero)
        end if
     end if ! id = 1,2
  end if
  ! ----------
  ! RIGHT
  ! ----------
  if(rightV /= periodic) then
     if(id == 1) then! horizontal velocity vh/uphi = 0
        call bcs(what, dirichlet, id, right, level, zero)
     else if(id == 2) then
        if(level > 0.0 .or. rightV == neumann .or. rightV == dirichlet) then
           call bcs(what, rightV, id, right, level, rightValueV)
           !warning: works only if rightValueV=0 for neumann
        else if(level == 0) then
           call bcs(what, rightV, id, right, level, zero)
        end if
     end if  ! id=1, 2
  end if
  
  nullify(pv)
  
  return
end subroutine bcsVel
! ======================================================= 
subroutine bcsSimple(what, level, id)
  ! zero value on top and bot
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level, id

  call bcs(what, dirichlet, id, top, level, zero)
  call bcs(what, dirichlet, id, bot, level, zero)

  if(leftV == periodic) then
     call periodicSides(what, level, id)
  else
     ! simple='reflecting' -> neumann
     call bcs(what, neumann, id, left, level, zero)
     call bcs(what, neumann, id, right, level, zero)
  end if

  return
end subroutine bcsSimple
! =======================================================
subroutine bcsSf(what, level)
  ! only implemented for id=4
  ! NOTE: for Cartesian geometry, coefficients
  !    topfreeslip, topdirichlet... are set to one in 
  !    geometry.F90
  use general
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level
  integer:: iex, iez, i, j
  integer, parameter:: id=4
  real(rk), pointer:: ptr(:,:) => null()
  real(rk), pointer:: pr(:) => null()
  real(rk):: d(1:2), val

  iex = iend(level, 1)
  iez = iend(level, 2)
  call getD(d, level) ! dx and dz
  call associatePtrScal(ptr, what, level, id)

  if(topV == freeslip) then  ! free-slip
     do i=i0-1, iex+1
        ptr(i, iez+1) = 0.0 
        ptr(i, iez+2) = -ptr(i, iez) * grid(level)%topfreeslip
     end do ! i
  else if(topV == dirichlet) then  ! rigid
     if(what == sf) then
        val = topValueV
     else
        val = zero
     end if
     do i=i0-1, iex+1
        ptr(i, iez+1) = 0.0
        ptr(i, iez+2) = -2.0 * val * d(2) + &
             grid(level)%topdirichlet(2) * ptr(i, iez)
        ptr(i, iez+2) = grid(level)%topdirichlet(1) * ptr(i, iez+2)
     end do
  else if(topV == "file") then
     if(what == error .or. what == res) then
        do i=i0-1, iex+1
           ptr(i, iez+1) = 0.0
           ptr(i, iez+2) = grid(level)%topdirichlet(1) * &
                grid(level)%topdirichlet(2) * ptr(i, iez)
        end do
     else
        do i=i0-1, iex+1
           ptr(i, iez+1) = 0.0
           ptr(i, iez+2) = -2.0 * vxtop(i) * d(2) + &
                grid(level)%topdirichlet(2) * ptr(i, iez)
           ptr(i, iez+2) = grid(level)%topdirichlet(1) * ptr(i, iez+2)
        end do
     end if
  end if

  if(botV == freeslip) then  ! free-slip
     do i=i0-1, iex+1
        ptr(i, i0) = 0.0 
        ptr(i, i0-1) = -ptr(i, i0+1) * grid(level)%botfreeslip
     end do ! i
  else if(botV == dirichlet) then  ! rigid
     if(what == sf) then
        val = botValueV
     else
        val = zero
     end if
     do i=i0-1, iex+1
        ptr(i, i0) = 0.0
        ptr(i, i0-1) = 2.0 * val * d(2) + &
             grid(level)%botdirichlet(2) * ptr(i, i0+1)
        ptr(i, i0-1) = grid(level)%botdirichlet(1) * ptr(i, i0-1)
     end do
  end if

 if(leftV == periodic) then ! rightV = periodic too then
     call periodicSides(what, level, id) ! 2 points here too
  else
     call associatePtrLine(pr, r, level)

     if(leftV == dirichlet) then
        if(what == sf) then
           val = leftValueV
        else
           val = zero
        end if
        do j=i0-1, iez+1
           ptr(i0, j) = 0.0
           ptr(i0-1, j) = ptr(i0+1, j) - 2.0 * val * d(1) * pr(j)
        end do
     else if(leftV == neumann) then
        do j=i0-1, iez+1
           ptr(i0, j) = 0.0
           ptr(i0-1, j) = -ptr(i0+1, j)
        end do
     end if

     if(rightV == dirichlet) then
        if(what == sf) then
           val = rightValueV
        else
           val = zero
        end if
        do j=i0-1, iez+1
           ptr(iex+1, j) = 0.0
           ptr(iex+2, j) = ptr(iex, j) + 2.0 * val * d(1) * pr(j)
        end do
     else if(rightV == neumann) then
        do j=i0-1, iez+1
           ptr(iex+1, j) = 0.0
           ptr(iex+2, j) = -ptr(iex, j)
        end do
     end if
   
     nullify(pr)
  end if ! not periodic

  nullify(ptr)

  return
end subroutine bcsSf
! =======================================================
subroutine applyBcs(what, level, id)
  use grids
  implicit none
  character(len=10), intent(in):: what
  integer, intent(in):: level, id

  if(what == "temp") then
     call bcsTemp(what, level, id)
  else if (what == "etaT" .or. what == "etaEff" .or. what == "eta")&
       & then
     call bcsEta(what, level, id)
  else if (what == "omega") then
     call bcsSimple(what, level, id)
  else if (what == "sf" .or. what == "corrSf" .or. what == "error"&
       & .or. what == "res") then
     call bcsSf(what, level) ! only implemented for id=4
  else if (what == "vel") then
     call bcsVel(what, level, id)
  else
     write(*,'(3a,i3)') "Boundary condition not defined for ", what, &
          ", id = ", id
     call endRun()
  end if

  return
end subroutine applyBcs
! =======================================================
