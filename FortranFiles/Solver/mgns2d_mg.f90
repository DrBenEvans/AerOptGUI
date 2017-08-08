!**********************************************
!*mgns2d_mg.f90                               *
!*                                            *
!*Solver file at multigrid abstraction level  *
!*                                            *
!*           mgns2d v 3.0                     *
!*                                            *
!*                                            *
!*Made by                                     *
!*Kaare A Sorensen,                           *
!*17.09.98-                                   *
!**********************************************
!*Description:                                *
!* This file contains the main driver sub-    *
!* routines for the multigrid solver, i.e.    *
!* it creates and sets up the grids and does  *
!* the multigrid solution procedure.          *
!**********************************************

module MultigridSolver
 use GridSolver
 use BoundarySolver
 use RestrictionOperator
 use ProlongationOperator
 use InputVariables

 type(GridSolverData),pointer :: grp(:) ! contains the single grid solvers
 type(InputVariablesData) :: ivd        ! contains the input variables 

 integer :: numberOfGrids
 real :: averageTimePerCycle
 integer :: time0(8),time1(8),time2(8),time3(8),time4(8)
 character*150 :: computationFileName,inputFileName,resultFileName,surfFileName
 character*150 :: startupFileName,residualFileName,wallInfoFileName 
 character*30 :: timeBuffer
 integer :: COMPINFILE,RESOUTFILE,INPINFILE,nameLength,allocateStatus,surfNameLength,RESOUTFILEt
 integer :: STARTFILE,RESIDUALOUTFILE,SURFOUTFILE,residualNameLength,resultNameLength
 integer :: WALLINFILE,SWITCHREG,SURFOUTFILE1
 real :: initialResidual,currentResidual,stopAtResidualReduction
 !added
 character*150 :: resultFileName2
 integer :: RESOUTFILE2
contains

!----------------------------------------------------------------------
 subroutine execute()
 ! main procedure
 IMPLICIT NONE
  call communicate()
  call multiGridCycles()
 end subroutine execute
!----------------------------------------------------------------------
 subroutine communicate()
 IMPLICIT NONE

 integer :: i
 integer :: NoNests
 character(len=8) :: FileStatus
 character(len=:), allocatable :: istr 
 
 write(*,*) ""
 write(*,'(A)') "*************************************************"
 write(*,'(A)') "   WELCOME TO MGNS2D - 2D MULTIGRID N-S SOLVER   "
 write(*,'(A)') "*************************************************"
 write(*,*) ""
 
 write(*,'(A)',advance="no") "Enter control filename: "
 read(*,'(A)') inputFileName
 nameLength = nameLen(inputFileName)
 if(nameLength>0) then 
  INPINFILE = 21
  open(INPINFILE,&
      file=inputFileName(1:nameLength),form='formatted',status='old')
  write(*,*) 'File '//inputFileName(1:nameLength) // ' opened ...'
 end if
 call readInputVariables(ivd,INPINFILE)
 if(nameLength>0) then 
  close(INPINFILE)
 end if

 write(*,'(A)',advance="no") "Enter computation filename: "
 read(*,'(A)') computationFileName
 nameLength = nameLen(computationFileName)
 COMPINFILE = 20
 open(COMPINFILE,&
    file=computationFileName(1:nameLength),form='unformatted',status='old')
 write(*,*) 'File '//computationFileName(1:nameLength) // ' opened ...'
 call readComputationData(COMPINFILE)
 close(COMPINFILE)

! write(*,'(A)',advance="no") "Enter wall information filename: "
! read(*,'(A)') wallInfoFileName
 wallInfoFileName = ""
 nameLength = nameLen(wallInfoFileName)
if(nameLength>0) then 
WALLINFILE = 27
 SWITCHREG = 28
 open(WALLINFILE,&
    file=wallInfoFileName(1:nameLength),form='formatted',status='old')
 write(*,*) 'File '//wallInfoFileName(1:nameLength) // ' opened ...'
 open(SWITCHREG,&
    file='switch.reg',form='formatted',status='old')
 write(*,*) 'File '//' switch.reg' // ' opened ...'
 call readWallInfoData(WALLINFILE,SWITCHREG)
 close(WALLINFILE)
 close(SWITCHREG)
endif

 STARTFILE = 0
 write(*,'(A)',advance="no") "Enter startup filename: "
 read(*,'(A)') startupFileName
 nameLength = nameLen(startupFileName)
 if(nameLength>0) then 
  STARTFILE = 22
  open(STARTFILE,&
      file=startupFileName(1:nameLength),form='formatted',status='old')
  write(*,*) 'File '//startupFileName(1:nameLength) // ' opened ...'
 end if
 call setUpInitialField(ivd,grp(1),STARTFILE)
 if(nameLength>0) then 
  close(STARTFILE)
 end if

 do i=2,numberOfGrids
  call setUpInitialField(ivd,grp(i),0)
 end do 

 RESOUTFILE=26
 RESOUTFILEt=226
 write(*,'(A)',advance="no") "Enter result filename: "
 read(*,'(A)') resultFileName

 RESIDUALOUTFILE=0
 write(*,'(A)',advance="no") "Enter residual filename: "
 read(*,'(A)') residualFileName
 residualNameLength = nameLen(residualFileName)
 if(residualNameLength>0) then 
  RESIDUALOUTFILE = 25
 end if

 SURFOUTFILE = 0
! write(*,'(A)',advance="no") "Enter surface filename: "
! read(*,'(A)') surfFileName
 surfFileName = ""
 surfNameLength = nameLen(surfFileName)
 if(surfNameLength>0) then
  SURFOUTFILE = 59
  SURFOUTFILE1 = 69
 end if

 end subroutine communicate 
!----------------------------------------------------------------------
 subroutine multiGridCycles()
! handles entire multigrid process  
 IMPLICIT NONE
 
! type (GridSolverData) :: grp
! type (InputVariablesData) :: ivd
 integer :: i,j 
 
 averageTimePerCycle = 0.0
 i = 0

 ! check what kind of iteration scheme to use

 initialResidual = 10000.0 
 currentResidual = 10000.0 
 stopAtResidualReduction = -1.0 
 if(ivd%numberOfMGIterations<0) then 
  stopAtResidualReduction = abs(dble(ivd%numberOfMGIterations))
 end if

 write(*,*) "Starting timestepping..."

 if(ivd%HighOrder) then
  do j=1,numberOfGrids
   call calculateHOCoefficients(grp(j),ivd)
  end do
 end if

 call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time0)
 if(stopAtResidualReduction<0) then 
  do i=1,ivd%numberOfMGIterations 
   do j=1,numberOfGrids   
    grp(j)%iterationNumber = i
   end do
   call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time1)
   call oneCycle(i)
   call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time2)
 
   ! write residuals
   call writeResiduals(grp(1),i)
!   call writeResiduals(grp(numberOfGrids),i)

   ! write to file if needed 

   if(mod(i,ivd%writeToFileInterval) == 0) then 
    call writeResultsToFile()
   end if 
  
   ! calculate average time
   time1 = time2 - time1
   averageTimePerCycle=averageTimePerCycle+&
       3600.0*time1(5)+60.0*time1(6)+1.0*time1(7)+0.001*time1(8)
  end do
 else
  do while((log(initialResidual)-log(currentResidual)&
                <log(10.0)*stopAtResidualReduction).AND.(i<ivd%maxitt))
   print * ,"Residual", (log(initialResidual)-log(currentResidual))/log(10.0)
   i = i+1
   do j=1,numberOfGrids   
    grp(j)%iterationNumber = i
   end do
   call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time1)
   print * , "Iteration", i
   call oneCycle(i)
   call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time2)
 
   ! write residuals
!   print * , "Write residuals"
   call writeResiduals(grp(1),i)
!!   call writeResiduals(grp(numberOfGrids),i)
   
   ! write to file if needed 

   if(mod(i,ivd%writeToFileInterval) == 0) then 
    call writeResultsToFile()
   end if 
  
   ! calculate average time
   time1 = time2 - time1
   averageTimePerCycle=averageTimePerCycle+&
       3600.0*time1(5)+60.0*time1(6)+1.0*time1(7)+0.001*time1(8)
  end do
  print * ,"Final residual", (log(initialResidual)-log(currentResidual))/log(10.0)
 end if

 ! always write results at end of calculations
 call writeResultsToFile()

 write(*,*) "Average time per cycle: ",averageTimePerCycle/dble(i)," seconds"
 end subroutine multiGridCycles
!-----------------------------------------------------------------------
 subroutine oneCycle(cycleNo)
! one multigrid cycle
 IMPLICIT NONE

 integer :: cycleNo

 integer :: i,j,numberOfSweeps,finestSweepLevel,localNumberOfGrids
 logical :: calculateTurbulence 

 if(cycleNo.le.ivd%numberOfTurbulenceSteps.and.ivd%ReynoldsNumber>0.0) then
  calculateTurbulence = .true.
 else
  calculateTurbulence = .false.
 end if

! call calculateHOCoefficients(grp(1),ivd)

 if(ivd%multigridScheme==1) then 
! V-cycle

 ! restriction part

  do i=1,numberOfGrids-1
  ! call grid solver
   call doTimeIterations(grp(i),ivd,ivd%numberOfRelaxationSteps,i==1,calculateTurbulence,cycleNo)
   call restrictGrid(grp(i),grp(i+1),ivd,i,calculateTurbulence)
!   call mapFineToCoarse(grp(i)%u,grp(i+1),grp(i+1)%u)
!   call mapTurbFineToCoarse(grp(i)%turbulentViscosity,grp(i+1),grp(i+1)%turbulentViscosity)
  end do
  
  call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo)
!  call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,.true.)

 ! prolongation part

  do i=numberOfGrids,2,-1
   call prolongateGrid(grp(i),grp(i-1),ivd)
!   call mapCoarseToFine(grp(i)%u,grp(i),grp(i-1),grp(i-1)%u)
!   call mapTurbCoarseToFine(grp(i)%turbulentViscosity,grp(i),grp(i-1),grp(i-1)%turbulentViscosity)
  end do

 else if(ivd%multigridScheme==2) then 
! W-cycle

  if(numberOfGrids<=2) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = 2*numberOfGrids - 4
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo)
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo)

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else if(i<=numberOfSweeps/2) then
    finestSweepLevel = finestSweepLevel - 1 
   else
    finestSweepLevel = finestSweepLevel + 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do

  end do
   ! prolongate back to finest grid again 

   do j=numberOfGrids,2,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do

 else if(ivd%multigridScheme==3) then  
 ! growing V-cycle

  if(numberOfGrids==1) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = numberOfGrids - 1 
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo)
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo)

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>=2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else
    finestSweepLevel = finestSweepLevel - 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do

  end do
 
 else if(ivd%multigridScheme==4) then 
 ! decaying V-cycle

  if(numberOfGrids==1) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = numberOfGrids - 1 
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo)
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo)

   ! set start and end points for sweep

   finestSweepLevel = finestSweepLevel + 1 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do
  end do

  ! prolongate back to finest grid again 

   do j=numberOfGrids,2,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do

 else if(ivd%multigridScheme==5) then 
! weighted W-cycle

  if(numberOfGrids<=2) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = 2*numberOfGrids - 4
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,j,j==finestSweepLevel,calculateTurbulence,cycleNo)
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,numberOfGrids,numberOfGrids==1,calculateTurbulence,cycleNo)

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else if(i<=numberOfSweeps/2) then
    finestSweepLevel = finestSweepLevel - 1 
   else
    finestSweepLevel = finestSweepLevel + 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do

  end do
   ! prolongate back to finest grid again 

   do j=numberOfGrids,2,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do
 else if(ivd%multigridScheme==7) then  
 ! growing V-cycle, every other grid
  localNumberOfGrids = numberOfGrids - 0.5*(numberOfGrids-1)
  if(numberOfGrids==1) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = numberOfGrids - 1 
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo)
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo)

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>=2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else
    finestSweepLevel = finestSweepLevel - 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
   end do

  end do
 
 else if(ivd%multigridScheme==4) then 
 else
  STOP "ERROR: Specified multigrid scheme not implemented"
 end if
 end subroutine oneCycle


!-----------------------------------------------------------------------
 subroutine  readWallInfoData(WALLINFILE,SWITCHREG)
! reads wall data from file    IMPLICIT NONE

  integer :: WALLINFILE,SWITCHREG
  integer :: i,j,k,allocateStatus,ibuff
  integer :: nonit,npswap
  integer :: ponit(100) 
  integer,pointer :: pswap(:,:)

 ! get number of grids to use

  read(WALLINFILE,*) numberOfWallPoints
  read(SWITCHREG,*) npswap

  write(*,*) "Number of wall points: ",numberOfWallPoints
  write(*,*) "Number of points swapped: ",npswap

  allocate(pswap(2,npswap),stat=allocateStatus)
  if (allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readWallInfoData"
  do i =1 , npswap
   read(SWITCHREG,*) pswap(1,i),pswap(2,i)
  end do
  do i = 1 , numberOfWallPoints
   read(WALLINFILE,*) nonit,(ponit(j),j=1,nonit)
   do j = 1 , nonit
    do k = 1 , npswap
     if(ponit(j).eq.pswap(1,k)) then
       ponit(j) = pswap(2,k)
     end if
     if(ponit(j).eq.pswap(2,k)) then
       ponit(j) = pswap(1,k)
     end if
    end do
   end do
   do j = 1 , nonit
    if(j.lt.nonit) grp(1)%AfterAndBefore(ponit(j),1) = ponit(j+1)
    if(j.gt.1) grp(1)%AfterAndBefore(ponit(j),2) = ponit(j-1)
   end do
 end do
 deallocate(pswap,stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: readWallInfoData couldn't deallocate"

! do i = 1 , grp(1)%numberOfNodes
!  write(99,*) i,grp(1)%AfterANdBefore(i,1),grp(1)%AfterANdBefore(i,2)
! end do
 end subroutine readWallInfoData


!-----------------------------------------------------------------------
 subroutine readComputationData(INFILE)
! reads data from file 
  IMPLICIT NONE

  integer :: INFILE

  integer :: i,j,k,allocateStatus,ibuff 
  real :: rbuff
  type(LinkedInteger),pointer :: currentInteger
  type(LinkedReal),pointer :: currentReal

 ! get number of grids to use
 
  read(INFILE) numberOfGrids

  if(ivd%numberOfGridsToUse>0) then 
   numberOfGrids = min(numberOfGrids,ivd%numberOfGridsToUse)
  end if
  write(*,*) "Number of grids: ",numberOfGrids

  allocate(grp(numberOfGrids),stat=allocateStatus)
  if (allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readComputationData"  
  
  allocate(grp(1)%brp,stat=allocateStatus)
  if (allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readComputationData"  
  grp(1)%gridNumber = 1
  call readBoundarySolverData(grp(1)%brp,INFILE)
  call readGridComputationData(.true.,grp(1),grp(1),grp(1),ivd,INFILE)
  nullify(grp(1)%inverseSideLengthArray) 
  do i=2,numberOfGrids
   allocate(grp(i)%brp,stat=allocateStatus)
   if (allocateStatus/=0)&
      STOP "ERROR: Not enough memory in readComputationData"  
   grp(i)%gridNumber = i 
   call readBoundarySolverData(grp(i)%brp,INFILE)
   call readGridComputationData(.false.,grp(1),grp(i-1),grp(i),ivd,INFILE)
  end do

 end subroutine readComputationData
!----------------------------------------------------------------------
 subroutine writeResiduals(grp,iterationNumber)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 integer :: iterationNumber

 integer :: i,j
 real :: res(4)
 real :: maxres(4),cres
 real :: sec,lift,drag,momentum,pressureLift,pressureDrag,frictionDrag

 do i=1,4
  res(i) = sqrt(sum((grp%u(:,i)-grp%uprev(:,i))&
             *(grp%u(:,i)-grp%uprev(:,i))))
 end do

 maxres = 0.
 do i=1,grp%numberOfNodes
  do j=1,4
   cres = abs(grp%u(i,j)-grp%uprev(i,j))
   if(cres>maxres(j)) maxres(j) = cres
  end do
 end do



! set residual monitor

! currentResidual = log(res(1)*grp%residualScalingFactor)/log(10.0)
 currentResidual = res(1)
 if(initialResidual==10000.0) then 
  initialResidual = res(1) 
 end if

! print * , "calculate lift and drag"
 call calculateLiftAndDrag(grp,ivd,lift,drag,momentum,pressureLift,pressureDrag,frictionDrag)

! write to screen
! write(*,501) iterationNumber,res(1)*grp%residualScalingFactor," (",maxres(1),") ",&
!                              res(2)*grp%residualScalingFactor," (",maxres(2),") ",&
!                              res(3)*grp%residualScalingFactor," (",maxres(3),") ",&
!                              res(4)*grp%residualScalingFactor," (",maxres(4),") ",&
!                              lift,drag 
! write(*,*) iterationNumber,res*grp%residualScalingFactor,lift,drag 


! write to file
 if(RESIDUALOUTFILE>0) then 
  if(iterationNumber==1) then 
   open(RESIDUALOUTFILE,&
      file=residualFileName(1:residualNameLength),form='formatted',status='unknown')
  else
   open(RESIDUALOUTFILE,&
      file=residualFileName(1:residualNameLength),form='formatted',status='unknown',&
                                                           position='append')
  end if
  if(ivd%useTimeResidual) then 
   sec = 3600.0*(time2(5)-time0(5))+60.0*(time2(6)-time0(6))+&
         1.0*(time2(7)-time0(7))+0.001*(time2(8)-time0(8))
!   write(RESIDUALOUTFILE,*) sec,log(res*grp%residualScalingFactor)/log(10.0)
   write(RESIDUALOUTFILE,'(I7,7E14.5)') sec,log(res(1)/initialResidual)/log(10.0),lift,drag,momentum,pressureLift,pressureDrag,frictionDrag
  else
!   write(RESIDUALOUTFILE,*) iterationNumber,log(res*grp%residualScalingFactor)/log(10.0)
   write(RESIDUALOUTFILE,'(I7,7E14.5)') iterationNumber,log(res(1)/initialResidual)/log(10.0),lift,drag,momentum,pressureLift,pressureDrag,frictionDrag
  end if
  close(RESIDUALOUTFILE)
 end if

!  open(15,&
!      file="increment.res",form='formatted',status='unknown')
!  do i=1,grp%numberOfNodes
!   write(15,*) i,grp%u(i,:)-grp%uprev(i,:),0,0 
!  end do
!  close(15)

501 format(I7,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E12.4,E12.4)

 end subroutine writeResiduals
!-----------------------------------------------------------------------
 subroutine writeResultsToFile()
 ! writes the solution fields of the finest grid to file
 IMPLICIT NONE

 integer :: i 
 real :: turbulenceCoefficient

 nameLength = nameLen(resultFileName)
 open(RESOUTFILE,&
     file=resultFileName(1:nameLength),form='formatted',status='replace')
 write(*,*) 'File '//resultFileName(1:nameLength) // ' opened ...'
! open(RESOUTFILEt,&
!     file='V_'//resultFileName(1:nameLength),form='formatted',status='replace')
! write(*,*) 'File V_'//resultFileName(1:nameLength) // ' opened ...'

! call writeGridResults(grp(1),RESOUTFILE,RESOUTFILEt,ivd)
 call writeGridResults(grp(1),RESOUTFILE,ivd)

 close(RESOUTFILE)
! close(RESOUTFILEt)

 if(SURFOUTFILE>0) then 
  nameLength = nameLen(surfFileName)
!  open(SURFOUTFILE,&
!      file=surfFileName(1:nameLength),form='formatted',status='replace')
!  write(*,*) 'File '//surfFileName(1:nameLength) // ' opened ...'
!  open(SURFOUTFILE1,&
!      file=surfFileName(1:nameLength)//'1',form='formatted',status='replace')
!  write(*,*) 'File '//surfFileName(1:nameLength) // ' opened ...'

!  call writeSurface(grp(1),SURFOUTFILE1,ivd)
!  call writeSurface2(grp(1),SURFOUTFILE,ivd)

!  close(SURFOUTFILE)
!  close(SURFOUTFILE1)
 end if

! open(15,&
!     file="increment.res",form='formatted',status='unknown')
! do i=1,grp(1)%numberOfNodes
!  write(15,'(I5,6E18.6)') i,grp(1)%u(i,:)-grp(1)%uprev(i,:),0.0,0.0 
! end do
! close(15)

! open(15,&
!     file="turbIncr.res",form='formatted',status='unknown')
! do i=1,grp(1)%numberOfNodes
!  write(15,'(I5,6E18.6)') i,1.0,grp(1)%u(i,5)-grp(1)%uprev(i,5),0.0,0.0,0.0,0.0 
! end do
! close(15)

 end subroutine writeResultsToFile
!-----------------------------------------------------------------------
! subroutine writeResultsToFile2()
 ! writes the solution fields at engine inlet of the finest grid to file
! IMPLICIT NONE

! integer :: i 
! real :: turbulenceCoefficient

! nameLength = nameLen(resultFileName2)
! open(RESOUTFILE2,&
!     file=resultFileName2(1:nameLength),form='formatted',status='replace')
! write(*,*) 'File '//resultFileName2(1:nameLength) // ' opened ...'

! call writeGridResults2(grp(1),RESOUTFILE2)

! close(RESOUTFILE2) 

! end subroutine writeResultsToFile2
 !----------------------------------------------------------------------
 integer function nameLen(fn)
 IMPLICIT NONE
 
 character*150 :: fn
 integer :: i

 do i=150,1,-1
  nameLen = i
  if(fn(i:i)/=' ') GOTO 77 ! EXIT
 end do
 nameLen = 0 
77  end function nameLen
!-----------------------------------------------------------------------
end module MultigridSolver



