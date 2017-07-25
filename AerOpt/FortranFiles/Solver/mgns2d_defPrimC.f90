!**************************************
!*mgns2d_def.f90                      *
!*                                    *
!*Module definition file for program  *
!*                                    *
!*         mgns2d v 3.0               *
!*                                    *
!*                                    *
!*Made by                             * 
!*Kaare A Sorensen,                   *
!*17.09.98-                           *
!**************************************
!*Description:                        *
!* This file includes the single grid *
!* solver together with prolongation  *
!* and restriction operators needed   *
!* for the multigrid algorithm.       *
!**************************************

module Toolbox
! various handy datastructures

 type LinkedInteger
  integer :: int  
  type(linkedInteger),pointer :: next
 end type linkedInteger

 type LinkedIntegerP
  type(LinkedInteger),pointer :: first
 end type LinkedIntegerP

 type LinkedIntegerPArray
  type(linkedIntegerPArray),pointer :: lipa(:)
 end type linkedIntegerPArray

 type LinkedReal
  real :: re 
  type(LinkedReal),pointer :: next 
 end type LinkedReal

 type LinkedRealP
  type(LinkedReal),pointer :: first 
 end type LinkedRealP

 type LinkedRealPArray
  type(LinkedRealP),pointer :: lrpa(:)
 end type LinkedRealPArray

 type RealArrayPointers
  real,pointer :: ar(:,:)
 end type realArrayPointers 

 type IntegerArrayPointers
  integer,pointer :: ar(:,:)
 end type IntegerArrayPointers

end module Toolbox

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module InputVariables
! Holds the input variables specified by the user, such as 
! inflow variables etc. 

 type InputVariablesData
  integer :: numberOfMGIterations
  integer :: maxitt ! added by Sean
  integer :: writeToFileInterval
  real :: ReynoldsNumber,MachNumber,PrandtlNumber,CFLNumber
  real :: turbulentPrandtlNumber
  real :: inflowTemperature
  real :: gamma ! you know - 1.4
!  real :: R ! real gas constant
  real :: alpha ! angle of attack
  real :: wallTemperature ! for isothermal wall
  real :: HartensCorrectionFactor ! used in roe matrix construction for
                                  ! dissipation terms
  logical :: useMatrixDissipation ! decides if matrix or scalar smoothing
                                  ! is to be applied
  real :: secondOrderDissipationFactor,fourthOrderDissipationFactor
  real :: coarseGridDissipationFactor ! for multigrid coarse grids
  logical :: useDissipationWeighting
  integer :: viscosityScheme
  integer :: coarseGridDissipationScheme
  integer :: numberOfRelaxationSteps ! number of iteration for each grid
  integer :: numberOfRSSteps ! number of residual smoothing steps
  real :: inflowField(5,2) ! inflow variables, 2d: 1-5 where 5 is pressure
  real :: ThomsonRelaxationFactor1,ThomsonRelaxationFactor2
  real :: normalComponentRelaxation
  real :: residualSmoothingFactor
  real :: prolongationSmoothingFactor
  integer :: numberOfPSSteps ! number of prolongation smoothing steps
  real :: vonKarmanConstant ! used in trubulence models
  real :: BLA ! used in Baldwin-Lomax turbulence model
  real :: BLAlpha
  real :: BLBeta
  real :: BLCw
  real :: turbulenceCFLFactor
  real :: turbulenceSmoothingFactor
  real :: secondOrderTurbVisc
  real :: maxIncrementFactor ! maximum percentage of unknown to be added
  real :: incrRelax1 ! increment relaxation parameters
  real :: incrRelax2 
  real :: wakeSectionDensity ! for use in B-L turbulence model
  real :: tripFactor
  real :: prolongationRelaxation
  real :: momentumPoint(2)
  integer :: numberOfWakeSections ! for use in B-L turbulence model
  integer :: tripNodes(10)
  integer :: numberOfTurbulenceSteps
  integer :: sizeOfSeparationField
  integer :: numberOfTriperations
  integer :: multigridScheme
  integer :: turbulenceModel
  integer :: numberOfTurbulenceSubsteps
  integer :: residualSmoothingIterations
  integer :: dissipationScheme
  integer :: numberOfGridsToUse
  logical :: useTimeResidual ! if true, plots residual vs. time instead of iterations
  integer :: boundaryTerm 

  integer :: prolongationScheme

  logical :: wallsAreAdiabatic
  logical :: HighOrder
  
  logical :: triggerTurbulence
  real :: triggerRadius
  real :: triggerValue
  integer :: numberOfTriggerSteps
  integer :: numberOfDissipationLayers
  real :: sixthOrderDissipationFactor
  
  ! for engine inlet boundary type
  real :: enginesFrontMassFlow
  real :: engineBCRelaxation
 
 end type InputVariablesData

contains

!-----------------------------------------------------------------------
 subroutine readInputVariables(ivd,INFILE)
 ! as the name says
 IMPLICIT NONE

 type(InputVariablesData) :: ivd
 integer :: INFILE

 integer :: allocateStatus
 real :: PI,alphaRad
 
 namelist /InputVariables/ ivd

  ! set default values 
 
  ivd%numberOfMGIterations = 100 
  ivd%numberOfRelaxationSteps = 5
  ivd%writeToFileInterval = 25
  ivd%gamma = 1.4
  ivd%ReynoldsNumber = 1.0E06 
  ivd%MachNumber = 0.5
  ivd%PrandtlNumber = 0.72 
  ivd%turbulentPrandtlNumber = 0.9
  ivd%CFLNumber = 1.0 
  ivd%inflowTemperature = 528.0
  ivd%wallTemperature = 528.0
  ivd%HartensCorrectionFactor = 0.0
  ivd%coarseGridDissipationScheme = 1
  ivd%useMatrixDissipation = .false.
  ivd%dissipationScheme = 1
  ivd%multigridScheme = 1
  ivd%turbulenceModel = 0
  ivd%secondOrderDissipationFactor = 0.4
  ivd%fourthOrderDissipationFactor = 0.2
  ivd%inflowField(1,1) = 1.0 ! density 
  ivd%inflowField(2,1) = 1.0 ! x-momentum
  ivd%inflowField(3,1) = 0.0 ! y-momentum
  ivd%inflowField(4,1) = 1.0 ! total energy
  ivd%inflowField(5,1) = ivd%inflowField(1,1)*(ivd%MachNumber**2)/ivd%gamma ! pressure 
  ivd%inflowField(:,2) = 0.0
  ivd%normalComponentRelaxation = 1.0
  ivd%residualSmoothingFactor = 0.0 
  ivd%prolongationSmoothingFactor = 0.0
  ivd%ThomsonRelaxationFactor1 = -0.2
  ivd%ThomsonRelaxationFactor2 = 2.0
  ivd%useTimeResidual = .false.
  ivd%numberOfGridsToUse = 0 ! this means all
  ivd%numberOfRSSteps = 1
  ivd%numberOfPSSteps = 1
  ivd%vonKarmanConstant = 0.41
  ivd%BLA = 26.0 
  ivd%BLAlpha = 0.3 
  ivd%BLBeta = 1.6 
  ivd%BLCw = 1.0
  ivd%tripFactor = 1.0
  ivd%numberOfTurbulenceSteps = 1000000
  ivd%numberOfTurbulenceSubsteps = 1
  ivd%turbulenceCFLFactor = 1.0
  ivd%numberOfTriperations = 0
  ivd%secondOrderTurbVisc = 0.0
  ivd%wakeSectionDensity = 0.25
  ivd%numberOfWakeSections = 20 
  ivd%turbulenceSmoothingFactor = 0.002
  ivd%maxIncrementFactor = 0.25
  ivd%incrRelax1 = -0.2
  ivd%incrRelax2 = 2.0
  ivd%useDissipationWeighting = .false.
  ivd%tripNodes = 0
  ivd%sizeOfSeparationField = 10
  ivd%prolongationRelaxation = 1.0
  ivd%momentumPoint(1) = 0.25
  ivd%momentumPoint(2) = 0.0
  ivd%viscosityScheme = 1
  ivd%boundaryTerm=1

  ivd%wallsAreAdiabatic = .true.
  ivd%HighOrder = .false.

  ivd%triggerTurbulence = .false.
  ivd%triggerRadius = 0.0
  ivd%triggerValue = 0.0
  ivd%numberOfTriggerSteps = 0 

  ivd%prolongationScheme = 1
  ivd%numberOfDissipationLayers=0
  ivd%sixthOrderDissipationFactor = 0.05
  
  ! for engine inlet boundary type
  ivd%enginesFrontMassFlow = 0.1
  ivd%engineBCRelaxation = 1.0

  ! added by Sean : max nb of iterations
  ivd%maxitt = 1000000

  ! read input variables file
  
  if(INFILE>0) then  
   read(INFILE,InputVariables) 
  end if 
  
  write(*,*) "Engine front mass flow : ", ivd%enginesFrontMassFlow

  PI = 4.0*atan(1.0)

  ivd%wallTemperature = (1.0/((ivd%gamma-1.0)*(ivd%MachNumber**2)))*(ivd%wallTemperature/ivd%inflowTemperature) 

  alphaRad = ivd%alpha*PI/180. 
  ivd%inflowField(1,1) = 1.0
  ivd%inflowField(2,1) = cos(alphaRad) 
  if (abs(ivd%inflowField(2,1)) < 1.d-15) then
    ivd%inflowField(2,1) = 0.D0
  endif
  ivd%inflowField(3,1) = sin(alphaRad) 
  if (abs(ivd%inflowField(3,1)) < 1.d-15) then
    ivd%inflowField(3,1) = 0.D0
  endif 

  write(*,*) "S: ",ivd%alpha,alphaRad,ivd%inflowField(2,1),ivd%inflowField(3,1)

  ivd%inflowField(4,1) = (1./((ivd%gamma-1.)*ivd%gamma*ivd%MachNumber**2))+0.5
  ivd%inflowfield(5,1) = 1./(ivd%gamma*ivd%MachNumber**2)

 end subroutine readInputVariables
!-----------------------------------------------------------------------
end module InputVariables

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module ProlongationOperator
use Toolbox

 type ProlongationOperatorData
  ! the prolongationArray is a chain array over which 
  ! nodes in the fine grid is merged into each coarse grid
  ! node
  type(LinkedIntegerP),pointer :: prolongationArray(:)
  integer :: numberOfNodes
 end type ProlongationOperatorData

contains

!-----------------------------------------------------------------------
 subroutine readProlongationOperatorData(pod,INFILE)
 IMPLICIT NONE
 
 type(ProlongationOperatorData) :: pod
 integer :: INFILE

 type(LinkedInteger),pointer :: currentInteger
 integer :: allocateStatus,ibuff,j

  allocate(pod%prolongationArray(pod%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readProlongationOperatorData "

  do j=1,pod%numberOfNodes
   nullify(pod%prolongationArray(j)%first)
   read(INFILE) ibuff
   do while(ibuff > 0) ! a zero means end of chain 
    allocate(currentInteger,stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    currentInteger%int=ibuff
    currentInteger%next=>pod%prolongationArray(j)%first
    pod%prolongationArray(j)%first=>currentInteger  
    read(INFILE) ibuff 
   end do
  end do
 end subroutine readProlongationOperatorData
!-----------------------------------------------------------------------
end module ProlongationOperator

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module RestrictionOperator
use Toolbox

 type RestrictionOperatorData
 ! the restrictionArray is a chain array of reals holding
 ! the coefficients to be multiplied to unknown field 
 ! in fine grid in mapping to coarse grid. Which nodes 
 ! in fine grid to be multiplied by these coefficients are
 ! given in the prolongationOperator module 
  type(LinkedRealP),pointer :: restrictionArray(:)
  integer :: numberOfNodes
 end type RestrictionOperatorData

contains

!-----------------------------------------------------------------------
 subroutine readRestrictionOperatorData(rod,INFILE)
 IMPLICIT NONE
 type(RestrictionOperatorData) :: rod 
 integer :: INFILE

 type(LinkedReal),pointer :: currentReal
  integer :: allocateStatus,j
  real :: rbuff

  ! restriction operator (fine to coarse)
  allocate(rod%restrictionArray(rod%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
          STOP "ERROR: Not enough memory in readRestrictionOperatorData"

  do j=1,rod%numberOfNodes
   nullify(rod%restrictionArray(j)%first)
   read(INFILE) rbuff
   do while(rbuff > 0) ! a zero means end of chain 
    allocate(currentReal,stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    currentReal%re=rbuff
    currentReal%next=>rod%restrictionArray(j)%first
    rod%restrictionArray(j)%first=>currentReal  
    read(INFILE) rbuff 
   end do
  end do

 end subroutine readRestrictionOperatorData
!-----------------------------------------------------------------------
end module RestrictionOperator

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module BoundarySolver
! does BC's etc

 type BoundarySolverData
  integer,pointer :: faceIndexArray(:,:)
  integer,pointer :: faceIndicatorArray(:)
  real,pointer :: faceWeightsArray(:,:)
  real,pointer :: faceHelpArray(:,:)
  real,pointer :: faceTangentArray(:,:)
  real,pointer :: faceWeightNorms(:)

  real,pointer :: sideLengthArray(:)

  real,pointer :: IOBoundaryCoordinates(:,:) ! holds coordinates for internal outflow boundaries
  real,pointer :: IOChainCoordinates(:,:) ! holds coordinates of nodes in interpolation chain 
  integer,pointer :: IOBoundaryRegister(:)
  integer,pointer :: IOIntervals(:,:) 
  integer :: numberOfIOBoundaries,numberOfIONodes

  integer :: numberOfBoundaryFaces,numberOfCoarseBoundaryFaces
  
  ! added for engine inlet boundary type (from 3D)
  integer :: numberOfEngineInletSides
  integer,pointer :: engineInletSideIndexes(:,:)
  real,pointer :: engineInletSideCoefficients(:,:)
  real :: engineInletNormals(2)
  real :: engineInletAreas
  integer,pointer :: nodeIndicatorRegister(:)
  integer,pointer :: nodeIndicatorArray(:)
  real,pointer :: nodeNormalArray(:,:)
  
 end type BoundarySolverData

contains

!-----------------------------------------------------------------------
 subroutine readBoundarySolverData(bsd,INFILE)
 IMPLICIT NONE
 type(BoundarySolverData) :: bsd
 integer :: INFILE

 integer :: allocateStatus,i,j,k,nodeNumber,numberOfIOBoundaryNodes


  read(INFILE) bsd%numberOfBoundaryFaces 
  write(*,*) ' Number of boundary sides=', bsd%numberOfBoundaryFaces
  allocate(bsd%faceIndexArray(bsd%numberOfBoundaryFaces,3),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  do i=1,bsd%numberOfBoundaryFaces
!   read(INFILE) bsd%faceIndexArray(i,1),bsd%faceIndexArray(i,2),&
!              bsd%faceIndexArray(i,3)
   read(INFILE) bsd%faceIndexArray(i,1:3)
  end do
  
  allocate(bsd%faceIndicatorArray(-20:bsd%numberOfBoundaryFaces),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  do i=-20,bsd%numberOfBoundaryFaces
   read(INFILE) bsd%faceIndicatorArray(i)
  end do

  allocate(bsd%faceWeightsArray(bsd%numberOfBoundaryFaces,2),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  allocate(bsd%faceWeightNorms(bsd%numberOfBoundaryFaces),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  do i=1,bsd%numberOfBoundaryFaces
   read(INFILE) bsd%faceWeightsArray(i,1),bsd%faceWeightsArray(i,2)
   bsd%faceWeightNorms(i) = sqrt(sum(bsd%faceWeightsArray(i,:)*bsd%faceWeightsArray(i,:)))
  end do

  allocate(bsd%faceTangentArray(bsd%numberOfBoundaryFaces,2),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  do i=1,bsd%numberOfBoundaryFaces
   read(INFILE) bsd%faceTangentArray(i,1),bsd%faceTangentArray(i,2)
  end do

  ! read internal outflow (IO) data
  print * , "reading IO info"
  numberOfIOBoundaryNodes = bsd%faceIndicatorArray(-7) - bsd%faceIndicatorArray(-8) + 1
  read(INFILE) bsd%numberOfIOBoundaries
  read(INFILE) bsd%numberOfIONodes
  if(bsd%numberOfIONodes>0) then 
   allocate(bsd%IOBoundaryRegister(bsd%numberOfIONodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
   allocate(bsd%IOChainCoordinates(bsd%numberOfIONodes,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  else
   nullify(bsd%IOBoundaryRegister)
   nullify(bsd%IOChainCoordinates)
  end if 
  if(numberOfIOBoundaryNodes>0) then 
   allocate(bsd%IOBoundaryCoordinates(bsd%numberOfBoundaryFaces,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  else
   nullify(bsd%IOBoundaryCoordinates)
  end if
  if(bsd%numberOfIOBoundaries>0) then 
   allocate(bsd%IOIntervals(bsd%numberOfIOBoundaries,4),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  else
   nullify(bsd%IOIntervals)
  end if
  j = 1  
  k = 1 
  write(*,*) "Number of IO nodes: ",bsd%numberOfIONodes,bsd%numberOfIOBoundaries
  do i=1,bsd%numberOfIONodes
   read(INFILE) nodeNumber
   bsd%IOBoundaryRegister(i) = nodeNumber
   write(*,*) "I: ",i,nodeNumber
   if(nodeNumber.le.bsd%numberOfBoundaryFaces) then 
    if(k==1) then 
     bsd%IOIntervals(j,1) = nodeNumber
     bsd%IOIntervals(j,3) = i
     k = 2
     write(*,*) "K1: ",j,nodeNumber,i
    else
     bsd%IOIntervals(j,2) = nodeNumber
     bsd%IOIntervals(j,4) = i
     j = j+1
     k = 1
     write(*,*) "K2: ",j,nodeNumber,i
    end if 
   end if
  end do 
  
  ! read engine inlet data
  read(INFILE) bsd%numberOfEngineInletSides
  write(*,*) 'Number of engine inlet sides:', bsd%numberOfEngineInletSides
  allocate(bsd%engineInletSideIndexes(bsd%numberOfEngineInletSides,2),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
  allocate(bsd%engineInletSideCoefficients(bsd%numberOfEngineInletSides,2),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"

  read(INFILE) bsd%engineInletSideIndexes(1:bsd%numberOfEngineInletSides,1:2)
  read(INFILE) bsd%engineInletSideCoefficients(1:bsd%numberOfEngineInletSides,1:2)



 end subroutine readBoundarySolverData
!-----------------------------------------------------------------------
 subroutine setIOCoordinates(bsd,coordinates)
 IMPLICIT NONE

 type(BoundarySolverData),pointer :: bsd
 real,pointer :: coordinates(:,:)

 integer :: i,j,k,it

  ! first set coordinates for chain coordinates
  do i=1,bsd%numberOfIONodes
   bsd%IOChainCoordinates(i,:) = coordinates(bsd%IOBoundaryRegister(i),:)
  end do 
  k = 0
!  do i=1,bsd%numberOfIOBoundaries
!   if(bsd%IOIntervals(i,1)<bsd%IOIntervals(i,2)) then 
!    it = 1
!   else
!    it = -1
!   end if
!   do j=bsd%IOIntervals(i,1)+it,bsd%IOIntervals(i,2)-it,it
!    k = k+1
!    bsd%IOBoundaryCoordinates(k,:) = coordinates(j,:) 
!   end do
!  end do 

 if(bsd%numberOfIONodes>0) then 
  bsd%IOBoundaryCoordinates(1:bsd%numberOfBoundaryfaces,:) =&
                  coordinates(1:bsd%numberOfBoundaryFaces,:)
 end if
 end subroutine setIOCoordinates
!-----------------------------------------------------------------------
end module BoundarySolver

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module GridSolver
use Toolbox
use RestrictionOperator
use ProlongationOperator
use BoundarySolver
use InputVariables 

! basically a one-grid solver with the added functionality of 
! intergrid mappings. 

type GridSolverData
 real,pointer :: sideWeightsArray(:,:)
 integer,pointer :: sideIndexArray(:,:)       
 real,pointer :: sideLengthArray(:,:) ! two first entries are dx and dy  
                                      ! components, the third is length
                                      ! of side
 real,pointer :: inverseSideLengthArray(:,:) ! used for agglomerated
                                             ! grids only
 real,pointer :: coordinates(:,:)
 real,pointer :: nodeHelpArray(:,:) ! node-based help array 
 real,pointer :: sideHelpArray(:,:) ! side-based help array
 real,pointer :: uprev(:,:) ! unknowns at previous time step
 real,pointer :: u(:,:)
 real,pointer :: dissipation(:,:)
 real,pointer :: laminarViscosity(:)
 real,pointer :: vorticity(:)
 real,pointer :: divergence(:)
 integer,pointer :: AfterAndBefore(:,:) ! nodes after and before the side nodes
 real,pointer :: hoc(:,:) ! high-order coefficients of the variable discretisations
 real,pointer :: hoc2(:,:) ! high-order coefficients for gradient interpolations
 real,pointer :: hoc3(:,:) ! high-order coefficients for Jacobian interpolations
 integer,pointer :: nodeConnectivityArray(:) ! number of sides per node 
 real,pointer :: rhs(:,:) ! right hand side in calculation
 real,pointer :: sourceTerm(:,:) ! used in FAS multigrid
 real,pointer :: p(:) ! pressure field
 real,pointer :: turbMGSourceTerm(:)
 real,pointer :: nodeVolume(:) ! the nodal volume (lumped mass matrix) 
 real,pointer :: localTimeSteps(:) ! time steps used
 real,pointer :: localInviscidTimeSteps(:) 
 real,pointer :: wallDistance(:) ! used in turbulence modeling
 real,pointer :: wallLength(:) ! used in calculating lift and drag
 real,pointer :: wallStress(:,:) ! used in lift and drag calculations
 real,pointer :: tripSourceArray(:)
 integer,pointer :: wallDistanceFaceArray(:) ! for turbulence modeling
 integer,pointer :: boundaryFaceNodeMappings(:) ! to translate indices of old
                                                ! to the ones in the new one
 real,pointer :: localFaceTangentArray(:,:) ! weighted average for coarse mesh
 logical,pointer :: isAdiabaticBoundary(:) ! for use in makeViscosity
 integer,pointer :: tripNodeFieldIndexes(:,:)
 integer :: numberOfSeparationFields,sizeOfSeparationFields
 real,pointer :: tripNodeFieldDistances(:,:) 
 real,pointer :: tripWallLength(:)

 integer :: numberOfPlotNodes
 integer,pointer :: plotNodePairs(:,:)
 real,pointer :: plotNodeDistance(:)

 type(RestrictionOperatorData),pointer :: rod
 type(ProlongationOperatorData),pointer :: pod 
 type(BoundarySolverData),pointer :: brp
 integer :: numberOfSides,numberOfNodes,numberOfBoundaryCVs
 real :: residualScalingFactor
 real,pointer :: wakeDataArray(:,:)
 real,pointer :: JamesonCoefficients(:) ! for timestepping
 integer :: gridNumber,numberOfTripNodes,iterationNumber
end type GridSolverData

contains

!-----------------------------------------------------------------------
 subroutine prepareForTimeIterations(ivd,grp)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! call setUpInitialField(ivd,grp)

 end subroutine prepareForTimeIterations
!-----------------------------------------------------------------------
 subroutine doTimeIterations(grp,ivd,relaxationSteps,calculateDissipation,calculateTurbulence,iterationNumber)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 integer :: relaxationSteps
 logical :: calculateDissipation,calculateTurbulence
 integer :: iterationNumber

 integer :: i 

! call calculateHOCoefficients(grp,ivd)
 ! time loop
 do i=1,relaxationSteps
  grp%uprev = grp%u
  call singleIteration(grp,ivd,calculateDissipation,calculateTurbulence,iterationNumber)
 end do

 end subroutine doTimeIterations
!-----------------------------------------------------------------------
 subroutine singleIteration(grp,ivd,calculateDissipation,calculateTurbulence,iterationNumber)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 logical :: calculateDissipation,calculateTurbulence
 integer :: iterationNumber

 integer :: i,j,k,ji,ind

 ind = 3260

 if(calculateDissipation) then 
  ! This section is not needed for coarse grids since
  ! timesteps and dissipation has been calculated in the 
  ! restrictGrid subroutine

  ! make pressure field 
  call makePressureField(grp,ivd,grp%u)
  grp%p = abs(grp%p)
  ! make local timesteps 
  call makeTimeSteps(grp,ivd)
  ! make artificial dissipation
  if(grp%gridNumber==1) then 
   if(ivd%dissipationScheme==1) then 
    call makeArtificialDissipation1(grp,ivd)
   else if(ivd%dissipationScheme==2) then 
    call makeArtificialDissipation2(grp,ivd)
   else
    STOP "ERROR: Selected dissipation scheme not implemented"
   end if
  else
   if(ivd%coarseGridDissipationScheme==1) then
    call makeArtificialDissipation1(grp,ivd)
   else if(ivd%coarseGridDissipationScheme==2) then
    call makeArtificialDissipation2(grp,ivd)
   else if(ivd%coarseGridDissipationScheme==3) then
    call makeArtificialDissipation3(grp,ivd)
   else
    STOP "ERROR: Selected dissipation scheme not implemented"
   end if
  end if
  if(ivd%ReynoldsNumber>0.0) then 
!   if(ivd%viscosityScheme==2.and.grp%gridNumber==1) then
   if(ivd%viscosityScheme==2) then
    call makeViscosityTerm2(grp,ivd,grp%dissipation,.true.)
   else 
    call makeViscosityTerm(grp,ivd,grp%dissipation,.true.)
   end if
   if(ivd%turbulenceModel==2) then 
    call makeBLTurbulence(grp,ivd)
   end if 
  end if
 end if
 do i=1,3  ! Jameson iterations 
  ! make pressure field 
  call makePressureField(grp,ivd,grp%u)
  ! make RHS
  call makeRHS(grp,ivd)
  grp%rhs(:,5) = 0.0
  if(ivd%turbulenceModel==1) then
  if(iterationNumber<ivd%numberOfTriggerSteps.and.grp%gridNumber==1) then
   call triggerTurbulenceField(grp,ivd)
  end if
  end if
!  if(ivd%turbulenceModel==1.and.grp%gridNumber<2) then ! KAS
  if(ivd%turbulenceModel==1.and.grp%gridNumber<4.and.calculateTurbulence) then
   call makeSARHS(grp,ivd)
   call makeSADiffusion(grp,ivd)
   call makeSASourceTerm(grp,ivd)
   if(grp%gridNumber==1) then
    call makeSATripTerm(grp,ivd)
!   else
!    grp%rhs(:,5) = grp%rhs(:,5) - grp%tripSourceArray(:)
   end if
   grp%rhs(:,5) = ivd%turbulenceCFLFactor*grp%rhs(:,5)
  end if
  ! fix RHS on outer boundaries
  !
  !OH commented this to check weak BC
  !call setRHSAtBoundary(grp,ivd) 
  ! multiply by timestep and invert mass matrix 
  call solveSystem(grp,ivd,i,grp%rhs)
  ! fix boundary conditions on increment field
  call setBCsOnIncrementField(grp,ivd,grp%rhs)
  ! smooth residual if wanted
  if(ivd%residualSmoothingFactor>0.0) then 
   do j=1,ivd%numberOfRSSteps
    call smoothResidual(grp,ivd,grp%rhs,ivd%residualSmoothingFactor)
   end do
   call setBCsOnIncrementField(grp,ivd,grp%rhs)
  end if
  ! nullify increment at outer boundary for coarse grid 
! NOR if(associated(grp%sourceTerm)) then
! NOR  call setOuterBoundaryConditions4(grp,ivd,grp%rhs)
! NOR  end if
  ! update unknown field
  call updateSolutionVector(grp%u,grp%uprev,grp%rhs,grp,ivd)
  ! fix boundary conditions on solution field
  call setBCsOnSolutionField(grp,ivd,grp%u)
  if(associated(grp%sourceTerm)) then 
   ! This means that grid is coarse and the solution
   ! field is added a source term
   grp%u = grp%u - grp%JamesonCoefficients(i)*grp%sourceTerm
  end if
 ! fix boundary conditions at outer boundary
  call setOuterBoundaryConditions3(grp,ivd,grp%u,grp%p,.true.)
 end do

 ! do turbulence iterations

if(.false.) then 
 if(ivd%turbulenceModel==1.and.grp%gridNumber<4.and.calculateTurbulence) then
  do i=1,ivd%numberOfTurbulenceSubsteps
   do j=1,3
    grp%rhs(:,5) = 0.0
    call makeSARHS(grp,ivd)
    call makeSADiffusion(grp,ivd)
    call makeSASourceTerm(grp,ivd)
    if(grp%gridNumber==1) then
     call makeSATripTerm(grp,ivd)
    end if
    grp%rhs(:,5) = ivd%turbulenceCFLFactor*grp%rhs(:,6)
    call solveTurbulenceSystem(grp,ivd,j,grp%rhs(:,5))
    grp%u(:,5) = grp%uprev(:,5) - grp%rhs(:,5)
    call setBCsOnIncrementField(grp,ivd,grp%rhs)
    if(grp%gridNumber>1) then
     grp%u(:,5) = grp%u(:,5) - grp%JamesonCoefficients(j)*grp%sourceTerm(:,5)
    end if
    call setBCsOnSolutionField(grp,ivd,grp%u)

    do k=1,grp%numberOfNodes
     grp%u(k,5) = max(grp%u(k,5),0.0)
    end do

    call setOuterBoundaryConditions3(grp,ivd,grp%u,grp%p,.true.)
   end do
  end do
 end if
end if

! do i=1,grp%numberOfNodes
!  write(798,*) i
!  write(798,'(7E15.7)') (grp%u(i,j),j=1,5),(grp%coordinates(i,j),j=1,2)
! end do

 end subroutine singleIteration 
!-----------------------------------------------------------------------
 subroutine updateSolutionVector(u,uprev,delu,grp,ivd)
 ! updates solution increment, relaxes if update is to large
 IMPLICIT NONE 
 
 real :: u(:,:),uprev(:,:),delu(:,:)
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j
 real :: increment,sign
 logical :: hasRelaxed,hasRelaxedTurb
 real :: maxincr(5)
 integer :: maxind(5)

 maxincr = -1.0
 hasRelaxed = .false.
 hasRelaxedTurb = .false.
 do i=1,grp%numberOfNodes
  ! relax density if necessary
  increment = delu(i,1)
!  if(abs(increment/uprev(i,1))<ivd%maxIncrementFactor) then 
   u(i,1) = uprev(i,1) - increment
!  else 
!   sign = increment/abs(increment)
!   u(i,1) = uprev(i,1) - sign*ivd%maxIncrementFactor*uprev(i,1)
!   hasRelaxed = .true.
!  end if
  
  ! velocities are not relaxed
  u(i,2) = uprev(i,2) - delu(i,2)
  u(i,3) = uprev(i,3) - delu(i,3)
  
  ! total energy not relaxed
  u(i,4) = uprev(i,4) - delu(i,4)

!  increment = delu(i,5)
!  if(abs(delu(i,5))>uprev(i,5)) then 
!   increment = delu(i,5)*uprev(i,5)/abs(delu(i,5))
!   hasRelaxedTurb = .true.
!  end if
  u(i,5) = uprev(i,5) - delu(i,5)
!  u(i,5) = uprev(i,5) - increment
  if(u(i,5)<0.0) u(i,5) = 0.0
  if(maxincr(5)<abs(u(i,5))) then 
   maxincr(5) = abs(u(i,5))
   maxind(5) = i
  end if

 end do
 if(hasRelaxed) then 
  write(*,*) "Warning: density has been relaxed in solution"
  maxincr = -1.0
  maxind = 0
  do i=1,grp%numberOfNodes
   do j=1,5
    if(abs(delu(i,j))>maxincr(j)) then 
     maxincr(j) = abs(delu(i,j))
     maxind(j) = i
    end if
   end do
  end do
  write(*,*) "max increment: ",maxincr
  write(*,*) "appearing at nodes: ",maxind
 end if
 end subroutine updateSolutionVector
!-----------------------------------------------------------------------
 subroutine solveSystem(grp,ivd,JamesonIteration,rhs)
 ! inverts mass matrix to get unknown from RHS
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 integer :: JamesonIteration
 real :: rhs(:,:)

 integer :: i
 real :: dt,dtJ
! multiply with time increments and divide by lumped mass

 if(JamesonIteration > 0) then 
  dtJ = ivd%CFLNumber*grp%JamesonCoefficients(JamesonIteration)
 else
  dtJ = ivd%CFLNumber
 end if

 do i=1,grp%numberOfNodes
  dt = dtJ*grp%localTimeSteps(i)
  rhs(i,:) = dt*rhs(i,:)/grp%nodeVolume(i)  
 end do
 end subroutine solveSystem 
!-----------------------------------------------------------------------
 subroutine solveTurbulenceSystem(grp,ivd,JamesonIteration,rhs)
 ! inverts mass matrix to get unknown from RHS
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 integer :: JamesonIteration
 real :: rhs(:)

 integer :: i
 real :: dt,dtJ
! multiply with time increments and divide by lumped mass

 if(JamesonIteration > 0) then
  dtJ = ivd%CFLNumber*grp%JamesonCoefficients(JamesonIteration)
 else
  dtJ = ivd%CFLNumber
 end if

 do i=1,grp%numberOfNodes
  dt = dtJ*grp%localTimeSteps(i)
  rhs(i) = dt*rhs(i)/grp%nodeVolume(i)
 end do
 end subroutine solveTurbulenceSystem
!-----------------------------------------------------------------------
 subroutine calculateHOCoefficients(grp,ivd)
 ! calculates the coefficients of the high-order
 ! discretisations of the variables 

 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,ind1,ind2,i0,i1,i2,i3,i4,i5,i6
 real :: y0,y1,y2,y3,y4,y5,y6,x0,x1,x2,x3,x4,x5,x6
 real :: h0,h1,h2,h3,h4,y,w1,w2,w3,w4,k

 k=1.0
 do i=1,grp%numberOfSides
  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  if(ivd%HighOrder.and.grp%gridNumber==1) then
    if(grp%AfterAndBefore(ind1,1).eq.ind2.or.grp%AfterAndBefore(ind1,2).eq.ind2)then
      if(grp%AfterAndBefore(ind1,1).eq.ind2) then
        i2 = ind1
        i3 = ind2
        i4 = grp%AfterAndBefore(i3,1)
        i1 = grp%AfterAndBefore(i2,2)
        if(i4.ne.0) then
          i5 = grp%AfterAndBefore(i4,1)
        else
          i5 = 0
        end if
        if(i1.ne.0) then
          i0 = grp%AfterAndBefore(i1,2)
        else
          i0 = 0
        end if
      else if(grp%AfterAndBefore(ind1,2).eq.ind2) then
        i2 = ind2
        i3 = ind1
        i4 = grp%AfterAndBefore(i3,1)
        i1 = grp%AfterAndBefore(i2,2)
        if(i4.ne.0) then
          i5 = grp%AfterAndBefore(i4,1)
        else
          i5 = 0
        end if
        if(i1.ne.0) then
          i0 = grp%AfterAndBefore(i1,2)
        else
          i0 = 0
        end if
      else
        print *,' should not be here'
      end if
      if(i1.eq.0) then
        i0 = 1 !required?
        i6 = grp%AfterAndBefore(i5,1)

        y2 = grp%coordinates(i2,2)
        y3 = grp%coordinates(i3,2)
        y4 = grp%coordinates(i4,2)
        y5 = grp%coordinates(i5,2)
        y6 = grp%coordinates(i6,2)

        x2 = grp%coordinates(i2,1)
        x3 = grp%coordinates(i3,1)
        x4 = grp%coordinates(i4,1)
        x5 = grp%coordinates(i5,1)
        x6 = grp%coordinates(i6,1)

        h0 = sqrt((x3-x2)**2+(y3-y2)**2)
        h1 = sqrt((x4-x3)**2+(y4-y3)**2)
        h2 = sqrt((x5-x4)**2+(y5-y4)**2)
        h3 = sqrt((x6-x5)**2+(y6-y5)**2)

        w1 = h0/2.0!CV widths
        w2 = (h0+h1)/2.0
        w3 = (h1+h2)/2.0
        w4 = (h2+h3)/2.0

        y = w1
        
        !debugging
        !if (k==1.0) then
        !  w1=0.05
        !  w2=0.1
        !  w3=0.1
        !  w4=0.1
        !  y=0.05
        !  k=0.0
        !endif


        !coefficients for equation(19) which assumes a point boundary value
        grp%hoc(i,1) =(2.0*w1*w3*w4+4.0*w1*w2*w4+w2*w3*w4+8.0*w1*w2*w3-4.0*y**3&
                      +4.0*w1**3+9.0*w1**2*w2+6.0*w1*w2**2+6.0*w1**2*w3+2.0*w1&
                      *w3**2+2.0*w2**2*w3+w2*w3**2+w2**3+3.0*w1**2*w4+w2**2*w4&
                      +9.0*y**2*w2+6.0*y**2*w3+12.0*y**2*w1-6.0*y*w2**2-2.0*&
                      y*w3**2+3.0*y**2*w4-12.0*y*w1**2-18.0*y*w1*w2-4.0*y*w2*w4-&
                      8.0*y*w2*w3-2.0*y*w3*w4-6.0*y*w1*w4-12.0*y*w1*w3)/(w2**3&
                      +2.0*w2**2*w3+w2**2*w4+6.0*w1*w2**2+w2*w3**2+w2*w3*w4+8.0&
                      *w1*w2*w3+9.0*w1**2*w2+4.0*w1*w2*w4+2.0*w1*w3**2+6.0*&
                      w1**2*w3+2.0*w1*w3*w4+3.0*w1**2*w4+4.0*w1**3)

        grp%hoc(i,2) = y*(12.0*w1*w3**3+6.0*w1**2*w4**2+32.0*w1**3*w3+16.0*&
                       w1**3*w4+12.0*w1*w2**3+30.0*w1**2*w2**2+4.0*w2*w3*w4&
                       **2+32.0*w1**3*w2+36.0*w1*w2*w3*w4+30.0*w1**2*w3**2+&
                       12.0*w1**4+12.0*w2**2*w3**2+8.0*w2*w3**3+2.0*w2**2*w4&
                       **2+4.0*w2**3*w4+8.0*w2**3*w3+18.0*w1*w3**2*w4+30.0*&
                       w1**2*w3*w4+6.0*w1*w3*w4**2+60.0*w1**2*w2*w3+36.0*w1&
                       *w2*w3**2+12.0*w2*w3**2*w4+6*w1*w2*w4**2+30.0*w1**2*&
                       w2*w4+18.0*w1*w2**2*w4+2.0*w3**4+12.0*w2**2*w3*w4+&
                       36.0*w1*w2**2*w3+2.0*w2**4+16.0*y**2*w1*w3+16.0*y**2&
                       *w1*w2+8.0*y**2*w2*w3+4.0*y**2*w2*w4+8.0*y**2*w1*w4-&
                       6.0*y*w1*w4**2-48.0*y*w1**2*w2-30.0*y*w1*w2**2-3.0*y&
                       *w2*w4**2-30.0*y*w1*w3**2-18.0*y*w2*w3**2-9.0*y*w3**2&
                       *w4-48.0*y*w1**2*w3-18.0*y*w2**2*w3-3.0*y*w3*w4**2&
                       -9.0*y*w2**2*w4-24.0*y*w1**2*w4+4.0*y**2*w3**2-6.0*y&
                       *w3**3+4.0*y**2*w3*w4+2.0*w3**2*w4**2+4.0*w3**3*w4+&
                       12.0*y**2*w1**2+4.0*y**2*w2**2-6.0*y*w2**3-24.0*y*w1&
                       **3-30.0*y*w1*w2*w4-18.0*y*w2*w3*w4-60.0*y*w1*w2*w3-&
                       30.0*y*w1*w3*w4)/(w3+w4)/(w2**3+2.0*w2**2*w3+w2**2.0*w4&
                       +6.0*w1*w2**2+w2*w3**2+w2*w3*w4+8.0*w1*w2*w3+9.0*w1**2&
                       *w2+4.0*w1*w2*w4+2.0*w1*w3**2+6.0*w1**2*w3+2.0*w1*w3*&
                       w4+3.0*w1**2*w4+4.0*w1**3)/w2/w3

        grp%hoc(i,3) = -y*(6.0*w1**2*w4**2+16.0*w1**3*w3+16.0*w1**3*w4+12.0*&
                       w1*w2**3+30.0*w1**2*w2**2+32.0*w1**3*w2+12.0*w1*w2*w3&
                       *w4+6.0*w1**2*w3**2+12.0*w1**4+2.0*w2**2*w3**2+2.0*w2&
                       **2*w4**2+4.0*w2**3*w4+4.0*w2**3*w3+12.0*w1**2*w3*w4+&
                       30.0*w1**2*w2*w3+6.0*w1*w2*w3**2+6.0*w1*w2*w4**2+30.0&
                       *w1**2*w2*w4+18.0*w1*w2**2*w4+4.0*w2**2*w3*w4+18.0*w1&
                       *w2**2*w3+2*w2**4+8.0*y**2*w1*w3+16.0*y**2*w1*w2+4.0*&
                       y**2*w2*w3+4.0*y**2*w2*w4+8.0*y**2*w1*w4-6.0*y*w1*w4&
                       **2-48.0*y*w1**2*w2-30*y*w1*w2**2-3.0*y*w2*w4**2-6.0*&
                       y*w1*w3**2-3.0*y*w2*w3**2-24.0*y*w1**2*w3-9.0*y*w2**2&
                       *w3-9.0*y*w2**2*w4-24.0*y*w1**2*w4+12.0*y**2*w1**2+4.0&
                       *y**2*w2**2-6.0*y*w2**3-24.0*y*w1**3-30.0*y*w1*w2*w4&
                       -6.0*y*w2*w3*w4-30.0*y*w1*w2*w3-12.0*y*w1*w3*w4)/&
                       (w2+w3)/(w2**3+2.0*w2**2*w3+w2**2*w4+6.0*w1*w2**2+w2*&
                       w3**2+w2*w3*w4+8.0*w1*w2*w3+9.0*w1**2*w2+4.0*w1*w2*w4&
                       +2.0*w1*w3**2+6.0*w1**2*w3+2*w1*w3*w4+3.0*w1**2*w4+4.0&
                       *w1**3)/w3/w4

        grp%hoc(i,4) = y*(16.0*w1**3*w3+12.0*w1*w2**3+30.0*w1**2*w2**2+32.0*&
                       w1**3*w2+6.0*w1**2*w3**2+12.0*w1**4+2.0*w2**2*w3**2+&
                       4.0*w2**3*w3+30.0*w1**2*w2*w3+6.0*w1*w2*w3**2+18.0*w1&
                       *w2**2*w3+2.0*w2**4+8.0*y**2*w1*w3+16.0*y**2*w1*w2+&
                       4.0*y**2*w2*w3-48.0*y*w1**2*w2-30.0*y*w1*w2**2-6.0*&
                       y*w1*w3**2-3.0*y*w2*w3**2-24.0*y*w1**2*w3-9.0*y*w2**2&
                       *w3+12.0*y**2*w1**2+4.0*y**2*w2**2-6.0*y*w2**3-24.0*y&
                       *w1**3-30.0*y*w1*w2*w3)/(w3+w4)/(w2+w3+w4)/(w2**3+2.0&
                       *w2**2*w3+w2**2*w4+6.0*w1*w2**2+w2*w3**2+w2*w3*w4&
                      +8.0*w1*w2*w3+9.0*w1**2*w2+4.0*w1*w2*w4+2.0*w1*w3**2+6.0*&
                      w1**2*w3+2.0*w1*w3*w4+3.0*w1**2*w4+4.0*w1**3)/w4
         
         !coefficients for equation(33) which assumes mean boundary values

         grp%hoc2(i,1) = -(-2.0*w1*w2*w4-w2*w3*w4-w1*w3*w4+8.0*y*w2*w3+&
                      4*y**3+8.0*y*w1*w3+12.0*y*w1*w2-4.0*w1*w2*w3+6.0*y*&
                      w2**2-9.0*y**2*w1-9.0*y**2*w2-6.0*y**2*w3+6.0*y*w1&
                      **2-w1**3-3.0*w1**2*w2-3.0*w1*w2**2-2.0*w1**2*w3-w1&
                      *w3**2-2.0*w2**2*w3-w2*w3**2-w2**3-w1**2*w4-w2**2*w4&
                      +2.0*y*w3*w4+4.0*y*w1*w4+4.0*y*w2*w4-3.0*y**2*w4+2.0&
                      *y*w3**2)/(w3+w2)/(w2+w3+w4)/w1/w2

        grp%hoc2(i,2) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-6.0*y**2*w3-3.0*&
                        y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+8.0*y*w1*w3+4.0*y*&
                        w1*w4+2.0*y*w3**2+4.0*y*w2*w3+2.0*y*w3*w4+2.0*y*w2&
                        **2+2.0*y*w2*w4-w1**3-2.0*w1**2*w2-2.0*w1**2*w3-w1&
                        **2*w4-w1*w3**2-2.0*w1*w2*w3-w1*w3*w4-w1*w2**2-w1*&
                        w2*w4)/(w3+w4)/(w1+w2)/w2/w3

        grp%hoc2(i,3) = -(4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3-3.0&
                        *y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+4.0*y&
                        *w1*w4+2.0*y*w2**2+2.0*y*w2*w3+2.0*y*w2*w4-w1**3-&
                        2.0*w1**2*w2-w1**2*w3-w1**2*w4-w1*w2**2-w1*w2*w3-&
                        w1*w2*w4)/(w3+w2)/(w1+w3+w2)/w3/w4

        grp%hoc2(i,4) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3+6.0&
                        *y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+2.0*y*w2**2+2.0*&
                        y*w2*w3-w1**3-2.0*w1**2*w2-w1**2*w3-w1*w2**2-w1*w2&
                        *w3)/(w3+w4)/(w2+w3+w4)/(w2+w3+w4+w1)/w4

         !coefficients for equation(20) which assume mean boundary values   
      y = 0.
         grp%hoc3(i,1) = -(-2.0*w1*w2*w4-w2*w3*w4-w1*w3*w4+8.0*y*w2*w3+&
                      4.0*y**3+8.0*y*w1*w3+12.0*y*w1*w2-4.0*w1*w2*w3+6.0*y*&
                      w2**2-9.0*y**2*w1-9.0*y**2*w2-6.0*y**2*w3+6.0*y*w1&
                      **2-w1**3-3.0*w1**2*w2-3.0*w1*w2**2-2.0*w1**2*w3-w1&
                      *w3**2-2.0*w2**2*w3-w2*w3**2-w2**3-w1**2*w4-w2**2*w4&
                      +2.0*y*w3*w4+4.0*y*w1*w4+4.0*y*w2*w4-3.0*y**2*w4+2.0&
                      *y*w3**2)/(w3+w2)/(w2+w3+w4)/w1/w2

        grp%hoc3(i,2) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-6.0*y**2*w3-3.0*&
                        y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+8.0*y*w1*w3+4.0*y*&
                        w1*w4+2.0*y*w3**2+4.0*y*w2*w3+2.0*y*w3*w4+2.0*y*w2&
                        **2+2.0*y*w2*w4-w1**3-2.0*w1**2*w2-2.0*w1**2*w3-w1&
                        **2*w4-w1*w3**2-2.0*w1*w2*w3-w1*w3*w4-w1*w2**2-w1*&
                        w2*w4)/(w3+w4)/(w1+w2)/w2/w3

        grp%hoc3(i,3) = -(4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3-3.0&
                        *y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+4.0*y&
                        *w1*w4+2.0*y*w2**2+2.0*y*w2*w3+2.0*y*w2*w4-w1**3-&
                        2.0*w1**2*w2-w1**2*w3-w1**2*w4-w1*w2**2-w1*w2*w3-&
                        w1*w2*w4)/(w3+w2)/(w1+w3+w2)/w3/w4

        grp%hoc3(i,4) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3+6.0&
                        *y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+2.0*y*w2**2+2.0*&
                        y*w2*w3-w1**3-2.0*w1**2*w2-w1**2*w3-w1*w2**2-w1*w2&
                        *w3)/(w3+w4)/(w2+w3+w4)/(w2+w3+w4+w1)/w4


      else if(i0.eq.0) then
        y1 = grp%coordinates(i1,2)
        y2 = grp%coordinates(i2,2)
        y3 = grp%coordinates(i3,2)
        y4 = grp%coordinates(i4,2)
        y5 = grp%coordinates(i5,2)

        x1 = grp%coordinates(i1,1)
        x2 = grp%coordinates(i2,1)
        x3 = grp%coordinates(i3,1)
        x4 = grp%coordinates(i4,1)
        x5 = grp%coordinates(i5,1)

        h0 = sqrt((x2-x1)**2+(y2-y1)**2)
        h1 = sqrt((x3-x2)**2+(y3-y2)**2)
        h2 = sqrt((x4-x3)**2+(y4-y3)**2)
        h3 = sqrt((x5-x4)**2+(y5-y4)**2)

        w1 = h0/2.0!CV widths
        w2 = (h0+h1)/2.0
        w3 = (h1+h2)/2.0
        w4 = (h2+h3)/2.0

        y = w1+w2
        
        !debugging
!        if (k==1.0) then
!          w1=0.05
!          w2=0.1
!          w3=0.1
!          w4=0.1
!          y=0.15
!          k=0.0
!         end if

        !coefficients for equation(19) which assumes a point boundary value
        grp%hoc(i,1) =(2.0*w1*w3*w4+4.0*w1*w2*w4+w2*w3*w4+8.0*w1*w2*w3-4.0*y**3&
                      +4.0*w1**3+9.0*w1**2*w2+6.0*w1*w2**2+6.0*w1**2*w3+2.0*w1&
                      *w3**2+2.0*w2**2*w3+w2*w3**2+w2**3+3.0*w1**2*w4+w2**2*w4&
                      +9.0*y**2*w2+6.0*y**2*w3+12.0*y**2*w1-6.0*y*w2**2-2.0*&
                      y*w3**2+3*y**2*w4-12.0*y*w1**2-18.0*y*w1*w2-4.0*y*w2*w4-&
                      8.0*y*w2*w3-2.0*y*w3*w4-6.0*y*w1*w4-12.0*y*w1*w3)/(w2**3&
                      +2.0*w2**2*w3+w2**2*w4+6.0*w1*w2**2+w2*w3**2+w2*w3*w4+8.0&
                      *w1*w2*w3+9.0*w1**2*w2+4.0*w1*w2*w4+2.0*w1*w3**2+6.0*&
                      w1**2*w3+2.0*w1*w3*w4+3.0*w1**2*w4+4.0*w1**3)

        grp%hoc(i,2) = y*(12.0*w1*w3**3+6.0*w1**2*w4**2+32.0*w1**3*w3+16.0*&
                       w1**3*w4+12.0*w1*w2**3+30.0*w1**2*w2**2+4.0*w2*w3*w4&
                       **2+32.0*w1**3*w2+36.0*w1*w2*w3*w4+30.0*w1**2*w3**2+&
                       12.0*w1**4+12.0*w2**2*w3**2+8.0*w2*w3**3+2.0*w2**2*w4&
                       **2+4.0*w2**3*w4+8.0*w2**3*w3+18.0*w1*w3**2*w4+30.0*&
                       w1**2*w3*w4+6.0*w1*w3*w4**2+60.0*w1**2*w2*w3+36.0*w1&
                       *w2*w3**2+12.0*w2*w3**2*w4+6*w1*w2*w4**2+30.0*w1**2*&
                       w2*w4+18.0*w1*w2**2*w4+2.0*w3**4+12.0*w2**2*w3*w4+&
                       36.0*w1*w2**2*w3+2.0*w2**4+16.0*y**2*w1*w3+16.0*y**2&
                       *w1*w2+8.0*y**2*w2*w3+4.0*y**2*w2*w4+8.0*y**2*w1*w4-&
                       6.0*y*w1*w4**2-48.0*y*w1**2*w2-30.0*y*w1*w2**2-3.0*y&
                       *w2*w4**2-30.0*y*w1*w3**2-18.0*y*w2*w3**2-9.0*y*w3**2&
                       *w4-48.0*y*w1**2*w3-18.0*y*w2**2*w3-3.0*y*w3*w4**2&
                       -9.0*y*w2**2*w4-24.0*y*w1**2*w4+4.0*y**2*w3**2-6.0*y&
                       *w3**3+4.0*y**2*w3*w4+2.0*w3**2*w4**2+4.0*w3**3*w4+&
                       12.0*y**2*w1**2+4.0*y**2*w2**2-6.0*y*w2**3-24.0*y*w1&
                       **3-30.0*y*w1*w2*w4-18.0*y*w2*w3*w4-60.0*y*w1*w2*w3-&
                       30.0*y*w1*w3*w4)/(w3+w4)/(w2**3+2*w2**2*w3+w2**2.0*w4&
                       +6.0*w1*w2**2+w2*w3**2+w2*w3*w4+8.0*w1*w2*w3+9.0*w1**2&
                       *w2+4.0*w1*w2*w4+2.0*w1*w3**2+6.0*w1**2*w3+2.0*w1*w3*&
                       w4+3.0*w1**2*w4+4.0*w1**3)/w2/w3

        grp%hoc(i,3) = -y*(6.0*w1**2*w4**2+16.0*w1**3*w3+16.0*w1**3*w4+12.0*&
                       w1*w2**3+30.0*w1**2*w2**2+32.0*w1**3*w2+12.0*w1*w2*w3&
                       *w4+6.0*w1**2*w3**2+12.0*w1**4+2.0*w2**2*w3**2+2.0*w2&
                       **2*w4**2+4.0*w2**3*w4+4.0*w2**3*w3+12.0*w1**2*w3*w4+&
                       30.0*w1**2*w2*w3+6.0*w1*w2*w3**2+6.0*w1*w2*w4**2+30.0&
                       *w1**2*w2*w4+18.0*w1*w2**2*w4+4.0*w2**2*w3*w4+18.0*w1&
                       *w2**2*w3+2*w2**4+8.0*y**2*w1*w3+16.0*y**2*w1*w2+4.0*&
                       y**2*w2*w3+4.0*y**2*w2*w4+8.0*y**2*w1*w4-6.0*y*w1*w4&
                       **2-48.0*y*w1**2*w2-30*y*w1*w2**2-3.0*y*w2*w4**2-6.0*&
                       y*w1*w3**2-3.0*y*w2*w3**2-24.0*y*w1**2*w3-9.0*y*w2**2&
                       *w3-9.0*y*w2**2*w4-24.0*y*w1**2*w4+12.0*y**2*w1**2+4.0&
                       *y**2*w2**2-6.0*y*w2**3-24.0*y*w1**3-30.0*y*w1*w2*w4&
                       -6.0*y*w2*w3*w4-30.0*y*w1*w2*w3-12.0*y*w1*w3*w4)/&
                       (w2+w3)/(w2**3+2.0*w2**2*w3+w2**2*w4+6.0*w1*w2**2+w2*&
                       w3**2+w2*w3*w4+8.0*w1*w2*w3+9.0*w1**2*w2+4.0*w1*w2*w4&
                       +2.0*w1*w3**2+6.0*w1**2*w3+2*w1*w3*w4+3.0*w1**2*w4+4.0&
                       *w1**3)/w3/w4

        grp%hoc(i,4) = y*(16.0*w1**3*w3+12.0*w1*w2**3+30.0*w1**2*w2**2+32.0*&
                       w1**3*w2+6.0*w1**2*w3**2+12.0*w1**4+2.0*w2**2*w3**2+&
                       4.0*w2**3*w3+30.0*w1**2*w2*w3+6.0*w1*w2*w3**2+18.0*w1&
                       *w2**2*w3+2.0*w2**4+8.0*y**2*w1*w3+16.0*y**2*w1*w2+&
                       4.0*y**2*w2*w3-48.0*y*w1**2*w2-30.0*y*w1*w2**2-6.0*&
                       y*w1*w3**2-3.0*y*w2*w3**2-24.0*y*w1**2*w3-9.0*y*w2**2&
                       *w3+12.0*y**2*w1**2+4.0*y**2*w2**2-6.0*y*w2**3-24.0*y&
                       *w1**3-30.0*y*w1*w2*w3)/(w3+w4)/(w2+w3+w4)/(w2**3+2.0&
                       *w2**2*w3+w2**2*w4+6.0*w1*w2**2+w2*w3**2+w2*w3*w4&
                      +8.0*w1*w2*w3+9.0*w1**2*w2+4.0*w1*w2*w4+2.0*w1*w3**2+6*&
                      w1**2*w3+2.0*w1*w3*w4+3.0*w1**2*w4+4.0*w1**3)/w4

         !coefficients for equation(33) which assumes mean boundary values

         grp%hoc2(i,1) = -(-2.0*w1*w2*w4-w2*w3*w4-w1*w3*w4+8.0*y*w2*w3+&
                      4*y**3+8.0*y*w1*w3+12.0*y*w1*w2-4.0*w1*w2*w3+6.0*y*&
                      w2**2-9.0*y**2*w1-9.0*y**2*w2-6.0*y**2*w3+6.0*y*w1&
                      **2-w1**3-3.0*w1**2*w2-3.0*w1*w2**2-2.0*w1**2*w3-w1&
                      *w3**2-2.0*w2**2*w3-w2*w3**2-w2**3-w1**2*w4-w2**2*w4&
                      +2.0*y*w3*w4+4.0*y*w1*w4+4.0*y*w2*w4-3.0*y**2*w4+2.0&
                      *y*w3**2)/(w3+w2)/(w2+w3+w4)/w1/w2

        grp%hoc2(i,2) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-6.0*y**2*w3-3.0*&
                        y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+8.0*y*w1*w3+4.0*y*&
                        w1*w4+2.0*y*w3**2+4.0*y*w2*w3+2.0*y*w3*w4+2.0*y*w2&
                        **2+2.0*y*w2*w4-w1**3-2.0*w1**2*w2-2.0*w1**2*w3-w1&
                        **2*w4-w1*w3**2-2.0*w1*w2*w3-w1*w3*w4-w1*w2**2-w1*&
                        w2*w4)/(w3+w4)/(w1+w2)/w2/w3

        grp%hoc2(i,3) = -(4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3-3.0&
                        *y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+4.0*y&
                        *w1*w4+2.0*y*w2**2+2.0*y*w2*w3+2.0*y*w2*w4-w1**3-&
                        2.0*w1**2*w2-w1**2*w3-w1**2*w4-w1*w2**2-w1*w2*w3-&
                        w1*w2*w4)/(w3+w2)/(w1+w3+w2)/w3/w4

        grp%hoc2(i,4) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3+6.0&
                        *y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+2.0*y*w2**2+2.0*&
                        y*w2*w3-w1**3-2.0*w1**2*w2-w1**2*w3-w1*w2**2-w1*w2&
                        *w3)/(w3+w4)/(w2+w3+w4)/(w2+w3+w4+w1)/w4

         !coefficients for equation(20) which assume mean boundary values 
      y=0.
         grp%hoc3(i,1) = -(-2.0*w1*w2*w4-w2*w3*w4-w1*w3*w4+8.0*y*w2*w3+&
                      4*y**3+8.0*y*w1*w3+12.0*y*w1*w2-4.0*w1*w2*w3+6.0*y*&
                      w2**2-9.0*y**2*w1-9.0*y**2*w2-6.0*y**2*w3+6.0*y*w1&
                      **2-w1**3-3.0*w1**2*w2-3.0*w1*w2**2-2.0*w1**2*w3-w1&
                      *w3**2-2.0*w2**2*w3-w2*w3**2-w2**3-w1**2*w4-w2**2*w4&
                      +2.0*y*w3*w4+4.0*y*w1*w4+4.0*y*w2*w4-3.0*y**2*w4+2.0&
                      *y*w3**2)/(w3+w2)/(w2+w3+w4)/w1/w2

        grp%hoc3(i,2) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-6.0*y**2*w3-3.0*&
                        y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+8.0*y*w1*w3+4.0*y*&
                        w1*w4+2.0*y*w3**2+4.0*y*w2*w3+2.0*y*w3*w4+2.0*y*w2&
                        **2+2.0*y*w2*w4-w1**3-2.0*w1**2*w2-2.0*w1**2*w3-w1&
                        **2*w4-w1*w3**2-2.0*w1*w2*w3-w1*w3*w4-w1*w2**2-w1*&
                        w2*w4)/(w3+w4)/(w1+w2)/w2/w3

        grp%hoc3(i,3) = -(4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3-3.0&
                        *y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+4.0*y&
                        *w1*w4+2.0*y*w2**2+2.0*y*w2*w3+2.0*y*w2*w4-w1**3-&
                        2.0*w1**2*w2-w1**2*w3-w1**2*w4-w1*w2**2-w1*w2*w3-&
                        w1*w2*w4)/(w3+w2)/(w1+w3+w2)/w3/w4

        grp%hoc3(i,4) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3+6.0&
                        *y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+2.0*y*w2**2+2.0*&
                        y*w2*w3-w1**3-2.0*w1**2*w2-w1**2*w3-w1*w2**2-w1*w2&
                        *w3)/(w3+w4)/(w2+w3+w4)/(w2+w3+w4+w1)/w4


      else if(i5.eq.0) then
      else 
        y0 = grp%coordinates(i0,2)
        y1 = grp%coordinates(i1,2)
        y2 = grp%coordinates(i2,2)
        y3 = grp%coordinates(i3,2)
        y4 = grp%coordinates(i4,2)
        y5 = grp%coordinates(i5,2)
        x0 = grp%coordinates(i0,1)
        x1 = grp%coordinates(i1,1)
        x2 = grp%coordinates(i2,1)
        x3 = grp%coordinates(i3,1)
        x4 = grp%coordinates(i4,1)
        x5 = grp%coordinates(i5,1)
        h0 = sqrt((x1-x0)**2+(y1-y0)**2)!distance between nodes
        h1 = sqrt((x2-x1)**2+(y2-y1)**2)
        h2 = sqrt((x3-x2)**2+(y3-y2)**2)
        h3 = sqrt((x4-x3)**2+(y4-y3)**2)
        h4 = sqrt((x5-x4)**2+(y5-y4)**2)

        w1 = (h0+h1)/2.0!CV widths
        w2 = (h1+h2)/2.0
        w3 = (h2+h3)/2.0
        w4 = (h3+h4)/2.0

        y = w1+w2
        !debugging
!        if (k==1.0) then
!          w1=0.1
!          w2=0.1
!          w3=0.1
!          w4=0.1
!          y=0.2
!          k=0.0
!         endif
        !coefficients for equation(10) 
         grp%hoc(i,1) = -(-2.0*w1*w2*w4-w2*w3*w4-w1*w3*w4+8.0*y*w2*w3+&
                      4*y**3+8.0*y*w1*w3+12.0*y*w1*w2-4.0*w1*w2*w3+6.0*y*&
                      w2**2-9.0*y**2*w1-9.0*y**2*w2-6.0*y**2*w3+6.0*y*w1&
                      **2-w1**3-3.0*w1**2*w2-3.0*w1*w2**2-2.0*w1**2*w3-w1&
                      *w3**2-2.0*w2**2*w3-w2*w3**2-w2**3-w1**2*w4-w2**2*w4&
                      +2.0*y*w3*w4+4.0*y*w1*w4+4.0*y*w2*w4-3.0*y**2*w4+2.0&
                      *y*w3**2)/(w3+w2)/(w2+w3+w4)/w1/w2

        grp%hoc(i,2) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-6.0*y**2*w3-3.0*&
                        y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+8.0*y*w1*w3+4.0*y*&
                        w1*w4+2.0*y*w3**2+4.0*y*w2*w3+2.0*y*w3*w4+2.0*y*w2&
                        **2+2.0*y*w2*w4-w1**3-2.0*w1**2*w2-2.0*w1**2*w3-w1&
                        **2*w4-w1*w3**2-2.0*w1*w2*w3-w1*w3*w4-w1*w2**2-w1*&
                        w2*w4)/(w3+w4)/(w1+w2)/w2/w3

        grp%hoc(i,3) = -(4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3-3.0&
                        *y**2*w4+6.0*y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+4.0*y&
                        *w1*w4+2.0*y*w2**2+2.0*y*w2*w3+2.0*y*w2*w4-w1**3-&
                        2.0*w1**2*w2-w1**2*w3-w1**2*w4-w1*w2**2-w1*w2*w3-&
                        w1*w2*w4)/(w3+w2)/(w1+w3+w2)/w3/w4

        grp%hoc(i,4) = (4.0*y**3-9.0*y**2*w1-6.0*y**2*w2-3.0*y**2*w3+6.0&
                        *y*w1**2+8.0*y*w1*w2+4.0*y*w1*w3+2.0*y*w2**2+2.0*&
                        y*w2*w3-w1**3-2.0*w1**2*w2-w1**2*w3-w1*w2**2-w1*w2&
                        *w3)/(w3+w4)/(w2+w3+w4)/(w2+w3+w4+w1)/w4

        !coefficients for equation(15)
        grp%hoc2(i,1) = grp%hoc(i,1) 
        grp%hoc2(i,2) = grp%hoc(i,2)
        grp%hoc2(i,3) = grp%hoc(i,3) 
        grp%hoc2(i,4) = grp%hoc(i,4)
 
      end if
     end if
    end if
  end do
 write(*,'(A)') "   High Order Coefficients Calculated"
 end subroutine calculateHOCoefficients
!-----------------------------------------------------------------------
 subroutine makeRHS(grp,ivd)
 ! makes right hand side by utilizing a side-based structure
 ! with coefficients created in preprocessor

 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,ind1,ind2,faceIndicator,i0,i1,i2,i3,i4,i5,im1

 real :: pr1,pr2,rinv1,rinv2,ru1,rv1,ru2,rv2,u1,v1,u2,v2,rh1,rh2,k
 real :: pr0,pr3,rinv0,rinv3,ru0,rv0,ru3,rv3,u0,v0,u3,v3,rh0,rh3,a
 real :: pr4,rinv4,ru4,rv4,u4,v4,rh4,ro0,ro1,ro2,ro3,ro4,ro5
 real :: pr5,rinv5,ru5,rv5,u5,v5,rh5,a0,a1,a2,a3,a4,a5,am
 real :: c0,c1,c2,c3,ym1,y0,y1,y2,y3,y4,y5,xm1,x0,x1,x2,x3,x4,x5
 real :: h0,h1,h2,h3,h4,w1,w2,w3,w4,ia0,ia1,ia2,ia3,ia4,amb
 real :: ut0,ut1,ut2,ut3,ut4,ut5,un0,un1,un2,un3,un4,un5,en0,en1,en2,en3,en4,en5
 real :: rom,utm,unm,uxm,vym,prm,enm,anx,any
 real :: f01x,f01y,f02x,f02y,f03x,f03y,f04x,f04y
 real :: f11x,f11y,f12x,f12y,f13x,f13y,f14x,f14y
 real :: f21x,f21y,f22x,f22y,f23x,f23y,f24x,f24y
 real :: f31x,f31y,f32x,f32y,f33x,f33y,f34x,f34y
 real :: f1mx,f1my,f2mx,f2my,f3mx,f3my,f4mx,f4my
 real :: rm,rum,rvm,rem,um,vm,em,pm,rhm,gammaMinusOne 
 real :: f1,f2,f3,f4,wx,wy,rkx,rky,rnorm,unorm1,unorm2,sc
 ! initialize

real :: rnx,rny,rho1,normalVelocity1
real :: rho2,normalVelocity2,vml,aml
real :: ro,uo,vo,eo,po,unr,rhor,uxr,vyr,wzr,epsr,presr,hr
real :: f11,f21,f31,f41,f51,f12,f22,f32,f42,f52,eps1,eps1l,eps2l
real :: fr11,fr21,fr31,fr41,fr51,fr12,fr22,fr32,fr42,fr52
real :: fl11,fl21,fl31,fl41,fl51,fl12,fl22,fl32,fl42,fl52
real :: di,d1,ui,vi,wi,hi,ci2,ci,af,ucp,rh1l,rh2l,h1l,h2l
real :: du1,du2,du3,du4,du5,rlam1,rlam2,rlam3,gam1,epslm
real :: s1,s2,al1x,al2x,cc1,cc2,f1l,f2l,f3l,f4l,f5l
real :: pr1l,rho1l,rinv1l,ru1l,rv1l,rw1l,u1l,v1l,w1l,normalVelocity1l
real :: pr2l,rho2l,rinv2l,ru2l,rv2l,rw2l,u2l,v2l,w2l,normalVelocity2l
integer :: ind

real :: r1,E1,r2,E2,pu1,pv1,pw1,pu2,pv2,pw2




 epslm = 1.0e-5
 eps1  = 1./(max(epslm,1.e-5))
 gammaMinusOne = ivd%gamma-1.0
 grp%rhs = grp%dissipation

 rewind (698)
!if (k.ne.2.0) then
!    k=1.0
!end if
 do i=1,grp%numberOfSides
  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  pr1 = grp%p(ind1)             ! pressure
  rinv1 = 1./grp%u(ind1,1)      ! rho inverse
  ru1 = grp%u(ind1,2)           ! x-momentum
  rv1 = grp%u(ind1,3)           ! y-momentum
  u1 = rinv1*ru1                ! x-velocity
  v1 = rinv1*rv1                ! y-velocity
  rh1 = grp%u(ind1,4) + pr1     ! enthalpy
   
  pr2 = grp%p(ind2)
  rinv2 = 1./grp%u(ind2,1)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  u2 = rinv2*ru2
  v2 = rinv2*rv2
  rh2 = grp%u(ind2,4) + pr2 

  f11x = ru1
  f12x = ru1*u1+pr1
  f13x = rv1*u1
  f14x = rh1*u1
  f11y = rv1
  f12y = ru1*v1
  f13y = rv1*v1+pr1
  f14y = rh1*v1

  f21x = ru2
  f22x = ru2*u2+pr2
  f23x = rv2*u2
  f24x = rh2*u2
  f21y = rv2
  f22y = ru2*v2
  f23y = rv2*v2+pr2
  f24y = rh2*v2

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)  
  anx = wx/sqrt(wx**2+wy**2)
  any = wy/sqrt(wx**2+wy**2)

  if(ivd%HighOrder.and.grp%gridNumber==1) then
  ! if(.false.) then
    if(grp%AfterAndBefore(ind1,1).eq.ind2.or.grp%AfterAndBefore(ind1,2).eq.ind2)then
      if(grp%AfterAndBefore(ind1,1).eq.ind2) then
        i1 = ind1
        i2 = ind2
        i3 = grp%AfterAndBefore(i2,1)
        i4 = grp%AfterAndBefore(i3,1)
        i0 = grp%AfterAndBefore(i1,2)
        im1 = grp%AfterAndBefore(i0,2)
      else if(grp%AfterAndBefore(ind1,2).eq.ind2) then
        i1 = ind2
        i2 = ind1
        i3 = grp%AfterAndBefore(i2,1)
        i4 = grp%AfterAndBefore(i3,1)
        i0 = grp%AfterAndBefore(i1,2)
        im1 = grp%AfterAndBefore(i0,2)
      else
        print *,' should not be here'
      end if
    
      if(i0.eq.0.and.im1.eq.0) then
        im1 = 1!required?
        i5 = grp%AfterAndBefore(i4,1)

        pr1   = grp%p(i1)        
        ro1   = grp%u(i1,1) 
        rinv1 = 1./grp%u(i1,1) 
        ru1   = grp%u(i1,2)      
        rv1   = grp%u(i1,3)      
        u1    = rinv1*ru1         
        v1    = rinv1*rv1         
        ut1   = u1*any-v1*anx
        un1   = u1*anx+v1*any
        en1   = grp%u(i1,4)*rinv1
        rh1   = grp%u(i1,4) + pr1
        a1    = grp%nodeVolume(i1)
 
        pr2   = grp%p(i2)        
        ro2   = grp%u(i2,1) 
        rinv2 = 1./grp%u(i2,1) 
        ru2   = grp%u(i2,2)      
        rv2   = grp%u(i2,3)      
        u2    = rinv2*ru2         
        v2    = rinv2*rv2         
        ut2   = u2*any-v2*anx
        un2   = u2*anx+v2*any
        en2   = grp%u(i2,4)*rinv2
        rh2   = grp%u(i2,4) + pr2 
        a2    = grp%nodeVolume(i2)

        pr3   = grp%p(i3)        
        ro3   = grp%u(i3,1) 
        rinv3 = 1./grp%u(i3,1) 
        ru3   = grp%u(i3,2)      
        rv3   = grp%u(i3,3)      
        u3    = rinv3*ru3         
        v3    = rinv3*rv3         
        ut3   = u3*any-v3*anx
        un3   = u3*anx+v3*any
        en3   = grp%u(i3,4)*rinv3
        rh3   = grp%u(i3,4) + pr3 
        a3    = grp%nodeVolume(i3)

        pr4   = grp%p(i4)        
        ro4   = grp%u(i4,1) 
        rinv4 = 1./grp%u(i4,1) 
        ru4   = grp%u(i4,2)      
        rv4   = grp%u(i4,3)      
        u4    = rinv4*ru4         
        v4    = rinv4*rv4         
        ut4   = u4*any-v4*anx
        un4   = u4*anx+v4*any
        en4   = grp%u(i4,4)*rinv4
        rh4   = grp%u(i4,4) + pr4 
        a4    = grp%nodeVolume(i4)

         !debugging
        !if (k==1.0) then
        !  a1=0.05
        !  a2=0.1
        !  a3=0.1
        !  a4=0.1
          !x^3-2*x
        !  ro1=-199.0/1000.0 !point value at x=0.1
        !  ro2=-783.0/2000.0 !mean values
        !  ro3=-2289.0/4000.0
        !  ro4=-147.0/200.0
        !endif

        !Approximate Jacobian

        !equation(14)
        ia1=a1
        ia2=a2+ia1
        ia3=a3+ia2
        ia4=a4+ia3

        
        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location
        am=c0*ia1+c1*ia2+c2*ia3+c3*ia4        !equation(15)

        c0 = grp%hoc3(i,1)
        c1 = grp%hoc3(i,2)
        c2 = grp%hoc3(i,3)
        c3 = grp%hoc3(i,4)
        !jacobian at boundary
        amb=c0*ia1+c1*ia2+c2*ia3+c3*ia4        !equation(20)

        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)

        !primitive functions 
    
        !Note:at boundary node point value is used

        !equation(8)
        ro1=ro1*amb
        !equation(9)
        ro2=ro2*a2
        ro3=ro2+ro3*a3
        ro4=ro3+ro4*a4
!same equations for following
        ut1=ut1*amb
        ut2=ut2*a2
        ut3=ut2+ut3*a3
        ut4=ut3+ut4*a4

        un1=un1*amb
        un2=un2*a2
        un3=un2+un3*a3
        un4=un3+un4*a4

        en1=en1*amb
        en2=en2*a2
        en3=en2+en3*a3
        en4=en3+en4*a4
 
       !interpolation 
        !equation(19)
        rom= c0*ro1+c1*ro2+c2*ro3+c3*ro4
        utm= c0*ut1+c1*ut2+c2*ut3+c3*ut4
        unm= c0*un1+c1*un2+c2*un3+c3*un4
        enm= c0*en1+c1*en2+c2*en3+c3*en4

        !divide through by Jacobian
        !equation(17)
        rom=rom/am
        utm=utm/am
        unm=unm/am
        enm=enm/am
        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',rom
        !  k=0.0
          !should be -0.296625
        !endif

        uxm= utm*any+unm*anx
        vym=-utm*anx+unm*any
        prm= (pr1+pr2)/2.

        f1mx = rom*uxm
        f2mx = rom*uxm**2+prm
        f3mx = rom*vym*uxm 
        f4mx = (rom*enm+prm)*uxm 
        f1my = rom*vym
        f2my = rom*vym*uxm
        f3my = rom*vym**2+prm
        f4my = (rom*enm+prm)*vym 
        f1 = wx*(2*f1mx)+wy*(2*f1my)
        f2 = wx*(2*f2mx)+wy*(2*f2my)
        f3 = wx*(2*f3mx)+wy*(2*f3my)
        f4 = wx*(2*f4mx)+wy*(2*f4my)
        !print*,1,i1

      else if(im1.eq.0) then
        pr0   = grp%p(i0)        
        ro0   = grp%u(i0,1) 
        rinv0 = 1./grp%u(i0,1) 
        ru0   = grp%u(i0,2)      
        rv0   = grp%u(i0,3)      
        u0    = rinv0*ru0         
        v0    = rinv0*rv0         
        ut0   = u0*any-v0*anx
        un0   = u0*anx+v0*any
        en0   = grp%u(i0,4)*rinv0
        rh0   = grp%u(i0,4) + pr0 
        a0    = grp%nodeVolume(i0)

        pr1   = grp%p(i1)        
        ro1   = grp%u(i1,1) 
        rinv1 = 1./grp%u(i1,1) 
        ru1   = grp%u(i1,2)      
        rv1   = grp%u(i1,3)      
        u1    = rinv1*ru1         
        v1    = rinv1*rv1         
        ut1   = u1*any-v1*anx
        un1   = u1*anx+v1*any
        en1   = grp%u(i1,4)*rinv1
        rh1   = grp%u(i1,4) + pr1 
        a1    = grp%nodeVolume(i1)

        pr2   = grp%p(i2)        
        ro2   = grp%u(i2,1) 
        rinv2 = 1./grp%u(i2,1) 
        ru2   = grp%u(i2,2)      
        rv2   = grp%u(i2,3)      
        u2    = rinv2*ru2         
        v2    = rinv2*rv2         
        ut2   = u2*any-v2*anx
        un2   = u2*anx+v2*any
        en2   = grp%u(i2,4)*rinv2
        rh2   = grp%u(i2,4) + pr2 
        a2    = grp%nodeVolume(i2)

        pr3   = grp%p(i3)        
        ro3   = grp%u(i3,1) 
        rinv3 = 1./grp%u(i3,1) 
        ru3   = grp%u(i3,2)      
        rv3   = grp%u(i3,3)      
        u3    = rinv3*ru3         
        v3    = rinv3*rv3         
        ut3   = u3*any-v3*anx
        un3   = u3*anx+v3*any
        en3   = grp%u(i3,4)*rinv3
        rh3   = grp%u(i3,4) + pr3 
        a3    = grp%nodeVolume(i3)

        !debugging
        !if (k==1.0) then
        !  a0=0.05
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
        !  ro0=199.0/1000.0
        !  ro1=783.0/2000.0
        !  ro2=2289.0/4000.0
        !  ro3=147.0/200.0
        !end if

        !Approximate Jacobian at interpolation face
        !equation(14)
        ia0=a0
        ia1=a1+ia0
        ia2=a2+ia1
        ia3=a3+ia2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !jacobian at interpolation location        
        am=c0*ia0+c1*ia1+c2*ia2+c3*ia3        !equation(15)

        c0 = grp%hoc3(i,1)
        c1 = grp%hoc3(i,2)
        c2 = grp%hoc3(i,3)
        c3 = grp%hoc3(i,4)
        !jacobian at boundary
        amb=c0*ia0+c1*ia1+c2*ia2+c3*ia3       !equation(20)

        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)

       !primitive functions

        !Note:at boundary node point value is used

        !equation(8)
        ro0=ro0*amb
        !equation(9)i
        ro1=ro1*a1
        ro2=ro1+ro2*a2
        ro3=ro2+ro3*a3
!same equations for following
        ut0=ut0*amb
        ut1=ut1*a1
        ut2=ut1+ut2*a2
        ut3=ut2+ut3*a3

        un0=un0*amb
        un1=un1*a1
        un2=un1+un2*a2
        un3=un2+un3*a3

        en0=en0*amb
        en1=en1*a1
        en2=en1+en2*a2
        en3=en2+en3*a3

        !equation(19)
        rom= c0*ro0+c1*ro1+c2*ro2+c3*ro3
        utm= c0*ut0+c1*ut1+c2*ut2+c3*ut3
        unm= c0*un0+c1*un1+c2*un2+c3*un3
        enm= c0*en0+c1*en1+c2*en2+c3*en3

        !a=sqrt((2.*wx)**2+(2.*wy)**2)
        !rom=rom/a
        !utm=utm/a
        !unm=unm/a
        !enm=enm/a

        !divide through by Jacobian
        !equation(17)
        rom=rom/am
        utm=utm/am
        unm=unm/am
        enm=enm/am
        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',enm
        !  k=0.0
        ! should be -0.484375
        !end if


        uxm= utm*any+unm*anx
        vym=-utm*anx+unm*any
        prm= (pr1+pr2)/2.

        f1mx = rom*uxm
        f2mx = rom*uxm**2+prm
        f3mx = rom*vym*uxm 
        f4mx = (rom*enm+prm)*uxm 
        f1my = rom*vym
        f2my = rom*vym*uxm
        f3my = rom*vym**2+prm
        f4my = (rom*enm+prm)*vym 
        f1 = wx*(2*f1mx)+wy*(2*f1my)
        f2 = wx*(2*f2mx)+wy*(2*f2my)
        f3 = wx*(2*f3mx)+wy*(2*f3my)
        f4 = wx*(2*f4mx)+wy*(2*f4my)
        !print*,1,i0

      else if(i4.eq.0) then
        f1 = wx*(f11x+f21x) + wy*(f11y+f21y)
        f2 = wx*(f12x+f22x) + wy*(f12y+f22y)
        f3 = wx*(f13x+f23x) + wy*(f13y+f23y)
        f4 = wx*(f14x+f24x) + wy*(f14y+f24y)
      else 
        pr0   = grp%p(i0)        
        ro0   = grp%u(i0,1) 
        rinv0 = 1./grp%u(i0,1) 
        ru0   = grp%u(i0,2)      
        rv0   = grp%u(i0,3)      
        u0    = rinv0*ru0         
        v0    = rinv0*rv0         
        ut0   = u0*any-v0*anx
        un0   = u0*anx+v0*any
        en0   = grp%u(i0,4)*rinv0
        rh0   = grp%u(i0,4) + pr0 
        a0    = grp%nodeVolume(i0)
        !print*,a0

        pr1   = grp%p(i1)        
        ro1   = grp%u(i1,1) 
        rinv1 = 1./grp%u(i1,1) 
        ru1   = grp%u(i1,2)      
        rv1   = grp%u(i1,3)      
        u1    = rinv1*ru1         
        v1    = rinv1*rv1         
        ut1   = u1*any-v1*anx
        un1   = u1*anx+v1*any
        en1   = grp%u(i1,4)*rinv1
        rh1   = grp%u(i1,4) + pr1 
        a1    = grp%nodeVolume(i1)

        pr2   = grp%p(i2)        
        ro2   = grp%u(i2,1) 
        rinv2 = 1./grp%u(i2,1) 
        ru2   = grp%u(i2,2)      
        rv2   = grp%u(i2,3)      
        u2    = rinv2*ru2         
        v2    = rinv2*rv2         
        ut2   = u2*any-v2*anx
        un2   = u2*anx+v2*any
        en2   = grp%u(i2,4)*rinv2
        rh2   = grp%u(i2,4) + pr2 
        a2    = grp%nodeVolume(i2)

        pr3   = grp%p(i3)        
        ro3   = grp%u(i3,1) 
        rinv3 = 1./grp%u(i3,1) 
        ru3   = grp%u(i3,2)      
        rv3   = grp%u(i3,3)      
        u3    = rinv3*ru3         
        v3    = rinv3*rv3         
        ut3   = u3*any-v3*anx
        un3   = u3*anx+v3*any
        en3   = grp%u(i3,4)*rinv3
        rh3   = grp%u(i3,4) + pr3 
        a3    = grp%nodeVolume(i3)

        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)

        !primitive functions 

        !debugging
        !if (k==1.0) then
        !  a0=0.1
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
          !ro0=-159.0/800.0
          !ro1=-783.0/2000.0
          !ro2=-2289.0/4000.0
          !ro3=-147.0/200.0
         !endif

        !equation(9)
        ro0=ro0*a0
        ro1=ro0+ro1*a1
        ro2=ro1+ro2*a2
        ro3=ro2+ro3*a3

        ut0=ut0*a0
        ut1=ut0+ut1*a1
        ut2=ut1+ut2*a2
        ut3=ut2+ut3*a3

        un0=un0*a0
        un1=un0+un1*a1
        un2=un1+un2*a2
        un3=un2+un3*a3

        en0=en0*a0
        en1=en0+en1*a1
        en2=en1+en2*a2
        en3=en2+en3*a3

        !equation(10)
        rom= c0*ro0+c1*ro1+c2*ro2+c3*ro3       
        utm= c0*ut0+c1*ut1+c2*ut2+c3*ut3
        unm= c0*un0+c1*un1+c2*un2+c3*un3
        enm= c0*en0+c1*en1+c2*en2+c3*en3

        !Approximate Jacobian
        !equation(14)
        ia0=a0
        ia1=a1+ia0
        ia2=a2+ia1
        ia3=a3+ia2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        
        am=c0*ia0+c1*ia1+c2*ia2+c3*ia3           !equation(15)
        !print*,am
        !a=sqrt((2.*wx)**2+(2.*wy)**2)
        !rom=rom/a
        !utm=utm/a
        !unm=unm/a
        !enm=enm/a

        !divide through by Jacobian
        !equation(17)
        rom=rom/am
        utm=utm/am
        unm=unm/am
        enm=enm/am
        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',enm
        !  k=0.0
          !should give -0.484375
        !endif


        !print*,rom

        uxm= utm*any+unm*anx
        vym=-utm*anx+unm*any
        prm= (pr1+pr2)/2.

        f1mx = rom*uxm
        f2mx = rom*uxm**2+prm
        f3mx = rom*vym*uxm 
        f4mx = (rom*enm+prm)*uxm
        f1my = rom*vym
        f2my = rom*vym*uxm
        f3my = rom*vym**2+prm
        f4my = (rom*enm+prm)*vym

        f1 = wx*(2*f1mx)+wy*(2*f1my)
        f2 = wx*(2*f2mx)+wy*(2*f2my)
        f3 = wx*(2*f3mx)+wy*(2*f3my)
        f4 = wx*(2*f4mx)+wy*(2*f4my)
      end if
    else
      f1 = wx*(f11x+f21x) + wy*(f11y+f21y)
      f2 = wx*(f12x+f22x) + wy*(f12y+f22y)
      f3 = wx*(f13x+f23x) + wy*(f13y+f23y)
      f4 = wx*(f14x+f24x) + wy*(f14y+f24y)
    end if
  else
! side weights from preprocessor 
 
    f1 = wx*(f11x+f21x) + wy*(f11y+f21y)
    f2 = wx*(f12x+f22x) + wy*(f12y+f22y)
    f3 = wx*(f13x+f23x) + wy*(f13y+f23y)
    f4 = wx*(f14x+f24x) + wy*(f14y+f24y)
  end if
! assemble right hand side 
  
! write(698,*) i,ind1,ind2
! write(698,'(4E15.7)') f1,f2,f3,f4

  grp%rhs(ind1,1) = grp%rhs(ind1,1) + f1
  grp%rhs(ind1,2) = grp%rhs(ind1,2) + f2 
  grp%rhs(ind1,3) = grp%rhs(ind1,3) + f3  
  grp%rhs(ind1,4) = grp%rhs(ind1,4) + f4 
  grp%rhs(ind2,1) = grp%rhs(ind2,1) - f1 
  grp%rhs(ind2,2) = grp%rhs(ind2,2) - f2
  grp%rhs(ind2,3) = grp%rhs(ind2,3) - f3 
  grp%rhs(ind2,4) = grp%rhs(ind2,4) - f4 
 end do 

!debugging, print inviscid fluxes
!if (k==1.0) then
!    open(14,file='Iflux.dat')
!    do i=1,grp%numberOfNodes
!        write(14,'(I5,4(E17.8))')i,grp%rhs(i,1),grp%rhs(i,2),grp%rhs(i,3),grp%rhs(i,4)
!    enddo
!    close(14)
!k=2.0
!endif


! face integral
do i=1,grp%brp%numberOfBoundaryFaces
 ind1  = grp%brp%faceIndexArray(i,1)
 ind2  = grp%brp%faceIndexArray(i,2)
 faceIndicator = grp%brp%faceIndexArray(i,3)
 wx  = grp%brp%faceWeightsArray(i,1)
 wy  = grp%brp%faceWeightsArray(i,2)
 rnorm=sqrt(wx*wx+wy*wy)
 rkx = wx/rnorm
 rky = wy/rnorm
 unorm1 = (grp%u(ind1,2)*rkx+grp%u(ind1,3)*rky)/grp%u(ind1,1)
 unorm2 = (grp%u(ind2,2)*rkx+grp%u(ind2,3)*rky)/grp%u(ind2,1)

 if(faceIndicator.ne.3) then
!if(faceIndicator.ne.3.or.(unorm1.ge.0.0.and.unorm2.ge.0.0)) then
!if(unorm1.ge.0.0.and.unorm2.ge.0.0) then 
!if(faceIndicator.ne.3.and.faceIndicator.ne.1) then 
  pr1 = grp%p(ind1) ! pressure
  rinv1  = 1./grp%u(ind1,1)
  ru1 = grp%u(ind1,2)
  rv1 = grp%u(ind1,3)
  rh1 = grp%u(ind1,4) + pr1
  u1  = ru1*rinv1
  v1  = rv1*rinv1

  pr2 = grp%p(ind2)
  rinv2  = 1./grp%u(ind2,1)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rh2 = grp%u(ind2,4) + pr2
  u2  = ru2*rinv2
  v2  = rv2*rinv2

  f1 = wx*(ru1+ru2)               + wy*(rv1+rv2) 
  f2 = wx*(ru1*u1+pr1+ru2*u2+pr2) + wy*(ru1*v1+ru2*v2)  
  f3 = wx*(rv1*u1+rv2*u2)         + wy*(rv1*v1+pr1+rv2*v2+pr2)   
  f4 = wx*(rh1*u1+rh2*u2)         + wy*(rh1*v1+rh2*v2)
  if(ivd%boundaryTerm==1) then 
   grp%rhs(ind1,1) = grp%rhs(ind1,1)+f1+2.*(wx*ru1          + wy*rv1         ) 
   grp%rhs(ind1,2) = grp%rhs(ind1,2)+f2+2.*(wx*(ru1*u1+pr1) + wy*(ru1*v1    )) 
   grp%rhs(ind1,3) = grp%rhs(ind1,3)+f3+2.*(wx*(rv1*u1    ) + wy*(rv1*v1+pr1)) 
   grp%rhs(ind1,4) = grp%rhs(ind1,4)+f4+2.*(wx*rh1*u1       + wy*rh1*v1      ) 

   grp%rhs(ind2,1) = grp%rhs(ind2,1)+f1+2.*(wx*ru2          + wy*rv2         ) 
   grp%rhs(ind2,2) = grp%rhs(ind2,2)+f2+2.*(wx*(ru2*u2+pr2) + wy*(ru2*v2    )) 
   grp%rhs(ind2,3) = grp%rhs(ind2,3)+f3+2.*(wx*(rv2*u2    ) + wy*(rv2*v2+pr2)) 
   grp%rhs(ind2,4) = grp%rhs(ind2,4)+f4+2.*(wx*rh2*u2       + wy*rh2*v2      )
  else
   grp%rhs(ind1,1) = grp%rhs(ind1,1)+4.*(wx*ru1          + wy*rv1         )
   grp%rhs(ind1,2) = grp%rhs(ind1,2)+4.*(wx*(ru1*u1+pr1) + wy*(ru1*v1    ))
   grp%rhs(ind1,3) = grp%rhs(ind1,3)+4.*(wx*(rv1*u1    ) + wy*(rv1*v1+pr1))
   grp%rhs(ind1,4) = grp%rhs(ind1,4)+4.*(wx*rh1*u1       + wy*rh1*v1      )

   grp%rhs(ind2,1) = grp%rhs(ind2,1)+4.*(wx*ru2          + wy*rv2         )
   grp%rhs(ind2,2) = grp%rhs(ind2,2)+4.*(wx*(ru2*u2+pr2) + wy*(ru2*v2    ))
   grp%rhs(ind2,3) = grp%rhs(ind2,3)+4.*(wx*(rv2*u2    ) + wy*(rv2*v2+pr2))
   grp%rhs(ind2,4) = grp%rhs(ind2,4)+4.*(wx*rh2*u2       + wy*rh2*v2      )
  end if
  else if(faceIndicator.eq.3) then ! faceIndicator >2
  rnorm=sqrt(wx*wx+wy*wy)
  if(rnorm>1.0e-10) then
   rnx = wx/rnorm
   rny = wy/rnorm

   ro    = ivd%inflowField(1,1)
   ro1   = 1./ro
   uo    = ivd%inflowField(2,1)
   vo    = ivd%inflowField(3,1)
   eo    = ivd%inflowField(4,1)
   po    = ivd%inflowField(5,1)

   pr1l = grp%p(ind1) ! pressure
   rho1l = grp%u(ind1,1)
   rinv1l  = 1./rho1l
   ru1l = grp%u(ind1,2)
   rv1l = grp%u(ind1,3)
   rh1l = grp%u(ind1,4) + pr1l
   u1l  = ru1l*rinv1l
   v1l  = rv1l*rinv1l
   eps1l = grp%u(ind1,4)
   h1l  = rh1l*rinv1l
   normalVelocity1l = u1l*rnx+v1l*rny

   fl11 = rho1l*normalVelocity1l
   fl21 = pr1l*rnx+rho1l*u1l*normalVelocity1l
   fl31 = pr1l*rny+rho1l*v1l*normalVelocity1l
   fl41 = rh1l*normalVelocity1l

   vml   = u1l*u1l+v1l*v1l
   aml   = sqrt((vml*rho1l)/(ivd%gamma*pr1))
 ! if( normalVelocity1l .gt. 0.0 .and. aml .gt. 1.0 ) then
 !  rhor  = rho1l
 !  uxr   = u1l
 !  vyr   = v1l
 !  epsr  = eps1l
 !  presr = pr1l
 !  hr    = h1l 
 ! else

    rhor  = ro
    uxr   = uo
    vyr   = vo
    epsr  = eo*ro
    presr = po
    hr    = eo + po*ro1
 ! end if
   unr   = uxr*rnx+vyr*rny
   fr11  = rhor*unr
   fr21  = rnx*presr+rhor*uxr*unr
   fr31  = rny*presr+rhor*vyr*unr
   fr41  = (epsr+presr)*unr

   di    = sqrt(rhor/rho1l)
   d1    = 1.0/(di+1.0)
   ui    = (di*uxr+u1l)*d1
   vi    = (di*vyr+v1l)*d1
   hi    = (di*hr+h1l)*d1
   ci2   = gammaMinusOne*(hi-0.5*(ui*ui+vi*vi))
   ci2   = max(ci2,1.0e-5)
   ci    = sqrt(ci2)
   af    = 0.5*(ui*ui+vi*vi)
   ucp   = ui*rnx+vi*rny

   du1   = rhor-rho1l
   du2   = rhor*uxr-rho1l*u1l
   du3   = rhor*vyr-rho1l*v1l
   du4   = epsr-eps1l

   rlam1 = abs(ucp+ci)
   rlam2 = abs(ucp-ci)
   rlam3 = abs(ucp)
   if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
   if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
   if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

   s1    = 0.5*(rlam1+rlam2)
   s2    = 0.5*(rlam1-rlam2)
   al1x  = gammaMinusOne*(af*du1-ui*du2-vi*du3+du4)
   al2x  = -ucp*du1+du2*rnx+du3*rny
   cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
   cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x
   f11   = 0.5*(fr11+fl11)-0.5*(rlam3*du1+cc1           )
   f21   = 0.5*(fr21+fl21)-0.5*(rlam3*du2+cc1*ui+cc2*rnx)
   f31   = 0.5*(fr31+fl31)-0.5*(rlam3*du3+cc1*vi+cc2*rny)
   f41   = 0.5*(fr41+fl41)-0.5*(rlam3*du4+cc1*hi+cc2*ucp)

   pr2l = grp%p(ind2) ! pressure
   rho2l = grp%u(ind2,1)
   rinv2l  = 1./rho2l
   ru2l = grp%u(ind2,2)
   rv2l = grp%u(ind2,3)
   rh2l = grp%u(ind2,4) + pr2l
   u2l = ru2l*rinv2l
   v2l  = rv2l*rinv2l
   h2l  = rh2l*rinv2l
   eps2l = grp%u(ind2,4)
   normalVelocity2l = u2l*rnx+v2l*rny

   fl12 = rho2l*normalVelocity2l
   fl22 = pr2l*rnx+rho2l*u2l*normalVelocity2l
   fl32 = pr2l*rny+rho2l*v2l*normalVelocity2l
   fl42 = rh2l*normalVelocity2l

   vml   = u2l*u2l+v2l*v2l
   aml   = sqrt((vml*rho2l)/(ivd%gamma*pr2))
 ! if( normalVelocity2l .gt. 0.0 .and. aml .gt. 1.0 ) then
 !  rhor  = rho2l
 !  uxr   = u2l
 !  vyr   = v2l
 !  epsr  = eps2l
 !  presr = pr2l
 !  hr    = h2l
 ! else

    rhor  = ro
    uxr   = uo
    vyr   = vo
    epsr  = eo*ro
    presr = po
    hr    = eo + po*ro1
 ! end if

   unr   = uxr*rnx+vyr*rny
   fr12  = rhor*unr
   fr22  = rnx*presr+rhor*uxr*unr
   fr32  = rny*presr+rhor*vyr*unr
   fr42  = (epsr+presr)*unr

   di    = sqrt(rhor/rho2l)
   d1    = 1.0/(di+1.0)
   ui    = (di*uxr+u2l)*d1
   vi    = (di*vyr+v2l)*d1
   hi    = (di*hr+h2l)*d1
   ci2   = gammaMinusOne*(hi-0.5*(ui*ui+vi*vi))
   ci2   = max(ci2,1.0e-5)
   ci    = sqrt(ci2)
   af    = 0.5*(ui*ui+vi*vi)
   ucp   = ui*rnx+vi*rny

   du1   = rhor-rho2l
   du2   = rhor*uxr-rho2l*u2l
   du3   = rhor*vyr-rho2l*v2l
   du4   = epsr-eps2l

   rlam1 = abs(ucp+ci)
   rlam2 = abs(ucp-ci)
   rlam3 = abs(ucp)
   if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
   if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
   if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

   s1    = 0.5*(rlam1+rlam2)
   s2    = 0.5*(rlam1-rlam2)
   al1x  = gammaMinusOne*(af*du1-ui*du2-vi*du3+du4)
   al2x  = -ucp*du1+du2*rnx+du3*rny
   cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
   cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x
   f12   = 0.5*(fr12+fl12)-0.5*(rlam3*du1+cc1           )
   f22   = 0.5*(fr22+fl22)-0.5*(rlam3*du2+cc1*ui+cc2*rnx)
   f32   = 0.5*(fr32+fl32)-0.5*(rlam3*du3+cc1*vi+cc2*rny)
   f42   = 0.5*(fr42+fl42)-0.5*(rlam3*du4+cc1*hi+cc2*ucp)


   if(ivd%boundaryTerm==1) then 
    grp%rhs(ind1,1) = grp%rhs(ind1,1)+rnorm*(3.*f11+f12)
    grp%rhs(ind1,2) = grp%rhs(ind1,2)+rnorm*(3.*f21+f22)
    grp%rhs(ind1,3) = grp%rhs(ind1,3)+rnorm*(3.*f31+f32)
    grp%rhs(ind1,4) = grp%rhs(ind1,4)+rnorm*(3.*f41+f42)

    grp%rhs(ind2,1) = grp%rhs(ind2,1)+rnorm*(f11+3.*f12)
    grp%rhs(ind2,2) = grp%rhs(ind2,2)+rnorm*(f21+3.*f22) 
    grp%rhs(ind2,3) = grp%rhs(ind2,3)+rnorm*(f31+3.*f32)
    grp%rhs(ind2,4) = grp%rhs(ind2,4)+rnorm*(f41+3.*f42)
   else
    grp%rhs(ind1,1) = grp%rhs(ind1,1)+rnorm*4.*f11
    grp%rhs(ind1,2) = grp%rhs(ind1,2)+rnorm*4.*f21
    grp%rhs(ind1,3) = grp%rhs(ind1,3)+rnorm*4.*f31
    grp%rhs(ind1,4) = grp%rhs(ind1,4)+rnorm*4.*f41

    grp%rhs(ind2,1) = grp%rhs(ind2,1)+rnorm*4.*f12
    grp%rhs(ind2,2) = grp%rhs(ind2,2)+rnorm*4.*f22
    grp%rhs(ind2,3) = grp%rhs(ind2,3)+rnorm*4.*f32
    grp%rhs(ind2,4) = grp%rhs(ind2,4)+rnorm*4.*f42
   end if
  end if
 end if
end do
end subroutine makeRHS
!-----------------------------------------------------------------------
subroutine setUpInitialField(ivd,grp,INFILE)
IMPLICIT NONE

type(InputVariablesData) :: ivd
type(GridSolverData) :: grp
integer :: INFILE
double precision :: xu(5)
integer :: i,j,buffi

if(INFILE>0) then 
! start up from given field
 write(*,*) "starting up from given field"
 do i=1,grp%numberOfNodes
  read(INFILE,*) buffi,grp%u(i,1:4),grp%p(i) 
  grp%u(i,2:4) = grp%u(i,1)*grp%u(i,2:4)
 end do

else
! start up from freestream
 write(*,*) "starting up from freestream"
 do i=1,grp%numberOfNodes
  grp%u(i,1) = ivd%inflowField(1,1)
  grp%u(i,2) = ivd%inflowField(2,1)
  grp%u(i,3) = ivd%inflowField(3,1)
  grp%u(i,4) = ivd%inflowField(4,1)
 end do

 if(ivd%turbulenceModel==1) then
  grp%u(:,5) = 0.1
 else
  grp%u(:,5) = 0.0
 end if

 call triggerTurbulenceField(grp,ivd)
end if
if(ivd%ReynoldsNumber==0) then 
 write(*,*) "setting viscosity to zero..."
 grp%laminarViscosity = 0.0
else
 grp%laminarViscosity = 1.0/ivd%ReynoldsNumber 
end if

call setBCsOnSolutionField(grp,ivd,grp%u)
end subroutine setUpInitialField
!-----------------------------------------------------------------------
subroutine triggerTurbulenceField(grp,ivd)
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

integer :: i,j,ind
real :: dt


 if(ivd%triggerTurbulence) then
  do i=1,grp%numberOfTripNodes
   do j=1,grp%tripNodeFieldIndexes(i,0)
    ind = grp%tripNodeFieldIndexes(i,j)
    dt = grp%tripNodeFieldDistances(i,j)
    if(dt<ivd%triggerRadius) then
     grp%u(ind,5) = ivd%triggerValue
    end if
   end do
  end do
 end if

end subroutine triggerTurbulenceField
!-----------------------------------------------------------------------
subroutine setRHSAtBoundary(grp,ivd)
IMPLICIT NONE

! nasty subroutine to fix right hand side at boundaries

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

integer :: i,i1,i2,faceIndicator
real :: wx,wy,rkx,rky,rnorm,pr1,r1,ru1,rv1,u1,v1,rek1,red1,pr2,r2,ru2,rv2,u2,v2
real :: rek2,red2,f1,f2,f3,f4,f11,f21,f31,f41,f12,f22,f32,f42,fl1l,fl2l,fl3l,fl4l
real :: ro,ro1,uo,vo,po,f11l,f12l,f13l,f14l,f21l,f22l,f23l,f24l,unorm1,unorm2
real :: rnx,rny,rhol,rho1,rh1,rh2,eo,ho,arf,vmo1,uxl,vyl,epsl,presl,hl,unl,rhor
real :: uxr,vyr,epsr,presr,hr,unr,fl1r,fl2r,fl3r,fl4r,ekr,edr,z

do i=1,grp%brp%numberOfBoundaryFaces
 i1 = grp%brp%faceIndexArray(i,1)
 i2 = grp%brp%faceIndexArray(i,2) 
 faceIndicator = grp%brp%faceIndexArray(i,3)
 wx  = grp%brp%faceWeightsArray(i,1)
 wy  = grp%brp%faceWeightsArray(i,2)
 rnorm=sqrt(wx*wx+wy*wy)
 rkx = wx/rnorm
 rky = wy/rnorm
 unorm1 = (grp%u(i1,2)*rkx+grp%u(i1,3)*rky)/grp%u(i1,1)
 unorm2 = (grp%u(i2,2)*rkx+grp%u(i2,3)*rky)/grp%u(i2,1)

! Inflow, solve a Riemann problem for freestream and boundary

 !if(faceIndicator.eq.3.and.(unorm1.lt.0.0.or.unorm2.lt.0.0)) then
if(faceIndicator.eq.3) then
  ro    = ivd%inflowField(1,1) 
  ro1   = 1./ro
  uo    = ivd%inflowField(2,1) 
  vo    = ivd%inflowField(3,1) 
  eo    = ivd%inflowField(4,1) 
  po    = ivd%inflowField(5,1)
  ho    = eo + po*ro1 

  rnx   = wx
  rny   = wy
  arf   = sqrt(rnx*rnx+rny*rny)
  vmo1  = 1./arf
  rnx   = rnx*vmo1
  rny   = rny*vmo1

  rhol  = grp%u(i1,1)
  rho1  = 1./rhol
  uxl   = grp%u(i1,2)*rho1
  vyl   = grp%u(i1,3)*rho1
  epsl  = grp%u(i1,4)
  presl = grp%p(i1)
  hl    = (epsl + presl)*rho1
  unl   = uxl*rnx+vyl*rny

  fl1l  = rhol*unl
  fl2l  = rnx*presl+rhol*uxl*unl
  fl3l  = rny*presl+rhol*vyl*unl
  fl4l  = (epsl+presl)*unl

  f11l  = fl1l
  f12l  = fl2l
  f13l  = fl3l
  f14l  = fl4l

  rhor  = ro
  uxr   = uo
  vyr   = vo
  epsr  = eo*ro
  presr = po
  hr    = eo + po*ro1
  unr   = uxr*rnx+vyr*rny

  fl1r  = rhor*unr
  fl2r  = rnx*presr+rhor*uxr*unr
  fl3r  = rny*presr+rhor*vyr*unr
  fl4r  = (epsr+presr)*unr

  z = 1.0 ! dummy
  call calculateRoeMatrix(i1,rnx,rny,rhor,rhol,uxr,uxl,vyr,vyl,hr,hl,&
            epsr,epsl,ivd%gamma,z,grp)

  f11   = 0.5*(fl1l+fl1r)-0.5*grp%nodeHelpArray(i1,1)
  f21   = 0.5*(fl2l+fl2r)-0.5*grp%nodeHelpArray(i1,2)
  f31   = 0.5*(fl3l+fl3r)-0.5*grp%nodeHelpArray(i1,3)
  f41   = 0.5*(fl4l+fl4r)-0.5*grp%nodeHelpArray(i1,4)

  rhol  = grp%u(i2,1)
  rho1  = 1./rhol
  uxl   = grp%u(i2,2)*rho1
  vyl   = grp%u(i2,3)*rho1
  epsl  = grp%u(i2,4)
  presl = grp%p(i2)
  hl    = (epsl + presl)*rho1
  unl   = uxl*rnx+vyl*rny

  fl1l  = rhol*unl
  fl2l  = rnx*presl+rhol*uxl*unl
  fl3l  = rny*presl+rhol*vyl*unl
  fl4l  = (epsl+presl)*unl

  f21l  = fl1l
  f22l  = fl2l
  f23l  = fl3l
  f24l  = fl4l

  rhor  = ro
  uxr   = uo
  vyr   = vo
  epsr  = eo*ro
  presr = po
  hr    = eo + po*ro1
  unr   = uxr*rnx+vyr*rny

  fl1r  = rhor*unr
  fl2r  = rnx*presr+rhor*uxr*unr
  fl3r  = rny*presr+rhor*vyr*unr
  fl4r  = (epsr+presr)*unr

  z = 1.0 ! dummy
  call calculateRoeMatrix(i2,rnx,rny,rhor,rhol,uxr,uxl,vyr,vyl,hr,hl,&
            epsr,epsl,ivd%gamma,z,grp)

  f12   = 0.5*(fl1l+fl1r)-0.5*grp%nodeHelpArray(i2,1)
  f22   = 0.5*(fl2l+fl2r)-0.5*grp%nodeHelpArray(i2,2)
  f32   = 0.5*(fl3l+fl3r)-0.5*grp%nodeHelpArray(i2,3)
  f42   = 0.5*(fl4l+fl4r)-0.5*grp%nodeHelpArray(i2,4)

  grp%rhs(i1,1) = grp%rhs(i1,1) + arf*(4.*f11 + 2.*f12 + f11l - f21l)
  grp%rhs(i1,2) = grp%rhs(i1,2) + arf*(4.*f21 + 2.*f22 + f12l - f22l)
  grp%rhs(i1,3) = grp%rhs(i1,3) + arf*(4.*f31 + 2.*f32 + f13l - f23l)
  grp%rhs(i1,4) = grp%rhs(i1,4) + arf*(4.*f41 + 2.*f42 + f14l - f24l)

  grp%rhs(i2,1) = grp%rhs(i2,1) + arf*(4.*f12 + 2.*f11 + f21l - f11l)
  grp%rhs(i2,2) = grp%rhs(i2,2) + arf*(4.*f22 + 2.*f21 + f22l - f12l)
  grp%rhs(i2,3) = grp%rhs(i2,3) + arf*(4.*f32 + 2.*f31 + f23l - f13l)
  grp%rhs(i2,4) = grp%rhs(i2,4) + arf*(4.*f42 + 2.*f41 + f24l - f14l)
 endif
end do
end subroutine setRHSAtBoundary
!---------------------------------------------------------------------------
subroutine setOuterBoundaryConditions(grp,ivd,u,p)

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:),p(:)

integer :: i,j,k,i1,i2,ies,ien,ip,nbfar,ind,ibf,ist,ib,n1 
real :: gamm1,gam11,c1,c2,c3,c4,c5,c6,p1,alpha,bt,co,cv,eo,ho,cs
real :: sn,fv,uv,vv,fc,dn,pr,uo,vo,ss,ul,uf,cf,fn,ft,sc,uc,vc,cc,ri,ri1,vi,ui,ei,pi
real :: ue,sum1

gamm1 = ivd%gamma - 1.
gam11 = 1./gamm1

nbfar = grp%brp%faceIndicatorArray(-1) ! number of farfield boundaries

grp%nodeHelpArray(:,9:14) = 0.0 ! initialize

! farfield boundaries
 ies = grp%brp%faceIndicatorArray(-16)
 ien = grp%brp%faceIndicatorArray(-17+2*nbfar)
 do k=ies,ien
  if(associated(grp%boundaryFaceNodeMappings)) then 
  ! agglomerated grid, need to map unknowns on boundary from 
  ! fine to coarse numbering systems
   ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(k))
  else
   ip = grp%brp%faceIndicatorArray(k) 
  end if
  grp%nodeHelpArray(ip,14) = 1.
 end do 
 do j=1,grp%numberOfSides
  i1  = grp%sideIndexArray(j,1)
  i2  = grp%sideIndexArray(j,2)
! if one side node is on farfield boundary and the other is not
! set value of boundary node equal to average of values connected to it
 sum1 = grp%nodeHelpArray(i1,14)+grp%nodeHelpArray(i2,14)
 if((sum1<1.5).and.(sum1>0.5)) then 
   grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + u(i2,1)
   grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + u(i2,2)
   grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + u(i2,3)
   grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + u(i2,4)
   grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + u(i1,1)
   grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + u(i1,2)
   grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + u(i1,3)
   grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + u(i1,4) 
   grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + 1.
   grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) + 1.
  end if
 end do

if(.false.) then 
do k=ies,ien
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated grid, need to map unknowns on boundary from 
 ! fine to coarse numbering systems
  ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(k))
 else
  ip = grp%brp%faceIndicatorArray(k) 
 end if 
 if(grp%nodeHelpArray(ip,13)>0.5) then 
  write(*,*) "OB1: ",ip,grp%nodeHelpArray(ip,12),&
                      grp%nodeHelpArray(ip,13),u(ip,4)
  u(ip,1) = grp%nodeHelpArray(ip,9)/grp%nodeHelpArray(ip,13)
  u(ip,2) = grp%nodeHelpArray(ip,10)/grp%nodeHelpArray(ip,13)
  u(ip,3) = grp%nodeHelpArray(ip,11)/grp%nodeHelpArray(ip,13)
  u(ip,4) = grp%nodeHelpArray(ip,12)/grp%nodeHelpArray(ip,13)
 end if
end do
end if 

ind = 17
! run through farfield boundaries
do ibf = 1,nbfar
 ist = grp%brp%faceIndicatorArray(-ind+1)
 ien = grp%brp%faceIndicatorArray(-ind+2)
 ind = ind - 2
 if(ien-ist>0) then 
  if(ibf.le.2) then ! only allow two outflow boundaries
   ri  = ivd%inflowField(1,ibf) 
   ri1 = 1./ri
   ui  = ivd%inflowField(2,ibf) 
   vi  = ivd%inflowField(3,ibf) 
   ei  = ivd%inflowField(4,ibf) 
   pi  = ivd%inflowField(5,ibf)
  end if
  do ib = ist,ien
   if(associated(grp%boundaryFaceNodeMappings)) then 
   ! agglomerated grid, need to map unknowns on boundary from 
   ! fine to coarse numbering systems
    n1 = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(ib))
   else
    n1 = grp%brp%faceIndicatorArray(ib) 
   end if 
   ip = grp%brp%faceIndicatorArray(ib)
   dn   =   -grp%brp%faceTangentArray(ip,1)*u(n1,3)+&
     grp%brp%faceTangentArray(ip,2)*u(n1,2)
   dn   = dn/u(n1,1)
   if( dn.lt. -0.05 ) then
!   if(.true.) then 
    ! inflow
    u(n1,1) = ri
    u(n1,2) = ri*ui
    u(n1,3) = ri*vi
    if(ivd%MachNumber.lt.1.0) then
     u(n1,4) = p(n1)/gamm1+0.5*ri*(ui*ui+vi*vi)
    else
     u(n1,4) = ri*ei
    end if
    if(ivd%MachNumber.ge.1.0) p(n1)= pi ! if supersonic inflow, set pressure to
                                        ! inflow variable
   else ! outflow
    dn  = 0.5*(u(n1,2)**2+u(n1,3)**2)/u(n1,1)
    if(ivd%MachNumber.lt.1.0) then
     p(n1) = pi
     u(n1,4) = pi / gamm1 + dn
    else
     p(n1) = ( ivd%gamma - 1.0 ) * ( u(n1,4) - dn )
    end if
   end if
  end do
 end if
end do
end subroutine setOuterBoundaryConditions
!-------------------------------------------------------------------------
subroutine setOuterBoundaryConditions3(grp,ivd,u,p,changePressure)

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:),p(:)
logical :: changePressure
 
integer :: ist,ien,ib,ip,coarseNodeNumber,i1,i2,j,k,ri
real :: gammaMinusOne,inpRho,inpU1,inpU2,inpE,inpPressure,nx,ny,normalVelocity
real :: rho,U1,U2,E,ww,s,t,xs,xc,x1,x2,ys,yc,y1,y2,sx,sy,rx,ry
logical :: intersectionIsValid,cantFind,doPlot,solutionFound
real :: A,B,C,radicand,v1,v2,t1,t2

gammaMinusOne = ivd%gamma-1.0
ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)
if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib =ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  grp%nodeHelpArray(coarseNodeNumber,10:12) = 0.0
 end do

 do ib =ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  ww = grp%brp%faceWeightNorms(ip)
  grp%nodeHelpArray(coarseNodeNumber,11) = grp%nodeHelpArray(coarseNodeNumber,11) +&
                                            ww*grp%brp%faceTangentArray(ip,2)
  grp%nodeHelpArray(coarseNodeNumber,12) = grp%nodeHelpArray(coarseNodeNumber,12) -&
                                            ww*grp%brp%faceTangentArray(ip,1)
  grp%nodeHelpArray(coarseNodeNumber,10) = grp%nodeHelpArray(coarseNodeNumber,10) + ww
 end do
end if

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  nx = grp%nodeHelpArray(coarseNodeNumber,11)/grp%nodeHelpArray(coarseNodeNumber,10) 
  ny = grp%nodeHelpArray(coarseNodeNumber,12)/grp%nodeHelpArray(coarseNodeNumber,10) 
 else
  coarseNodeNumber = ip
  nx = grp%brp%faceTangentArray(ip,2) 
  ny = -1.0*grp%brp%faceTangentArray(ip,1) 
 end if

if(grp%gridNumber.ne.1) then
 inpRho = ivd%inflowField(1,1)
 inpU1 = ivd%inflowField(2,1)
 inpU2 = ivd%inflowField(3,1)
 inpE = ivd%inflowField(4,1)
! inpPressure = gammaMinusOne*(inpE-0.5*(inpU1**2+inpU2**2)/inpRho)
 inpPressure = ivd%inflowField(5,1)
 rho = u(coarseNodeNumber,1) 
 U1 = u(coarseNodeNumber,2) 
 U2 = u(coarseNodeNumber,3) 
! E = u(coarseNodeNumber,4)
! pressure = gammaMinusOne*(E-0.5*(U1**2+U2**2)/rho)
 normalVelocity = u(coarseNodeNumber,2)*nx + u(coarseNodeNumber,3)*ny
! normalVelocity = inpU1*nx + inpU2*ny
! if(normalVelocity < 0.05.or.(ip==1.or.ip==304)) then 
! if(normalVelocity < 0.05.or.ip==1) thena ! this was used
 if(normalVelocity < 0.001) then 
 !if(.true.) then  ! CHANGE 
 ! inflow
  if(ivd%MachNumber>1.0) then 
   ! supersonic inflow
   u(coarseNodeNumber,1) = inpRho 
   u(coarseNodeNumber,2) = inpU1
   u(coarseNodeNumber,3) = inpU2
   u(coarseNodeNumber,4) = inpE
   if(changePressure) then 
    grp%p(coarseNodeNumber) = inpPressure 
   end if
   if(ivd%turbulenceModel==1) then
    u(coarseNodeNumber,5) = 0.1
   else
    u(coarseNodeNumber,5) = 0.0
    end if
  else
  ! subsonic inflow  - assume inviscid at farfield
   u(coarseNodeNumber,1) = inpRho
   u(coarseNodeNumber,2) = inpU1
   u(coarseNodeNumber,3) = inpU2
!   u(coarseNodeNumber,4) = inpE 
   u(coarseNodeNumber,4) = p(coarseNodeNumber)/gammaMinusOne&
                           +0.5*inpRho*(inpU1*inpU1+inpU2*inpU2)
   if(ivd%turbulenceModel==1) then 
    u(coarseNodeNumber,5) = 0.1
   else
    u(coarseNodeNumber,5) = 0.0
   end if
!   if(changePressure) then 
!    grp%p(coarseNodeNumber) = gammaMinusOne*(inpE-0.5*(U1**2+U2**2)/rho) 
!   end if
  end if
 else
 ! outflow 
  if(ivd%MachNumber<1.0) then 
  ! subsonic outflow
   !u(coarseNodeNumber,4) = inpPressure/gammaMinusOne + 0.5*(U1**2+U2**2)/rho
   !u(coarseNodeNumber,5) = 0.1
   u(coarseNodeNumber,1) = inpRho
   u(coarseNodeNumber,2) = inpU1
   u(coarseNodeNumber,3) = inpU2
   u(coarseNodeNumber,4) = p(coarseNodeNumber)/gammaMinusOne&
                           +0.5*inpRho*(inpU1*inpU1+inpU2*inpU2)
   if(ivd%turbulenceModel==1) then 
    u(coarseNodeNumber,5) = 0.1
   else
    u(coarseNodeNumber,5) = 0.0
   end if
  else
   if(ivd%turbulenceModel==1) then 
    u(coarseNodeNumber,5) = 0.1
   else
    u(coarseNodeNumber,5) = 0.0
   end if
  end if
! for supersonic outflow
 end if
else
 if(ivd%turbulenceModel==1) then 
  u(coarseNodeNumber,5) = 0.1
 else
  u(coarseNodeNumber,5) = 0.0
 end if 
end if
end do

! IO BC's


ist = grp%brp%faceIndicatorArray(-8)
ien = grp%brp%faceIndicatorArray(-7)


if(.not.associated(grp%boundaryFaceNodeMappings)) then 
do ib=ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 rho = u(ip,1) 
! U1 = u(ip,2) 
! U2 = u(ip,3) 
  U1 = 1.0
  U2 = 0.0

! localMachNumber = sqrt((U1*U1+U2*U2)/(ivd%gamma*p(ip))) 
 nx = grp%brp%faceTangentArray(ip,2) 
 ny = -1.0*grp%brp%faceTangentArray(ip,1) 

  doPlot = .false.
 ! first check if boundary really is outflow
 if(U1*nx+U2*ny<0) then  
  write(*,*) "Warning: Node ",ip," is specified as IO but is inflow" 
!  u(ip,1) = ivd%inflowField(1,1)
!  u(ip,2) = ivd%inflowField(2,1)
!  u(ip,3) = ivd%inflowField(3,1)
!  u(ip,4) = ivd%inflowField(4,1)
  doPlot = .true.
!  U1 = nx
!  U2 = ny
  U1 = -U1
  U2 = -U2
 end if 

 cantFind = .false.
55 continue

!  do j=1,grp%brp%numberOfIOBoundaries
!   if((grp%brp%IOIntervals(j,1)<ip.and.grp%brp%IOIntervals(j,2)>ip).or.&
!      (grp%brp%IOIntervals(j,1)>ip.and.grp%brp%IOIntervals(j,2)<ip)) exit
!  end do
  ! j now contains which chain ip should be interpolated to

  ! the object is now to find where the line emanating from ip in the upstream 
  ! direction crosses the IO chain line. The flow variables at this point (xc,yc)
  ! are copied to IO boundary. 

  sx = -U1
  sy = -U2
  xs = grp%brp%IOBoundaryCoordinates(ip,1)
  ys = grp%brp%IOBoundaryCoordinates(ip,2)
  intersectionIsValid = .false.
  do j=1,grp%brp%numberOfIOBoundaries
   do k=grp%brp%IOIntervals(j,3),grp%brp%IOIntervals(j,4)-1
    x1 = grp%brp%IOChainCoordinates(k,1)
    y1 = grp%brp%IOChainCoordinates(k,2)
    x2 = grp%brp%IOChainCoordinates(k+1,1)
    y2 = grp%brp%IOChainCoordinates(k+1,2)
    rx = x2-x1
    ry = y2-y1
    if(abs(ry*sx-rx*sy)>1.0e-16) then
     s = ((ys-y1)*rx-(xs-x1)*ry)/(ry*sx-rx*sy)
     xc = sx*s + xs
     yc = sy*s + ys
     if(s>0.0) then 
      if(x1.le.xc) then 
       if(x2.ge.xc) then 
        if(y1.le.yc) then 
         if(y2.ge.yc) then 
         ! BINGO - have found a valid intersection
          intersectionIsValid = .true.
          exit
         end if
        else if(y1.ge.yc) then 
         if(y2.le.yc) then 
         ! BINGO - have found a valid intersection
          intersectionIsValid = .true.
          exit
         end if
        end if
       end if
      else if(x1.ge.xc) then 
       if(x2.le.xc) then 
        if(y1.le.yc) then 
         if(y2.ge.yc) then 
         ! BINGO - have found a valid intersection
          intersectionIsValid = .true.
          exit
         end if
        else if(y1.ge.yc) then 
         if(y2.le.yc) then 
         ! BINGO - have found a valid intersection
          intersectionIsValid = .true.
          exit
         end if
        end if
       end if
      end if 
     end if 
    end if
   end do 
   if(intersectionIsValid) exit
  end do
  if(intersectionIsValid) then 
  ! find value to translate
   if(abs(rx)>abs(ry)) then 
    ri = 1
   else
    ri = 2
   end if 
   if(ri==1) then 
    t = (sx*s+xs-x1)/rx 
   else
    t = (sy*s+ys-y1)/ry
   end if
   if(abs(x2-x1)>abs(y2-y1)) then 
    ri = 1
   else
    ri = 2
   end if
   if(ri==1) then 
    t = rx*t/(x2-x1)
   else
    t = ry*t/(y2-y1)
   end if 
   ! for any variable Z we now have that Z(xc,yc)=(1-t)*Z(x1,y1)+t*Z(x2,y2)
   i1 = grp%brp%IOBoundaryRegister(k)
   i2 = grp%brp%IOBoundaryRegister(k+1)
   grp%u(ip,:) = (1.0-t)*grp%u(i1,:)+t*grp%u(i2,:)
   p(ip) = (1.0-t)*p(i1)+t*p(i2)
  else
   write(*,*) "Warning: Couldn't find a valid intersection for IO boundary at node ",ip
   if(.not.cantFind) then 
    cantFind = .true.
    U1 = nx
    U2 = ny
    goto 55
   end if
  end if
! end if
end do
end if

end subroutine setOuterBoundaryConditions3
!-------------------------------------------------------------------------
subroutine setOuterBoundaryConditions4(grp,ivd,delu)
! nullifies coarse mesh increments at farfield

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: delu(:,:)

integer :: ist,ien,ib,ip

ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)

do ib =ist,ien
 if(associated(grp%boundaryFaceNodeMappings)) then 
  ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(ib))
 else
  ip = grp%brp%faceIndicatorArray(ib)
 end if
 delu(ip,:) = 0.0
end do 

ist = grp%brp%faceIndicatorArray(-8)
ien = grp%brp%faceIndicatorArray(-7)

do ib =ist,ien
 if(associated(grp%boundaryFaceNodeMappings)) then 
  ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(ib))
 else
  ip = grp%brp%faceIndicatorArray(ib)
 end if
 delu(ip,:) = 0.0
end do 

ist = grp%brp%faceIndicatorArray(-20)
ien = grp%brp%faceIndicatorArray(-19)

do ib =ist,ien
 if(associated(grp%boundaryFaceNodeMappings)) then 
  ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(ib))
 else
  ip = grp%brp%faceIndicatorArray(ib)
 end if
 delu(ip,:) = 0.0
end do 

end subroutine setOuterBoundaryConditions4
!-------------------------------------------------------------------------
subroutine nullifyIOChainIncrement(grp,delu)
IMPLICIT NONE

type(GridSolverData) :: grp
real :: delu(:,:)

integer :: ip,ib

do ib=1,grp%brp%numberOfIONodes
 ip = grp%brp%IOBoundaryRegister(ib)
 delu(ip,:) = 0.0
end do

end subroutine nullifyIOChainIncrement
!-------------------------------------------------------------------------
subroutine setOuterBoundaryConditions2(grp,ivd,u,p)

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:),p(:)

integer :: ist,ien,ib,ip,coarseNodeNumber
real :: gammaMinusOne

gammaMinusOne = ivd%gamma-1.0
ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 u(coarseNodeNumber,1) = ivd%inflowField(1,1)
 u(coarseNodeNumber,2) = ivd%inflowField(2,1)
 u(coarseNodeNumber,3) = ivd%inflowField(3,1)
 u(coarseNodeNumber,4) = ivd%inflowField(4,1)
 grp%p(coarseNodeNumber) = gammaMinusOne*(u(coarseNodeNumber,4)&
                -0.5*(u(coarseNodeNumber,2)**2+u(coarseNodeNumber,3)**2)&
                /u(coarseNodeNumber,1)) 
end do

end subroutine setOuterBoundaryConditions2
!-------------------------------------------------------------------------
subroutine nullifyOuterBoundary(grp,u)

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:)

integer :: ist,ien,ib,ip,coarseNodeNumber
ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 u(coarseNodeNumber,1) = 0.0 
 u(coarseNodeNumber,2) = 0.0 
 u(coarseNodeNumber,3) = 0.0 
 u(coarseNodeNumber,4) = 0.0 
end do
end subroutine nullifyOuterBoundary
!-------------------------------------------------------------------------
subroutine setOldBCs(grp,ivd,u,delu)
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:),delu(:,:)

integer :: ist,ien,ip,ib
real :: tx,ty,nx,ny,un

ist = grp%brp%faceIndicatorArray(-10)
ien = grp%brp%faceIndicatorArray(-9)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 delu(ip,2) = u(ip,2)
 delu(ip,3) = u(ip,3)
end do

do ib=ien+1,grp%brp%numberOfBoundaryFaces
 ip = grp%brp%faceIndicatorArray(ib)
 delu(ip,2) = u(ip,2)
 delu(ip,3) = u(ip,3)
end do
end subroutine setOldBCs
!-------------------------------------------------------------------------
subroutine setBCsOnSolutionField(grp,ivd,u)
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:)

integer :: ist,ien,ib,ip,coarseNodeNumber,nbf,ind1,ind2,i
real :: tx,ty,nx,ny,normalVelocity,tangentialVelocity,ww,vectorNorm
real :: storeMassFactor,engineMassFlows,area,ptfunc
real :: wx(2),u1(2),u2(2),massflowfunction,engineMassFactor
real :: press1,press2,press,temp1,temp2,temp
real :: normal(2),intE,factV(2),factor

!!!!! Inviscid wall !!!!!!!!!!!!
ist = grp%brp%faceIndicatorArray(-20)
ien = grp%brp%faceIndicatorArray(-19)

if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib = ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  grp%nodeHelpArray(coarseNodeNumber,11:12) = 0.0
 end do
end if

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 tx = -1.0*grp%brp%faceTangentArray(ip,1) 
 ty = -1.0*grp%brp%faceTangentArray(ip,2) 
 nx = ty
 ny = -1.0*tx
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  ww = grp%brp%faceWeightNorms(ip)
 else
  coarseNodeNumber = ip
 end if
 normalVelocity = u(coarseNodeNumber,2)*nx + u(coarseNodeNumber,3)*ny
 normalVelocity = (1.0-ivd%normalComponentRelaxation)*normalVelocity
 normalVelocity = 0.0
 tangentialVelocity = u(coarseNodeNumber,2)*tx + u(coarseNodeNumber,3)*ty  
 if(associated(grp%boundaryFaceNodeMappings)) then 
  grp%nodeHelpArray(coarseNodeNumber,11) = grp%nodeHelpArray(coarseNodeNumber,11) + ww*tx 
  grp%nodeHelpArray(coarseNodeNumber,12) = grp%nodeHelpArray(coarseNodeNumber,12) + ww*ty 
 else
 ! fine grid
  u(coarseNodeNumber,2) = tangentialVelocity*tx+normalVelocity*nx
  u(coarseNodeNumber,3) = tangentialVelocity*ty+normalVelocity*ny
 end if
end do

if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib = ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  tangentialVelocity = u(coarseNodeNumber,2)*grp%nodeHelpArray(coarseNodeNumber,11) +&
                       u(coarseNodeNumber,3)*grp%nodeHelpArray(coarseNodeNumber,12)
  vectorNorm = grp%nodeHelpArray(coarseNodeNumber,11)**2+grp%nodeHelpArray(coarseNodeNumber,12)**2
  tangentialVelocity = tangentialVelocity/vectorNorm
  u(coarseNodeNumber,2) = tangentialVelocity*grp%nodeHelpArray(coarseNodeNumber,11)
  u(coarseNodeNumber,3) = tangentialVelocity*grp%nodeHelpArray(coarseNodeNumber,12)
 end do
end if


!!! Isothermal viscous wall   !!!!!!!!!!!!!!!!!

ist = grp%brp%faceIndicatorArray(-12)
ien = grp%brp%faceIndicatorArray(-11)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 u(coarseNodeNumber,2:3) = 0.0
 u(coarseNodeNumber,4) = u(coarseNodeNumber,1)*ivd%wallTemperature/ivd%gamma 
 u(coarseNodeNumber,5) = 0.0
 grp%p(coarseNodeNumber) = (ivd%gamma/(ivd%gamma-1.0))*u(coarseNodeNumber,4)/u(coarseNodeNumber,1)
end do


!!!!! Adiabatic viscous wall  !!!!!!!!!!!!

ist = grp%brp%faceIndicatorArray(-10)
ien = grp%brp%faceIndicatorArray(-9)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 u(coarseNodeNumber,2:3) = 0.0

 u(coarseNodeNumber,5) = 0.0

 if(.not.ivd%wallsAreAdiabatic) then 
  u(coarseNodeNumber,4) = u(coarseNodeNumber,1)*ivd%wallTemperature/ivd%gamma
  grp%p(coarseNodeNumber) = (ivd%gamma/(ivd%gamma-1.0))*u(coarseNodeNumber,4)/u(coarseNodeNumber,1)
 end if

end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! engine inlet - enforce the given mass flow function (mass flow * sqrt(temp) over total pressure)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! mass flow function iteration
  storeMassFactor = 99999999999.

  engineMassFlows = 0.
  engineMassFactor = 0.
  area = 0.
  ptfunc = 0.
  massflowfunction = 0.

  if(grp%gridNumber==1) then ! 1 grid
    do i=1,grp%brp%numberOfEngineInletSides
      ind1 = grp%brp%engineInletSideIndexes(i,1)
      ind2 = grp%brp%engineInletSideIndexes(i,2)
      wx = grp%brp%engineInletSideCoefficients(i,:)
      u1 = u(ind1,2:3)
      u2 = u(ind2,2:3)
      engineMassFlows= engineMassFlows + sum(u1*wx) + sum(u2*wx)
      press1 = grp%p(ind1) + 0.5*sum(u(ind1,2:3)*u(ind1,2:3))/u(ind1,1)
      press2 = grp%p(ind2) + 0.5*sum(u(ind2,2:3)*u(ind2,2:3))/u(ind2,1)
      press = 0.5*(press1 + press2)
      temp1 = (ivd%gamma/u(ind1,1))*(u(ind1,4)-(0.5*sum(u(ind1,2:3)*u(ind1,2:3)))/u(ind1,1))
      temp2 = (ivd%gamma/u(ind2,1))*(u(ind2,4)-(0.5*sum(u(ind2,2:3)*u(ind2,2:3)))/u(ind2,1))
      temp = 0.5*(temp1+temp2)
      ptfunc = ptfunc + (sqrt(temp)/press)*sqrt(sum(wx*wx))
      area = area + sqrt(sum(wx*wx))
    end do
   else ! multigrid
    do i=1,grp%brp%numberOfEngineInletSides
      ind1 = grp%brp%engineInletSideIndexes(i,1)
      ind2 = grp%brp%engineInletSideIndexes(i,2)
      wx = grp%brp%engineInletSideCoefficients(i,:)
      engineMassFlows = engineMassFlows + sqrt(sum(wx*wx))
    end do
  end if

  if(engineMassFlows>0.0) then
      massflowfunction = engineMassFlows*(ptfunc/area)
      engineMassFactor = ivd%enginesFrontMassFlow/massflowfunction
! if exceeding the max possible mass flow factor for this duct cross section print warning and leave the loop
      if(engineMassFactor.gt.storeMassFactor)then
        print*,'!!!WARNING!!! You are trying to exceed the maximum achievable mass flow function for engine 1'
        goto 9999
      endif
        storeMassFactor = engineMassFactor
  end if

  if(grp%gridNumber==1) then
    if(ivd%enginesFrontMassFlow>0.0) then
      write(*,*) "Engine mass flow function (calculated): ", massflowfunction
!      write(*,*) "Engine mass flow function (prescribed): ", ivd%enginesFrontMassFlow
!      write(*,*) "Engine mass factor: ",engineMassFactor
!      write(*,*) "Engine inlet length: ", area
    end if
  end if
    
  if(grp%gridNumber==1) then
    ist = grp%brp%faceIndicatorArray(-6)
    ien = grp%brp%faceIndicatorArray(-5)

    do ib =ist,ien
      ip = grp%brp%faceIndicatorArray(ib)
      factV = grp%brp%engineInletNormals*ivd%enginesFrontMassFlow/grp%brp%engineInletAreas
      intE = grp%u(ip,4) - 0.5*sum(grp%u(ip,2:3)*grp%u(ip,2:3))/grp%u(ip,1)
      u(ip,2:3) = u(ip,2:3)*engineMassFactor
      u(ip,4) = intE + 0.5*sum(u(ip,2:3)*u(ip,2:3))/u(ip,1)
      if(ivd%turbulenceModel > 0) then
        u(ip,5) = 0.0
      end if
    end do  

! recompute the pressure field

    call makePressureField(grp,ivd,grp%u)
    
    if(ivd%enginesFrontMassFlow>0.0)then
      factor = abs(massflowfunction - ivd%enginesFrontMassFlow)/ivd%enginesFrontMassFlow
     else
      factor = 0.0
    endif

    if(factor.lt.0.01)goto 9999

  end if
9999 continue

! make sure there are no negative velocities in engine inlet

  if(grp%gridNumber==1) then
    ist = grp%brp%faceIndicatorArray(-6)
    ien = grp%brp%faceIndicatorArray(-5)

    do ib =ist,ien
      ip = grp%brp%faceIndicatorArray(ib)
      normal = grp%brp%engineInletNormals
      normalVelocity = sum(normal*u(ip,2:3))/u(ip,1)
      if(normalVelocity<0.00) then
        write(*,*) "WARNING: Normal velocity for engine inlet node ",ip," is negative ",normalVelocity,normal
        u(ip,2) = u(ip,2)-u(ip,1)*(normalVelocity-0.01)*normal(1)
        u(ip,3) = u(ip,3)-u(ip,1)*(normalVelocity-0.01)*normal(2)
     end if
    end do

  end if
  
  
 ien = grp%brp%faceIndicatorArray(-5)

!
! trailing edges
!

! DO NOT VECTORIZE:
if(ivd%ReynoldsNumber>0.0) then 
do ib=ien+1,grp%brp%numberOfBoundaryFaces
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 u(coarseNodeNumber,2:3) = 0.0
 u(coarseNodeNumber,5) = 0.0
end do
end if
end subroutine setBCsOnSolutionField 
!-------------------------------------------------------------------------
subroutine setBCsOnIncrementField(grp,ivd,rhs)
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: rhs(:,:)

integer :: i,ist,ien,ib,ip,coarseNodeNumber
real :: tx,ty,tangentialVelocityIncrement,ww,vectorNorm,Tw

!
! inviscid wall
!

ist = grp%brp%faceIndicatorArray(-20)
ien = grp%brp%faceIndicatorArray(-19)
if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib = ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  grp%nodeHelpArray(coarseNodeNumber,11:12) = 0.0
 end do
end if

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 tx = -1.0*grp%brp%faceTangentArray(ip,1) 
 ty = -1.0*grp%brp%faceTangentArray(ip,2) 
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  ww = grp%brp%faceWeightNorms(ip)
 else
  coarseNodeNumber = ip
 end if
 tangentialVelocityIncrement = rhs(coarseNodeNumber,2)*tx + rhs(coarseNodeNumber,3)*ty  
 if(associated(grp%boundaryFaceNodeMappings)) then 
  grp%nodeHelpArray(coarseNodeNumber,11) = grp%nodeHelpArray(coarseNodeNumber,11) + ww*tx 
  grp%nodeHelpArray(coarseNodeNumber,12) = grp%nodeHelpArray(coarseNodeNumber,12) + ww*ty 
 else
 ! fine grid
  rhs(coarseNodeNumber,2) = tangentialVelocityIncrement*tx
  rhs(coarseNodeNumber,3) = tangentialVelocityIncrement*ty
 end if
end do
  
if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib = ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  tangentialVelocityIncrement = rhs(coarseNodeNumber,2)*grp%nodeHelpArray(coarseNodeNumber,11) +&
                                rhs(coarseNodeNumber,3)*grp%nodeHelpArray(coarseNodeNumber,12)
  vectorNorm = grp%nodeHelpArray(coarseNodeNumber,11)**2+grp%nodeHelpArray(coarseNodeNumber,12)**2
  tangentialVelocityIncrement = tangentialVelocityIncrement/vectorNorm
  rhs(coarseNodeNumber,2) = tangentialVelocityIncrement*grp%nodeHelpArray(coarseNodeNumber,11)
  rhs(coarseNodeNumber,3) = tangentialVelocityIncrement*grp%nodeHelpArray(coarseNodeNumber,12)
 end do
end if

!
! isothermal viscid wall  
!

ist = grp%brp%faceIndicatorArray(-12)
ien = grp%brp%faceIndicatorArray(-11)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if

 rhs(coarseNodeNumber,2:5) = 0.0
end do

!
! adiabatic viscid wall  
!

ist = grp%brp%faceIndicatorArray(-10)
ien = grp%brp%faceIndicatorArray(-9)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 rhs(coarseNodeNumber,2:3) = 0.0

 rhs(coarseNodeNumber,5) = 0.0
end do

ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
! rhs(coarseNodeNumber,1:4) = 0.0
 else
  coarseNodeNumber = ip
! if(ivd%MachNumber<1) then
!   rhs(coarseNodeNumber,2:4) = 0.0
! end if
 end if
! if(ivd%MachNumber<1) then
!   rhs(coarseNodeNumber,2:4) = 0.0
! end if
end do


! relax engine inlets
! engine inflow 

ist = grp%brp%faceIndicatorArray(-6)
ien = grp%brp%faceIndicatorArray(-5)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
  if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 rhs(coarseNodeNumber,:) = ivd%engineBCRelaxation*rhs(coarseNodeNumber,:)
end do


!
! trailing edges
!
ien = grp%brp%faceIndicatorArray(-5)

! DO NOT VECTORIZE:
if(ivd%ReynoldsNumber>0.0) then 
do ib=ien+1,grp%brp%numberOfBoundaryFaces
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 rhs(coarseNodeNumber,2:3) = 0.0
 rhs(coarseNodeNumber,5) = 0.0
end do
end if
end subroutine setBCsOnIncrementField
!-------------------------------------------------------------------------
subroutine setBCsOnCoarseIncrementField(grp,ivd,delu)
! nullifies coarse mesh increments at inner boundary

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: delu(:,:)

integer :: ist,ien,ib,ip

ist = grp%brp%faceIndicatorArray(-12)
ien = grp%brp%faceIndicatorArray(-11)

do ib =ist,ien
 ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(ib))
 delu(ip,:) = 0.0
end do 

ist = grp%brp%faceIndicatorArray(-10)
ien = grp%brp%faceIndicatorArray(-9)

do ib =ist,ien
 ip = grp%boundaryFaceNodeMappings(grp%brp%faceIndicatorArray(ib))
 delu(ip,:) = 0.0
end do 
end subroutine setBCsOnCoarseIncrementField
!-------------------------------------------------------------------------
subroutine setTemperatureForAdiabaticWall(grp,dTx,dTy)
! sets zero gradient on temperature as weak BC
IMPLICIT NONE

type(GridSolverData) :: grp
real :: dTx(:)
real :: dTy(:)
integer :: ib,ist,ien,ip,coarseNodeNumber
real :: tangentialTemperatureGradient,vectorNorm,tx,ty,ww

ist = grp%brp%faceIndicatorArray(-10)
ien = grp%brp%faceIndicatorArray(-9)

if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib = ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  grp%nodeHelpArray(coarseNodeNumber,11:12) = 0.0
 end do
end if

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 tx = -1.0*grp%brp%faceTangentArray(ip,1) 
 ty = -1.0*grp%brp%faceTangentArray(ip,2) 
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  ww = grp%brp%faceWeightNorms(ip)
 else
  coarseNodeNumber = ip
 end if
 if(associated(grp%boundaryFaceNodeMappings)) then 
  grp%nodeHelpArray(coarseNodeNumber,11) = grp%nodeHelpArray(coarseNodeNumber,11) + ww*tx 
  grp%nodeHelpArray(coarseNodeNumber,12) = grp%nodeHelpArray(coarseNodeNumber,12) + ww*ty 
 else
 ! fine grid
  tangentialTemperatureGradient = dTx(coarseNodeNumber)*tx + dTy(coarseNodeNumber)*ty  
  dTx(coarseNodeNumber) = tangentialTemperatureGradient*tx
  dTy(coarseNodeNumber) = tangentialTemperatureGradient*ty
 end if
end do

if(associated(grp%boundaryFaceNodeMappings)) then 
 do ib = ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
  tangentialTemperatureGradient = dTx(coarseNodeNumber)*grp%nodeHelpArray(coarseNodeNumber,11)+&
                                  dTy(coarseNodeNumber)*grp%nodeHelpArray(coarseNodeNumber,12)
  vectorNorm = grp%nodeHelpArray(coarseNodeNumber,11)**2+grp%nodeHelpArray(coarseNodeNumber,12)**2
  tangentialTemperatureGradient = tangentialTemperatureGradient/vectorNorm
  dTx(coarseNodeNumber) = tangentialTemperatureGradient*grp%nodeHelpArray(coarseNodeNumber,11)
  dTy(coarseNodeNumber) = tangentialTemperatureGradient*grp%nodeHelpArray(coarseNodeNumber,12)
 end do
end if

! trailing edges

! DO NOT VECTORIZE:
do ib=ien+1,grp%brp%numberOfBoundaryFaces
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 dTx(coarseNodeNumber) = 0.0
 dTy(coarseNodeNumber) = 0.0
end do

end subroutine setTemperatureForAdiabaticWall
!-------------------------------------------------------------------------
subroutine makeViscosityTerm(grp,ivd,rhs,recalculateTurbulence)
! makes the viscosity fluxes 
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: rhs(:,:)
logical :: recalculateTurbulence

integer :: i1,i2,i,ip1,ip2,jm1,j0,j1,j2,j3,j4,j5,ic1,ic2,ic3,im2,im1,ip3,ip4
integer :: ind1,ind2,ind1m2,ind2m2,ind1m1,ind2m1,ind1p1,ind2p1,ind1p2,ind2p2,ko
real :: r0,r3,u0,u3,v0,v3,t3,rm,um,vm,tm,r4,u4,v4,t4,r5,u5,v5,t5
real :: c0,c1,c2,c3,h0,h1,h2,h3,h4,w1,w2,w3,w4
real :: xm1,x0,x1,x2,x3,x4,x5,x6,ym1,y0,y1,y2,y3,y4,y5,y6
real :: r1,r2,u1,u2,v1,v2,t1,t2,wx,wy,ux,vx,tx,uy,vy,ty,oneOverReynoldsNumber
real :: utm,unm,ut0,ut1,ut2,ut3,ut4,ut5,un0,un1,un2,un3,un4,un5
real :: f12x,f13x,f14x,f12y,f13y,f14y,f22x,f23x,f24x,f22y,f23y,f24y
real :: f02x,f03x,f04x,f02y,f03y,f04y,f32x,f33x,f34x,f32y,f33y,f34y
real :: f52x,f53x,f54x,f52y,f53y,f54y,f42x,f43x,f44x,f42y,f43y,f44y
real :: f12mx,f13mx,f14mx,f12my,f13my,f14my
real :: f22mx,f23mx,f24mx,f22my,f23my,f24my
real :: f2mx,f3mx,f4mx,f2my,f3my,f4my,a0,a1,a2,a3,a4,am
real :: ft02,ft12,ft22,ft32,ft42,ft52,ia0,ia1,ia2,ia3,ia4,amb
real :: ft03,ft13,ft23,ft33,ft43,ft53
real :: ft04,ft14,ft24,ft34,ft44,ft54
real :: fn02,fn12,fn22,fn32,fn42,fn52
real :: fn03,fn13,fn23,fn33,fn43,fn53
real :: fn04,fn14,fn24,fn34,fn44,fn54
real :: f2tm,f3tm,f4tm,f2nm,f3nm,f4nm
real :: mu,k,T,T0,f2,f3,f4,rx,ry,vt,t11,t21,t31,t41,t51,t12,t22,t32,Ksi
real :: t10, t20, t30, t40, t50, u10, u20, t13, t23, t33, t43, t53, u13, u23
real :: t14, t24, t34, t44, t54, u14, u24!, f42x,f43x,f44x,f42y,f43y,f44y
real :: t1m, t2m, t3m, t4m, t5m, u1m, u2m, um1, um2, vm1, vm2, tm1, tm2
real :: t42,t52,dudx,dvdx,dTdx,dudy,dvdy,dTdy,dt1,dt2,dt3,dt4,dt5,oneOverPrandtlNumber
real :: dudyFD,dvdyFD,dTdyFD,dy,anx,any
real :: twoThirds,fourThirds,gammaOverGammaMinusOne,u11,u12,u21,u22,dut1,dut2
real :: wallTangent(2),wallNormal(2),tempConv,dut11,dut12,dut21,dut22,inflowTemp
real :: adTemp1,adTemp2,turbViscosityCoefficient,turbDiffusionCoefficient
real :: viscosityCoefficient,diffusionCoefficient,oneOverTurbPrandtlNumber
integer :: wallInd

grp%nodeHelpArray(1:grp%numberOfNodes,1:6) = 0.0

gammaOverGammaMinusOne = ivd%gamma/(ivd%gamma-1.0)
tempConv = (ivd%gamma-1.0)*ivd%MachNumber**2
inflowTemp = ivd%inflowTemperature

rewind (599)
!ko=1.0
do i=1,grp%numberOfSides
! indexes of nodes in side
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
 wx = grp%sideWeightsArray(i,1)  
 wy = grp%sideWeightsArray(i,2)  
 anx = wx/sqrt(wx**2+wy**2)
 any = wy/sqrt(wx**2+wy**2)
 if(ivd%HighOrder.and.grp%gridNumber==1) then
   if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
      if(grp%AfterAndBefore(i1,1).eq.i2) then
        j1 = i1
        j2 = i2
        j3 = grp%AfterAndBefore(j2,1)
        j4 = grp%AfterAndBefore(j3,1)
        j0 = grp%AfterAndBefore(j1,2)
        jm1 = grp%AfterAndBefore(j0,2)
      else if(grp%AfterAndBefore(i1,2).eq.i2) then
        j1 = i2
        j2 = i1
        j3 = grp%AfterAndBefore(j2,1)
        j4 = grp%AfterAndBefore(j3,1)
        j0 = grp%AfterAndBefore(j1,2)
        jm1 = grp%AfterAndBefore(j0,2)
      else
        print *,' should not be here'
      end if
      if(j0.eq.0) then
        !debugging
        !if (j1==450) then
        !    print*,'nodes',j1,j2,j3,j4
        !ko=1.0
        !endif
        jm1=1!required?
        j5 = grp%AfterAndBefore(j4,1)
        r1 = grp%u(j1,1)            
        u1 = grp%u(j1,2)/r1        
        v1 = grp%u(j1,3)/r1       
        t1 = gammaOverGammaMinusOne*grp%p(j1)/r1
        ut1= u1*any-v1*anx
        un1= u1*anx+v1*any
        a1 = grp%nodeVolume(j1)

        r2 = grp%u(j2,1)            
        u2 = grp%u(j2,2)/r2        
        v2 = grp%u(j2,3)/r2       
        t2 = gammaOverGammaMinusOne*grp%p(j2)/r2
        ut2= u2*any-v2*anx
        un2= u2*anx+v2*any 
        a2 = grp%nodeVolume(j2)

        r3 = grp%u(j3,1)            
        u3 = grp%u(j3,2)/r3        
        v3 = grp%u(j3,3)/r3       
        t3 = gammaOverGammaMinusOne*grp%p(j3)/r3
        ut3= u3*any-v3*anx
        un3= u3*anx+v3*any
        a3 = grp%nodeVolume(j3)

        r4 = grp%u(j4,1)            
        u4 = grp%u(j4,2)/r4        
        v4 = grp%u(j4,3)/r4       
        t4 = gammaOverGammaMinusOne*grp%p(j4)/r4
        ut4= u4*any-v4*anx
        un4= u4*anx+v4*any  
        a4 = grp%nodeVolume(j4)   
       !debugging
        !if (ko==1.0) then
          !a1=0.05
          !a2=0.1
          !a3=0.1
          !a4=0.1
          !x^3-2*x
         ! r1=-199.0/1000.0 !point value at x=0.1
         ! r2=-783.0/2000.0 !mean values
         ! r3=-2289.0/4000.0
         ! r4=-147.0/200.0
        !endif

      !Approximate Jacobian
        !equation(14)
        ia1=a1
        ia2=a2+ia1
        ia3=a3+ia2
        ia4=a4+ia3

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !jacobian at interpolation location 
        !equation(15)       
        am=c0*ia1+c1*ia2+c2*ia3+c3*ia4

        c0 = grp%hoc3(i,1)
        c1 = grp%hoc3(i,2)
        c2 = grp%hoc3(i,3)
        c3 = grp%hoc3(i,4)
        !jacobian at boundary  
        !equation(20)     
        amb=c0*ia1+c1*ia2+c2*ia3+c3*ia4
        !debugging
        !if (ko==1.0) then
        !  print*,'jacobian at boundary',amb
        !  ko=0.0
        !endif
        !print*,'jacobian',j1,j2,amb
        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)
        !Primitive function
         !equation(8)
        r1 = r1*amb
        !equation(9) 
        r2 = r2*a2
        r3 = r2+r3*a3
        r4 = r3+r4*a4

        ut1 = ut1*amb
        ut2 = ut2*a2
        ut3 = ut2+ut3*a3
        ut4 = ut3+ut4*a4

        un1 = un1*amb
        un2 = un2*a2
        un3 = un2+un3*a3
        un4 = un3+un4*a4

        t1 = t1*amb
        t2 = t2*a2
        t3 = t2+t3*a3
        t4 = t3+t4*a4

        !equation(19)
        rm = c0*r1+c1*r2+c2*r3+c3*r4
        utm= c0*ut1+c1*ut2+c2*ut3+c3*ut4
        unm= c0*un1+c1*un2+c2*un3+c3*un4
        tm = c0*t1+c1*t2+c2*t3+c3*t4

  
        !divide through by Jacobian
        !equation(17) 
        rm=rm/am
        utm=utm/am
        unm=unm/am
        tm=tm/am
        !debugging
        !if (ko==1.0) then
        !  print*,'interpolated value',tm
        !  ko=0.0
          !should be -0.296625
        !endif


        um = utm*any+unm*anx
        vm =-utm*anx+unm*any
        ux = (um+um)*wx
        vx = (vm+vm)*wx
        tx = (tm+tm)*wx
        uy = (um+um)*wy
        vy = (vm+vm)*wy
        ty = (tm+tm)*wy
!
!       rm = (3*r1+6*r2-r3)/8.
!print*,um,vm
      !print*,i1
      else if(jm1.eq.0) then
        r0 = grp%u(j0,1)            
        u0 = grp%u(j0,2)/r0        
        v0 = grp%u(j0,3)/r0       
        t0 = gammaOverGammaMinusOne*grp%p(j0)/r0
        ut0= u0*any-v0*anx
        un0= u0*anx+v0*any
        a0 = grp%nodeVolume(j0)

        r1 = grp%u(j1,1)            
        u1 = grp%u(j1,2)/r1        
        v1 = grp%u(j1,3)/r1       
        t1 = gammaOverGammaMinusOne*grp%p(j1)/r1
        ut1= u1*any-v1*anx
        un1= u1*anx+v1*any
        a1 = grp%nodeVolume(j1)

        r2 = grp%u(j2,1)            
        u2 = grp%u(j2,2)/r2        
        v2 = grp%u(j2,3)/r2       
        t2 = gammaOverGammaMinusOne*grp%p(j2)/r2
        ut2= u2*any-v2*anx
        un2= u2*anx+v2*any
        a2 = grp%nodeVolume(j2)

        r3 = grp%u(j3,1)            
        u3 = grp%u(j3,2)/r3        
        v3 = grp%u(j3,3)/r3       
        t3 = gammaOverGammaMinusOne*grp%p(j3)/r3
        ut3= u3*any-v3*anx
        un3= u3*anx+v3*any
        a3 = grp%nodeVolume(j3)
        
       !debugging
        !if (ko==1.0) then
        !  a0=0.05
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
        !  r0=199.0/1000.0
        !  r1=783.0/2000.0
        !  r2=2289.0/4000.0
        !  r3=147.0/200.0
        !end if

       !Approximate Jacobian
        !equation(14) 
        ia0=a0
        ia1=a1+ia0
        ia2=a2+ia1
        ia3=a3+ia2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location  
        !equation(15)         
        am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

        c0 = grp%hoc3(i,1)
        c1 = grp%hoc3(i,2)
        c2 = grp%hoc3(i,3)
        c3 = grp%hoc3(i,4)
        !Jacobian at boundary   
        !equation(20)      
        amb=c0*ia0+c1*ia1+c2*ia2+c3*ia3

        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)

        !Primitive function
        !equation(8) 
        r0 = r0*amb
        !equation(9) 
        r1 = r1*a1
        r2 = r1+r2*a2
        r3 = r2+r3*a3

        ut0 = ut0*amb
        ut1 = ut1*a1
        ut2 = ut1+ut2*a2
        ut3 = ut2+ut3*a3

        un0 = un0*amb
        un1 = un1*a1
        un2 = un1+un2*a2
        un3 = un2+un3*a3

        t0 = t0*amb
        t1 = t1*a1
        t2 = t1+t2*a2
        t3 = t2+t3*a3

        !equation(10) 
        rm = c0*r0+c1*r1+c2*r2+c3*r3
        utm= c0*ut0+c1*ut1+c2*ut2+c3*ut3
        unm= c0*un0+c1*un1+c2*un2+c3*un3
        tm = c0*t0+c1*t1+c2*t2+c3*t3

 
        !a=sqrt((2.*wx)**2+(2.*wy)**2)
        !rom=rom/a
        !utm=utm/a
        !unm=unm/a
        !enm=enm/a

        !divide through by Jacobian
        !equation(17) 
        rm=rm/am
        utm=utm/am
        unm=unm/am
        tm=tm/am
        !print*, 'tang,norm',utm, unm
       !debugging
        !if (ko==1.0) then
        !  print*,'interpolated value',tm
        !  ko=0.0
        ! should be 0.484375
        !end if

        um = utm*any+unm*anx
        vm =-utm*anx+unm*any
        !print*,'u,v', um, unm

!       rm = (3*r1+6*r2-r3)/8.
!       um = (3*u1+6*u2-u3)/8.
!       vm = (3*v1+6*v2-v3)/8.
!       tm = (3*t1+6*t2-t3)/8.
!
        ux = (um+um)*wx
        vx = (vm+vm)*wx
        tx = (tm+tm)*wx
        uy = (um+um)*wy
        vy = (vm+vm)*wy
        ty = (tm+tm)*wy
       !print*,2,j0
      else if(j4.eq.0) then

        r1 = grp%u(i1,1)              
        u1 = grp%u(i1,2)/r1          
        v1 = grp%u(i1,3)/r1         
        t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
        r2 = grp%u(i2,1)
        u2 = grp%u(i2,2)/r2
        v2 = grp%u(i2,3)/r2
        t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
        ux = (u1+u2)*wx
        vx = (v1+v2)*wx
        tx = (t1+t2)*wx
        uy = (u1+u2)*wy
        vy = (v1+v2)*wy
        ty = (t1+t2)*wy

      else
        r0 = grp%u(j0,1)            
        u0 = grp%u(j0,2)/r0        
        v0 = grp%u(j0,3)/r0       
        t0 = gammaOverGammaMinusOne*grp%p(j0)/r0
        ut0= u0*any-v0*anx
        un0= u0*anx+v0*any
        a0 = grp%nodeVolume(j0)

        r1 = grp%u(j1,1)            
        u1 = grp%u(j1,2)/r1        
        v1 = grp%u(j1,3)/r1       
        t1 = gammaOverGammaMinusOne*grp%p(j1)/r1
        ut1= u1*any-v1*anx
        un1= u1*anx+v1*any
        a1 = grp%nodeVolume(j1)

        r2 = grp%u(j2,1)            
        u2 = grp%u(j2,2)/r2        
        v2 = grp%u(j2,3)/r2       
        t2 = gammaOverGammaMinusOne*grp%p(j2)/r2
        ut2= u2*any-v2*anx
        un2= u2*anx+v2*any
        a2 = grp%nodeVolume(j2)

        r3 = grp%u(j3,1)            
        u3 = grp%u(j3,2)/r3        
        v3 = grp%u(j3,3)/r3       
        t3 = gammaOverGammaMinusOne*grp%p(j3)/r3
        ut3= u3*any-v3*anx
        un3= u3*anx+v3*any
        a3 = grp%nodeVolume(j3)


        c0 = grp%hoc(i,1)!!!!could also be hoc2 
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)

!Primitive function
       !debugging
!        if (ko==1.0) then
!          a0=0.1
!          a1=0.1
!          a2=0.1
!          a3=0.1
!          r0=-159.0/800.0
!          r1=-783.0/2000.0
!          r2=-2289.0/4000.0
!          r3=-147.0/200.0
!         endif

        !equation(9) 
        r0 = r0*a0
        r1 = r0+r1*a1
        r2 = r1+r2*a2
        r3 = r2+r3*a3

        ut0 = ut0*a0
        ut1 = ut0+ut1*a1
        ut2 = ut1+ut2*a2
        ut3 = ut2+ut3*a3

        un0 = un0*a0
        un1 = un0+un1*a1
        un2 = un1+un2*a2
        un3 = un2+un3*a3

        t0 = t0*a0
        t1 = t0+t1*a1
        t2 = t1+t2*a2
        t3 = t2+t3*a3

        !equation(10) 
        rm = c0*r0+c1*r1+c2*r2+c3*r3
        utm= c0*ut0+c1*ut1+c2*ut2+c3*ut3
        unm= c0*un0+c1*un1+c2*un2+c3*un3
        tm = c0*t0+c1*t1+c2*t2+c3*t3
!
!       rm = (-r0+9*r1+9*r2-r3)/16.
!       umt= (-ut0+9*ut1+9*ut2-ut3)/16.
!       vmt= (-vt0+9*vt1+9*vt2-vt3)/16.
!       tm = (-t0+9*t1+9*t2-t3)/16.

       !Approximate Jacobian
        !equation(14) 
        a0=a0
        a1=a1+a0
        a2=a2+a1
        a3=a3+a2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !equation(15)        
        am=c0*a0+c1*a1+c2*a2+c3*a3

        !divide through by Jacobian (if primitive function is based on area)
        !equation(17) 
        rm=rm/am
        utm=utm/am
        unm=unm/am
        tm=tm/am
        !debugging
        !if (ko==1.0) then
        !  print*,'interpolated value',tm
        !  ko=0.0
          !should give -0.484375
        !endif
        
        um =  utm*any+unm*anx
        vm = -utm*anx+unm*any
        ux = (um+um)*wx
        vx = (vm+vm)*wx
        tx = (tm+tm)*wx
        uy = (um+um)*wy
        vy = (vm+vm)*wy
        ty = (tm+tm)*wy
 !debugging
      !print*,i1
      !if (i1==10) then
      !print*,'test',ux,wx,uy,wy
      !endif
      end if

   else

      r1 = grp%u(i1,1)
      u1 = grp%u(i1,2)/r1
      v1 = grp%u(i1,3)/r1
      t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
      r2 = grp%u(i2,1)
      u2 = grp%u(i2,2)/r2
      v2 = grp%u(i2,3)/r2
      t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
      ux = (u1+u2)*wx
      vx = (v1+v2)*wx
      tx = (t1+t2)*wx
      uy = (u1+u2)*wy
      vy = (v1+v2)*wy
      ty = (t1+t2)*wy

   end if

 else

      r1 = grp%u(i1,1)
      u1 = grp%u(i1,2)/r1
      v1 = grp%u(i1,3)/r1
      t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
      r2 = grp%u(i2,1)
      u2 = grp%u(i2,2)/r2
      v2 = grp%u(i2,3)/r2
      t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
      ux = (u1+u2)*wx
      vx = (v1+v2)*wx
      tx = (t1+t2)*wx
      uy = (u1+u2)*wy
      vy = (v1+v2)*wy
      ty = (t1+t2)*wy

 end if


 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + ux
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + vx
 grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + tx
 grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + uy
 grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + vy
 grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + ty

 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - ux
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - vx
 grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - tx
 grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - uy
 grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - vy
 grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - ty            
!debugging
!if (i1==10) then
!print*,ux,uy,vx,tx,vy,ty
!endif
!if (i2==10) then
!print*,ux,uy,vx,tx,vy,ty
!endif

end do

do i=1,grp%brp%numberOfBoundaryFaces
 i1 = grp%brp%faceIndexArray(i,1)
 i2 = grp%brp%faceIndexArray(i,2)
 wx   = grp%brp%faceWeightsArray(i,1)
 wy   = grp%brp%faceWeightsArray(i,2) 

 r1   = grp%u(i1,1)
 u1   = grp%u(i1,2)/r1
 v1   = grp%u(i1,3)/r1
 t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
 r2   = grp%u(i2,1)
 u2   = grp%u(i2,2)/r2
 v2   = grp%u(i2,3)/r2
 t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
 ux   = u1+u2
 vx   = v1+v2
 tx   = t1+t2
 uy   = u1+u2
 vy   = v1+v2
 ty   = t1+t2

! adding boundary face contributions        

 if(ivd%boundaryTerm==1) then 
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*u1+ux)*wx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*v1+vx)*wx
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*t1+tx)*wx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*u1+uy)*wy
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*v1+vy)*wy
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*t1+ty)*wy
       
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*u2+ux)*wx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*v2+vx)*wx
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*t2+tx)*wx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*u2+uy)*wy
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*v2+vy)*wy
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*t2+ty)*wy     
 else
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + 4.*u1*wx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + 4.*v1*wx
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + 4.*t1*wx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + 4.*u1*wy
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + 4.*v1*wy
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + 4.*t1*wy
  
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + 4.*u2*wx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + 4.*v2*wx
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + 4.*t2*wx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + 4.*u2*wy
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + 4.*v2*wy
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + 4.*t2*wy
 end if 
end do

!debugging, print gradient
!open(12,file='gradient.dat')
!do i=1,grp%numberOfNodes
!    write(12,'(I5,2(E17.8))')i,grp%nodeHelpArray(i,1),grp%nodeHelpArray(i,4)
!enddo
!close(12)
    
twoThirds = 2./3.
fourThirds = 4./3.
T0 = 1./((ivd%gamma-1)*ivd%MachNumber**2) 
oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
oneOverPrandtlNumber = 1./ivd%PrandtlNumber
oneOverTurbPrandtlNumber = 1./ivd%turbulentPrandtlNumber 
! divide by volume to remove mass matrix on LHS
do i=1,grp%numberOfNodes
 vt = 1./grp%nodeVolume(i)
 dudx = grp%nodeHelpArray(i,1)*vt
 dvdx = grp%nodeHelpArray(i,2)*vt
 dTdx = grp%nodeHelpArray(i,3)*vt
 dudy = grp%nodeHelpArray(i,4)*vt
 dvdy = grp%nodeHelpArray(i,5)*vt
 dTdy = grp%nodeHelpArray(i,6)*vt


! create the viscosity matrix  

 T = gammaOverGammaMinusOne*grp%p(i)/grp%u(i,1)
 if(T<0) then 
  write(*,*) "T negative ",T," at node ",i,grp%u(i,1),grp%p(i)
  T = abs(T)
 end if

 ! Sunderlands law of viscosity 
 mu = oneOverReynoldsNumber*(((tempConv*T)**1.5)*((inflowTemp+198.6)/(tempConv*inflowTemp*T+198.6))) 
! mu = oneOverReynoldsNumber  
 k = oneOverPrandtlNumber*mu
 grp%laminarViscosity(i) = mu
 grp%vorticity(i) = sqrt((dudy-dvdx)**2) ! vorticity array
 grp%divergence(i) = dudx + dvdy

 ! add turbulence effects
 if(ivd%turbulenceModel==1) then ! Spalart-Allmaras
  turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i,5)

  turbViscosityCoefficient = grp%u(i,1)*turbViscosityCoefficient  ! KAS

  Ksi = turbViscosityCoefficient/mu
  turbViscosityCoefficient =&
  turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)

  mu = mu + turbViscosityCoefficient
  k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber
 end if

! stress tensor
 grp%nodeHelpArray(i,1) = mu*(fourThirds*dudx-twoThirds*dvdy) 
 grp%nodeHelpArray(i,2) = mu*(dudy+dvdx) 
 grp%nodeHelpArray(i,3) = mu*(fourThirds*dvdy-twoThirds*dudx) 
 grp%nodeHelpArray(i,4) = k*dTdx 
 grp%nodeHelpArray(i,5) = k*dTdy 
end do


! for adiabatic wall, set temperature gradient equal to zero
 call setTemperatureForAdiabaticWall(grp,grp%nodeHelpArray(:,4),grp%nodeHelpArray(:,5))

! make stress tensor along boundary

if(.not.associated(grp%sourceTerm)) then 
 do i=1,grp%brp%numberOfBoundaryFaces
  wallTangent(1) = grp%brp%faceTangentArray(i,1)
  wallTangent(2) = grp%brp%faceTangentArray(i,2)
  wallNormal(1) = wallTangent(2)
  wallNormal(2) = -1.0*wallTangent(1)
  grp%wallStress(i,1) = (grp%nodeHelpArray(i,1)*wallNormal(1)&
                    +grp%nodeHelpArray(i,2)*wallNormal(2))
  grp%wallStress(i,2) = (grp%nodeHelpArray(i,2)*wallNormal(1)&
                    +grp%nodeHelpArray(i,3)*wallNormal(2)) 
 end do
end if

grp%nodeHelpArray(1:grp%numberOfNodes,9:11) = 0.0

do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
 wx = grp%sideWeightsArray(i,1)
 wy = grp%sideWeightsArray(i,2)  
 anx = wx/sqrt(wx**2+wy**2)
 any = wy/sqrt(wx**2+wy**2)

 if(ivd%HighOrder.and.grp%gridNumber==1) then
 ! if(.false.) then
   if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
      if(grp%AfterAndBefore(i1,1).eq.i2) then
        j1 = i1
        j2 = i2
        j3 = grp%AfterAndBefore(j2,1)
        j4 = grp%AfterAndBefore(j3,1)
        j0 = grp%AfterAndBefore(j1,2)
        jm1 = grp%AfterAndBefore(j0,2)
      else if(grp%AfterAndBefore(i1,2).eq.i2) then
        j1 = i2
        j2 = i1
        j3 = grp%AfterAndBefore(j2,1)
        j4 = grp%AfterAndBefore(j3,1)
        j0 = grp%AfterAndBefore(j1,2)
        jm1 = grp%AfterAndBefore(j0,2)
      else
        print *,' should not be here'
      end if
      if(j0.eq.0) then
        j5 = grp%AfterAndBefore(j4,1)

        t11= grp%nodeHelpArray(j1,1)
        t21= grp%nodeHelpArray(j1,2)
        t31= grp%nodeHelpArray(j1,3)
        t41= grp%nodeHelpArray(j1,4)
        t51= grp%nodeHelpArray(j1,5)
        u11= grp%u(j1,2)/grp%u(j1,1)
        u21= grp%u(j1,3)/grp%u(j1,1) 
        a1 = grp%nodeVolume(j1)

        t12= grp%nodeHelpArray(j2,1)
        t22= grp%nodeHelpArray(j2,2)
        t32= grp%nodeHelpArray(j2,3)
        t42= grp%nodeHelpArray(j2,4)
        t52= grp%nodeHelpArray(j2,5)
        u12= grp%u(j2,2)/grp%u(j2,1)
        u22= grp%u(j2,3)/grp%u(j2,1) 
        a2 = grp%nodeVolume(j2)

        t13= grp%nodeHelpArray(j3,1)
        t23= grp%nodeHelpArray(j3,2)
        t33= grp%nodeHelpArray(j3,3)
        t43= grp%nodeHelpArray(j3,4)
        t53= grp%nodeHelpArray(j3,5)
        u13= grp%u(j3,2)/grp%u(j3,1)
        u23= grp%u(j3,3)/grp%u(j3,1) 
        a3 = grp%nodeVolume(j3)

        t14= grp%nodeHelpArray(j4,1)
        t24= grp%nodeHelpArray(j4,2)
        t34= grp%nodeHelpArray(j4,3)
        t44= grp%nodeHelpArray(j4,4)
        t54= grp%nodeHelpArray(j4,5)
        u14= grp%u(j4,2)/grp%u(j4,1)
        u24= grp%u(j4,3)/grp%u(j4,1) 
        a4 = grp%nodeVolume(j4)

        f12x = t11
        f13x = t21 
        f14x = t41+u11*t11+u21*t21
        f12y = t21
        f13y = t31
        f14y = t51+u11*t21+u21*t31
        ft12= f12x*any-f12y*anx
        fn12= f12x*anx+f12y*any
        ft13= f13x*any-f13y*anx
        fn13= f13x*anx+f13y*any
        ft14= f14x*any-f14y*anx
        fn14= f14x*anx+f14y*any
        f22x = t12
        f23x = t22 
        f24x = t42+u12*t12+u22*t22
        f22y = t22
        f23y = t32
        f24y = t52+u12*t22+u22*t32
        ft22= f22x*any-f22y*anx
        fn22= f22x*anx+f22y*any
        ft23= f23x*any-f23y*anx
        fn23= f23x*anx+f23y*any
        ft24= f24x*any-f24y*anx
        fn24= f24x*anx+f24y*any
        f32x = t13
        f33x = t23 
        f34x = t43+u13*t13+u23*t23
        f32y = t23
        f33y = t33
        f34y = t53+u13*t23+u23*t33
        ft32= f32x*any-f32y*anx
        fn32= f32x*anx+f32y*any
        ft33= f33x*any-f33y*anx
        fn33= f33x*anx+f33y*any
        ft34= f34x*any-f34y*anx
        fn34= f34x*anx+f34y*any
        f42x = t14
        f43x = t24 
        f44x = t44+u14*t14+u24*t24
        f42y = t24
        f43y = t34
        f44y = t54+u14*t24+u24*t34
        ft42= f42x*any-f42y*anx
        fn42= f42x*anx+f42y*any
        ft43= f43x*any-f43y*anx
        fn43= f43x*anx+f43y*any
        ft44= f44x*any-f44y*anx
        fn44= f44x*anx+f44y*any

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)

        !primitive functions
    
        !debugging
!        if (ko==1.0) then
!          a1=0.05
!          a2=0.1
!          a3=0.1
!          a4=0.1
          !x^3-2*x mean gradients
!          ft12=-1.9525 
!          ft22=-1.8775 
!          ft32=-1.7275
!          ft42=-1.5175
!        endif

        !equation(32) 
        ft12 = ft12*a1
        ft22 = ft12+ft22*a2
        ft32 = ft22+ft32*a3
        ft42 = ft32+ft42*a4

        fn12 = fn12*a1
        fn22 = fn12+fn22*a2
        fn32 = fn22+fn32*a3
        fn42 = fn32+fn42*a4

        ft13 = ft13*a1
        ft23 = ft13+ft23*a2
        ft33 = ft23+ft33*a3
        ft43 = ft33+ft43*a4

        fn13 = fn13*a1
        fn23 = fn13+fn23*a2
        fn33 = fn23+fn33*a3
        fn43 = fn33+fn43*a4

        ft14 = ft14*a1
        ft24 = ft14+ft24*a2
        ft34 = ft24+ft34*a3
        ft44 = ft34+ft44*a4

        fn14 = fn14*a1
        fn24 = fn14+fn24*a2
        fn34 = fn24+fn34*a3
        fn44 = fn34+fn44*a4

        !equation(33) 
        f2tm= c0*ft12+c1*ft22+c2*ft32+c3*ft42
        f2nm= c0*fn12+c1*fn22+c2*fn32+c3*fn42
        f3tm= c0*ft13+c1*ft23+c2*ft33+c3*ft43
        f3nm= c0*fn13+c1*fn23+c2*fn33+c3*fn43
        f4tm= c0*ft14+c1*ft24+c2*ft34+c3*ft44
        f4nm= c0*fn14+c1*fn24+c2*fn34+c3*fn44

        !Approximate Jacobian
        !equation(14) 
        a1= a1
        a2= a2+a1
        a3= a3+a2
        a4= a4+a3
 
        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)

        !equation(20)        
        am=c0*a1+c1*a2+c2*a3+c3*a4

        !divide through by Jacobian
        !equation(37) 
        f2tm = f2tm/am
        f2nm = f2nm/am
        f3tm = f3tm/am
        f3nm = f3nm/am
        f4tm = f4tm/am
        f4nm = f4nm/am
       !debugging
        !if (ko==1.0) then
        !  print*,'interpolated value',f4nm
        !  ko=0.0
          !should be -1.9325
        !endif

        f2mx=( f2tm*any+f2nm*anx)*2
        f2my=(-f2tm*anx+f2nm*any)*2
        f3mx=( f3tm*any+f3nm*anx)*2
        f3my=(-f3tm*anx+f3nm*any)*2
        f4mx=( f4tm*any+f4nm*anx)*2
        f4my=(-f4tm*anx+f4nm*any)*2

      else if(jm1.eq.0) then
        !j5 = grp%AfterAndBefore(j4,1)
        t10= grp%nodeHelpArray(j0,1)
        t20= grp%nodeHelpArray(j0,2)
        t30= grp%nodeHelpArray(j0,3)
        t40= grp%nodeHelpArray(j0,4)
        t50= grp%nodeHelpArray(j0,5)
        u10= grp%u(j0,2)/grp%u(j0,1)
        u20= grp%u(j0,3)/grp%u(j0,1) 
        a0 = grp%nodeVolume(j0)

        t11= grp%nodeHelpArray(j1,1)
        t21= grp%nodeHelpArray(j1,2)
        t31= grp%nodeHelpArray(j1,3)
        t41= grp%nodeHelpArray(j1,4)
        t51= grp%nodeHelpArray(j1,5)
        u11= grp%u(j1,2)/grp%u(j1,1)
        u21= grp%u(j1,3)/grp%u(j1,1) 
        a1 = grp%nodeVolume(j1)

        t12= grp%nodeHelpArray(j2,1)
        t22= grp%nodeHelpArray(j2,2)
        t32= grp%nodeHelpArray(j2,3)
        t42= grp%nodeHelpArray(j2,4)
        t52= grp%nodeHelpArray(j2,5)
        u12= grp%u(j2,2)/grp%u(j2,1)
        u22= grp%u(j2,3)/grp%u(j2,1) 
        a2 = grp%nodeVolume(j2)

        t13= grp%nodeHelpArray(j3,1)
        t23= grp%nodeHelpArray(j3,2)
        t33= grp%nodeHelpArray(j3,3)
        t43= grp%nodeHelpArray(j3,4)
        t53= grp%nodeHelpArray(j3,5)
        u13= grp%u(j3,2)/grp%u(j3,1)
        u23= grp%u(j3,3)/grp%u(j3,1) 
        a3 = grp%nodeVolume(j3)

        f02x = t10
        f03x = t20 
        f04x = t40+u10*t10+u20*t20
        f02y = t20
        f03y = t30
        f04y = t50+u10*t20+u20*t30
        ft02= f02x*any-f02y*anx
        fn02= f02x*anx+f02y*any
        ft03= f03x*any-f03y*anx
        fn03= f03x*anx+f03y*any
        ft04= f04x*any-f04y*anx
        fn04= f04x*anx+f04y*any
        f12x = t11
        f13x = t21 
        f14x = t41+u11*t11+u21*t21
        f12y = t21
        f13y = t31
        f14y = t51+u11*t21+u21*t31
        ft12= f12x*any-f12y*anx
        fn12= f12x*anx+f12y*any
        ft13= f13x*any-f13y*anx
        fn13= f13x*anx+f13y*any
        ft14= f14x*any-f14y*anx
        fn14= f14x*anx+f14y*any
        f22x = t12
        f23x = t22 
        f24x = t42+u12*t12+u22*t22
        f22y = t22
        f23y = t32
        f24y = t52+u12*t22+u22*t32
        ft22= f22x*any-f22y*anx
        fn22= f22x*anx+f22y*any
        ft23= f23x*any-f23y*anx
        fn23= f23x*anx+f23y*any
        ft24= f24x*any-f24y*anx
        fn24= f24x*anx+f24y*any
        f32x = t13
        f33x = t23 
        f34x = t43+u13*t13+u23*t23
        f32y = t23
        f33y = t33
        f34y = t53+u13*t23+u23*t33
        ft32= f32x*any-f32y*anx
        fn32= f32x*anx+f32y*any
        ft33= f33x*any-f33y*anx
        fn33= f33x*anx+f33y*any
        ft34= f34x*any-f34y*anx
        fn34= f34x*anx+f34y*any

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)

                !primitive functions
        !debugging
!        if (ko==1.0) then
!          a0=0.05
!          a1=0.1
!          a2=0.1
!          a3=0.1
          !x^3-2*x mean gradients
!          ft02=-1.9525 
!          ft12=-1.8775 
!          ft22=-1.7275
!          ft32=-1.5175
!        endif

        !equation(32) 
        ft02 = ft02*a0
        ft12 = ft02+ft12*a1
        ft22 = ft12+ft22*a2
        ft32 = ft22+ft32*a3

        fn02 = fn02*a0
        fn12 = fn02+fn12*a1
        fn22 = fn12+fn22*a2
        fn32 = fn22+fn32*a3

        ft03 = ft03*a0
        ft13 = ft03+ft13*a1
        ft23 = ft13+ft23*a2
        ft33 = ft23+ft33*a3

        fn03 = fn03*a0
        fn13 = fn03+fn13*a1
        fn23 = fn13+fn23*a2
        fn33 = fn23+fn33*a3

        ft04 = ft04*a0
        ft14 = ft04+ft14*a1
        ft24 = ft14+ft24*a2
        ft34 = ft24+ft34*a3

        fn04 = fn04*a0
        fn14 = fn04+fn14*a1
        fn24 = fn14+fn24*a2
        fn34 = fn24+fn34*a3

        !equation(33) 
        f2tm= c0*ft02+c1*ft12+c2*ft22+c3*ft32
        f2nm= c0*fn02+c1*fn12+c2*fn22+c3*fn32
        f3tm= c0*ft03+c1*ft13+c2*ft23+c3*ft33
        f3nm= c0*fn03+c1*fn13+c2*fn23+c3*fn33
        f4tm= c0*ft04+c1*ft14+c2*ft24+c3*ft34
        f4nm= c0*fn04+c1*fn14+c2*fn24+c3*fn34

       !Approximate Jacobian
        !equation(14) 
        a0=a0
        a1=a1+a0
        a2=a2+a1
        a3=a3+a2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !equation(20)         
        am=c0*a0+c1*a1+c2*a2+c3*a3

        !divide through by Jacobian
        !equation(37) 
        f2tm = f2tm/am
        f2nm = f2nm/am
        f3tm = f3tm/am
        f3nm = f3nm/am
        f4tm = f4tm/am
        f4nm = f4nm/am

        !debugging
        !if (ko==1.0) then
        !  print*,'interpolated value',f4nm
        !  ko=0.0
          !should be -1.8125
        !endif

        f2mx=( f2tm*any+f2nm*anx)*2
        f2my=(-f2tm*anx+f2nm*any)*2
        f3mx=( f3tm*any+f3nm*anx)*2
        f3my=(-f3tm*anx+f3nm*any)*2
        f4mx=( f4tm*any+f4nm*anx)*2
        f4my=(-f4tm*anx+f4nm*any)*2

      else if(j4.eq.0) then

        t11 = grp%nodeHelpArray(i1,1)
        t21 = grp%nodeHelpArray(i1,2)
        t31 = grp%nodeHelpArray(i1,3)
        t41 = grp%nodeHelpArray(i1,4)
        t51 = grp%nodeHelpArray(i1,5)
        u11 = grp%u(i1,2)/grp%u(i1,1)
        u21 = grp%u(i1,3)/grp%u(i1,1) 
 
        t12 = grp%nodeHelpArray(i2,1)
        t22 = grp%nodeHelpArray(i2,2)
        t32 = grp%nodeHelpArray(i2,3)
        t42 = grp%nodeHelpArray(i2,4)
        t52 = grp%nodeHelpArray(i2,5)
        u12 = grp%u(i2,2)/grp%u(i2,1)
        u22 = grp%u(i2,3)/grp%u(i2,1)

        f12x = t11
        f13x = t21 
        f14x = t41+u11*t11+u21*t21
        f12y = t21
        f13y = t31
        f14y = t51+u11*t21+u21*t31
        f22x = t12
        f23x = t22 
        f24x = t42+u12*t12+u22*t22
        f22y = t22
        f23y = t32
        f24y = t52+u12*t22+u22*t32
        f2mx = f12x+f22x
        f3mx = f13x+f23x
        f4mx = f14x+f24x
        f2my = f12y+f22y
        f3my = f13y+f23y
        f4my = f14y+f24y
 
      else
        t10= grp%nodeHelpArray(j0,1)
        t20= grp%nodeHelpArray(j0,2)
        t30= grp%nodeHelpArray(j0,3)
        t40= grp%nodeHelpArray(j0,4)
        t50= grp%nodeHelpArray(j0,5)
        u10= grp%u(j0,2)/grp%u(j0,1)
        u20= grp%u(j0,3)/grp%u(j0,1) 
        a0 = grp%nodeVolume(j0)

        t11= grp%nodeHelpArray(j1,1)
        t21= grp%nodeHelpArray(j1,2)
        t31= grp%nodeHelpArray(j1,3)
        t41= grp%nodeHelpArray(j1,4)
        t51= grp%nodeHelpArray(j1,5)
        u11= grp%u(j1,2)/grp%u(j1,1)
        u21= grp%u(j1,3)/grp%u(j1,1) 
        a1 = grp%nodeVolume(j1)

        t12= grp%nodeHelpArray(j2,1)
        t22= grp%nodeHelpArray(j2,2)
        t32= grp%nodeHelpArray(j2,3)
        t42= grp%nodeHelpArray(j2,4)
        t52= grp%nodeHelpArray(j2,5)
        u12= grp%u(j2,2)/grp%u(j2,1)
        u22= grp%u(j2,3)/grp%u(j2,1) 
        a2 = grp%nodeVolume(j2)

        t13= grp%nodeHelpArray(j3,1)
        t23= grp%nodeHelpArray(j3,2)
        t33= grp%nodeHelpArray(j3,3)
        t43= grp%nodeHelpArray(j3,4)
        t53= grp%nodeHelpArray(j3,5)
        u13= grp%u(j3,2)/grp%u(j3,1)
        u23= grp%u(j3,3)/grp%u(j3,1) 
        a3 = grp%nodeVolume(j3)

        f02x = t10
        f03x = t20 
        f04x = t40+u10*t10+u20*t20
        f02y = t20
        f03y = t30
        f04y = t50+u10*t20+u20*t30
        ft02= f02x*any-f02y*anx
        fn02= f02x*anx+f02y*any
        ft03= f03x*any-f03y*anx
        fn03= f03x*anx+f03y*any
        ft04= f04x*any-f04y*anx
        fn04= f04x*anx+f04y*any
        f12x = t11
        f13x = t21 
        f14x = t41+u11*t11+u21*t21
        f12y = t21
        f13y = t31
        f14y = t51+u11*t21+u21*t31
        ft12= f12x*any-f12y*anx
        fn12= f12x*anx+f12y*any
        ft13= f13x*any-f13y*anx
        fn13= f13x*anx+f13y*any
        ft14= f14x*any-f14y*anx
        fn14= f14x*anx+f14y*any
        f22x = t12
        f23x = t22 
        f24x = t42+u12*t12+u22*t22
        f22y = t22
        f23y = t32
        f24y = t52+u12*t22+u22*t32
        ft22= f22x*any-f22y*anx
        fn22= f22x*anx+f22y*any
        ft23= f23x*any-f23y*anx
        fn23= f23x*anx+f23y*any
        ft24= f24x*any-f24y*anx
        fn24= f24x*anx+f24y*any
        f32x = t13
        f33x = t23 
        f34x = t43+u13*t13+u23*t23
        f32y = t23
        f33y = t33
        f34y = t53+u13*t23+u23*t33
        ft32= f32x*any-f32y*anx
        fn32= f32x*anx+f32y*any
        ft33= f33x*any-f33y*anx
        fn33= f33x*anx+f33y*any
        ft34= f34x*any-f34y*anx
        fn34= f34x*anx+f34y*any

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
               !primitive functions
        !debugging
!        if (ko==1.0) then
!          a0=0.1
!          a1=0.1
!          a2=0.1
!          a3=0.1
          !x^3-2*x mean gradients
!          ft02=-1.9675 
!          ft12=-1.8775 
!          ft22=-1.7275
!          ft32=-1.5175
!        endif

        !equation(32) 
        ft02 = ft02*a0
        ft12 = ft02+ft12*a1
        ft22 = ft12+ft22*a2
        ft32 = ft22+ft32*a3

        fn02 = fn02*a0
        fn12 = fn02+fn12*a1
        fn22 = fn12+fn22*a2
        fn32 = fn22+fn32*a3

        ft03 = ft03*a0
        ft13 = ft03+ft13*a1
        ft23 = ft13+ft23*a2
        ft33 = ft23+ft33*a3

        fn03 = fn03*a0
        fn13 = fn03+fn13*a1
        fn23 = fn13+fn23*a2
        fn33 = fn23+fn33*a3

        ft04 = ft04*a0
        ft14 = ft04+ft14*a1
        ft24 = ft14+ft24*a2
        ft34 = ft24+ft34*a3

        fn04 = fn04*a0
        fn14 = fn04+fn14*a1
        fn24 = fn14+fn24*a2
        fn34 = fn24+fn34*a3

        !equation(33) 
        f2tm= c0*ft02+c1*ft12+c2*ft22+c3*ft32
        f2nm= c0*fn02+c1*fn12+c2*fn22+c3*fn32
        f3tm= c0*ft03+c1*ft13+c2*ft23+c3*ft33
        f3nm= c0*fn03+c1*fn13+c2*fn23+c3*fn33
        f4tm= c0*ft04+c1*ft14+c2*ft24+c3*ft34
        f4nm= c0*fn04+c1*fn14+c2*fn24+c3*fn34
 
        !Approximate Jacobian
        !equation(14) 
        a0=a0
        a1=a1+a0
        a2=a2+a1
        a3=a3+a2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !equation(15)         
        am=c0*a0+c1*a1+c2*a2+c3*a3

        !divide through by Jacobian
        !equation(37) 
        f2tm = f2tm/am
        f2nm = f2nm/am
        f3tm = f3tm/am
        f3nm = f3nm/am
        f4tm = f4tm/am
        f4nm = f4nm/am

       !debugging
        !if (ko==1.0) then
        !  print*,'interpolated value',f4nm
        !  ko=0.0
          !should be -1.8125
        !endif

        f2mx=( f2tm*any+f2nm*anx)*2
        f2my=(-f2tm*anx+f2nm*any)*2
        f3mx=( f3tm*any+f3nm*anx)*2
        f3my=(-f3tm*anx+f3nm*any)*2
        f4mx=( f4tm*any+f4nm*anx)*2
        f4my=(-f4tm*anx+f4nm*any)*2

      end if

   else

       t11 = grp%nodeHelpArray(i1,1)
       t21 = grp%nodeHelpArray(i1,2)
       t31 = grp%nodeHelpArray(i1,3)
       t41 = grp%nodeHelpArray(i1,4)
       t51 = grp%nodeHelpArray(i1,5)
       u11 = grp%u(i1,2)/grp%u(i1,1)
       u21 = grp%u(i1,3)/grp%u(i1,1) 

       t12 = grp%nodeHelpArray(i2,1)
       t22 = grp%nodeHelpArray(i2,2)
       t32 = grp%nodeHelpArray(i2,3)
       t42 = grp%nodeHelpArray(i2,4)
       t52 = grp%nodeHelpArray(i2,5)
       u12 = grp%u(i2,2)/grp%u(i2,1)
       u22 = grp%u(i2,3)/grp%u(i2,1)

       f12x = t11
       f13x = t21 
       f14x = t41+u11*t11+u21*t21
       f12y = t21
       f13y = t31
       f14y = t51+u11*t21+u21*t31
       f22x = t12
       f23x = t22 
       f24x = t42+u12*t12+u22*t22
       f22y = t22
       f23y = t32
       f24y = t52+u12*t22+u22*t32
       f2mx = f12x+f22x
       f3mx = f13x+f23x
       f4mx = f14x+f24x
       f2my = f12y+f22y
       f3my = f13y+f23y
       f4my = f14y+f24y

   end if

 else

     t11 = grp%nodeHelpArray(i1,1)
     t21 = grp%nodeHelpArray(i1,2)
     t31 = grp%nodeHelpArray(i1,3)
     t41 = grp%nodeHelpArray(i1,4)
     t51 = grp%nodeHelpArray(i1,5)
     u11 = grp%u(i1,2)/grp%u(i1,1)
     u21 = grp%u(i1,3)/grp%u(i1,1) 

     t12 = grp%nodeHelpArray(i2,1)
     t22 = grp%nodeHelpArray(i2,2)
     t32 = grp%nodeHelpArray(i2,3)
     t42 = grp%nodeHelpArray(i2,4)
     t52 = grp%nodeHelpArray(i2,5)
     u12 = grp%u(i2,2)/grp%u(i2,1)
     u22 = grp%u(i2,3)/grp%u(i2,1)

     f12x = t11
     f13x = t21 
     f14x = t41+u11*t11+u21*t21
     f12y = t21
     f13y = t31
     f14y = t51+u11*t21+u21*t31
     f22x = t12
     f23x = t22 
     f24x = t42+u12*t12+u22*t22
     f22y = t22
     f23y = t32
     f24y = t52+u12*t22+u22*t32
     f2mx = f12x+f22x
     f3mx = f13x+f23x
     f4mx = f14x+f24x
     f2my = f12y+f22y
     f3my = f13y+f23y
     f4my = f14y+f24y

 end if

!write(599,*) i,i1,i2
!write(599,'(6F15.7)') f2mx,f2my,f3mx,f3my,f4mx,f4my

 f2 = wx*f2mx+ wy*f2my
 f3 = wx*f3mx+ wy*f3my
 f4 = wx*f4mx+ wy*f4my


 grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + f2 
 grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + f3  
 grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + f4 

 grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - f2 
 grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - f3  
 grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - f4 
end do

do i=1,grp%brp%numberOfBoundaryFaces
 i1 = grp%brp%faceIndexArray(i,1)
 i2 = grp%brp%faceIndexArray(i,2)
 wx   = grp%brp%faceWeightsArray(i,1)
 wy   = grp%brp%faceWeightsArray(i,2) 
 
 t11 = grp%nodeHelpArray(i1,1)
 t21 = grp%nodeHelpArray(i1,2)
 t31 = grp%nodeHelpArray(i1,3)
 t41 = grp%nodeHelpArray(i1,4)
 t51 = grp%nodeHelpArray(i1,5)
 u11 = grp%u(i1,2)/grp%u(i1,1)
 u21 = grp%u(i1,3)/grp%u(i1,1) 

 t12 = grp%nodeHelpArray(i2,1)
 t22 = grp%nodeHelpArray(i2,2)
 t32 = grp%nodeHelpArray(i2,3)
 t42 = grp%nodeHelpArray(i2,4)
 t52 = grp%nodeHelpArray(i2,5)
 u12 = grp%u(i2,2)/grp%u(i2,1)
 u22 = grp%u(i2,3)/grp%u(i2,1)

 dt1 = t11 + t12 
 dt2 = t21 + t22 
 dt3 = t31 + t32 
 dt4 = t41 + t42
 dt5 = t51 + t52

 dut11 = u11*t11 + u21*t21
 dut21 = u11*t21 + u21*t31
 dut12 = u12*t12 + u22*t22
 dut22 = u12*t22 + u22*t32

 dut1 = dut11 + dut12 
 dut2 = dut21 + dut22  

 f2 = wx*dt1 + wy*dt2
 f3 = wx*dt2 + wy*dt3
 f4 = wx*(dt4+dut1) + wy*(dt5+dut2)
 
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  ip1 = grp%boundaryFaceNodeMappings(i1)
  ip2 = grp%boundaryFaceNodeMappings(i2)
 else
  ip1 = i1
  ip2 = i2
 end if


 if(ivd%boundaryTerm==1) then 
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + f2 + 2.0*(wx*t11 + wy*t21) 
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + f3 + 2.0*(wx*t21 + wy*t31) 
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + f4 + 2.0*(wx*(t41+dut11) + wy*(t51+dut21)) 

  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + f2 + 2.0*(wx*t12 + wy*t22) 
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + f3 + 2.0*(wx*t22 + wy*t32) 
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + f4 + 2.0*(wx*(t42+dut12) + wy*(t52+dut22)) 
 else
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + 4.0*(wx*t11 + wy*t21)
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + 4.0*(wx*t21 + wy*t31)
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + 4.0*(wx*(t41+dut11) + wy*(t51+dut21))

  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + 4.0*(wx*t12 + wy*t22)
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) +  4.0*(wx*t22 + wy*t32)
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) +  4.0*(wx*(t42+dut12) + wy*(t52+dut22))
 end if
end do

!debugging, print viscid fluxes
!open(13,file='Vflux.dat')
!do i=1,grp%numberOfNodes
!    write(13,'(I5,3(E17.8))')i,grp%nodeHelpArray(i,9),grp%nodeHelpArray(i,10),grp%nodeHelpArray(i,11)
!enddo
!close(13)
! add viscosity to dissipation vector

!grp%dissipation = 0.0
do i=1,grp%numberOfNodes
 rhs(i,2:4) = rhs(i,2:4) - grp%nodeHelpArray(i,9:11)
end do
end subroutine makeViscosityTerm
!-----------------------------------------------------------------------
subroutine makeViscosityTerm2(grp,ivd,rhs,recalculateTurbulence)
! makes the viscosity fluxes
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: rhs(:,:)
logical :: recalculateTurbulence

integer :: i1,i2,i,ip1,ip2
real :: r1,r2,u1,u2,v1,v2,t1,t2,wx,wy,ux,vx,tx,uy,vy,ty,oneOverReynoldsNumber
real :: mu,k,T,T0,f2,f3,f4,rx,ry,vt,t11,t21,t31,t41,t51,t12,t22,t32,Ksi
real :: t42,t52,dudx,dvdx,dTdx,dudy,dvdy,dTdy,dt1,dt2,dt3,dt4,dt5,oneOverPrandtlNumber
real :: twoThirds,fourThirds,gammaOverGammaMinusOne,u11,u12,u21,u22,dut1,dut2
real :: wallTangent(2),wallNormal(2),tempConv,dut11,dut12,dut21,dut22,inflowTemp
real :: adTemp1,adTemp2,turbViscosityCoefficient,turbDiffusionCoefficient
real :: viscosityCoefficient,diffusionCoefficient,oneOverTurbPrandtlNumber
integer :: wallInd

double precision :: buff,fact

double precision :: dudx1,dudx2,dudy1,dudy2,dvdx1,dvdx2,dvdy1,dvdy2,dTdx1,dTdx2
double precision :: dTdy1,dTdy2,mmu1,mmu2,k1,k2,dudr,dvdr
double precision :: dTdr,rhoinv1,rhoinv2,mu1,mu2,mv1,mv2,kT1,kT2,r(2),inverseSideLength
double precision :: du,dv,dT,c1,c2,c3

grp%nodeHelpArray(1:grp%numberOfNodes,1:6) = 0.0

gammaOverGammaMinusOne = ivd%gamma/(ivd%gamma-1.0)
tempConv = (ivd%gamma-1.0)*ivd%MachNumber**2
inflowTemp = ivd%inflowTemperature


do i=1,grp%numberOfSides
! indexes of nodes in side
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor
 wx = grp%sideWeightsArray(i,1)
 wy = grp%sideWeightsArray(i,2)

 r1   = grp%u(i1,1)                                  ! density
 u1   = grp%u(i1,2)/r1                               ! x - velocity
 v1   = grp%u(i1,3)/r1                               ! y - velocity
 t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
 r2   = grp%u(i2,1)
 u2   = grp%u(i2,2)/r2
 v2   = grp%u(i2,3)/r2
 t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
 ux   = (u1+u2)*wx
 vx   = (v1+v2)*wx
 tx   = (t1+t2)*wx
 uy   = (u1+u2)*wy
 vy   = (v1+v2)*wy
 ty   = (t1+t2)*wy

 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + ux
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + vx
 grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + tx
 grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + uy
 grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + vy
 grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + ty

 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - ux
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - vx
 grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - tx
 grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - uy
 grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - vy
 grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - ty
end do

if(ivd%boundaryTerm==1) then 
 do i=1,grp%brp%numberOfBoundaryFaces
  i1 = grp%brp%faceIndexArray(i,1)
  i2 = grp%brp%faceIndexArray(i,2)
  wx   = grp%brp%faceWeightsArray(i,1)
  wy   = grp%brp%faceWeightsArray(i,2)

  r1   = grp%u(i1,1)
  u1   = grp%u(i1,2)/r1
  v1   = grp%u(i1,3)/r1
  t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
  r2   = grp%u(i2,1)
  u2   = grp%u(i2,2)/r2
  v2   = grp%u(i2,3)/r2
  t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
  ux   = u1+u2
  vx   = v1+v2
  tx   = t1+t2
  uy   = u1+u2
  vy   = v1+v2
  ty   = t1+t2

! adding boundary face contributions

  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*u1+ux)*wx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*v1+vx)*wx
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*t1+tx)*wx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*u1+uy)*wy
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*v1+vy)*wy
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*t1+ty)*wy

  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*u2+ux)*wx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*v2+vx)*wx
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*t2+tx)*wx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*u2+uy)*wy
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*v2+vy)*wy
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*t2+ty)*wy
 end do
else
! traditional FV boundary term
 do i=1,grp%brp%numberOfBoundaryFaces
  i1 = grp%brp%faceIndexArray(i,1)
  i2 = grp%brp%faceIndexArray(i,2)
  wx   = grp%brp%faceWeightsArray(i,1)
  wy   = grp%brp%faceWeightsArray(i,2)

  r1   = grp%u(i1,1)
  u1   = grp%u(i1,2)/r1
  v1   = grp%u(i1,3)/r1
  t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
  r2   = grp%u(i2,1)
  u2   = grp%u(i2,2)/r2
  v2   = grp%u(i2,3)/r2
  t2 = gammaOverGammaMinusOne*grp%p(i2)/r2

! adding boundary face contributions

  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + 4.*u1*wx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + 4.*v1*wx
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + 4.*t1*wx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + 4.*u1*wy
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + 4.*v1*wy
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + 4.*t1*wy

  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + 4.*u2*wx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + 4.*v2*wx
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + 4.*t2*wx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + 4.*u2*wy
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + 4.*v2*wy
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + 4.*t2*wy
 end do
end if

twoThirds = 2./3.
fourThirds = 4./3.
T0 = 1./((ivd%gamma-1)*ivd%MachNumber**2)
oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
oneOverPrandtlNumber = 1./ivd%PrandtlNumber
oneOverTurbPrandtlNumber = 1./ivd%turbulentPrandtlNumber
! divide by volume to remove mass matrix on LHS
do i=1,grp%numberOfNodes
 vt = 1./grp%nodeVolume(i)
 grp%nodeHelpArray(i,1) = grp%nodeHelpArray(i,1)*vt
 grp%nodeHelpArray(i,2) = grp%nodeHelpArray(i,2)*vt
 grp%nodeHelpArray(i,3) = grp%nodeHelpArray(i,3)*vt
 grp%nodeHelpArray(i,4) = grp%nodeHelpArray(i,4)*vt
 grp%nodeHelpArray(i,5) = grp%nodeHelpArray(i,5)*vt
 grp%nodeHelpArray(i,6) = grp%nodeHelpArray(i,6)*vt

 dudx = grp%nodeHelpArray(i,1)
 dvdx = grp%nodeHelpArray(i,2)
 dTdx = grp%nodeHelpArray(i,3)
 dudy = grp%nodeHelpArray(i,4)
 dvdy = grp%nodeHelpArray(i,5)
 dTdy = grp%nodeHelpArray(i,6)

! create the viscosity matrix

 T = gammaOverGammaMinusOne*grp%p(i)/grp%u(i,1)
 if(T<0) then
  write(*,*) "T negative ",T," at node ",i,grp%u(i,1),grp%p(i)
  T = abs(T)
 end if


 ! Sunderlands law of viscosity
 mu = oneOverReynoldsNumber*(((tempConv*T)**1.5)*((inflowTemp+198.6)/(tempConv*inflowTemp*T+198.6)))
! mu = oneOverReynoldsNumber
 k = oneOverPrandtlNumber*mu
 grp%laminarViscosity(i) = mu
 grp%vorticity(i) = sqrt((dudy-dvdx)**2) ! vorticity array
 grp%divergence(i) = dudx + dvdy

 ! add turbulence effects
 if(ivd%turbulenceModel==1) then ! Spalart-Allmaras
  turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i,5)

  turbViscosityCoefficient = grp%u(i,1)*turbViscosityCoefficient  ! KAS

  Ksi = turbViscosityCoefficient/mu
  turbViscosityCoefficient =&
  turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)

  mu = mu + turbViscosityCoefficient
  k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber
 end if

 grp%nodeHelpArray(i,1:2) = grp%nodeHelpArray(i,1:2)*mu
 grp%nodeHelpArray(i,4:5) = grp%nodeHelpArray(i,4:5)*mu
 grp%nodeHelpArray(i,3) = grp%nodeHelpArray(i,3)*k
 grp%nodeHelpArray(i,6) = grp%nodeHelpArray(i,6)*k

 grp%nodeHelpArray(i,7) = mu
 grp%nodeHelpArray(i,8) = k
end do


! for adiabatic wall, set temperature gradient equal to zero
! call setTemperatureForAdiabaticWall(grp,grp%nodeHelpArray(:,4:5))
 call setTemperatureForAdiabaticWall(grp,grp%nodeHelpArray(:,3),grp%nodeHelpArray(:,6))

! make stress tensor along boundary

if(.not.associated(grp%sourceTerm)) then
 do i=1,grp%brp%numberOfBoundaryFaces
  wallTangent(1) = grp%brp%faceTangentArray(i,1)
  wallTangent(2) = grp%brp%faceTangentArray(i,2)
  wallNormal(1) = wallTangent(2)
  wallNormal(2) = -1.0*wallTangent(1)
 
  c1 = fourThirds*grp%nodeHelpArray(i,1)-twoThirds*grp%nodeHelpArray(i,5)
  c2 = grp%nodeHelpArray(i,4)+grp%nodeHelpArray(i,2)
  c3 = fourThirds*grp%nodeHelpArray(i,5)-twoThirds*grp%nodeHelpArray(i,1)

  grp%wallStress(i,1) = c1*wallNormal(1)+c2*wallNormal(2)
  grp%wallStress(i,2) = c2*wallNormal(1)+c3*wallNormal(2)
 end do
end if

grp%nodeHelpArray(1:grp%numberOfNodes,9:11) = 0.0

buff = 0.0


! KKK

 do i=1,grp%brp%numberOfBoundaryFaces
  grp%nodeHelpArray(i,3) = 0.0  
  grp%nodeHelpArray(i,6) = 0.0  
 end do

! KKK


do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! if(grp%gridNumber==1) then 
  inverseSideLength = 1.0/grp%sideLengthArray(i,3)
! else
!  inverseSideLength = 1.0/grp%sideLengthArray(i,3)
! end if


! side weights from preprocessor
 wx = grp%sideWeightsArray(i,1)
 wy = grp%sideWeightsArray(i,2)


 mmu1 = grp%nodeHelpArray(i1,7)
 mmu2 = grp%nodeHelpArray(i2,7)

 k1 = grp%nodeHelpArray(i1,8)
 k2 = grp%nodeHelpArray(i2,8)


 rhoinv1 = 1.0/grp%u(i1,1)
 u1 = grp%u(i1,2)*rhoinv1
 v1 = grp%u(i1,3)*rhoinv1
 mu1 = mmu1*u1
 mv1 = mmu1*v1
 T1 = gammaOverGammaMinusOne*grp%p(i1)*rhoinv1
 kT1 = k1*T1

 rhoinv2 = 1.0/grp%u(i2,1)
 u2 = grp%u(i2,2)*rhoinv2
 v2 = grp%u(i2,3)*rhoinv2
 mu2 = mmu2*u2
 mv2 = mmu2*v2
 T2 = gammaOverGammaMinusOne*grp%p(i2)*rhoinv2
 kT2 = k2*T2

! if(grp%gridNumber==1) then 
!  r(1) = grp%coordinates(i2,1)-grp%coordinates(i1,1)
!  r(2) = grp%coordinates(i2,2)-grp%coordinates(i1,2)
! else
  r(1) = grp%sideLengthArray(i,1)
  r(2) = grp%sideLengthArray(i,2)
! end if
 r = r*inverseSideLength

 dudx1 = grp%nodeHelpArray(i1,1)
 dudy1 = grp%nodeHelpArray(i1,4)
 dudx2 = grp%nodeHelpArray(i2,1)
 dudy2 = grp%nodeHelpArray(i2,4)
 dvdx1 = grp%nodeHelpArray(i1,2)
 dvdy1 = grp%nodeHelpArray(i1,5)
 dvdx2 = grp%nodeHelpArray(i2,2)
 dvdy2 = grp%nodeHelpArray(i2,5)
 dTdx1 = grp%nodeHelpArray(i1,3)
 dTdy1 = grp%nodeHelpArray(i1,6)
 dTdx2 = grp%nodeHelpArray(i2,3)
 dTdy2 = grp%nodeHelpArray(i2,6)

 dudx = dudx1 + dudx2
 dudy = dudy1 + dudy2
 dvdx = dvdx1 + dvdx2
 dvdy = dvdy1 + dvdy2
 dTdx = dTdx1 + dTdx2
 dTdy = dTdy1 + dTdy2

! remove derivative component along side

 dudr = dudx*r(1)+dudy*r(2)
 dvdr = dvdx*r(1)+dvdy*r(2)
 dTdr = dTdx*r(1)+dTdy*r(2)

 dudx = dudx - dudr*r(1)
 dudy = dudy - dudr*r(2)
 dvdx = dvdx - dvdr*r(1)
 dvdy = dvdy - dvdr*r(2)
! dTdx = dTdx - dTdr*r(1)
! dTdy = dTdy - dTdr*r(2)

! add finite difference along edge

 du = (mmu1+mmu2)*(u2-u1)*inverseSideLength  ! factor of 2 since wx,wy have factor of 0.5
 dv = (mmu1+mmu2)*(v2-v1)*inverseSideLength
 dT = (k1+k2)*(T2-T1)*inverseSideLength

 dudx = dudx + du*r(1)
 dudy = dudy + du*r(2)
 dvdx = dvdx + dv*r(1)
 dvdy = dvdy + dv*r(2)
 dTdx = dTdx + dT*r(1)
 dTdy = dTdy + dT*r(2)

! stress tensor

 dt1 = fourThirds*dudx-twoThirds*dvdy
 dt2 = dudy+dvdx
 dt3 = fourThirds*dvdy-twoThirds*dudx
 dt4 = dTdx
 dt5 = dTdy

! if(i1<grp%brp%numberOfBoundaryFaces.or.i2<grp%brp%numberOfBoundaryFaces) then 
!  dt4 = 0.0
!  dt5 = 0.0
! end if


! dissipation terms

 dut1 = 0.5*(u1+u2)*dt1+0.5*(v1+v2)*dt2 
 dut2 = 0.5*(u1+u2)*dt2+0.5*(v1+v2)*dt3

 f2 = wx*dt1 + wy*dt2
 f3 = wx*dt2 + wy*dt3 
 f4 = wx*(dt4+dut1) + wy*(dt5+dut2)

 grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + f2
 grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + f3
 grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + f4

 grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - f2
 grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - f3
 grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - f4

end do


if(ivd%boundaryTerm==1) then 
 do i=1,grp%brp%numberOfBoundaryFaces
  i1 = grp%brp%faceIndexArray(i,1)
  i2 = grp%brp%faceIndexArray(i,2)


  STOP "ERROR: Wrong boundary term"

  inverseSideLength = 1.0/grp%sideLengthArray(i,3)


  wx   = grp%brp%faceWeightsArray(i,1)
  wy   = grp%brp%faceWeightsArray(i,2)
 
  mmu1 = grp%nodeHelpArray(i1,7)
  mmu2 = grp%nodeHelpArray(i2,7)
 
  k1 = grp%nodeHelpArray(i1,8)
  k2 = grp%nodeHelpArray(i2,8)
 
 
  rhoinv1 = 1.0/grp%u(i1,1)
  u1 = grp%u(i1,2)*rhoinv1
  v1 = grp%u(i1,3)*rhoinv1
  mu1 = mmu1*u1
  mv1 = mmu1*v1
  T1 = gammaOverGammaMinusOne*grp%p(i1)*rhoinv1
  kT1 = k1*T1
 
  rhoinv2 = 1.0/grp%u(i2,1)
  u2 = grp%u(i2,2)*rhoinv2
  v2 = grp%u(i2,3)*rhoinv2
  mu2 = mmu2*u2
  mv2 = mmu2*v2
  T2 = gammaOverGammaMinusOne*grp%p(i2)*rhoinv2
  kT2 = k2*T2
 
  r = grp%coordinates(i2,:)-grp%coordinates(i1,:)
  r = r*inverseSideLength
 
  dudx1 = grp%nodeHelpArray(i1,1)
  dudy1 = grp%nodeHelpArray(i1,4)
  dudx2 = grp%nodeHelpArray(i2,1)
  dudy2 = grp%nodeHelpArray(i2,4)
  dvdx1 = grp%nodeHelpArray(i1,2)
  dvdy1 = grp%nodeHelpArray(i1,5)
  dvdx2 = grp%nodeHelpArray(i2,2)
  dvdy2 = grp%nodeHelpArray(i2,5)
  dTdx1 = grp%nodeHelpArray(i1,3)
  dTdy1 = grp%nodeHelpArray(i1,6)
  dTdx2 = grp%nodeHelpArray(i2,3)
  dTdy2 = grp%nodeHelpArray(i2,6)

  dudx = dudx1 + dudx2
  dudy = dudy1 + dudy2
  dvdx = dvdx1 + dvdx2
  dvdy = dvdy1 + dvdy2
  dTdx = dTdx1 + dTdx2
  dTdy = dTdy1 + dTdy2
 
 ! remove derivative component along side
 
  dudr = dudx*r(1)+dudy*r(2)
  dvdr = dvdx*r(1)+dvdy*r(2)
  dTdr = dTdx*r(1)+dTdy*r(2)
 
  dudx = dudx - dudr*r(1)
  dudy = dudy - dudr*r(2)
  dvdx = dvdx - dvdr*r(1)
  dvdy = dvdy - dvdr*r(2)
!  dTdx = dTdx - dTdr*r(1)
!  dTdy = dTdy - dTdr*r(2)
 
   
 
 ! add finite difference along edge
 
  du = (mmu1+mmu2)*(u2-u1)*inverseSideLength  
  dv = (mmu1+mmu2)*(v2-v1)*inverseSideLength
  dT = (k1+k2)*(T2-T1)*inverseSideLength
 
  dudx = dudx + du*r(1)
  dudy = dudy + du*r(2)
  dvdx = dvdx + dv*r(1)
  dvdy = dvdy + dv*r(2)
  dTdx = dTdx + dT*r(1)
  dTdy = dTdy + dT*r(2)
 
 ! stress tensor
 
  t11 = fourThirds*dudx1-twoThirds*dvdy1
  t21 = dudy1+dvdx1
  t31 = fourThirds*dvdy1-twoThirds*dudx1
  t41 = dTdx1
  t51 = dTdy1
 
  t12 = fourThirds*dudx2-twoThirds*dvdy2
  t22 = dudy2+dvdx2
  t32 = fourThirds*dvdy2-twoThirds*dudx2
  t42 = dTdx2
  t52 = dTdy2
 
  dt1 = fourThirds*dudx-twoThirds*dvdy
  dt2 = dudy+dvdx
  dt3 = fourThirds*dvdy-twoThirds*dudx
  dt4 = dTdx
  dt5 = dTdy
 
  ! isentropic boundary 
 
  dt4 = 0.0
  dt5 = 0.0
 
 ! dissipation terms
 
  dut11 = u1*t11+v1*t21
  dut12 = u2*t12+v2*t22
  dut21 = u1*t21+v1*t31
  dut22 = u2*t22+v2*t32
 
  dut1 = 0.5*(u1+u2)*dt1+0.5*(v1+v2)*dt2
  dut2 = 0.5*(u1+u2)*dt2+0.5*(v1+v2)*dt3
 
  f2 = wx*dt1 + wy*dt2 
  f3 = wx*dt2 + wy*dt3 
  f4 = wx*(dt4+dut1) + wy*(dt5+dut2) 
 
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + f2 + 2.0*(wx*t11+wy*t21)
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + f3 + 2.0*(wx*t21+wy*t41)
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + f4 + 2.0*(wx*(t41+dut11) + wy*(t51+dut21))
 
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + f2 + 2.0*(wx*t12+wy*t22)
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + f3 + 2.0*(wx*t22+wy*t42)
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + f4 + 2.0*(wx*(t42+dut12) + wy*(t52+dut22))
 end do
else
 do i=1,grp%brp%numberOfBoundaryFaces
  i1 = grp%brp%faceIndexArray(i,1)
  i2 = grp%brp%faceIndexArray(i,2)

  wx   = grp%brp%faceWeightsArray(i,1)
  wy   = grp%brp%faceWeightsArray(i,2)

  mmu1 = grp%nodeHelpArray(i1,7)
  mmu2 = grp%nodeHelpArray(i2,7)

  k1 = grp%nodeHelpArray(i1,8)
  k2 = grp%nodeHelpArray(i2,8)


  rhoinv1 = 1.0/grp%u(i1,1)
  u1 = grp%u(i1,2)*rhoinv1
  v1 = grp%u(i1,3)*rhoinv1
  mu1 = mmu1*u1
  mv1 = mmu1*v1
  kT1 = k1*gammaOverGammaMinusOne*grp%p(i1)*rhoinv1

  rhoinv2 = 1.0/grp%u(i2,1)
  u2 = grp%u(i2,2)*rhoinv2
  v2 = grp%u(i2,3)*rhoinv2
  mu2 = mmu2*u2
  mv2 = mmu2*v2
  kT2 = k2*gammaOverGammaMinusOne*grp%p(i2)*rhoinv2

  dudx1 = grp%nodeHelpArray(i1,1)
  dudy1 = grp%nodeHelpArray(i1,4)
  dudx2 = grp%nodeHelpArray(i2,1)
  dudy2 = grp%nodeHelpArray(i2,4)
  dvdx1 = grp%nodeHelpArray(i1,2)
  dvdy1 = grp%nodeHelpArray(i1,5)
  dvdx2 = grp%nodeHelpArray(i2,2)
  dvdy2 = grp%nodeHelpArray(i2,5)
  dTdx1 = grp%nodeHelpArray(i1,3)
  dTdy1 = grp%nodeHelpArray(i1,6)
  dTdx2 = grp%nodeHelpArray(i2,3)
  dTdy2 = grp%nodeHelpArray(i2,6)

 ! stress tensor

  t11 = fourThirds*dudx1-twoThirds*dvdy1
  t21 = dudy1+dvdx1
  t31 = fourThirds*dvdy1-twoThirds*dudx1
  t41 = dTdx1
  t51 = dTdy1

  t12 = fourThirds*dudx2-twoThirds*dvdy2
  t22 = dudy2+dvdx2
  t32 = fourThirds*dvdy2-twoThirds*dudx2
  t42 = dTdx2
  t52 = dTdy2

  ! isentropic boundary

  t41 = 0.0
  t51 = 0.0
  t42 = 0.0
  t52 = 0.0

 ! dissipation terms

  dut11 = u1*t11+v1*t21
  dut12 = u2*t12+v2*t22
  dut21 = u1*t21+v1*t31
  dut22 = u2*t22+v2*t32

  dut1 = 0.5*(u1+u2)*dt1+0.5*(v1+v2)*dt2
  dut2 = 0.5*(u1+u2)*dt2+0.5*(v1+v2)*dt3

  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + 4.0*(wx*t11+wy*t21)
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + 4.0*(wx*t21+wy*t41)
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + 4.0*(wx*(t41+dut11) + wy*(t51+dut21))

  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + 4.0*(wx*t12+wy*t22)
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + 4.0*(wx*t22+wy*t42)
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + 4.0*(wx*(t42+dut12) + wy*(t52+dut22))
 end do
end if

! add viscosity to dissipation vector

!grp%dissipation = 0.0
if(grp%gridNumber>1) then 
 fact = 0.5**(grp%gridNumber-1)
 grp%nodeHelpArray(:,9:11) = fact*grp%nodeHelpArray(:,9:11)
end if
do i=1,grp%numberOfNodes
 rhs(i,2:4) = rhs(i,2:4) - grp%nodeHelpArray(i,9:11)
end do
end subroutine makeViscosityTerm2
!-----------------------------------------------------------------------
subroutine makeSARHS(grp,ivd)
IMPLICIT NONE
! Spalart-Allmaras turbulence model

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 real :: turbVisc1,turbVisc2,ru1,rv1,u1,v1,rinv1,rinv2,k
 real :: ru2,rv2,u2,v2,wx,wy,fconv,ax,ay,az,anx,any
 real :: turbVisc0,turbVisc3,turbVisc4,turbViscm
 real :: c0,c1,c2,c3,ym1,y0,y1,y2,y3,y4,y5,xm1,x0,x1,x2,x3,x4,x5
 real :: rinv0,ru0,rv0,u0,v0,h0,h1,h2,h3,h4
 real :: rinv3,ru3,rv3,u3,v3,w1,w2,w3,w4
 real :: rinv4,ru4,rv4,u4,v4,a0,a1,a2,a3,a4,ia0,ia1,ia2,ia3,ia4,amb,am
 real :: ut0,ut1,ut2,ut3,ut4,un0,un1,un2,un3,un4
 real :: unm,utm,uxm,vym
 integer :: i,i1,i2,ind,j1,j2,j3,j0,j4,jm1,j5

 ind = 39498
k=1.0
 do i=1,grp%numberOfSides
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)
 
! side weights from preprocessor

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  anx = wx/sqrt(wx**2+wy**2)
  any = wy/sqrt(wx**2+wy**2)

  rinv1 = 1./grp%u(i1,1)                 ! rho inverse
  turbVisc1 = grp%u(i1,5)
  ru1 = grp%u(i1,2)                      ! x-momentum
  rv1 = grp%u(i1,3)                      ! y-momentum
  u1 = rinv1*ru1                           ! x-velocity
  v1 = rinv1*rv1                           ! y-velocity
  ut1 = u1*any-v1*anx
  un1 = u1*anx+v1*any

  rinv2 = 1./grp%u(i2,1)
  turbVisc2 = grp%u(i2,5)
  ru2 = grp%u(i2,2)
  rv2 = grp%u(i2,3)
  u2 = rinv2*ru2
  v2 = rinv2*rv2
  ut2 = u2*any-v2*anx
  un2 = u2*anx+v2*any

  ! convection terms

  if(.false.) then!ivd%HighOrder.and.grp%gridNumber==1) then
    if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
       if(grp%AfterAndBefore(i1,1).eq.i2) then
         j1 = i1
         j2 = i2
         j3 = grp%AfterAndBefore(j2,1)
         j4 = grp%AfterAndBefore(j3,1)
         j0 = grp%AfterAndBefore(j1,2)
         jm1 = grp%AfterAndBefore(j0,2)
       else if(grp%AfterAndBefore(i1,2).eq.i2) then
         j1 = i2
         j2 = i1
         j3 = grp%AfterAndBefore(j2,1)
         j4 = grp%AfterAndBefore(j3,1)
         j0 = grp%AfterAndBefore(j1,2)
         jm1 = grp%AfterAndBefore(j0,2)
       else
         print *,' should not be here'
       end if
       if(j0.eq.0) then
         !jm1 = 1
         !j5 = grp%AfterAndBefore(j4,1)

         rinv1 = 1./grp%u(j1,1)
         turbVisc1 = grp%u(j1,5)
         ru1 = grp%u(j1,2)
         rv1 = grp%u(j1,3)
         u1 = rinv1*ru1
         v1 = rinv1*rv1
         ut1 = u1*any-v1*anx
         un1 = u1*anx+v1*any
         a1 = grp%nodeVolume(j1)

         rinv2 = 1./grp%u(j2,1)
         turbVisc2 = grp%u(j2,5)
         ru2 = grp%u(j2,2)
         rv2 = grp%u(j2,3)
         u2 = rinv2*ru2
         v2 = rinv2*rv2
         ut2 = u2*any-v2*anx
         un2 = u2*anx+v2*any
         a2 = grp%nodeVolume(j2)

         rinv3 = 1./grp%u(j3,1)
         turbVisc3 = grp%u(j3,5)
         ru3 = grp%u(j3,2)
         rv3 = grp%u(j3,3)
         u3 = rinv3*ru3
         v3 = rinv3*rv3
         ut3 = u3*any-v3*anx
         un3 = u3*anx+v3*any
         a3 = grp%nodeVolume(j3)

         rinv4 = 1./grp%u(j4,1)
         turbVisc4 = grp%u(j4,5)
         ru4 = grp%u(j4,2)
         rv4 = grp%u(j4,3)
         u4 = rinv4*ru4
         v4 = rinv4*rv4
         ut4 = u4*any-v4*anx
         un4 = u4*anx+v4*any
         a4 = grp%nodeVolume(j4)
         !debugging
        !if (k==1.0) then
        !  a1=0.05
        !  a2=0.1
        !  a3=0.1
        !  a4=0.1
          !x^3-2*x
        !  turbVisc1=-199.0/1000.0 !point value at x=0.1
        !  turbVisc2=-783.0/2000.0 !mean values
        !  turbVisc3=-2289.0/4000.0
        !  turbVisc4=-147.0/200.0
        !endif


         !Approximate Jacobian
         ia1=a1
         ia2=a2+ia1
         ia3=a3+ia2
         ia4=a4+ia3

         c0 = grp%hoc2(i,1)
         c1 = grp%hoc2(i,2)
         c2 = grp%hoc2(i,3)
         c3 = grp%hoc2(i,4)
         !Jacobian at interpolation location        
         am=c0*ia1+c1*ia2+c2*ia3+c3*ia4

         c0 = grp%hoc3(i,1)
         c1 = grp%hoc3(i,2)
         c2 = grp%hoc3(i,3)
         c3 = grp%hoc3(i,4)
         !jacobian at boundary        
         amb=c0*ia1+c1*ia2+c2*ia3+c3*ia4

         c0 = grp%hoc(i,1)
         c1 = grp%hoc(i,2)
         c2 = grp%hoc(i,3)
         c3 = grp%hoc(i,4)

        !Primitive function

         ut1 = ut1*amb
         ut2 = ut2*a2
         ut3 = ut2+ut3*a3
         ut4 = ut3+ut4*a4

         un1 = un1*amb
         un2 = un2*a2
         un3 = un2+un3*a3
         un4 = un3+un4*a4

         turbVisc1 = turbVisc1*amb
         turbVisc2 = turbVisc2*a2
         turbVisc3 = turbVisc2+turbVisc3*a3
         turbVisc4 = turbVisc3+turbVisc4*a4

         turbViscm = c0*turbVisc1+c1*turbVisc2+c2*turbVisc3+c3*turbVisc4
         utm = c0*ut1+c1*ut2+c2*ut3+c3*ut4
         unm = c0*un1+c1*un2+c2*un3+c3*un4

        !divide through by Jacobian
         turbViscm = turbViscm/am
         utm = utm/am
         unm = unm/am

        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',turbViscm
        !  k=0.0
          !should be -0.296625
        !endif

         uxm = utm*any+unm*anx
         vym =-utm*anx+unm*any

         fconv = wx*(2.*uxm*turbViscm)+wy*(2.*vym*turbViscm)

       else if(jm1.eq.0) then
         rinv0 = 1./grp%u(j0,1)
         turbVisc0 = grp%u(j0,5)
         ru0 = grp%u(j0,2)
         rv0 = grp%u(j0,3)
         u0 = rinv0*ru0
         v0 = rinv0*rv0
         ut0 = u0*any-v0*anx
         un0 = u0*anx+v0*any
         a0 = grp%nodeVolume(j0)

         rinv1 = 1./grp%u(j1,1)
         turbVisc1 = grp%u(j1,5)
         ru1 = grp%u(j1,2)
         rv1 = grp%u(j1,3)
         u1 = rinv1*ru1
         v1 = rinv1*rv1
         ut1 = u1*any-v1*anx
         un1 = u1*anx+v1*any
         a1 = grp%nodeVolume(j1)
 
         rinv2 = 1./grp%u(j2,1)
         turbVisc2 = grp%u(j2,5)
         ru2 = grp%u(j2,2)
         rv2 = grp%u(j2,3)
         u2 = rinv2*ru2
         v2 = rinv2*rv2
         ut2 = u2*any-v2*anx
         un2 = u2*anx+v2*any
         a2 = grp%nodeVolume(j2)

         rinv3 = 1./grp%u(j3,1)
         turbVisc3 = grp%u(j3,5)
         ru3 = grp%u(j3,2)
         rv3 = grp%u(j3,3)
         u3 = rinv3*ru3
         v3 = rinv3*rv3
         ut3 = u3*any-v3*anx
         un3 = u3*anx+v3*any
         a3 = grp%nodeVolume(j3)

         !debugging
        !if (k==1.0) then
        !  a0=0.05
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
          !x^3-2*x
        !  un0=-199.0/1000.0 !point value at x=0.1
        !  un1=-783.0/2000.0 !mean values
        !  un2=-2289.0/4000.0
        !  un3=-147.0/200.0
        !endif

        !Approximate Jacobian
        ia0=a0
        ia1=a1+ia0
        ia2=a2+ia1
        ia3=a3+ia2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location          
        am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

        c0 = grp%hoc3(i,1)
        c1 = grp%hoc3(i,2)
        c2 = grp%hoc3(i,3)
        c3 = grp%hoc3(i,4)
        !Jacobian at boundary        
        amb=c0*ia0+c1*ia1+c2*ia2+c3*ia3

        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)

        !Primitive function

        ut0 = ut0*amb
        ut1 = ut1*a1
        ut2 = ut1+ut2*a2
        ut3 = ut2+ut3*a3

        un0 = un0*amb
        un1 = un1*a1
        un2 = un1+un2*a2
        un3 = un2+un3*a3

        turbVisc0 = turbVisc0*amb
        turbVisc1 = turbVisc1*a1
        turbVisc2 = turbVisc1+turbVisc2*a2
        turbVisc3 = turbVisc2+turbVisc3*a3

         turbViscm = c0*turbVisc0+c1*turbVisc1+c2*turbVisc2+c3*turbVisc3
         utm = c0*ut0+c1*ut1+c2*ut2+c3*ut3
         unm = c0*un0+c1*un1+c2*un2+c3*un3

        !divide through by Jacobian
         turbViscm = turbViscm/am
         utm = utm/am
         unm = unm/am

        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',unm
        !  k=0.0
          !should be -0.484375
        !endif

         uxm = utm*any+unm*anx
         vym =-utm*anx+unm*any

         fconv = wx*(2.*uxm*turbViscm)+wy*(2.*vym*turbViscm)

       else if(j4.eq.0) then

         fconv = wx*(u1*turbVisc1+u2*turbVisc2)+wy*(v1*turbVisc1+v2*turbVisc2)

       else

         rinv0 = 1./grp%u(j0,1)
         turbVisc0 = grp%u(j0,5)
         ru0 = grp%u(j0,2)
         rv0 = grp%u(j0,3)
         u0 = rinv0*ru0
         v0 = rinv0*rv0
         ut0 = u0*any-v0*anx
         un0 = u0*anx+v0*any
         a0 = grp%nodeVolume(j0)

         rinv1 = 1./grp%u(j1,1)
         turbVisc1 = grp%u(j1,5)
         ru1 = grp%u(j1,2)
         rv1 = grp%u(j1,3)
         u1 = rinv1*ru1
         v1 = rinv1*rv1
         ut1 = u1*any-v1*anx
         un1 = u1*anx+v1*any
         a1 = grp%nodeVolume(j1)

         rinv2 = 1./grp%u(j2,1)
         turbVisc2 = grp%u(j2,5)
         ru2 = grp%u(j2,2)
         rv2 = grp%u(j2,3)
         u2 = rinv2*ru2
         v2 = rinv2*rv2
         ut2 = u2*any-v2*anx
         un2 = u2*anx+v2*any
         a2 = grp%nodeVolume(j2)

         rinv3 = 1./grp%u(j3,1)
         turbVisc3 = grp%u(j3,5)
         ru3 = grp%u(j3,2)
         rv3 = grp%u(j3,3)
         u3 = rinv3*ru3
         v3 = rinv3*rv3
         ut3 = u3*any-v3*anx
         un3 = u3*anx+v3*any
         a3 = grp%nodeVolume(j3)


       !debugging
        !if (k==1.0) then
        !  a0=0.1
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
        !  turbVisc0=-159.0/800.0
        !  turbVisc1=-783.0/2000.0
        !  turbVisc2=-2289.0/4000.0
        !  turbVisc3=-147.0/200.0
        ! endif

        !Approximate Jacobian
        ia0=a0
        ia1=a1+ia0
        ia2=a2+ia1
        ia3=a3+ia2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location          
        am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

         c0 = grp%hoc(i,1)
         c1 = grp%hoc(i,2)
         c2 = grp%hoc(i,3)
         c3 = grp%hoc(i,4)

        !Primitive function
 
        ut0 = ut0*a0
        ut1 = ut0+ut1*a1
        ut2 = ut1+ut2*a2
        ut3 = ut2+ut3*a3

        un0 = un0*a0
        un1 = un0+un1*a1
        un2 = un1+un2*a2
        un3 = un2+un3*a3

        turbVisc0 = turbVisc0*a0
        turbVisc1 = turbVisc0+turbVisc1*a1
        turbVisc2 = turbVisc1+turbVisc2*a2
        turbVisc3 = turbVisc2+turbVisc3*a3

         turbViscm = c0*turbVisc0+c1*turbVisc1+c2*turbVisc2+c3*turbVisc3
         utm = c0*ut0+c1*ut1+c2*ut2+c3*ut3
         unm = c0*un0+c1*un1+c2*un2+c3*un3

        !divide through by Jacobian
         turbViscm = turbViscm/am
         utm = utm/am
         unm = unm/am

       !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',turbViscm
        !  k=0.0
          !should give -0.484375
        !endif

         uxm = utm*any+unm*anx
         vym =-utm*anx+unm*any

         fconv = wx*(2.*uxm*turbViscm)+wy*(2.*vym*turbViscm)

       end if

    else

       fconv = wx*(u1*turbVisc1+u2*turbVisc2)+wy*(v1*turbVisc1+v2*turbVisc2)

    end if

  else

       fconv = wx*(u1*turbVisc1+u2*turbVisc2)+wy*(v1*turbVisc1+v2*turbVisc2)

  end if

! conservative:
!  fconv = wx*(u1*turbVisc1 + u2*turbVisc2) + wy*(v1*turbVisc1 + v2*turbVisc2) 
  if(abs(turbvisc1-turbvisc2)>1.0e-20) then
   ax = (u2*turbVisc2 - u1*turbVisc1)/(turbvisc1-turbvisc2)
   ay = (v2*turbVisc2 - v1*turbVisc1)/(turbvisc1-turbvisc2)
   fconv = fconv + abs(wx*ax+wy*ay)*(turbvisc1-turbvisc2)
  end if

! nonconservative:

!  fconv = wx*(u1*turbVisc1 + u1*turbVisc2) + wy*(v1*turbVisc1 + v1*turbVisc2)
!  if(abs(turbvisc1-turbvisc2)>1.0e-20) then
!   ax = (u1*turbVisc2 - u1*turbVisc1)/(turbvisc1-turbvisc2)
!   ay = (v1*turbVisc2 - v1*turbVisc1)/(turbvisc1-turbvisc2)
!   fconv = fconv + abs(wx*ax+wy*ay)*(turbvisc1-turbvisc2)
!  end if



  grp%rhs(i1,5) = grp%rhs(i1,5) + fconv
  grp%rhs(i2,5) = grp%rhs(i2,5) - fconv

 end do
end subroutine makeSARHS
!-----------------------------------------------------------------------
subroutine makeSADiffusion(grp,ivd)
IMPLICIT NONE
! Spalart-Allmaras turbulence model

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,ind1,ind2,ind
 integer :: i1,i2,j1,j2,j3,j4,jm1,j0,j5
 real :: ru1,ru2,rv1,rv2,rw1,rw2,fact,vt,visc1,visc2,fdiff1,fdiff2,fdiff3,k
 real :: turbVisc1,turbVisc2,tvx1,tvx2,tvy1,tvy2,tvz1,tvz2,dnudx,dnudy
 real :: turbVisc0,turbVisc3,turbVisc4,turbViscm,a0,a1,a2,a3,a4,ia0,ia1,ia2,ia3,ia4,amb,am
 real :: c0,c1,c2,c3,ym1,y0,y1,y2,y3,y4,y5,xm1,x0,x1,x2,x3,x4,x5,h0,h1,h2,h3,h4,w1,w2,w3,w4
 real :: visc0,visc3,visc4,viscdnudxm,viscdnudym,fdiffm
 real :: sigma,cb2,wx,wy,wz,dnudx1,dnudx2,dnudy1,dnudy2
 real :: dnudx0,dnudx3,dnudx4,dnudy0,dnudy3,dnudy4
 real :: vecProd,oneOverReynoldsNumber
 real :: viscdnudx0,viscdnudx1,viscdnudx2,viscdnudx3,viscdnudx4
 real :: viscdnudy0,viscdnudy1,viscdnudy2,viscdnudy3,viscdnudy4

 ind = 408

 sigma = 2./3.
 cb2 = 0.622

 oneOverReynoldsNumber = 1./ivd%ReynoldsNumber

! diffusion terms
k=1.0
 grp%nodeHelpArray(1:grp%numberOfNodes,9:14) = 0.0
 do i=1,grp%numberOfSides
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  if(.false.)then!ivd%HighOrder.and.grp%gridNumber==1) then
    if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
       if(grp%AfterAndBefore(i1,1).eq.i2) then
         j1 = i1
         j2 = i2
         j3 = grp%AfterAndBefore(j2,1)
         j4 = grp%AfterAndBefore(j3,1)
         j0 = grp%AfterAndBefore(j1,2)
         jm1 = grp%AfterAndBefore(j0,2)
       else if(grp%AfterAndBefore(i1,2).eq.i2) then
         j1 = i2
         j2 = i1
         j3 = grp%AfterAndBefore(j2,1)
         j4 = grp%AfterAndBefore(j3,1)
         j0 = grp%AfterAndBefore(j1,2)
         jm1 = grp%AfterAndBefore(j0,2)
       else
         print *,' should not be here'
       end if
       if(j0.eq.0) then
         !jm1=1!required?
         !j5 = grp%AfterAndBefore(j4,1)
         turbVisc1 = grp%u(j1,5)
         a1 = grp%nodeVolume(j1)              
         turbVisc2 = grp%u(j2,5) 
         a2 = grp%nodeVolume(j2)            
         turbVisc3 = grp%u(j3,5)
         a3 = grp%nodeVolume(j3)            
         turbVisc4 = grp%u(j4,5) 
         a4 = grp%nodeVolume(j4)
         !debugging
        !if (k==1.0) then
        !  a1=0.05
        !  a2=0.1
        !  a3=0.1
        !  a4=0.1
          !x^3-2*x
        !  turbVisc1=-199.0/1000.0 !point value at x=0.1
        !  turbVisc2=-783.0/2000.0 !mean values
        !  turbVisc3=-2289.0/4000.0
        !  turbVisc4=-147.0/200.0
        !endif

         !Approximate Jacobian
        ia1=a1
        ia2=a2+ia1
        ia3=a3+ia2
        ia4=a4+ia3

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location        
        am=c0*ia1+c1*ia2+c2*ia3+c3*ia4

        c0 = grp%hoc3(i,1)
        c1 = grp%hoc3(i,2)
        c2 = grp%hoc3(i,3)
        c3 = grp%hoc3(i,4)
        !Jacobian at boundary        
        amb=c0*ia1+c1*ia2+c2*ia3+c3*ia4

        c0 = grp%hoc(i,1)
        c1 = grp%hoc(i,2)
        c2 = grp%hoc(i,3)
        c3 = grp%hoc(i,4)
        
         !primitive functions

         turbVisc1 = turbVisc1*amb
         turbVisc2 = turbVisc2*a2
         turbVisc3 = turbVisc2+turbVisc3*a3
         turbVisc4 = turbVisc3+turbVisc4*a4
 
         turbViscm = grp%hoc(i,1)*turbVisc1+grp%hoc(i,2)*turbVisc2+ &
                     grp%hoc(i,3)*turbVisc3+grp%hoc(i,4)*turbVisc4

         !divide through by Jacobian
         turbViscm = turbViscm/am
        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',turbViscm
        !  k=0.0
          !should be -0.296625
        !endif
      
         fdiff1 = (turbViscm+turbViscm)*wx
         fdiff2 = (turbViscm+turbViscm)*wy

       else if(jm1.eq.0) then
         turbVisc0 = grp%u(j0,5) 
         a0 = grp%nodeVolume(j0)
         turbVisc1 = grp%u(j1,5)  
         a1 = grp%nodeVolume(j1)          
         turbVisc2 = grp%u(j2,5)  
         a2 = grp%nodeVolume(j2)          
         turbVisc3 = grp%u(j3,5) 
         a3 = grp%nodeVolume(j3)          

         !debugging
        !if (k==1.0) then
        !  a0=0.05
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
          !x^3-2*x
        !  turbVisc0=-199.0/1000.0 !point value at x=0.1
        !  turbVisc1=-783.0/2000.0 !mean values
        !  turbVisc2=-2289.0/4000.0
        !  turbVisc3=-147.0/200.0
        !endif
    
        !Approximate Jacobian
         ia0=a0
         ia1=a1+ia0
         ia2=a2+ia1
         ia3=a3+ia2

         c0 = grp%hoc2(i,1)
         c1 = grp%hoc2(i,2)
         c2 = grp%hoc2(i,3)
         c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location          
         am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

         c0 = grp%hoc3(i,1)
         c1 = grp%hoc3(i,2)
         c2 = grp%hoc3(i,3)
         c3 = grp%hoc3(i,4)
        !Jacobian at boundary        
         amb=c0*ia0+c1*ia1+c2*ia2+c3*ia3

         c0 = grp%hoc(i,1)
         c1 = grp%hoc(i,2)
         c2 = grp%hoc(i,3)
         c3 = grp%hoc(i,4)
       
         !primitive functions

         turbVisc0 = turbVisc0*amb
         turbVisc1 = turbVisc1*a2
         turbVisc2 = turbVisc1+turbVisc2*a2
         turbVisc3 = turbVisc2+turbVisc3*a3
 
         turbViscm = grp%hoc(i,1)*turbVisc0+grp%hoc(i,2)*turbVisc1+ &
                     grp%hoc(i,3)*turbVisc2+grp%hoc(i,4)*turbVisc3

        !divide through by Jacobian
         turbViscm = turbViscm/am
        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',turbViscm
        !  k=0.0
          !should be -0.484375
        !endif
         
         fdiff1 = (turbViscm+turbViscm)*wx
         fdiff2 = (turbViscm+turbViscm)*wy


       else if(j4.eq.0) then

         turbVisc1 = grp%u(i1,5)              
         turbVisc2 = grp%u(i2,5)
         fdiff1 = (turbVisc1+turbVisc2)*wx
         fdiff2 = (turbVisc1+turbVisc2)*wy

       else

         turbVisc0 = grp%u(j0,5) 
         a0 = grp%nodeVolume(j0)           
         turbVisc1 = grp%u(j1,5)
         a1 = grp%nodeVolume(j1)             
         turbVisc2 = grp%u(j2,5)  
         a2 = grp%nodeVolume(j2)         
         turbVisc3 = grp%u(j3,5)   
         a3 = grp%nodeVolume(j3) 
 
       !debugging
        !if (k==1.0) then
        !  a0=0.1
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
        !  turbVisc0=-159.0/800.0
        !  turbVisc1=-783.0/2000.0
        !  turbVisc2=-2289.0/4000.0
        !  turbVisc3=-147.0/200.0
        ! endif

        !Approximate Jacobian
         ia0=a0
         ia1=a1+ia0
         ia2=a2+ia1
         ia3=a3+ia2

         c0 = grp%hoc2(i,1)
         c1 = grp%hoc2(i,2)
         c2 = grp%hoc2(i,3)
         c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location          
         am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

         c0 = grp%hoc(i,1)
         c1 = grp%hoc(i,2)
         c2 = grp%hoc(i,3)
         c3 = grp%hoc(i,4)
 
         !primitive functions        
         turbVisc0 = turbVisc0*a0
         turbVisc1 = turbVisc0+turbVisc1*a1
         turbVisc2 = turbVisc1+turbVisc2*a2
         turbVisc3 = turbVisc2+turbVisc3*a3

         turbViscm = grp%hoc(i,1)*turbVisc0+grp%hoc(i,2)*turbVisc1+ &
                     grp%hoc(i,3)*turbVisc2+grp%hoc(i,4)*turbVisc3

        !divide through by Jacobian
         turbViscm = turbViscm/am

       !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',turbViscm
        !  k=0.0
          !should give -0.484375
        !endif
        
         fdiff1 = (turbViscm+turbViscm)*wx
         fdiff2 = (turbViscm+turbViscm)*wy

       end if

    else

       turbVisc1 = grp%u(i1,5)
       turbVisc2 = grp%u(i2,5)
       fdiff1 = wx*(turbVisc1+turbVisc2)
       fdiff2 = wy*(turbVisc1+turbVisc2)

    end if

  else

       turbVisc1 = grp%u(i1,5)
       turbVisc2 = grp%u(i2,5)
       fdiff1 = wx*(turbVisc1+turbVisc2)
       fdiff2 = wy*(turbVisc1+turbVisc2)

  end if

  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + fdiff1
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - fdiff1
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + fdiff2
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - fdiff2
 end do

 ! boundary faces
 do i=1,grp%brp%numberOfBoundaryFaces
  ind1 = grp%brp%faceIndexArray(i,1)
  ind2 = grp%brp%faceIndexArray(i,2)
  wx   = grp%brp%faceWeightsArray(i,1)
  wy   = grp%brp%faceWeightsArray(i,2)

! KAS  turbVisc1 = grp%u(ind1,5)/grp%u(ind1,1)
! KAS  turbVisc2 = grp%u(ind2,5)/grp%u(ind2,1)

  turbVisc1 = grp%u(ind1,5)
  turbVisc2 = grp%u(ind2,5)

  if(ivd%boundaryTerm==1) then 
  ! adding boundary face contributions

   tvx1 = (3.*turbVisc1 + turbVisc2)*wx
   tvx2 = (turbVisc1 + 3.*turbVisc2)*wx
   tvy1 = (3.*turbVisc1 + turbVisc2)*wy
   tvy2 = (turbVisc1 + 3.*turbVisc2)*wy

   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + tvx1
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + tvx2
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + tvy1
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + tvy2
  else
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + 4.*turbVisc1*wx
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + 4.*turbVisc2*wx 
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + 4.*turbVisc1*wy
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + 4.*turbVisc2*wy
  end if
 end do

 ! get turbulent viscosity gradient
 fact = oneOverReynoldsNumber*cb2/sigma
 do i=1,grp%numberOfNodes
  vt = 1./grp%nodeVolume(i)
  dnudx = grp%nodeHelpArray(i,9)*vt
  dnudy = grp%nodeHelpArray(i,10)*vt
  grp%nodeHelpArray(i,12) = dnudx
  grp%nodeHelpArray(i,13) = dnudy

  vecProd = grp%nodeHelpArray(i,9)*dnudx + grp%nodeHelpArray(i,10)*dnudy 

  ! add vector producs
! KAS  grp%rhs(i,5) = grp%rhs(i,5) - fact*grp%u(i,1)*vecProd
  grp%rhs(i,5) = grp%rhs(i,5) - fact*vecProd
 end do

 grp%nodeHelpArray(:,9) = 0.0

! fact = oneOverReynoldsNumber*(1.0+cb2)/sigma
 fact = oneOverReynoldsNumber/sigma
 do i=1,grp%numberOfSides
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)

  if(.false.) then !(ivd%HighOrder.and.grp%gridNumber==1) then
    if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
       if(grp%AfterAndBefore(i1,1).eq.i2) then
         j1 = i1
         j2 = i2
         j3 = grp%AfterAndBefore(j2,1)
         j4 = grp%AfterAndBefore(j3,1)
         j0 = grp%AfterAndBefore(j1,2)
         jm1 = grp%AfterAndBefore(j0,2)
       else if(grp%AfterAndBefore(i1,2).eq.i2) then
         j1 = i2
         j2 = i1
         j3 = grp%AfterAndBefore(j2,1)
         j4 = grp%AfterAndBefore(j3,1)
         j0 = grp%AfterAndBefore(j1,2)
         jm1 = grp%AfterAndBefore(j0,2)
       else
         print *,' should not be here'
       end if
       if(j0.eq.0) then
         !jm1=1
         !j5 = grp%AfterAndBefore(j4,1)
         visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(j1)/grp%u(j1,1)+grp%u(j1,5)
         visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(j2)/grp%u(j2,1)+grp%u(j2,5)
         visc3 = ivd%ReynoldsNumber*grp%laminarViscosity(j3)/grp%u(j3,1)+grp%u(j3,5)
         visc4 = ivd%ReynoldsNumber*grp%laminarViscosity(j4)/grp%u(j4,1)+grp%u(j4,5)

         dnudx1 = grp%nodeHelpArray(j1,12)
         dnudy1 = grp%nodeHelpArray(j1,13)
         dnudx2 = grp%nodeHelpArray(j2,12)
         dnudy2 = grp%nodeHelpArray(j2,13)
         dnudx3 = grp%nodeHelpArray(j3,12)
         dnudy3 = grp%nodeHelpArray(j3,13)
         dnudx4 = grp%nodeHelpArray(j4,12)
         dnudy4 = grp%nodeHelpArray(j4,13)
 
         a1 = grp%nodeVolume(j1)   
         a2 = grp%nodeVolume(j2)
         a3 = grp%nodeVolume(j3)  
         a4 = grp%nodeVolume(j4)     

         viscdnudx1 = visc1*dnudx1
         viscdnudx2 = visc2*dnudx2
         viscdnudx3 = visc3*dnudx3
         viscdnudx4 = visc4*dnudx4

         viscdnudy1 = visc1*dnudy1
         viscdnudy2 = visc2*dnudy2
         viscdnudy3 = visc3*dnudy3
         viscdnudy4 = visc4*dnudy4

         !debugging
        !if (k==1.0) then
        !  a1=0.05
        !  a2=0.1
        !  a3=0.1
        !  a4=0.1
          !x^3-2*x
        !  viscdnudy1=-1.9525
        !  viscdnudy2=-1.8775 !mean values
        !  viscdnudy3=-1.7275
        !  viscdnudy4=-1.5175
        !endif

        !Approximate Jacobian
        ia1=a1
        ia2=a2+ia1
        ia3=a3+ia2
        ia4=a4+ia3

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location        
        am=c0*ia1+c1*ia2+c2*ia3+c3*ia4

        !c0 = grp%hoc3(i,1)
        !c1 = grp%hoc3(i,2)
        !c2 = grp%hoc3(i,3)
        !c3 = grp%hoc3(i,4)
        !Jacobian at boundary        
        !amb=c0*ia1+c1*ia2+c2*ia3+c3*ia4

        !c0 = grp%hoc(i,1)
        !c1 = grp%hoc(i,2)
        !c2 = grp%hoc(i,3)
        !c3 = grp%hoc(i,4)

        !note:all values assumed to be mean values

        !primitive functions
         viscdnudx1 = viscdnudx1*a1
         viscdnudx2 = viscdnudx1+viscdnudx2*a2
         viscdnudx3 = viscdnudx2+viscdnudx3*a3
         viscdnudx4 = viscdnudx3+viscdnudx4*a4

         viscdnudy1 = viscdnudy1*a1
         viscdnudy2 = viscdnudy1+viscdnudy2*a2
         viscdnudy3 = viscdnudy2+viscdnudy3*a3
         viscdnudy4 = viscdnudy3+viscdnudy4*a4

         viscdnudxm = c0*viscdnudx1+c1*viscdnudx2+ &
                      c2*viscdnudx3+c3*viscdnudx4
         viscdnudym = c0*viscdnudy1+c1*viscdnudy2+ &
                      c2*viscdnudy3+c3*viscdnudy4

         !divide through by Jacobian
         viscdnudxm = viscdnudxm/am
         viscdnudym = viscdnudym/am

        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',viscdnudym
        !  print*,'check'
        !  k=0.0
          !should be -1.9325
        !endif

         fdiffm = wx*(2.*viscdnudxm)+wy*(2.*viscdnudym)

       else if(jm1.eq.0) then
         visc0 = ivd%ReynoldsNumber*grp%laminarViscosity(j0)/grp%u(j0,1)+grp%u(j0,5)
         visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(j1)/grp%u(j1,1)+grp%u(j1,5)
         visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(j2)/grp%u(j2,1)+grp%u(j2,5)
         visc3 = ivd%ReynoldsNumber*grp%laminarViscosity(j3)/grp%u(j3,1)+grp%u(j3,5)
         visc4 = ivd%ReynoldsNumber*grp%laminarViscosity(j4)/grp%u(j4,1)+grp%u(j4,5)
         dnudx0 = grp%nodeHelpArray(j0,12)
         dnudy0 = grp%nodeHelpArray(j0,13)
         dnudx1 = grp%nodeHelpArray(j1,12)
         dnudy1 = grp%nodeHelpArray(j1,13)
         dnudx2 = grp%nodeHelpArray(j2,12)
         dnudy2 = grp%nodeHelpArray(j2,13)
         dnudx3 = grp%nodeHelpArray(j3,12)
         dnudy3 = grp%nodeHelpArray(j3,13)

         a0 = grp%nodeVolume(j0)   
         a1 = grp%nodeVolume(j1)
         a2 = grp%nodeVolume(j2)  
         a3 = grp%nodeVolume(j3)

         viscdnudx0 = visc0*dnudx0
         viscdnudx1 = visc1*dnudx1
         viscdnudx2 = visc2*dnudx2
         viscdnudx3 = visc3*dnudx3

         viscdnudy0 = visc0*dnudy0
         viscdnudy1 = visc1*dnudy1
         viscdnudy2 = visc2*dnudy2
         viscdnudy3 = visc3*dnudy3

         c0 = grp%hoc2(i,1)
         c1 = grp%hoc2(i,2)
         c2 = grp%hoc2(i,3)
         c3 = grp%hoc2(i,4)

        !debugging
        !if (k==1.0) then
        !  a0=0.05
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
          !x^3-2*x mean gradients
        !  viscdnudy0=-1.9525 
        !  viscdnudy1=-1.8775 
        !  viscdnudy2=-1.7275
        !  viscdnudy3=-1.5175
        !endif

        !Approximate Jacobian
        ia0=a0
        ia1=a1+ia0
        ia2=a2+ia1
        ia3=a3+ia2

        c0 = grp%hoc2(i,1)
        c1 = grp%hoc2(i,2)
        c2 = grp%hoc2(i,3)
        c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location        
        am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

        !c0 = grp%hoc3(i,1)
        !c1 = grp%hoc3(i,2)
        !c2 = grp%hoc3(i,3)
        !c3 = grp%hoc3(i,4)
        !Jacobian at boundary        
        !amb=c0*ia0+c1*ia1+c2*ia2+c3*ia3

        !c0 = grp%hoc(i,1)
        !c1 = grp%hoc(i,2)
        !c2 = grp%hoc(i,3)
        !c3 = grp%hoc(i,4)

         !primitive functions
         viscdnudx0 = viscdnudx0*a0
         viscdnudx1 = viscdnudx0+viscdnudx1*a1
         viscdnudx2 = viscdnudx1+viscdnudx2*a2
         viscdnudx3 = viscdnudx2+viscdnudx3*a3

         viscdnudy0 = viscdnudy0*a0
         viscdnudy1 = viscdnudy0+viscdnudy1*a1
         viscdnudy2 = viscdnudy1+viscdnudy2*a2
         viscdnudy3 = viscdnudy2+viscdnudy3*a3

         viscdnudxm = c0*viscdnudx0+c1*viscdnudx1+ &
                      c2*viscdnudx2+c3*viscdnudx3
         viscdnudym = c0*viscdnudy0+c1*viscdnudy1+ &
                      c2*viscdnudy2+c3*viscdnudy3

         !divide through by Jacobian
         viscdnudxm = viscdnudxm/am
         viscdnudym = viscdnudym/am

        !debugging
        !if (k==1.0) then
        !  print*,'interpolated value',viscdnudym
        !  k=0.0
          !should be -1.8125
        !endif

         fdiffm = wx*(2.*viscdnudxm)+wy*(2.*viscdnudym)


       else if(j4.eq.0) then

         dnudx1 = grp%nodeHelpArray(j1,12)
         dnudy1 = grp%nodeHelpArray(j1,13)
         dnudx2 = grp%nodeHelpArray(j2,12)
         dnudy2 = grp%nodeHelpArray(j2,13)
         visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(j1)/grp%u(j1,1)+grp%u(j1,5)
         visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(j2)/grp%u(j2,1)+grp%u(j2,5)
         
         fdiffm = wx*(visc1*dnudx1+visc2*dnudx2)+wy*(visc1*dnudy1+visc2*dnudy2)

       else

         visc0 = ivd%ReynoldsNumber*grp%laminarViscosity(j0)/grp%u(j0,1)+grp%u(j0,5)
         visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(j1)/grp%u(j1,1)+grp%u(j1,5)
         visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(j2)/grp%u(j2,1)+grp%u(j2,5)
         visc3 = ivd%ReynoldsNumber*grp%laminarViscosity(j3)/grp%u(j3,1)+grp%u(j3,5)

         dnudx0 = grp%nodeHelpArray(j0,12)
         dnudy0 = grp%nodeHelpArray(j0,13)
         dnudx1 = grp%nodeHelpArray(j1,12)
         dnudy1 = grp%nodeHelpArray(j1,13)
         dnudx2 = grp%nodeHelpArray(j2,12)
         dnudy2 = grp%nodeHelpArray(j2,13)
         dnudx3 = grp%nodeHelpArray(j3,12)
         dnudy3 = grp%nodeHelpArray(j3,13)

         a0 = grp%nodeVolume(j0)   
         a1 = grp%nodeVolume(j1)
         a2 = grp%nodeVolume(j2)  
         a3 = grp%nodeVolume(j3)

         viscdnudx0 = visc0*dnudx0
         viscdnudx1 = visc1*dnudx1
         viscdnudx2 = visc2*dnudx2
         viscdnudx3 = visc3*dnudx3

         viscdnudy0 = visc0*dnudy0
         viscdnudy1 = visc1*dnudy1
         viscdnudy2 = visc2*dnudy2
         viscdnudy3 = visc3*dnudy3

         c0 = grp%hoc2(i,1)
         c1 = grp%hoc2(i,2)
         c2 = grp%hoc2(i,3)
         c3 = grp%hoc2(i,4)
               !primitive functions
        !debugging
        !if (k==1.0) then
        !  a0=0.1
        !  a1=0.1
        !  a2=0.1
        !  a3=0.1
          !x^3-2*x mean gradients
        !  viscdnudy0=-1.9675 
        !  viscdnudy1=-1.8775 
        !  viscdnudy2=-1.7275
        !  viscdnudy3=-1.5175
        !endif

        !Approximate Jacobian
         ia0=a0
         ia1=a1+ia0
         ia2=a2+ia1
         ia3=a3+ia2

         c0 = grp%hoc2(i,1)
         c1 = grp%hoc2(i,2)
         c2 = grp%hoc2(i,3)
         c3 = grp%hoc2(i,4)
        !Jacobian at interpolation location          
         am=c0*ia0+c1*ia1+c2*ia2+c3*ia3

         c0 = grp%hoc(i,1)
         c1 = grp%hoc(i,2)
         c2 = grp%hoc(i,3)
         c3 = grp%hoc(i,4)

         !primitive functions
         viscdnudx0 = viscdnudx0*a0
         viscdnudx1 = viscdnudx0+viscdnudx1*a1
         viscdnudx2 = viscdnudx1+viscdnudx2*a2
         viscdnudx3 = viscdnudx2+viscdnudx3*a3

         viscdnudy0 = viscdnudy0*a0
         viscdnudy1 = viscdnudy0+viscdnudy1*a1
         viscdnudy2 = viscdnudy1+viscdnudy2*a2
         viscdnudy3 = viscdnudy2+viscdnudy3*a3

         viscdnudxm = c0*viscdnudx0+c1*viscdnudx1+ &
                      c2*viscdnudx2+c3*viscdnudx3
         viscdnudym = c0*viscdnudy0+c1*viscdnudy1+ &
                      c2*viscdnudy2+c3*viscdnudy3

         !divide through by Jacobian
         viscdnudxm = viscdnudxm/am
         viscdnudym = viscdnudym/am
       !debugging
        !if (k==1.0) then
         ! print*,'interpolated value',viscdnudym
         ! k=0.0
          !should be -1.8125
        !endif

         fdiffm = wx*(2.*viscdnudxm)+wy*(2.*viscdnudym)
        
       end if

    else

       dnudx1 = grp%nodeHelpArray(i1,12)
       dnudy1 = grp%nodeHelpArray(i1,13)
       dnudx2 = grp%nodeHelpArray(i2,12)
       dnudy2 = grp%nodeHelpArray(i2,13)
       visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(i1)/grp%u(i1,1)+grp%u(i1,5)
       visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(i2)/grp%u(i2,1)+grp%u(i2,5)
         
       fdiffm = wx*(visc1*dnudx1+visc2*dnudx2)+wy*(visc1*dnudy1+visc2*dnudy2)

    end if

  else

       dnudx1 = grp%nodeHelpArray(i1,12)
       dnudy1 = grp%nodeHelpArray(i1,13)
       dnudx2 = grp%nodeHelpArray(i2,12)
       dnudy2 = grp%nodeHelpArray(i2,13)
       visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(i1)/grp%u(i1,1)+grp%u(i1,5)
       visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(i2)/grp%u(i2,1)+grp%u(i2,5)
         
       fdiffm = wx*(visc1*dnudx1+visc2*dnudx2)+wy*(visc1*dnudy1+visc2*dnudy2)

  end if

!KAS  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)+grp%u(ind1,5)
!KAS  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)+grp%u(ind2,5)

  ! add second derivatives
  grp%rhs(i1,5) = grp%rhs(i1,5) - fact*fdiffm
  grp%rhs(i2,5) = grp%rhs(i2,5) + fact*fdiffm
 end do

! write(*,*) "c2: ",grp%rhs(ind,5)
 ! boundary faces
 do i=1,grp%brp%numberOfBoundaryFaces
  ind1 = grp%brp%faceIndexArray(i,1)
  ind2 = grp%brp%faceIndexArray(i,2)
  wx   = grp%brp%faceWeightsArray(i,1)
  wy   = grp%brp%faceWeightsArray(i,2)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)

!KAS  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)+grp%u(ind1,5)
!KAS  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)+grp%u(ind2,5)

  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)/grp%u(ind1,1)+grp%u(ind1,5)
  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)/grp%u(ind2,1)+grp%u(ind2,5)

  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2)

  ! adding boundary face contributions

  if(ivd%boundaryTerm==1) then 
   grp%rhs(ind1,5) = grp%rhs(ind1,5) - fact*(3.*fdiff1+fdiff2)
   grp%rhs(ind2,5) = grp%rhs(ind2,5) - fact*(fdiff1+3.*fdiff2)
  else
   grp%rhs(ind1,5) = grp%rhs(ind1,5) - fact*4.*fdiff1
   grp%rhs(ind2,5) = grp%rhs(ind2,5) - fact*4.*fdiff2
  end if
 end do

 end subroutine makeSADiffusion
!-----------------------------------------------------------------------
 subroutine makeSASourceTerm(grp,ivd)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind
 double precision :: vt,mucurl,Ksi,ftt,fv1,fv2,fv3,Scurl,d,dt,KsiInThird,r
 double precision :: dSquared,muOverd,fw,deltau,deltauSquared,tripVorticity,tripDeltax
 double precision :: ft1,ft2,gt,cb1,cb2,kappa,kappaSquared,sigma,cw1,cw2,cw3,cv1,cv2
 double precision :: ct1,ct2,ct3,ct4,cv1Tripled,oneOverReynoldsNumber
 double precision :: radicand,g,cw3InSixth,oneOverSix

 ind = 39498

 oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
 cb1 = 0.1355
 cb2 = 0.622
 kappa = 0.41
 kappaSquared = kappa*kappa
 sigma = 2./3.
 cw1 = (cb1/kappaSquared) + (1.0+cb2)/sigma
 cw2 = 0.3
 cw3 = 2.0
 cv1 = 7.1
 cv2 = 5.0
 ct1 = 1.0
 ct2 = 2.0
! ct3 = 1.2 ! new value
 ct3 = 1.1 ! old value
! ct4 = 0.5 ! new value
 ct4 = 2.0 ! old value
 cv1Tripled = cv1**3
 cw3InSixth =cw3**6
 oneOverSix = 1./6.

 do i=1,grp%numberOfNodes
  vt = grp%nodeVolume(i)
! mass matrices are lumped for source terms
  mucurl = grp%u(i,5)
  mucurl = max(mucurl,0.0)
  if(grp%laminarViscosity(i)<0.1/ivd%ReynoldsNumber) then
   write(*,*) "L: ",i,grp%laminarViscosity(i)
   grp%laminarViscosity(i) = 0.1*ivd%ReynoldsNumber
  end if
! KAS  Ksi = oneOverReynoldsNumber*mucurl/grp%laminarViscosity(i)
  Ksi = grp%u(i,1)*oneOverReynoldsNumber*mucurl/grp%laminarViscosity(i)
  Ksi = max(Ksi,0.001)
  ft2 = ct3*exp(-ct4*Ksi*Ksi)
  KsiInThird = Ksi**3
  d = grp%wallDistance(i)
  d = max(d,1.0e-10)

! KAS  muOverd = mucurl/(d*grp%u(i,1))
  muOverd = mucurl/d

  dSquared = d*d
  fv1 = KsiInThird/(KsiInThird+cv1Tripled)
!  fv2 = 1.0/((1+Ksi/cv2)**3) ! new formulation
  fv2 = 1.0-(Ksi/(1.0+Ksi*fv1)) ! old formulation
!  fv3 = (1.0+Ksi*fv1)*(1.0-fv2)/Ksi ! new formulation
  fv3 = 1.0                        ! old formulation
! KASWWWW  Scurl = fv3*grp%vorticity(i) + oneOverReynoldsNumber*mucurl*fv2/(kappaSquared*dSquared*grp%u(i,1))

  Scurl = fv3*grp%vorticity(i) + oneOverReynoldsNumber*mucurl*fv2/(kappaSquared*dSquared)

!  Scurl = max(Scurl,1.0e-20)
  if(abs(Scurl)>1.0e-15) then 
! KAS   r = oneOverReynoldsNumber*mucurl/(Scurl*grp%u(i,1)*kappaSquared*dSquared)
   r = oneOverReynoldsNumber*mucurl/(Scurl*kappaSquared*dSquared)
   if(r>100.0) r = 100.0
  else
   r = 100.0
  end if

  g = r + cw2*(r**6-r)
  radicand = (1.0+cw3InSixth)/(g**6+cw3InSixth)
  if(radicand>1.0e-15) then
   fw = g*(exp(oneOverSix*log(radicand)))
  else
   fw = exp(oneOverSix*log(1.0+cw3InSixth))
  end if

! add source terms

  grp%rhs(i,5) = grp%rhs(i,5) - vt*cb1*(1.0-ft2)*Scurl*mucurl
! KAS  grp%rhs(i,5) = grp%rhs(i,5) + oneOverReynoldsNumber*vt*grp%u(i,1)*(cw1*fw-(cb1*ft2/kappaSquared))*muOverd*muOverd
  grp%rhs(i,5) = grp%rhs(i,5) + oneOverReynoldsNumber*vt*(cw1*fw-(cb1*ft2/kappaSquared))*muOverd*muOverd
  grp%rhs(i,5) = grp%rhs(i,5) - vt*mucurl*grp%divergence(i) 
 end do

 end subroutine makeSASourceTerm
!-----------------------------------------------------------------------
 subroutine makeSATripTerm(grp,ivd)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind
 real :: vt,tripVorticity,tripDeltax,deltau,deltausquared,d,dt,gt,ft1,ct1,ct2,ReynoldsNumber

 ReynoldsNumber = ivd%ReynoldsNumber

 ct1 = 1.0
 ct2 = 2.0

! add trip sources


  grp%tripSourceArray = 0.0
  do i=1,grp%numberOfTripNodes
   do j=1,grp%tripNodeFieldIndexes(i,0)
    ind = grp%tripNodeFieldIndexes(i,j)
    vt = grp%nodeVolume(ind)
    deltau = sqrt(sum(grp%u(ind,2:3)*grp%u(ind,2:3)))/grp%u(ind,1)
    deltausquared = deltau*deltau
    tripVorticity = grp%vorticity(ivd%tripNodes(i))
    tripDeltax = grp%tripWallLength(i)
    d = grp%wallDistance(ind)
    dt = grp%tripNodeFieldDistances(i,j)
    gt = min(0.1,deltau/(tripVorticity*tripDeltax))
    if(deltausquared>1.0e-20) then
! KAS      ft1 = grp%u(ind,1)*ct1*gt*exp(-ct2*((tripVorticity**2)/deltausquared)*(d**2+(gt**2)*(dt**2)))
     ft1 = ct1*gt*exp(-ct2*((tripVorticity**2)/deltausquared)*(d**2+(gt**2)*(dt**2)))
    else
     if(d>1.0e-20) then
      ft1 = 0.0
     else
      ft1 = grp%u(ind,1)*ct1*gt
     end if
    end if

    grp%rhs(ind,5) = grp%rhs(ind,5) - ivd%tripFactor*ft1*ReynoldsNumber*vt*deltausquared

    grp%tripSourceArray(ind) = ivd%tripFactor*ft1*ReynoldsNumber*vt*deltausquared    
   end do
  end do


 end subroutine makeSATripTerm
!-----------------------------------------------------------------------
subroutine makeBLTurbulence(grp,ivd)
! Baldwin-Lomax turbulence model
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

integer :: i,istart,wallIndex,ind
real :: l,wallTangent(2),wallNormal(2),wallShearStress
real :: yplus,Gamma,uNorm,ReynoldsSquareRoot
real :: innerTurbulenceCoefficient,outerTurbulenceCoefficient,turbulenceCoefficient
real :: GammaMax,yMax,ymin,F,wallDistance,uMax,Fw,delta
real :: oneOverTurbPrandtlNumber,turbDiffusionCoefficient,wakeGamma

ind = 11115

istart = 1

ReynoldsSquareRoot = sqrt(ivd%ReynoldsNumber)

grp%nodeHelpArray(i,10) = 0.0 
grp%nodeHelpArray(:,11:12) = -1.0
grp%nodeHelpArray(:,13) = 0.0
grp%nodeHelpArray(:,14) = 1.0e12
grp%wakeDataArray(:,1:2) = -1.0
grp%wakeDataArray(:,3) = 0.0
grp%wakeDataArray(:,4) = 1.0e12


do i=istart,grp%numberOfNodes
 ! calculate Gamma function
 wallIndex = grp%wallDistanceFaceArray(i)
 if(wallIndex>0) then  ! use wall as base
  if(.not.associated(grp%sourceTerm)) then
   wallTangent(1) = grp%brp%faceTangentArray(wallIndex,1)
   wallTangent(2) = grp%brp%faceTangentArray(wallIndex,2)
  else
   wallIndex = grp%boundaryFaceNodeMappings(wallIndex)
   wallTangent(1) = grp%localFaceTangentArray(wallIndex,1)
   wallTangent(2) = grp%localFaceTangentArray(wallIndex,2)
  end if
  wallNormal(1) = wallTangent(2)
  wallNormal(2) = -1.0*wallTangent(1)
  wallShearStress = grp%wallStress(wallIndex,1)*wallTangent(1)+grp%wallStress(wallIndex,2)*wallTangent(2)
  wallShearStress = abs(wallShearStress)
  yplus = (grp%u(i,1)*grp%wallDistance(i)/grp%laminarViscosity(i))&
       *sqrt(wallShearStress/grp%u(i,1))
  grp%nodeHelpArray(i,10) = grp%wallDistance(i)*(1-exp(-1.0*(yplus/ivd%BLA)))
 else ! use wake as base
  grp%nodeHelpArray(i,10) = grp%wallDistance(i)
 end if

 if(i==ind) write(*,*) "d: ",grp%nodeHelpArray(i,10),grp%wallDistance(i),grp%laminarViscosity(i),wallShearStress,wallIndex
 if(i==ind) write(*,*) "d2: ",wallShearStress,yplus,(1-exp(-1.0*(yplus/ivd%BLA)))

 Gamma = grp%nodeHelpArray(i,10)*grp%vorticity(i)
 uNorm = sqrt(grp%u(i,2)**2+grp%u(i,3)**2)/grp%u(i,1)
 if(grp%wallDistance(i)>1.0e-12) then
 if(wallIndex>0) then
  if(Gamma>grp%nodeHelpArray(wallIndex,11)) then
   grp%nodeHelpArray(wallIndex,11) = Gamma
   grp%nodeHelpArray(wallIndex,12) = grp%wallDistance(i)
 !  uMax = grp%u(i,2)
  end if
  if(uNorm>grp%nodeHelpArray(wallIndex,13)) grp%nodeHelpArray(wallIndex,13) = uNorm
 else
  if(Gamma>grp%wakeDataArray(-wallIndex,1)) then
   grp%wakeDataArray(-wallIndex,1) = Gamma
   grp%wakeDataArray(-wallIndex,2) = grp%wallDistance(i)
  end if
  if(uNorm>grp%wakeDataArray(-wallIndex,3)) grp%wakeDataArray(-wallIndex,3) = uNorm
 end if
 end if
end do

 
yMin = -1.0
do i=istart,grp%numberOfNodes
 wallIndex = grp%wallDistanceFaceArray(i)
 wallDistance = grp%wallDistance(i)

 if(wallIndex>0) then
  gammaMax = grp%nodeHelpArray(wallIndex,11)
  ymax = grp%nodeHelpArray(wallIndex,12)
  umax = grp%nodeHelpArray(wallIndex,13)
 else
  gammaMax = grp%wakeDataArray(-wallIndex,1)
  ymax = grp%wakeDataArray(-wallIndex,2)
  umax = grp%wakeDataArray(-wallIndex,3)
 end if
 ! calculate inner turbulent viscosity coefficient
 l = ivd%vonKarmanConstant*grp%nodeHelpArray(i,10)
 innerTurbulenceCoefficient = l*l*grp%vorticity(i)
 grp%nodeHelpArray(i,9) = innerTurbulenceCoefficient

 if(i==ind) write(*,*) "e: ",innerTurbulenceCoefficient,l,grp%vorticity(i),ivd%vonKarmanConstant,grp%nodeHelpArray(i,10) 
 ! calculate outer turbulent viscosity coefficient
 F = 1.0/(1.0+5.5*((ivd%BLAlpha*wallDistance/yMax)**6)) ! Klebanoff intermittency
 Fw = min(yMax*gammaMax,ivd%BLCw*ymax*(uMax**2)/gammaMax) ! wake factor
 outerTurbulenceCoefficient = 0.0168*ivd%BLBeta*F*Fw
 grp%nodeHelpArray(i,10) = outerTurbulenceCoefficient

! delta = 0.37*grp%coor(i,1)*exp(-0.2*log(grp%coor(i,1)*ivd%ReynoldsNumber))
! F = 1./(1. + (wallDistance/delta)**6)
! outerTurbulenceCoefficient = 0.0168*delta*F/8.
! grp%nodeHelpArray(i,10) = outerTurbulenceCoefficient
! if(i==1012) then
! write(*,*) i,delta,grp%coor(i,1),wallDistance,F,outerTurbulenceCoefficient
! end if
 if(wallIndex>0) then
  if(grp%nodeHelpArray(wallIndex,14)>wallDistance.and.&
     innerTurbulenceCoefficient>outerTurbulenceCoefficient) then
   grp%nodeHelpArray(wallIndex,14) = wallDistance
  end if
 else
  if(grp%wakeDataArray(-wallIndex,4)>wallDistance.and.&
    innerTurbulenceCoefficient>outerTurbulenceCoefficient) then
   grp%wakeDataArray(-wallIndex,4) = wallDistance
  end if
 end if
 if(i==ind) write(*,*) "f: ",outerTurbulenceCoefficient,F,Fw,ivd%BLBeta  
 if(i==ind) write(*,*) "f2: ",yMax*gammaMax,ivd%BLCw*ymax*(uMax**2)/gammaMax
 if(i==ind) write(*,*) "f3: ",yMax,gammaMax,uMax  
end do

 write(*,*) "a: ",grp%nodeHelpArray(ind,9:10)
 write(*,*) "b: ",grp%u(ind,5),ivd%ReynoldsNumber

do i=istart,grp%numberOfNodes
 wallIndex = grp%wallDistanceFaceArray(i)
 wallDistance = grp%wallDistance(i)
 if(wallIndex>0) then
  if(wallDistance<grp%nodeHelpArray(wallIndex,14)) then
   turbulenceCoefficient = grp%nodeHelpArray(i,9) ! inner turbulence coefficient
  else
   turbulenceCoefficient = grp%nodeHelpArray(i,10) ! outer turbulence coefficient
  end if
 else
  if(wallDistance<grp%wakeDataArray(-wallIndex,4)) then
   turbulenceCoefficient = grp%nodeHelpArray(i,9) ! inner turbulence coefficient
  else
   turbulenceCoefficient = grp%nodeHelpArray(i,10) ! outer turbulence coefficient
  end if
 end if

 if(grp%coordinates(i,1)>0.03) then  
 grp%uprev(i,5) = ivd%ReynoldsNumber*turbulenceCoefficient*grp%u(i,1) ! multiply with density
 grp%u(i,5) = grp%uprev(i,5)
 else
 grp%uprev(i,5) =  0.0
 grp%u(i,5) = grp%uprev(i,5)
 end if
end do
grp%rhs(:,5) = 0.0

 write(*,*) "c: ",grp%u(ind,5),grp%u(ind,1)
end subroutine makeBLTurbulence
!-----------------------------------------------------------------------
subroutine addTurbulenceEffects1(grp,ivd,stressTensor,viscosityCoefficient,diffusionCoefficient,vorticity) 
! applies turbulence effects using Baldwin-Lomax model
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: stressTensor(:,:)
real :: viscosityCoefficient(:),diffusionCoefficient(:),vorticity(:)

real :: yplus,l,wallShearStress,wallNormal(2),wallTangent(2)
real :: innerTurbulenceCoefficient,outerTurbulenceCoefficient,turbulenceCoefficient
real :: Gamma,GammaMax,yMax,ymin,F,wallDistance,uMax,uNorm,Fw,delta
real :: oneOverTurbPrandtlNumber,turbDiffusionCoefficient,wakeGamma
integer :: i,wallIndex,istart

oneOverTurbPrandtlNumber = 1./ivd%turbulentPrandtlNumber
gammaMax = -1.0
yMax = -1.0
uMax = 0.0
!istart = grp%numberOfBoundaryCVs + 1
istart = 1

!grp%nodeHelpArray(1:istart,11:12) = -1.0
!grp%nodeHelpArray(1:istart,13) = 0.0
!grp%nodeHelpArray(1:istart,14) = 1.0e12 
grp%nodeHelpArray(:,11:12) = -1.0
grp%nodeHelpArray(:,13) = 0.0
grp%nodeHelpArray(:,14) = 1.0e12 
grp%wakeDataArray(:,1:2) = -1.0 
grp%wakeDataArray(:,3) = 0.0 
grp%wakeDataArray(:,4) = 1.0e12 

do i=istart,grp%numberOfNodes
 ! calculate Gamma function
 wallIndex = grp%wallDistanceFaceArray(i)
 if(wallIndex>0) then  ! use wall as base
  if(.not.associated(grp%sourceTerm)) then 
   wallTangent(1) = grp%brp%faceTangentArray(wallIndex,1)
   wallTangent(2) = grp%brp%faceTangentArray(wallIndex,2)
  else
   wallIndex = grp%boundaryFaceNodeMappings(wallIndex) 
   wallTangent(1) = grp%localFaceTangentArray(wallIndex,1)
   wallTangent(2) = grp%localFaceTangentArray(wallIndex,2)
  end if
  wallNormal(1) = wallTangent(2)
  wallNormal(2) = -1.0*wallTangent(1)
  wallShearStress = (stressTensor(wallIndex,1)*wallNormal(1)&
                    +stressTensor(wallIndex,2)*wallNormal(2))*wallTangent(1)+&
                    (stressTensor(wallIndex,2)*wallNormal(1)&
                    +stressTensor(wallIndex,3)*wallNormal(2))*wallTangent(2) 
  wallShearStress = abs(wallShearStress)
  yplus = (grp%u(i,1)*grp%wallDistance(i)/viscosityCoefficient(i))*sqrt(wallShearStress/grp%u(i,1)) 
  grp%nodeHelpArray(i,10) = grp%wallDistance(i)*(1-exp(-1.0*(yplus/ivd%BLA)))
 else ! use wake as base
  grp%nodeHelpArray(i,10) = grp%wallDistance(i)
 end if
 Gamma = grp%nodeHelpArray(i,10)*vorticity(i)
 uNorm = sqrt(grp%u(i,2)**2+grp%u(i,3)**2)/grp%u(i,1)
 if(grp%wallDistance(i)>1.0e-12) then 
 if(wallIndex>0) then 
  if(Gamma>grp%nodeHelpArray(wallIndex,11)) then 
   grp%nodeHelpArray(wallIndex,11) = Gamma 
   grp%nodeHelpArray(wallIndex,12) = grp%wallDistance(i)
 !  uMax = grp%u(i,2)
  end if
  if(uNorm>grp%nodeHelpArray(wallIndex,13)) grp%nodeHelpArray(wallIndex,13) = uNorm
 else
  if(Gamma>grp%wakeDataArray(-wallIndex,1)) then 
   grp%wakeDataArray(-wallIndex,1) = Gamma
   grp%wakeDataArray(-wallIndex,2) = grp%wallDistance(i)
  end if
  if(uNorm>grp%wakeDataArray(-wallIndex,3)) grp%wakeDataArray(-wallIndex,3) = uNorm
 end if
 end if
end do

!if(.not.associated(grp%sourceTerm)) then 
! open(77,file="turbdata.plt",form='formatted',status='unknown')
!end if

yMin = -1.0
do i=istart,grp%numberOfNodes
 wallIndex = grp%wallDistanceFaceArray(i)
 wallDistance = grp%wallDistance(i)

 if(wallIndex>0) then 
  gammaMax = grp%nodeHelpArray(wallIndex,11)
  ymax = grp%nodeHelpArray(wallIndex,12)
  umax = grp%nodeHelpArray(wallIndex,13)
 else
  gammaMax = grp%wakeDataArray(-wallIndex,1) 
  ymax = grp%wakeDataArray(-wallIndex,2) 
  umax = grp%wakeDataArray(-wallIndex,3) 
 end if
 ! calculate inner turbulent viscosity coefficient
 l = ivd%vonKarmanConstant*grp%nodeHelpArray(i,10) 
 innerTurbulenceCoefficient = l*l*vorticity(i)  
 grp%nodeHelpArray(i,9) = innerTurbulenceCoefficient

 ! calculate outer turbulent viscosity coefficient
!if(.false.) then 
 F = 1.0/(1.0+5.5*((ivd%BLAlpha*wallDistance/yMax)**6)) ! Klebanoff intermittency
 Fw = min(yMax*gammaMax,ivd%BLCw*ymax*(uMax**2)/gammaMax) ! wake factor
 outerTurbulenceCoefficient = 0.0168*ivd%BLBeta*F*Fw
 grp%nodeHelpArray(i,10) = outerTurbulenceCoefficient
!end if

! delta = 0.37*grp%coor(i,1)*exp(-0.2*log(grp%coor(i,1)*ivd%ReynoldsNumber))
! F = 1./(1. + (wallDistance/delta)**6) 
! outerTurbulenceCoefficient = 0.0168*delta*F/8.
! grp%nodeHelpArray(i,10) = outerTurbulenceCoefficient
! if(i==1012) then 
! write(*,*) i,delta,grp%coor(i,1),wallDistance,F,outerTurbulenceCoefficient 
! end if
 if(wallIndex>0) then 
  if(grp%nodeHelpArray(wallIndex,14)>wallDistance.and.&
     innerTurbulenceCoefficient>outerTurbulenceCoefficient) then 
   grp%nodeHelpArray(wallIndex,14) = wallDistance 
  end if
 else
  if(grp%wakeDataArray(-wallIndex,4)>wallDistance.and.&
    innerTurbulenceCoefficient>outerTurbulenceCoefficient) then
   grp%wakeDataArray(-wallIndex,4) = wallDistance 
  end if
 end if
end do

!if(yMin<0) STOP "ERROR: Something's wrong in addTurbulenceEffects - 1"

do i=istart,grp%numberOfNodes
 wallIndex = grp%wallDistanceFaceArray(i)
 wallDistance = grp%wallDistance(i)
 if(wallIndex>0) then 
  if(wallDistance<grp%nodeHelpArray(wallIndex,14)) then 
   turbulenceCoefficient = grp%nodeHelpArray(i,9) ! inner turbulence coefficient
  else
   turbulenceCoefficient = grp%nodeHelpArray(i,10) ! outer turbulence coefficient
  end if
 else
  if(wallDistance<grp%wakeDataArray(-wallIndex,4)) then 
   turbulenceCoefficient = grp%nodeHelpArray(i,9) ! inner turbulence coefficient
  else
   turbulenceCoefficient = grp%nodeHelpArray(i,10) ! outer turbulence coefficient
  end if
 end if
 
 turbulenceCoefficient = turbulenceCoefficient*grp%u(i,1) ! multiply with density

! if(wallIndex<75) turbulenceCoefficient = 0.0

end do

end subroutine addTurbulenceEffects1
!-----------------------------------------------------------------------
subroutine makeArtificialDissipation1(grp,ivd)
 ! makes artificial dissipation and puts it in RHS vector
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! a bunch of help variables
 integer :: i,i1,i2,ist,ien,coarseNodeNumber,ip,ib
 real :: r1,r2,ru1,ru2,u1,u2,h1,h2,p1,p2,wx,wy,rx,ux,vx,tx,ry,uy,vy,ty 
 real :: dr,du,dv,dh,dl,vt,epslm,eps1,gam1,rnx,rny,al,vmo1,rhol,rho1
 real :: uxl,vyl,rhor,uxr,vyr,epsl,presl,hl,epsr,presr,hr,di,d1
 real :: ui,vi,hi,rlam1,rlam2,rlam3,s1,s2,al1x,al2x,cc14,cc24,cc12,cc22
 real :: v1,t1,dr2,du2,dv2,de2,dcf,def,dr4,du4,dv4,de4,d2,d4,d41,d42
 real :: v2,t2,de,dx,dy,df,dp,sp,d0,ci,ci2,af,ucp,u,v,c1,c2,weightNorm,dist
 real :: secondOrderDis,fourthOrderDis

 ! assumes that rhs is initiated

 ! Linearly preserving artificial dissipation by Crumpton et al
 ! Tenth int. conf. on num. meth. for lam. and turb. flow 

 ! initiate help array to zero
 ! vector is split as follows:
 ! 1-8:    dU_i/dx_j
 ! 9-12:   L_j(U)
 ! 13-14:  L_j(x)
 ! 15:     scaling coefficients 


 grp%dissipation = 0.0 
 grp%nodeHelpArray = 0.0



 do i=1,grp%numberOfSides
! indexes of nodes in side
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)  

  r1   = grp%u(i1,1)
  u1   = grp%u(i1,2)
  v1   = grp%u(i1,3)
  t1   = grp%u(i1,4)
  r2   = grp%u(i2,1)
  u2   = grp%u(i2,2)
  v2   = grp%u(i2,3)
  t2   = grp%u(i2,4)
  rx   = (r1+r2)*wx
  ux   = (u1+u2)*wx
  vx   = (v1+v2)*wx
  tx   = (t1+t2)*wx
  ry   = (r1+r2)*wy
  uy   = (u1+u2)*wy
  vy   = (v1+v2)*wy
  ty   = (t1+t2)*wy

! create dU_i/dx_j
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + rx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + ux
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + vx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + tx
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + ry
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + uy
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + vy
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + ty

  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - rx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - ux
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - vx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - tx
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - ry
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - uy
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - vy
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - ty            

! make L(U) and L_j(x)
      
  dr = grp%u(i1,1) - grp%u(i2,1)
  du = grp%u(i1,2) - grp%u(i2,2)
  dv = grp%u(i1,3) - grp%u(i2,3)
  de = grp%u(i1,4) - grp%u(i2,4)

  dx = grp%sideLengthArray(i,1)
  dy = grp%sideLengthArray(i,2)
  if(ivd%useDissipationWeighting) then
   if(associated(grp%inverseSideLengthArray))then 
   ! an agglomerated grid
    df = grp%inverseSideLengthArray(i,3)
    dl = grp%sideLengthArray(i,3)
   else
   ! the fine grid
    dl = grp%sideLengthArray(i,3) 
    if(dl.gt.0) then
     df = 1.0/dl    
    else
     write(*,*) 'ERROR: Grid points coincide'
     stop
    end if
   end if
  else
   df = 1.0
  end if

! set help variables for L_j(U)
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) - dr*df
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) - du*df
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) - dv*df
  grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) - de*df
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + dr*df
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + du*df
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + dv*df
  grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + de*df

  grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + df
  grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) + df 

! for linearly preserving factor, L_j(x)
  grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) - dx*df
  grp%nodeHelpArray(i1,14) = grp%nodeHelpArray(i1,14) - dy*df
  grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) + dx*df
  grp%nodeHelpArray(i2,14) = grp%nodeHelpArray(i2,14) + dy*df     
 end do

! scale if dissipation weighting is used

 if(ivd%useDissipationWeighting) then
  do i=1,grp%numberOfNodes
   grp%nodeHelpArray(i,15) = 3.0*grp%nodeHelpArray(i,15)
  end do
 end if

 grp%nodeHelpArray(1:grp%brp%numberOfBoundaryFaces,13:14) = 0.0 ! remove this?

do i=1,grp%brp%numberOfBoundaryFaces
 i1 = grp%brp%faceIndexArray(i,1)
 i2 = grp%brp%faceIndexArray(i,2)
 wx   = grp%brp%faceWeightsArray(i,1)
 wy   = grp%brp%faceWeightsArray(i,2) 
 r1   = grp%u(i1,1)
 u1   = grp%u(i1,2)
 v1   = grp%u(i1,3)
 t1   = grp%u(i1,4)
 r2   = grp%u(i2,1)
 u2   = grp%u(i2,2)
 v2   = grp%u(i2,3)
 t2   = grp%u(i2,4)
 rx   = r1+r2
 ux   = u1+u2
 vx   = v1+v2
 tx   = t1+t2
 ry   = r1+r2
 uy   = u1+u2
 vy   = v1+v2
 ty   = t1+t2

! adding boundary face contributions        

 if(ivd%boundaryTerm==1) then 
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*r1+rx)*wx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*u1+ux)*wx
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*v1+vx)*wx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*t1+tx)*wx
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*r1+ry)*wy
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*u1+uy)*wy
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*v1+vy)*wy
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*t1+ty)*wy
         
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*r2+rx)*wx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*u2+ux)*wx
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*v2+vx)*wx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*t2+tx)*wx
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*r2+ry)*wy
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*u2+uy)*wy
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*v2+vy)*wy
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*t2+ty)*wy     
 else
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + 4.*r1*wx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + 4.*u1*wx
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + 4.*v1*wx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + 4.*t1*wx
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + 4.*r1*wy
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + 4.*u1*wy
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + 4.*v1*wy
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + 4.*t1*wy

  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + 4.*r2*wx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + 4.*u2*wx
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + 4.*v2*wx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + 4.*t2*wx
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + 4.*r2*wy
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + 4.*u2*wy
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + 4.*v2*wy
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + 4.*t2*wy
 end if 
end do

! divide by volume to remove mass matrix on LHS
do i=1,grp%numberOfNodes
 vt = 1./grp%nodeVolume(i)
 grp%nodeHelpArray(i,1) = grp%nodeHelpArray(i,1)*vt
 grp%nodeHelpArray(i,2) = grp%nodeHelpArray(i,2)*vt
 grp%nodeHelpArray(i,3) = grp%nodeHelpArray(i,3)*vt
 grp%nodeHelpArray(i,4) = grp%nodeHelpArray(i,4)*vt
 grp%nodeHelpArray(i,5) = grp%nodeHelpArray(i,5)*vt
 grp%nodeHelpArray(i,6) = grp%nodeHelpArray(i,6)*vt
 grp%nodeHelpArray(i,7) = grp%nodeHelpArray(i,7)*vt
 grp%nodeHelpArray(i,8) = grp%nodeHelpArray(i,8)*vt
end do

! make L^hat*grp%nodeHelpArray(10,x) 
!               = (L(U)- dU_i/dx_j L_j(x))*grp%nodeHelpArray(10,x)

 do i=1,grp%numberOfNodes                        
  grp%nodeHelpArray(i,9) = grp%nodeHelpArray(i,9)&
    +grp%nodeHelpArray(i,1)*grp%nodeHelpArray(i,13)&
    +grp%nodeHelpArray(i,5)*grp%nodeHelpArray(i,14)       
  grp%nodeHelpArray(i,10) = grp%nodeHelpArray(i,10)&
    +grp%nodeHelpArray(i,2)*grp%nodeHelpArray(i,13)&
    +grp%nodeHelpArray(i,6)*grp%nodeHelpArray(i,14) 
  grp%nodeHelpArray(i,11) = grp%nodeHelpArray(i,11)& 
    +grp%nodeHelpArray(i,3)*grp%nodeHelpArray(i,13)&
    +grp%nodeHelpArray(i,7)*grp%nodeHelpArray(i,14)
  grp%nodeHelpArray(i,12) = grp%nodeHelpArray(i,12)&  
    +grp%nodeHelpArray(i,4)*grp%nodeHelpArray(i,13)&
    +grp%nodeHelpArray(i,8)*grp%nodeHelpArray(i,14)           
 end do

! remove fourth order viscosity at farfield boundary

ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 grp%nodeHelpArray(coarseNodeNumber,9:12) = 0.0
end do

! make pressure switch

! no need for gradient matrix anymore
do i = 1,grp%numberOfNodes
 grp%nodeHelpArray(i,1) = 0.0
 grp%nodeHelpArray(i,2) = 0.0
end do

do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
 dp = grp%p(i1) - grp%p(i2)
 sp = grp%p(i1) + grp%p(i2)
! set help variables
 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) - dp
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + sp
 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dp
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + sp
end do

! pressure switch, see Hirsh vol2 p.280
do i=1,grp%numberOfNodes
 grp%nodeHelpArray(i,16) = 12.*abs(grp%nodeHelpArray(i,1)/grp%nodeHelpArray(i,2))
end do

if(.not.ivd%useDissipationWeighting) then 
 do i=1,grp%numberOfNodes
  grp%nodeHelpArray(i,9:12) = grp%nodeHelpArray(i,9:12)/grp%nodeConnectivityArray(i)
 end do
end if

! assemble smoother

secondOrderDis = ivd%secondOrderDissipationFactor
if(grp%gridNumber==1) then
 fourthOrderDis = ivd%fourthOrderDissipationFactor
else
 fourthOrderDis = ivd%coarseGridDissipationFactor
end if


do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
! making L^hat differences
 dr4 = grp%nodeHelpArray(i1,9) - grp%nodeHelpArray(i2,9)
 du4 = grp%nodeHelpArray(i1,10) - grp%nodeHelpArray(i2,10)
 dv4 = grp%nodeHelpArray(i1,11) - grp%nodeHelpArray(i2,11)
 de4 = grp%nodeHelpArray(i1,12) - grp%nodeHelpArray(i2,12) 
! making second derivative differences
 dr2 = grp%u(i1,1) - grp%u(i2,1)
 du2 = grp%u(i1,2) - grp%u(i2,2)
 dv2 = grp%u(i1,3) - grp%u(i2,3)
 de2 = grp%u(i1,4) - grp%u(i2,4)        

 if(.not.ivd%useMatrixDissipation) then
 ! scalar dissipation coefficient
  d0  = amin1(grp%nodeVolume(i1)/grp%localTimeSteps(i1),&
              grp%nodeVolume(i2)/grp%localTimeSteps(i2))/&
     (grp%nodeConnectivityArray(i1)+grp%nodeConnectivityArray(i2))
  d2  = d0*secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,16),grp%nodeHelpArray(i2,16))
  d4  = d0*fourthOrderDis
  d4  = dim(d4,d2)
  if(ivd%useDissipationWeighting) then 
   d41  = d4/grp%nodeHelpArray(i1,15)
   d42  = d4/grp%nodeHelpArray(i2,15)
  else
   d41  = d4
   d42  = d4
  end if
 else
! Roe matrix dissipation
! (based on FLITE3d
!  Computational Dynamics Research Ltd.                 
!  Innovation Centre,
!  Swansea)   

  epslm = ivd%HartensCorrectionFactor 

  eps1 = 1./(max(epslm,1.e-5))
  gam1 = ivd%gamma-1.
  rnx   = grp%sideWeightsArray(i,1)
  rny   = grp%sideWeightsArray(i,2)
  al    = sqrt(rnx*rnx+rny*rny) ! CAN USE sideLengthsArray here 
  vmo1  = 1./al
  rnx   = rnx*vmo1
  rny   = rny*vmo1 

! *** node 1
  rhol  = grp%u(i1,1)
  rho1  = 1./rhol
  uxl   = grp%u(i1,2)*rho1
  vyl   = grp%u(i1,3)*rho1
  epsl  = grp%u(i1,4)
  presl = grp%p(i1)
  hl    = (epsl + presl)*rho1

! *** node 2
  rhor  = grp%u(i2,1)
  rho1  = 1./rhor
  uxr   = grp%u(i2,2)*rho1
  vyr   = grp%u(i2,3)*rho1
  epsr  = grp%u(i2,4)
  presr = grp%p(i2)
  hr    = (epsr + presr)*rho1

! *** Roe's averaging
  if(rhor/rhol<0) then 
   write(*,*) "NEGATIVE: ",i1,i2,rhor,rhol
  end if
  di    = sqrt(rhor/rhol)
  d1    = 1.0/(di+1.0)
  ui    = (di*uxr+uxl)*d1
  vi    = (di*vyr+vyl)*d1
  hi    = (di*hr+hl)*d1
  ci2   = gam1*(hi-0.5*(ui*ui+vi*vi))
  if(ci2<0) then 
   write(*,*) "NEGATIVE2: ",i1,i2,hi,ui,vi
  end if
  ci    = sqrt(ci2)
  af    = 0.5*(ui*ui+vi*vi)
  ucp   = ui*rnx+vi*rny

! *** eigenvalues

  rlam1 = abs(ucp+ci)
  rlam2 = abs(ucp-ci)
  rlam3 = abs(ucp)

! *** Harten's correction

  if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
  if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
  if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

! *** dissipation terms "a la Turkel"

  s1    = 0.5*(rlam1+rlam2)
  s2    = 0.5*(rlam1-rlam2)
  al1x  = gam1*(af*dr4-ui*du4-vi*dv4+de4)
  al2x  = -ucp*dr4+du4*rnx
  cc14   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc24   = (s2*al1x/ci)+(s1-rlam3)*al2x
     
  al1x  = gam1*(af*dr2-ui*du2-vi*dv2+de2)
  al2x  = -ucp*dr2+du2*rnx   
  cc12   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc22   = (s2*al1x/ci)+(s1-rlam3)*al2x     

! *** [A]*(Uj-Ui)

  dcf   = al
  dr4    = dcf*(rlam3*dr4+cc14           )
  du4    = dcf*(rlam3*du4+cc14*ui+cc24*rnx)
  dv4    = dcf*(rlam3*dv4+cc14*vi+cc24*rny)
  de4    = dcf*(rlam3*de4+cc14*hi+cc24*ucp)

  dr2    = dcf*(rlam3*dr2+cc12           )
  du2    = dcf*(rlam3*du2+cc12*ui+cc22*rnx)
  dv2    = dcf*(rlam3*dv2+cc12*vi+cc22*rny)
  de2    = dcf*(rlam3*de2+cc12*hi+cc22*ucp)

! set dissipation coefficients
  d2  = secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,16),grp%nodeHelpArray(i2,16))
  d4  = fourthOrderDis
  d4  = dim(d4,d2)
  d41  = d4
  d42  = d4
 endif

! update dissipation matrix 
 grp%dissipation(i1,1) = grp%dissipation(i1,1) - dr4*d41+dr2*d2
 grp%dissipation(i1,2) = grp%dissipation(i1,2) - du4*d41+du2*d2
 grp%dissipation(i1,3) = grp%dissipation(i1,3) - dv4*d41+dv2*d2
 grp%dissipation(i1,4) = grp%dissipation(i1,4) - de4*d41+de2*d2
 grp%dissipation(i2,1) = grp%dissipation(i2,1) + dr4*d42-dr2*d2
 grp%dissipation(i2,2) = grp%dissipation(i2,2) + du4*d42-du2*d2
 grp%dissipation(i2,3) = grp%dissipation(i2,3) + dv4*d42-dv2*d2
 grp%dissipation(i2,4) = grp%dissipation(i2,4) + de4*d42-de2*d2    
end do
end subroutine makeArtificialDissipation1
!---------------------------------------------------------------------------
subroutine makeArtificialDissipation2(grp,ivd)
 ! makes artificial dissipation and puts it in RHS vector
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! a bunch of help variables
 integer :: i,i1,i2,ist,ien,ib,coarseNodeNumber,ip,il,ic
 real :: r1,r2,ru1,ru2,u1,u2,h1,h2,p1,p2,wx,wy,rx,ux,vx,tx,ry,uy,vy,ty 
 real :: dr,du,dv,dh,dl,vt,epslm,eps1,gam1,rnx,rny,al,vmo1,rhol,rho1
 real :: uxl,vyl,rhor,uxr,vyr,epsl,presl,hl,epsr,presr,hr,di,d1
 real :: ui,vi,hi,rlam1,rlam2,rlam3,s1,s2,al1x,al2x,cc14,cc24,cc12,cc22
 real :: v1,t1,dr2,du2,dv2,de2,dcf,def,dr4,du4,dv4,de4,d2,d4,d41,d42
 real :: v2,t2,de,dx,dy,df,dp,sp,d0,ci,ci2,af,ucp
 real :: secondOrderDis,fourthOrderDis

 real :: dt1,dt2

 real :: relDist

 ! assumes that rhs is initiated

 ! Jameson d2-d4 dissipation 

 grp%dissipation = 0.0 
 grp%nodeHelpArray = 0.0


 do i=1,grp%numberOfSides
! indexes of nodes in side
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)

! make L(U) 
       
! if(ivd%HighOrder.and.grp%gridNumber==1.and.grp%AfterAndBefore(i1,1).ne.0) then
!   if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
!     dr = grp%u(i1,1) - grp%u(i2,1)
!     du = grp%u(i1,2) - grp%u(i2,2)
!     dv = grp%u(i1,3) - grp%u(i2,3)
!     de = grp%u(i1,4) - grp%u(i2,4) + grp%p(i1) - grp%p(i2)
!     dx = -1.0*grp%sideLengthArray(i,1) 
!     dy = -1.0*grp%sideLengthArray(i,2) 
!     if(ivd%useDissipationWeighting) then
!       if(associated(grp%inverseSideLengthArray))then
!      ! an agglomerated grid
!         df = grp%inverseSideLengthArray(i,3)
!         dl = grp%sideLengthArray(i,3)
!       else
!      ! the fine grid
!         dl = grp%sideLengthArray(i,3)
!         if(dl.gt.0) then
!           df = 1./dl
!         else
!           write(*,*) 'ERROR: Grid points coincide'
!           stop
!         end if
!       end if
!     else
!       df = 1.0
!     end if
!   else
!      dr = 0.
!      du = 0.
!      dv = 0.
!      de = 0.
!   end if
! else
    dr = grp%u(i1,1) - grp%u(i2,1)
    du = grp%u(i1,2) - grp%u(i2,2)
    dv = grp%u(i1,3) - grp%u(i2,3)
    de = grp%u(i1,4) - grp%u(i2,4) + grp%p(i1) - grp%p(i2)
    dx = -1.0*grp%sideLengthArray(i,1) 
    dy = -1.0*grp%sideLengthArray(i,2) 
    if(ivd%useDissipationWeighting) then
      if(associated(grp%inverseSideLengthArray))then
   ! an agglomerated grid
        df = grp%inverseSideLengthArray(i,3)
        dl = grp%sideLengthArray(i,3)
      else
   ! the fine grid
        dl = grp%sideLengthArray(i,3)
        if(dl.gt.0) then
          df = 1./dl
        else
          write(*,*) 'ERROR: Grid points coincide'
          stop
        end if
      end if
    else
      df = 1.0
    end if
! end if


  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) - dr*df
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) - du*df
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) - dv*df
  grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) - de*df
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + dr*df
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + du*df
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + dv*df
  grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + de*df

  grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + df
  grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) + df 
end do

! grp%nodeHelpArray(1:grp%brp%numberOfNodes,9:12) = 0.0  ! include this?

! remove fourth order viscosity at farfield boundary
ist = grp%brp%faceIndicatorArray(-16)
ien = grp%brp%faceIndicatorArray(-15)

do ib =ist,ien
 ip = grp%brp%faceIndicatorArray(ib)
 if(associated(grp%boundaryFaceNodeMappings)) then 
 ! agglomerated mesh
  coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
 else
  coarseNodeNumber = ip
 end if
 grp%nodeHelpArray(coarseNodeNumber,9:12) = 0.0
end do

!
! *** set dissipation to zero near wall  (OH)
!
do i=1,grp%numberOfNodes
  grp%nodeHelpArray(i,3) = 1.0
end do

  ic = 0
if(grp%gridNumber==1.and.ivd%numberOfDissipationLayers.ne.0) then

  ist = grp%brp%faceIndicatorArray(-12)
  ien = grp%brp%faceIndicatorArray(-9)

  do ib =ist,ien
    ip = grp%brp%faceIndicatorArray(ib)
    grp%nodeHelpArray(ip,3) = 0
    il = 0
 30 ip = grp%AfterAndBefore(ip,1)
   if(ip.ne.0.and.il.lt.ivd%numberOfDissipationLayers) then
     il = il + 1
     df = il/ivd%numberOfDissipationLayers
     grp%nodeHelpArray(ip,3) = df
     ic = ic + 1
     goto 30
   end if
  end do
end if
!
! *** end of OH
!

do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
 dp = grp%p(i1) - grp%p(i2)
 sp = grp%p(i1) + grp%p(i2)
! set help variables
 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) - dp
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + sp
 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dp
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + sp
end do

! pressure switch, see Hirsch vol2 p.280

do i=1,grp%numberOfNodes
 grp%nodeHelpArray(i,16) = 12.*abs(grp%nodeHelpArray(i,1)/grp%nodeHelpArray(i,2))
end do

if(.not.ivd%useDissipationWeighting) then 
 do i=1,grp%numberOfNodes
  grp%nodeHelpArray(i,9:12) = grp%nodeHelpArray(i,9:12)/grp%nodeConnectivityArray(i)
 end do
end if
!
! fourth order dissipation
!
 do i=1,grp%numberOfSides
! indexes of nodes in side
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)
      
  dr = grp%nodeHelpArray(i1,9) - grp%nodeHelpArray(i2,9)
  du = grp%nodeHelpArray(i1,10) - grp%nodeHelpArray(i2,10)
  dv = grp%nodeHelpArray(i1,11) - grp%nodeHelpArray(i2,11)
  de = grp%nodeHelpArray(i1,12) - grp%nodeHelpArray(i2,12)

  df = 1.0

  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) - dr*df
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) - du*df
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) - dv*df
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) - de*df
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + dr*df
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + du*df
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + dv*df
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + de*df

end do

 do i=1,grp%numberOfNodes
  grp%nodeHelpArray(i,5:8) = grp%nodeHelpArray(i,5:8)/grp%nodeConnectivityArray(i)
 end do

! assemble smoother

secondOrderDis = ivd%secondOrderDissipationFactor
if(grp%gridNumber==1) then
 fourthOrderDis = ivd%fourthOrderDissipationFactor
else
 fourthOrderDis = ivd%coarseGridDissipationFactor
end if

do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
!D6
 !if(.false.)then 
 if(ivd%HighOrder.and.grp%gridNumber==1) then
   if(grp%AfterAndBefore(i1,1).eq.i2.or.grp%AfterAndBefore(i1,2).eq.i2)then
     dr4 = grp%nodeHelpArray(i2,5) - grp%nodeHelpArray(i1,5)
     du4 = grp%nodeHelpArray(i2,6) - grp%nodeHelpArray(i1,6)
     dv4 = grp%nodeHelpArray(i2,7) - grp%nodeHelpArray(i1,7)
     de4 = grp%nodeHelpArray(i2,8) - grp%nodeHelpArray(i1,8) 
     fourthOrderDis = ivd%sixthOrderDissipationFactor
!    dr4 = grp%nodeHelpArray(i1,9) - grp%nodeHelpArray(i2,9)
!    du4 = grp%nodeHelpArray(i1,10) - grp%nodeHelpArray(i2,10)
!    dv4 = grp%nodeHelpArray(i1,11) - grp%nodeHelpArray(i2,11)
!    de4 = grp%nodeHelpArray(i1,12) - grp%nodeHelpArray(i2,12) 
   else
! making L differences
     dr4 = grp%nodeHelpArray(i1,9) - grp%nodeHelpArray(i2,9)
     du4 = grp%nodeHelpArray(i1,10) - grp%nodeHelpArray(i2,10)
     dv4 = grp%nodeHelpArray(i1,11) - grp%nodeHelpArray(i2,11)
     de4 = grp%nodeHelpArray(i1,12) - grp%nodeHelpArray(i2,12) 
     fourthOrderDis = ivd%fourthOrderDissipationFactor
   end if
 else
   dr4 = grp%nodeHelpArray(i1,9) - grp%nodeHelpArray(i2,9)
   du4 = grp%nodeHelpArray(i1,10) - grp%nodeHelpArray(i2,10)
   dv4 = grp%nodeHelpArray(i1,11) - grp%nodeHelpArray(i2,11)
   de4 = grp%nodeHelpArray(i1,12) - grp%nodeHelpArray(i2,12) 
 end if
! making second derivative differences
 dr2 = grp%u(i1,1) - grp%u(i2,1)
 du2 = grp%u(i1,2) - grp%u(i2,2)
 dv2 = grp%u(i1,3) - grp%u(i2,3)
 de2 = grp%u(i1,4) - grp%u(i2,4)        
 if(.not.ivd%useMatrixDissipation.or.grp%gridNumber>1) then
 ! scalar dissipation coefficient

  dt1 = grp%localTimeSteps(i1)
  dt2 = grp%localTimeSteps(i2)
  d0  = amin1(grp%nodeVolume(i1)/dt1,&
              grp%nodeVolume(i2)/dt2)/&
     (grp%nodeConnectivityArray(i1)+grp%nodeConnectivityArray(i2))
  d2  = d0*secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,16),grp%nodeHelpArray(i2,16))
  d4  = d0*fourthOrderDis
  d4  = dim(d4,d2)
  if(ivd%useDissipationWeighting) then
   d41  = d4/grp%nodeHelpArray(i1,15)
   d42  = d4/grp%nodeHelpArray(i2,15)
  else
   d41  = d4*grp%nodeHelpArray(i1,3)
   d42  = d4*grp%nodeHelpArray(i2,3)
  end if
 else
! Roe matrix dissipation
! (based on FLITE3d
!  Computational Dynamics Research Ltd.                 
!  Innovation Centre,
!  Swansea)   
  epslm = ivd%HartensCorrectionFactor 

  eps1 = 1./(max(epslm,1.e-5))
  gam1 = ivd%gamma-1.
  rnx   = grp%sideWeightsArray(i,1)
  rny   = grp%sideWeightsArray(i,2)
  al    = sqrt(rnx*rnx+rny*rny)  
  vmo1  = 1./al
  rnx   = rnx*vmo1
  rny   = rny*vmo1 

! *** node 1
  rhol  = grp%u(i1,1)
  rho1  = 1./rhol
  uxl   = grp%u(i1,2)*rho1
  vyl   = grp%u(i1,3)*rho1
  epsl  = grp%u(i1,4)
  presl = grp%p(i1)
  hl    = (epsl + presl)*rho1

! *** node 2
  rhor  = grp%u(i2,1)
  rho1  = 1./rhor
  uxr   = grp%u(i2,2)*rho1
  vyr   = grp%u(i2,3)*rho1
  epsr  = grp%u(i2,4)
  presr = grp%p(i2)
  hr    = (epsr + presr)*rho1

! *** Roe's averaging
  di    = sqrt(rhor/rhol)
  d1    = 1.0/(di+1.0)
  ui    = (di*uxr+uxl)*d1
  vi    = (di*vyr+vyl)*d1
  hi    = (di*hr+hl)*d1
  ci2   = gam1*(hi-0.5*(ui*ui+vi*vi))
  ci    = sqrt(ci2)
  af    = 0.5*(ui*ui+vi*vi)
  ucp   = ui*rnx+vi*rny

! *** eigenvalues

  rlam1 = abs(ucp+ci)
  rlam2 = abs(ucp-ci)
  rlam3 = abs(ucp)

! *** Harten's correction

  if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
  if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
  if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

! *** dissipation terms "a la Turkel"

 
  s1    = 0.5*(rlam1+rlam2)
  s2    = 0.5*(rlam1-rlam2)
  al1x  = gam1*(af*dr4-ui*du4-vi*dv4+de4)
  al2x  = -ucp*dr4+du4*rnx
  cc14   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc24   = (s2*al1x/ci)+(s1-rlam3)*al2x
     
  al1x  = gam1*(af*dr2-ui*du2-vi*dv2+de2)
  al2x  = -ucp*dr2+du2*rnx   
  cc12   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc22   = (s2*al1x/ci)+(s1-rlam3)*al2x     
 
! *** [A]*(Uj-Ui)

  dcf   = al
  dr4    = dcf*(rlam3*dr4+cc14           )
  du4    = dcf*(rlam3*du4+cc14*ui+cc24*rnx)
  dv4    = dcf*(rlam3*dv4+cc14*vi+cc24*rny)
  de4    = dcf*(rlam3*de4+cc14*hi+cc24*ucp)

  dr2    = dcf*(rlam3*dr2+cc12           )
  du2    = dcf*(rlam3*du2+cc12*ui+cc22*rnx)
  dv2    = dcf*(rlam3*dv2+cc12*vi+cc22*rny)
  de2    = dcf*(rlam3*de2+cc12*hi+cc22*ucp)

! set dissipation coefficients


  d2  = secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,16),grp%nodeHelpArray(i2,16))
  d4  = fourthOrderDis
  d4  = dim(d4,d2)
  d41  = d4
  d42  = d4

!  d41  = d4*grp%nodeConnectivityArray(i1)/grp%nodeHelpArray(i1,15)
!  d42  = d4*grp%nodeConnectivityArray(i2)/grp%nodeHelpArray(i2,15)                    


 endif

! update dissipation matrix 
 grp%dissipation(i1,1) = grp%dissipation(i1,1) - dr4*d41+dr2*d2
 grp%dissipation(i1,2) = grp%dissipation(i1,2) - du4*d41+du2*d2
 grp%dissipation(i1,3) = grp%dissipation(i1,3) - dv4*d41+dv2*d2
 grp%dissipation(i1,4) = grp%dissipation(i1,4) - de4*d41+de2*d2
 grp%dissipation(i2,1) = grp%dissipation(i2,1) + dr4*d42-dr2*d2
 grp%dissipation(i2,2) = grp%dissipation(i2,2) + du4*d42-du2*d2
 grp%dissipation(i2,3) = grp%dissipation(i2,3) + dv4*d42-dv2*d2
 grp%dissipation(i2,4) = grp%dissipation(i2,4) + de4*d42-de2*d2    
end do

  rewind(799)
  ist = grp%brp%faceIndicatorArray(-12)
  ien = grp%brp%faceIndicatorArray(-9)

  do ib =ist,ien
    ip = grp%brp%faceIndicatorArray(ib)
    il = 0
!   write(799,'(2i6,4E15.7)') il,ip,grp%u(ip,1:4)
!   write(799,'(2i6,4E15.7)') 1,1,grp%nodeHelpArray(ip,9:12)
!   write(799,'(2i6,4E15.7)') 2,2,grp%dissipation(ip,1:4)
!   write(799,'(2i6,4E15.7)') 3,3,grp%nodeHelpArray(ip,5:8)
 32 ip = grp%AfterAndBefore(ip,1)
   if(ip.ne.0) then
     il = il + 1
!   write(799,'(2i6,4E15.7)') il,ip,grp%u(ip,1:4)
!   write(799,'(2i6,4E15.7)') 1,1,grp%nodeHelpArray(ip,9:12)
!   write(799,'(2i6,4E15.7)') 2,2,grp%dissipation(ip,1:4)
!   write(799,'(2i6,4E15.7)') 3,3,grp%nodeHelpArray(ip,5:8)
     goto 32
   end if
  end do

end subroutine makeArtificialDissipation2
!--------------------------------------------------------------------------
subroutine makeArtificialDissipation3(grp,ivd)
 ! makes artificial dissipation and puts it in RHS vector
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! a bunch of help variables
 integer :: i,i1,i2
 real :: r1,r2,ru1,ru2,u1,u2,h1,h2,p1,p2,wx,wy,rx,ux,vx,tx,ry,uy,vy,ty 
 real :: dr,du,dv,dh,dl,vt,epslm,eps1,gam1,rnx,rny,al,vmo1,rhol,rho1
 real :: uxl,vyl,rhor,uxr,vyr,epsl,presl,hl,epsr,presr,hr,di,d1
 real :: ui,vi,hi,rlam1,rlam2,rlam3,s1,s2,al1x,al2x,cc14,cc24,cc12,cc22
 real :: v1,t1,dr2,du2,dv2,de2,dcf,def,dr4,du4,dv4,de4,d2,d4,d41,d42
 real :: v2,t2,de,dx,dy,df,dp,sp,d0,ci,ci2,af,ucp

 ! assumes that rhs is initiated

 ! d2 dissipation only 
grp%dissipation = 0.0 

! assemble smoother

do i=1,grp%numberOfSides
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! making second derivative differences
 dr2 = grp%u(i1,1) - grp%u(i2,1)
 du2 = grp%u(i1,2) - grp%u(i2,2)
 dv2 = grp%u(i1,3) - grp%u(i2,3)
 de2 = grp%u(i1,4) - grp%u(i2,4)        
 ! scalar dissipation coefficient
 d0  = amin1(grp%nodeVolume(i1)/grp%localTimeSteps(i1),&
             grp%nodeVolume(i2)/grp%localTimeSteps(i2))/&
    (grp%nodeConnectivityArray(i1)+grp%nodeConnectivityArray(i2))
! d0 = d0/grp%sideLengthArray(i,3)
! d0 = 1.0
 d2  = d0*ivd%coarseGridDissipationFactor
 ! update dissipation matrix 
 grp%dissipation(i1,1) = grp%dissipation(i1,1)+dr2*d2
 grp%dissipation(i1,2) = grp%dissipation(i1,2)+du2*d2
 grp%dissipation(i1,3) = grp%dissipation(i1,3)+dv2*d2
 grp%dissipation(i1,4) = grp%dissipation(i1,4)+de2*d2
 grp%dissipation(i2,1) = grp%dissipation(i2,1)-dr2*d2
 grp%dissipation(i2,2) = grp%dissipation(i2,2)-du2*d2
 grp%dissipation(i2,3) = grp%dissipation(i2,3)-dv2*d2
 grp%dissipation(i2,4) = grp%dissipation(i2,4)-de2*d2    
end do
end subroutine makeArtificialDissipation3
!-----------------------------------------------------------------------
 subroutine restrictGrid(finegrp,coarsegrp,ivd,gridnum,calculateTurbulence)
 ! Does a Full Approximation Storage (FAS) multigrid fine to coarse
 ! transformation
 IMPLICIT NONE

 type(GridSolverData) :: finegrp,coarsegrp
 type(InputVariablesData) :: ivd
 integer :: gridnum
 logical :: calculateTurbulence

 integer :: i,ind
 real :: fact

 ind = 394

! call mapScalarFineToCoarse(finegrp%tripSourceArray,coarsegrp,coarsegrp%tripSourceArray)

 ! create L_H I_h^H u_h

 call mapFineToCoarse(finegrp%u,coarsegrp,coarsegrp%u)
 call makePressureField(coarsegrp,ivd,coarsegrp%u)
 call setBCsOnSolutionField(coarsegrp,ivd,coarsegrp%u)
 call setOuterBoundaryConditions3(coarsegrp,ivd,coarsegrp%u,coarsegrp%p,.true.)
 coarsegrp%uprev = coarsegrp%u
 call makePressureField(coarsegrp,ivd,coarsegrp%u)
 coarsegrp%p = abs(coarsegrp%p)
 call makeTimeSteps(coarsegrp,ivd)
 if(ivd%coarseGridDissipationScheme==1) then
  call makeArtificialDissipation1(coarsegrp,ivd)
 else if(ivd%coarseGridDissipationScheme==2) then
  call makeArtificialDissipation2(coarsegrp,ivd)
 else 
  call makeArtificialDissipation3(coarsegrp,ivd)
 end if

 if(ivd%ReynoldsNumber>0.0) then
!  if(ivd%viscosityScheme==2.and.coarsegrp%gridNumber==1) then 
  if(ivd%viscosityScheme==2) then 
   call makeViscosityTerm2(coarsegrp,ivd,coarsegrp%dissipation,.true.)
  else
   call makeViscosityTerm(coarsegrp,ivd,coarsegrp%dissipation,.true.)
  end if
  if(ivd%turbulenceModel==2) then
   call makeBLTurbulence(coarsegrp,ivd)
  end if
 end if
 call makeRHS(coarsegrp,ivd)
 coarsegrp%rhs(:,5) = 0.0
! if(ivd%turbulenceModel==1.and.coarsegrp%gridNumber<2.and.calculateTurbulence) then  ! KAS
 if(ivd%turbulenceModel==1.and.coarsegrp%gridNumber<4.and.calculateTurbulence) then
  call makeSARHS(coarsegrp,ivd)
  call makeSADiffusion(coarsegrp,ivd)
  call makeSASourceTerm(coarsegrp,ivd)
!  coarsegrp%rhs(:,5) = coarsegrp%rhs(:,5) - coarsegrp%tripSourceArray(:)
  coarsegrp%rhs(:,5) = ivd%turbulenceCFLFactor*coarsegrp%rhs(:,5)
 end if
 !call setRHSAtBoundary(coarsegrp,ivd)
 call solveSystem(coarsegrp,ivd,0,coarsegrp%rhs)
 call setBCsOnIncrementField(coarsegrp,ivd,coarsegrp%rhs) 
 call nullifyOuterBoundary(coarsegrp,coarsegrp%rhs)

 coarsegrp%u = coarsegrp%uprev - coarsegrp%rhs

 do i=1,coarsegrp%numberOfNodes
  if(coarsegrp%u(i,5)<0.0) coarsegrp%u(i,5) = 0.0
 end do

 call setBCsOnSolutionField(coarsegrp,ivd,coarsegrp%u)

 call setOuterBoundaryConditions3(coarsegrp,ivd,coarsegrp%u,coarsegrp%p,.true.)

 coarsegrp%sourceTerm = coarsegrp%u - coarsegrp%uprev
 coarsegrp%u = coarsegrp%uprev


 ! create I_h^H L_h u_h
 call makePressureField(finegrp,ivd,finegrp%u)
 finegrp%p = abs(finegrp%p)

 if(finegrp%gridNumber>1) then 
  if (ivd%coarseGridDissipationScheme==1) then
  call makeArtificialDissipation1(finegrp,ivd)
  else if(ivd%coarseGridDissipationScheme==2) then
   call makeArtificialDissipation2(finegrp,ivd)
  else  ! NOR
   call makeArtificialDissipation3(finegrp,ivd)
  end if
 else
  if (ivd%dissipationScheme==1) then
   call makeArtificialDissipation1(finegrp,ivd)
  else if(ivd%dissipationScheme==2) then 
   call makeArtificialDissipation2(finegrp,ivd)
  else  ! NOR
   call makeArtificialDissipation3(finegrp,ivd)
  end if
 end if


!  if (ivd%dissipationScheme==1.and.finegrp%gridNumber==1) then
!   call makeArtificialDissipation1(finegrp,ivd)
!  else 
!   call makeArtificialDissipation2(finegrp,ivd)
!  end if

 if(ivd%ReynoldsNumber>0.0) then
!  if(ivd%viscosityScheme==2.and.finegrp%gridNumber==1) then 
  if(ivd%viscosityScheme==2) then 
   call makeViscosityTerm2(finegrp,ivd,finegrp%dissipation,.true.)
  else
   call makeViscosityTerm(finegrp,ivd,finegrp%dissipation,.true.)
  end if
  if(ivd%turbulenceModel==2) then
   call makeBLTurbulence(finegrp,ivd)
  end if
 end if

 call makeRHS(finegrp,ivd)
 finegrp%rhs(:,5) = 0.0
! if(ivd%turbulenceModel==1.and.finegrp%gridNumber<2.and.calculateTurbulence) then ! KAS 
 if(ivd%turbulenceModel==1.and.finegrp%gridNumber<4.and.calculateTurbulence) then
  call makeSARHS(finegrp,ivd)
  call makeSADiffusion(finegrp,ivd)
  call makeSASourceTerm(finegrp,ivd)
  if(finegrp%gridNumber==1) then
   call makeSATripTerm(finegrp,ivd)
!  else
!   finegrp%rhs(:,5) = finegrp%rhs(:,5) - finegrp%tripSourceArray(:)
  end if
  finegrp%rhs(:,5) = ivd%turbulenceCFLFactor*finegrp%rhs(:,5)
 end if
 !call setRHSAtBoundary(finegrp,ivd)
 call solveSystem(finegrp,ivd,0,finegrp%rhs)
 if(associated(finegrp%sourceTerm)) then
  call setBCsOnIncrementField(finegrp,ivd,finegrp%rhs) ! CHANGED 01.09.99
  call nullifyOuterBoundary(finegrp,finegrp%rhs)
 else
  call setBCsOnIncrementField(finegrp,ivd,finegrp%rhs)
 end if
 finegrp%rhs = finegrp%uprev - finegrp%rhs

 do i=1,finegrp%numberOfNodes
  if(finegrp%rhs(i,5)<0.0) finegrp%rhs(i,5) = 0.0
 end do

 call setBCsOnSolutionField(finegrp,ivd,finegrp%rhs)
 if(finegrp%gridNumber>1) then
  finegrp%rhs = finegrp%rhs - finegrp%sourceTerm
 end if
 call setOuterBoundaryConditions3(finegrp,ivd,finegrp%rhs,finegrp%p,.true.)
 finegrp%rhs = finegrp%rhs - finegrp%uprev
 call mapFineToCoarse(finegrp%rhs,coarsegrp,coarsegrp%rhs)


 coarsegrp%sourceTerm = coarsegrp%sourceTerm-coarsegrp%rhs

 ! As a smart move, initialize coarse grid vector using second part of source term
 coarsegrp%u = coarsegrp%u + coarsegrp%rhs

 do i=1,coarsegrp%numberOfNodes
  coarsegrp%u(i,5) = max(coarsegrp%u(i,5),0.0)
 end do

 end subroutine restrictGrid
!-----------------------------------------------------------------------
 subroutine prolongateGrid(coarsegrp,finegrp,ivd) 
 ! Does a Full Approximation Storage (FAS) multigrid coarse to fine
 ! transformation
 IMPLICIT NONE

 type(GridSolverData) :: coarsegrp,finegrp
 type(InputVariablesData) :: ivd

 integer :: i

! character*30 :: timeBuffer
! integer :: time1(8),time2(8)

 ! first get the fine grid approximation mapped to the coarse grid 
 call mapFineToCoarse(finegrp%u,coarsegrp,coarsegrp%uprev)
 call setBCsOnSolutionField(coarsegrp,ivd,coarsegrp%uprev)


 ! now create the coarse grid correction v_H = u_H-I_h^H u_h
 coarsegrp%rhs = coarsegrp%u  - coarsegrp%uprev
! NOR call setOuterBoundaryConditions4(coarsegrp,ivd,coarsegrp%rhs)

! call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time1)
 if(ivd%prolongationScheme==1) then 
  call mapCoarseToFine(coarsegrp%rhs,coarsegrp,finegrp,finegrp%rhs)
 else if(ivd%prolongationScheme==2) then 
  call mapCoarseToFine2(coarsegrp%rhs,coarsegrp,finegrp,finegrp%rhs)
 else if(ivd%prolongationScheme==3) then 
  call mapCoarseToFine3(coarsegrp%rhs,coarsegrp,finegrp,finegrp%rhs)
 else
  STOP "ERROR: Undefined prolongation scheme"
 end if

! call nullifyIOChainIncrement(finegrp,finegrp%rhs)

! call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time2)
!  write(*,*) "pro2: ",time2(7)+0.001*time2(8)-time1(7)-0.001*time1(8)
 if(ivd%prolongationSmoothingFactor>0.0) then 
  do i=1,ivd%numberOfPSSteps
   call smoothResidual(finegrp,ivd,finegrp%rhs,ivd%prolongationSmoothingFactor)
  end do
 end if


! if(finegrp%gridNumber==1) then 
!  open(12,FILE="Corr.res")
!  do i=1,finegrp%numberOfNodes
!   write(12,'(I7,5E17.7)') i,finegrp%rhs(i,1:5)
!  end do
!  close(12)
!  PAUSE
! end if


! NOR call setOuterBoundaryConditions4(finegrp,ivd,finegrp%rhs)
 finegrp%u = finegrp%u+ivd%prolongationRelaxation*finegrp%rhs
 call setBCsOnSolutionField(finegrp,ivd,finegrp%u)
 call setOuterBoundaryConditions3(finegrp,ivd,finegrp%u,finegrp%p,.true.)
 end subroutine prolongateGrid 
!-----------------------------------------------------------------------
 subroutine mapFineToCoarse(fineu,coarsegrp,u)
 ! maps fine grid unknown fields to coarse grid
 IMPLICIT NONE 

 real :: fineu(:,:)
 type(GridSolverData) :: coarsegrp
 real :: u(:,:) 
 
 integer :: i,j
 type(LinkedReal),pointer :: currentReal
 type(LinkedInteger),pointer :: currentInteger

 ! initialize coarse grid

 u = 0.0
 ! run through all nodes in coarse grid (mappings are always done
 ! from the coarse grid) and add the coefficients found in 
 ! restriction array multiplied by field values of nodes found
 ! in prolongation array (which here is just used as a register)
 ! to fields in coarse grid 
  do i=1,coarsegrp%numberOfNodes 
   currentInteger=>coarsegrp%pod%prolongationArray(i)%first
   currentReal=>coarsegrp%rod%restrictionArray(i)%first
   do while(associated(currentReal))
    u(i,:) = u(i,:)&
                +currentReal%re*fineu(currentInteger%int,:)
    currentInteger=>currentInteger%next 
    currentReal=>currentReal%next
   end do  
  end do
 end subroutine mapFineToCoarse
!-----------------------------------------------------------------------
 subroutine mapCoarseToFine(ucoarse,coarsegrp,finegrp,u)
 ! maps coarse grid fields to fine grid
 ! here simple injection is used
 IMPLICIT NONE

 type(GridSolverData) :: coarsegrp,finegrp
 real :: ucoarse(:,:),u(:,:)

 integer :: i
 type(LinkedInteger),pointer :: currentInteger


 do i=1,coarsegrp%numberOfNodes
  currentInteger=>coarsegrp%pod%prolongationArray(i)%first
  do while(associated(currentInteger))
   u(currentInteger%int,:) = ucoarse(i,:)
   currentInteger=>currentInteger%next
  end do
 end do

 end subroutine mapCoarseToFine
!-----------------------------------------------------------------------
 subroutine mapCoarseToFine2(ucoarse,coarsegrp,finegrp,u)
 ! maps coarse grid fields to fine grid
 ! using linear interpolation
 IMPLICIT NONE

 type(GridSolverData) :: coarsegrp,finegrp
 real :: ucoarse(:,:),u(:,:)

 integer :: i,j
 type(LinkedInteger),pointer :: currentInteger

 integer :: i1,i2
 real :: u1(5),u2(5),ux(5),uy(5),x0(2),r(2),wx,wy
 real :: linInt

 integer :: ind

 ind = 83

 coarsegrp%nodeHelpArray(1:coarsegrp%numberOfNodes,1:10) = 0.0

 do i=1,coarsegrp%numberOfSides
 ! indexes of nodes in side
  i1 = coarsegrp%sideIndexArray(i,1)
  i2 = coarsegrp%sideIndexArray(i,2)

 ! side weights from preprocessor
  wx = coarsegrp%sideWeightsArray(i,1)
  wy = coarsegrp%sideWeightsArray(i,2)

  u1 = ucoarse(i1,:)
  u2 = ucoarse(i2,:)

  ux = (u1+u2)*wx
  uy = (u1+u2)*wy
 
  coarsegrp%nodeHelpArray(i1,1:5) = coarsegrp%nodeHelpArray(i1,1:5) + ux
  coarsegrp%nodeHelpArray(i1,6:10) = coarsegrp%nodeHelpArray(i1,6:10) + uy
  coarsegrp%nodeHelpArray(i2,1:5) = coarsegrp%nodeHelpArray(i2,1:5) - ux
  coarsegrp%nodeHelpArray(i2,6:10) = coarsegrp%nodeHelpArray(i2,6:10) - uy
 end do

 do i=1,coarsegrp%brp%numberOfBoundaryFaces
  i1 = coarsegrp%brp%faceIndexArray(i,1)
  i2 = coarsegrp%brp%faceIndexArray(i,2)
  wx   = coarsegrp%brp%faceWeightsArray(i,1)
  wy   = coarsegrp%brp%faceWeightsArray(i,2)
 
  u1 = ucoarse(i1,:)
  u2 = ucoarse(i2,:)
 
  coarsegrp%nodeHelpArray(i1,1:5) = coarsegrp%nodeHelpArray(i1,1:5) + 4.*u1*wx
  coarsegrp%nodeHelpArray(i1,6:10) = coarsegrp%nodeHelpArray(i1,6:10) + 4.*u1*wy
  coarsegrp%nodeHelpArray(i2,1:5) = coarsegrp%nodeHelpArray(i2,1:5) + 4.*u2*wx
  coarsegrp%nodeHelpArray(i2,6:10) = coarsegrp%nodeHelpArray(i2,6:10) + 4.*u2*wy
 end do

 do i=1,coarsegrp%numberOfNodes
  coarsegrp%nodeHelpArray(i,1:10) = coarsegrp%nodeHelpArray(i,1:10)/coarsegrp%nodeVolume(i)
 end do



! do i=coarsegrp%brp%numberOfCoarseBoundaryFaces+1,coarsegrp%numberOfNodes
 do i=1,coarsegrp%numberOfNodes
  currentInteger=>coarsegrp%pod%prolongationArray(i)%first
  x0 = coarsegrp%coordinates(i,:)
  do while(associated(currentInteger))
   r = finegrp%coordinates(currentInteger%int,:) - x0
   do j=1,5
    linInt = ucoarse(i,j) + coarsegrp%nodeHelpArray(i,j)*r(1)+coarsegrp%nodeHelpArray(i,j+5)*r(2)
!    if(linInt*ucoarse(i,j)>0.0) then 
    if(abs(linInt)<1.0*abs(ucoarse(i,j))) then 
!    if(.false.) then 
     u(currentInteger%int,j) = linInt
    else
     u(currentInteger%int,j) = ucoarse(i,j)
    end if
   end do

!   u(currentInteger%int,1) = ucoarse(i,1) + coarsegrp%nodeHelpArray(i,1)*r(1)+coarsegrp%nodeHelpArray(i,6)*r(2)
!   u(currentInteger%int,2) = ucoarse(i,2) + coarsegrp%nodeHelpArray(i,2)*r(1)+coarsegrp%nodeHelpArray(i,7)*r(2)
!   u(currentInteger%int,3) = ucoarse(i,3) + coarsegrp%nodeHelpArray(i,3)*r(1)+coarsegrp%nodeHelpArray(i,8)*r(2)
!   u(currentInteger%int,4) = ucoarse(i,4) + coarsegrp%nodeHelpArray(i,4)*r(1)+coarsegrp%nodeHelpArray(i,9)*r(2)
!   u(currentInteger%int,5) = ucoarse(i,5) + coarsegrp%nodeHelpArray(i,5)*r(1)+coarsegrp%nodeHelpArray(i,10)*r(2)
   currentInteger=>currentInteger%next
  end do
 end do
 

 end subroutine mapCoarseToFine2
!-----------------------------------------------------------------------
 subroutine mapCoarseToFine3(ucoarse,coarsegrp,finegrp,u)
 ! maps coarse grid fields to fine grid
 ! using linear interpolation
 IMPLICIT NONE

 type(GridSolverData) :: coarsegrp,finegrp
 real :: ucoarse(:,:),u(:,:)

 integer :: i,j
 type(LinkedInteger),pointer :: currentInteger

 integer :: i1,i2,ind1,ind2
 real :: u1(5),u2(5),ux(5),uy(5),x0(2),r(2),wx,wy
 real :: linInt,mapping(5)

 integer :: ind,numberOfLimits

 coarsegrp%nodeHelpArray(1:coarsegrp%numberOfNodes,1:10) = 0.0

 do i=1,coarsegrp%numberOfSides
 ! indexes of nodes in side
  i1 = coarsegrp%sideIndexArray(i,1)
  i2 = coarsegrp%sideIndexArray(i,2)

 ! side weights from preprocessor
  wx = coarsegrp%sideWeightsArray(i,1)
  wy = coarsegrp%sideWeightsArray(i,2)

  u1 = ucoarse(i1,:)
  u2 = ucoarse(i2,:)

  ux = (u1+u2)*wx
  uy = (u1+u2)*wy

  coarsegrp%nodeHelpArray(i1,1:5) = coarsegrp%nodeHelpArray(i1,1:5) + ux
  coarsegrp%nodeHelpArray(i1,6:10) = coarsegrp%nodeHelpArray(i1,6:10) + uy
  coarsegrp%nodeHelpArray(i2,1:5) = coarsegrp%nodeHelpArray(i2,1:5) - ux
  coarsegrp%nodeHelpArray(i2,6:10) = coarsegrp%nodeHelpArray(i2,6:10) - uy
 end do

 do i=1,coarsegrp%brp%numberOfBoundaryFaces
  i1 = coarsegrp%brp%faceIndexArray(i,1)
  i2 = coarsegrp%brp%faceIndexArray(i,2)
  wx   = coarsegrp%brp%faceWeightsArray(i,1)
  wy   = coarsegrp%brp%faceWeightsArray(i,2)

  u1 = ucoarse(i1,:)
  u2 = ucoarse(i2,:)

  coarsegrp%nodeHelpArray(i1,1:5) = coarsegrp%nodeHelpArray(i1,1:5) + 4.*u1*wx
  coarsegrp%nodeHelpArray(i1,6:10) = coarsegrp%nodeHelpArray(i1,6:10) + 4.*u1*wy
  coarsegrp%nodeHelpArray(i2,1:5) = coarsegrp%nodeHelpArray(i2,1:5) + 4.*u2*wx
  coarsegrp%nodeHelpArray(i2,6:10) = coarsegrp%nodeHelpArray(i2,6:10) + 4.*u2*wy
 end do

 do i=1,coarsegrp%numberOfNodes
  coarsegrp%nodeHelpArray(i,1:10) = coarsegrp%nodeHelpArray(i,1:10)/coarsegrp%nodeVolume(i)
 end do



 coarsegrp%nodeHelpArray(1:coarsegrp%numberOfNodes,11:15) = 0.0

 do i=1,coarsegrp%numberOfNodes
  currentInteger=>coarsegrp%pod%prolongationArray(i)%first
  x0 = coarsegrp%coordinates(i,:)
  do while(associated(currentInteger))
   r = finegrp%coordinates(currentInteger%int,:) - x0
   coarsegrp%nodeHelpArray(currentInteger%int,11) =  coarsegrp%nodeHelpArray(i,1)*r(1)+coarsegrp%nodeHelpArray(i,6)*r(2)
   coarsegrp%nodeHelpArray(currentInteger%int,12) =  coarsegrp%nodeHelpArray(i,2)*r(1)+coarsegrp%nodeHelpArray(i,7)*r(2)
   coarsegrp%nodeHelpArray(currentInteger%int,13) =  coarsegrp%nodeHelpArray(i,3)*r(1)+coarsegrp%nodeHelpArray(i,8)*r(2)
   coarsegrp%nodeHelpArray(currentInteger%int,14) =  coarsegrp%nodeHelpArray(i,4)*r(1)+coarsegrp%nodeHelpArray(i,9)*r(2)
   coarsegrp%nodeHelpArray(currentInteger%int,15) =  coarsegrp%nodeHelpArray(i,5)*r(1)+coarsegrp%nodeHelpArray(i,10)*r(2)
   currentInteger=>currentInteger%next
  end do
 end do

 ! find maximum neighbour values
 coarsegrp%nodeHelpArray(1:coarsegrp%numberOfNodes,1:5) = ucoarse(1:coarsegrp%numberOfNodes,1:5)
 coarsegrp%nodeHelpArray(1:coarsegrp%numberOfNodes,6:10) = ucoarse(1:coarsegrp%numberOfNodes,1:5)
! write(*,*) "a0: ",ucoarse(225,2)
 do i=1,coarsegrp%numberOfSides
  ind1 = coarsegrp%sideIndexArray(i,1)
  ind2 = coarsegrp%sideIndexArray(i,2)
!  if(ind1==225) write(*,*) "a: ",ucoarse(ind1,2),coarsegrp%nodeHelpArray(ind1,2),coarsegrp%nodeHelpArray(ind1,7)
!  if(ind2==225) write(*,*) "c: ",ucoarse(ind2,2),coarsegrp%nodeHelpArray(ind2,2),coarsegrp%nodeHelpArray(ind2,7)
  do j=1,5
   if(ucoarse(ind2,j)>coarsegrp%nodeHelpArray(ind1,j+5)) then 
    coarsegrp%nodeHelpArray(ind1,j+5) = ucoarse(ind2,j) 
   end if
   if(ucoarse(ind1,j)>coarsegrp%nodeHelpArray(ind2,j+5)) then
    coarsegrp%nodeHelpArray(ind2,j+5) = ucoarse(ind1,j)
   end if
   if(ucoarse(ind2,j)<coarsegrp%nodeHelpArray(ind1,j)) then
    coarsegrp%nodeHelpArray(ind1,j) = ucoarse(ind2,j)
   end if
   if(ucoarse(ind1,j)<coarsegrp%nodeHelpArray(ind2,j)) then
    coarsegrp%nodeHelpArray(ind2,j) = ucoarse(ind1,j)
   end if
  end do
!  if(ind1==225) write(*,*) "b: ",ucoarse(ind1,2),coarsegrp%nodeHelpArray(ind1,2),coarsegrp%nodeHelpArray(ind1,7)
!  if(ind2==225) write(*,*) "d: ",ucoarse(ind2,2),coarsegrp%nodeHelpArray(ind2,2),coarsegrp%nodeHelpArray(ind2,7)
 end do

 ! do mapping

 numberOfLimits = 0
 do i=1,coarsegrp%numberOfNodes
  do j=1,5
   mapping(j) = ucoarse(i,j) + coarsegrp%nodeHelpArray(i,j+10)
   if(mapping(j)<coarsegrp%nodeHelpArray(i,j)) then
!    write(*,*) "r: ",i,j,mapping(j),coarsegrp%nodeHelpArray(i,j),coarsegrp%nodeHelpArray(i,j+5)
!    mapping(j) = coarsegrp%nodeHelpArray(i,j)
!    mapping(j) = ucoarse(i,j)
    mapping(j) = ucoarse(i,j) + 0.25*(coarsegrp%nodeHelpArray(i,j)-ucoarse(i,j))
    numberOfLimits = numberOfLimits + 1
!    pause
   else if(mapping(j)>coarsegrp%nodeHelpArray(i,j+5)) then
!    write(*,*) "s: ",i,j,mapping(j),coarsegrp%nodeHelpArray(i,j),coarsegrp%nodeHelpArray(i,j+5)
!    mapping = coarsegrp%nodeHelpArray(i,j+5)
!    mapping(j) = ucoarse(i,j)
    mapping(j) = ucoarse(i,j) + 0.25*(coarsegrp%nodeHelpArray(i,j)-ucoarse(i,j))
    numberOfLimits = numberOfLimits + 1
   end if
  end do
  currentInteger=>coarsegrp%pod%prolongationArray(i)%first
  do while(associated(currentInteger))
   u(currentInteger%int,:) = mapping
   currentInteger=>currentInteger%next
  end do
 end do

 write(*,*) "Number of limits: ",numberOfLimits,5*finegrp%numberOfNodes

 end subroutine mapCoarseToFine3
!-----------------------------------------------------------------------
 subroutine mapScalarFineToCoarse(fineu,coarsegrp,u)
 ! maps fine grid unknown fields to coarse grid
 IMPLICIT NONE 

 real :: fineu(:)
 type(GridSolverData) :: coarsegrp
 real :: u(:) 
 
 integer :: i,j
 type(LinkedReal),pointer :: currentReal
 type(LinkedInteger),pointer :: currentInteger

 ! initialize coarse grid

 u = 0.0
 ! run through all nodes in coarse grid (mappings are always done
 ! from the coarse grid) and add the coefficients found in 
 ! restriction array multiplied by field values of nodes found
 ! in prolongation array (which here is just used as a register)
 ! to fields in coarse grid 
  do i=1,coarsegrp%numberOfNodes 
   currentInteger=>coarsegrp%pod%prolongationArray(i)%first
   currentReal=>coarsegrp%rod%restrictionArray(i)%first
   do while(associated(currentReal))
    u(i) = u(i)&
                +currentReal%re*fineu(currentInteger%int)
    currentInteger=>currentInteger%next 
    currentReal=>currentReal%next
   end do  
  end do
 end subroutine mapScalarFineToCoarse
!-----------------------------------------------------------------------
 subroutine mapScalarCoarseToFine(ucoarse,coarsegrp,finegrp,u)
 IMPLICIT NONE 

 type(GridSolverData) :: coarsegrp,finegrp
 real :: ucoarse(:),u(:)

 integer :: i
 type(LinkedInteger),pointer :: currentInteger


 do i=1,coarsegrp%numberOfNodes
  currentInteger=>coarsegrp%pod%prolongationArray(i)%first
  do while(associated(currentInteger))
   u(currentInteger%int) = ucoarse(i)
   currentInteger=>currentInteger%next
  end do
 end do

 end subroutine mapScalarCoarseToFine
!-----------------------------------------------------------------------
 subroutine makeTimeSteps(grp,ivd)
 ! make local time steps
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,i1,i2
 real :: c1,c2,c3,c7,wx,wy,al,ro1,ro2,cc1,cc2,vn1,vn2,t1,am1,t2,am2
 real :: vis1,vis2,dv1,dv2,vnm,ccm,dvm 
 real :: cvis1,cvis2,Ksi,mu,turbViscosityCoefficient,oneOverReynoldsNumber

 ! initialize
 grp%localTimeSteps = 0.0
 grp%nodeHelpArray(:,1) = 0.0

 c1 = ivd%gamma
 c2 = c1 - 1.0
 c3 = 0.
 if(ivd%ReynoldsNumber.gt.1.0e-6) c3 = c2**1.5*ivd%MachNumber**3/ivd%ReynoldsNumber
 c7 = ivd%MachNumber**2*c2*ivd%inflowTemperature

 if(ivd%ReynoldsNumber>0.0) then  
  oneOverReynoldsNumber = 1.0/ivd%ReynoldsNumber
 else
  oneOverReynoldsNumber = 0.0
 end if
! run over sides 

 do i=1,grp%numberOfSides       
  i1  = grp%sideIndexArray(i,1) 
  i2  = grp%sideIndexArray(i,2) 
  wx  = grp%sideWeightsArray(i,1)
  wy  = grp%sideWeightsArray(i,2)
  al  = wx*wx + wy*wy 
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  cc1 = c1*grp%p(i1)*ro1
  cc2 = c1*grp%p(i2)*ro2
  cc1 = max(cc1,0.0)
  cc2 = max(cc2,0.0)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) )*ro1 
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) )*ro2 
  grp%localTimeSteps(i1) = grp%localTimeSteps(i1) + abs(vn1) + sqrt(cc1*al)
  grp%localTimeSteps(i2) = grp%localTimeSteps(i2) + abs(vn2) + sqrt(cc2*al)

  t1  = c1*grp%p(i1)/(grp%u(i1,1)*c2)
  t2  = c1*grp%p(i2)/(grp%u(i2,1)*c2)
  if(ivd%ReynoldsNumber>0.0) then 
   am1 = c3*t1**1.5*(ivd%inflowTemperature+198.6)/(t1*c7+198.6)
   am2 = c3*t2**1.5*(ivd%inflowTemperature+198.6)/(t2*c7+198.6)
  else
   am1 = 0.0
   am2 = 0.0 
  end if
 
  if(ivd%turbulenceModel>0) then 
   if(ivd%turbulenceModel==1) then 
    mu = grp%laminarViscosity(i1)
    turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i1,5)
    Ksi = turbViscosityCoefficient/mu
    turbViscosityCoefficient =&
    turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)
    mu = mu + turbViscosityCoefficient

    cvis1 = abs(mu)

    mu = grp%laminarViscosity(i2)
    turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i2,5)
    Ksi = turbViscosityCoefficient/mu
    turbViscosityCoefficient =&
    turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)
    mu = mu + turbViscosityCoefficient

    cvis2 = abs(mu)
   end if
  else if(ivd%ReynoldsNumber>0.0) then 
   cvis1 = abs(grp%laminarViscosity(i1))
   cvis2 = abs(grp%laminarViscosity(i2))
  else
   cvis1 = 0.0
   cvis2 = 0.0
  end if

!  vis1= amax1(am1,abs(grp%laminarViscosity(i1)))
!  vis2= amax1(am2,abs(grp%laminarViscosity(i2)))
  vis1= amax1(am1,cvis1)
  vis2= amax1(am2,cvis2)

  dv1 = 2.0*vis1*al/grp%u(i1,1)
  dv2 = 2.0*vis2*al/grp%u(i2,1)
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv2
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv1
 end do

! run over faces 

 do i=1,grp%brp%numberOfBoundaryFaces  
  i1  = grp%brp%faceIndexArray(i,1)
  i2  = grp%brp%faceIndexArray(i,2)
  wx  = grp%brp%faceWeightsArray(i,1)
  wy  = grp%brp%faceWeightsArray(i,2)
  al  = wx*wx + wy*wy 
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  cc1 = ivd%gamma*grp%p(i1)*ro1
  cc2 = ivd%gamma*grp%p(i2)*ro2
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) )*ro1 
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) )*ro2 
  vn1 = abs(vn1)
  vn2 = abs(vn2)
  vnm = vn1 + vn2 
  cc1 = sqrt(cc1*al)
  cc2 = sqrt(cc2*al)
  ccm = cc1 + cc2 


  grp%localTimeSteps(i1) = grp%localTimeSteps(i1) + vn1 + vnm + cc1 + ccm 
  grp%localTimeSteps(i2) = grp%localTimeSteps(i2) + vn2 + vnm + cc2 + ccm 

 
  t1  = c1*grp%p(i1)/(grp%u(i1,1)*c2)
  t2  = c1*grp%p(i2)/(grp%u(i2,1)*c2)
  if(ivd%ReynoldsNumber>0.0) then
   am1 = c3*t1**1.5*(ivd%inflowTemperature+198.6)/(t1*c7+198.6)
   am2 = c3*t2**1.5*(ivd%inflowTemperature+198.6)/(t2*c7+198.6)
  else
   am1 = 0.0
   am2 = 0.0
  end if

  if(ivd%turbulenceModel>0) then
   if(ivd%turbulenceModel==1) then
    mu = grp%laminarViscosity(i1)
    turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i1,5)
    Ksi = turbViscosityCoefficient/mu
    turbViscosityCoefficient =&
    turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)
    mu = mu + turbViscosityCoefficient

    cvis1 = abs(mu)

    mu = grp%laminarViscosity(i2)
    turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i2,5)
    Ksi = turbViscosityCoefficient/mu
    turbViscosityCoefficient =&
    turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)
    mu = mu + turbViscosityCoefficient

    cvis2 = abs(mu)
   end if
  else if(ivd%ReynoldsNumber>0.0) then
   cvis1 = abs(grp%laminarViscosity(i1))
   cvis2 = abs(grp%laminarViscosity(i2))
  else
   cvis1 = 0.0
   cvis2 = 0.0
  end if


!  vis1= amax1(am1,abs(grp%laminarViscosity(i1)))
!  vis2= amax1(am2,abs(grp%laminarViscosity(i2)))
  vis1= amax1(am1,cvis1)  
  vis2= amax1(am2,cvis2)  
  
  dv1 = 2.0*vis1*al/grp%u(i1,1)
  dv2 = 2.0*vis2*al/grp%u(i2,1)
  dvm = dv1 + dv2
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv1 + dvm
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv2 + dvm
 end do


 do i = 1,grp%numberOfNodes
  grp%localInviscidTimeSteps(i) = 0.5/&
    grp%localTimeSteps(i)/grp%nodeVolume(i)
  grp%localTimeSteps(i) = 0.5/&
    (grp%localTimeSteps(i)/grp%nodeVolume(i)&
     +4.0*grp%nodeHelpArray(i,1)/(grp%nodeVolume(i)**2))
 end do

 end subroutine makeTimeSteps 
!-----------------------------------------------------------------------
 subroutine makePressureField(grp,ivd,u)
 ! makes pressure field from the ideal gas equation
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real,pointer :: u(:,:)
 integer :: i
 real :: gammaMinusOne,twoThirds


 gammaMinusOne = ivd%gamma-1.0
 grp%p = gammaMinusOne*(u(:,4)&
                 -0.5*(u(:,2)**2+u(:,3)**2)/u(:,1)) 
 end subroutine makePressureField
!-----------------------------------------------------------------------
 subroutine smoothResidual(grp,ivd,rhs,sfactor)  
 ! smooths residual to increase allowable CFL number

  IMPLICIT NONE

  type(GridSolverData) :: grp
  type(InputVariablesData) :: ivd
  real :: rhs(:,:)
  real :: sfactor

  integer :: i,i1,i2,j
  real :: dw

  grp%nodeHelpArray(1:grp%numberOfNodes,1) = 0.0 ! initialize 

! make help variables
  do i=1,4 ! number of PDE's 
   do j=1,grp%numberOfSides
    i1 = grp%sideIndexArray(j,1)
    i2 = grp%sideIndexArray(j,2)
    dw = rhs(i1,i) - rhs(i2,i)
    grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) - dw
    grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dw
   end do

   do j=1,grp%numberOfNodes
    rhs(j,i) = (rhs(j,i)+sfactor*grp%nodeConnectivityArray(j)*rhs(j,i)+&
       sfactor*grp%nodeHelpArray(j,1))/(1.+grp%nodeConnectivityArray(j)*sfactor)
   end do
  end do
 end subroutine smoothResidual
!-------------------------------------------------------------------------
 subroutine findNodeDistanceFromWall(grp,ivd)
 ! finds the shortest distance from the nodes to a wall boundary
 ! this is for the moment only done on the finest grid but can 
 ! easily be implemented for agglomerated grid using weihghted averages
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,k,indicator,faceIndex,allocateStatus,nodeOnFace,ist,ien
 real :: minimumDistance,dist1,dist2,dist3,xp(2),x1(2),x2(2),tmin
 logical :: gridIsViscid

 write(*,*) "Finding node distance from wall..."

 allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
 allocate(grp%wallDistanceFaceArray(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"

 grp%wallDistance = 0.0
 grp%wallDistanceFaceArray = 0

 do i=1,grp%numberOfNodes
  xp = grp%coordinates(i,:)
  faceIndex = 0
  if(xp(1)>1.0) then  ! C
   if(xp(1)<4.0) then 
    minimumDistance = xp(2)*xp(2) ! distance from wake ! C 
    if(xp(2)>0) then
     faceIndex = anint(-ivd%wakeSectionDensity*(xp(1)-1.0)*ivd%numberOfWakeSections)
     if(faceIndex<-ivd%numberOfWakeSections) faceIndex = -ivd%numberOfWakeSections
    else
     faceIndex = anint(-ivd%wakeSectionDensity*(xp(1)-1.0)*ivd%numberOfWakeSections-ivd%numberOfWakeSections-1)
     if(faceIndex<-2*ivd%numberOfWakeSections-1) faceIndex = -2*ivd%numberOfWakeSections-1
    end if
   else
    xp(1) = grp%coordinates(i,1) - 4.0
    xp(2) = grp%coordinates(i,2)
    minimumDistance = sqrt(sum(xp*xp))
    if(xp(2)>0) then 
     faceIndex = -ivd%numberOfWakeSections
    else
     faceIndex = -2*ivd%numberOfWakeSections-1
    end if
   end if
  else ! C
   minimumDistance = 1.0e10
  end if ! C


  if(ivd%turbulenceModel==1) then 
   minimumDistance = 1.0e10  ! changed
  end if

  ist = grp%brp%faceIndicatorArray(-12)
  ien = grp%brp%faceIndicatorArray(-9)
  do j=ist,ien 
   x1 = grp%coordinates(grp%brp%faceIndicatorArray(j),:)
   x2 = grp%coordinates(grp%brp%faceIndicatorArray(j+1),:)
   tmin =((x2(1)-x1(1))*(xp(1)-x1(1))+(x2(2)-x1(2))*(xp(2)-x1(2)))/sum((x2-x1)*(x2-x1)) 
   if(tmin>=0.0.and.tmin<=1.0) then 
    ! there is an internal minimum
    dist1 = sum((xp-x1-(x2-x1)*tmin)*(xp-x1-(x2-x1)*tmin)) 
    if(tmin>0.5) then 
     nodeOnFace = 2
    else
     nodeOnFace = 1
    end if
   else
    dist2 = sum((xp-x1)*(xp-x1))
    dist3 = sum((xp-x2)*(xp-x2))
    dist1 = min(dist2,dist3)
    if(dist1==dist3) then 
     nodeOnFace = 2
    else
     nodeOnFace = 1
    end if
   end if    
   if(minimumDistance>dist1) then 
    minimumDistance = dist1
    faceIndex = grp%brp%faceIndicatorArray(j+nodeOnFace-1)
   end if
  end do
  grp%wallDistance(i) = sqrt(minimumDistance)
  grp%wallDistanceFaceArray(i) = faceIndex
 end do
 end subroutine findNodeDistanceFromWall
!-----------------------------------------------------------------------
 subroutine findWallDistanceForCoarseMesh(grp,hgrp)
  IMPLICIT NONE
  type(GridSolverData) :: hgrp,grp

  type(LinkedInteger),pointer :: currentNode
  type(LinkedReal),pointer :: currentReal
  integer :: i,allocateStatus

  allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
  allocate(grp%wallDistanceFaceArray(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
  grp%wallDistance = 0.0 
  do i=1,grp%numberOfNodes
   currentNode => grp%pod%prolongationArray(i)%first
   currentReal => grp%rod%restrictionArray(i)%first
   grp%wallDistanceFaceArray(i) = hgrp%wallDistanceFaceArray(currentNode%int)
   do while(associated(currentNode)) 
    grp%wallDistance(i) = grp%wallDistance(i)+&
                            hgrp%wallDistance(currentNode%int)*currentReal%re
    currentNode => currentNode%next
    currentReal => currentReal%next
   end do
  end do  

 end subroutine findWallDistanceForCoarseMesh
!-----------------------------------------------------------------------
 subroutine readGridComputationData(isInitial,finegrp,hgrp,grp,ivd,INFILE)
! reads data from file, created by preprocessor 
  IMPLICIT NONE

  logical :: isInitial
  type(GridSolverData) :: finegrp,hgrp,grp
  type(InputVariablesData) :: ivd
  integer :: INFILE

  integer :: i,ind1,ind2,j,k,allocateStatus,ip,coarseNodeNumber,ist,ien,ib
  real :: dx,dy,tx,ty,ww,wx(2)
  real :: x(2)
 
  if(isInitial) then ! finest mesh
   nullify(grp%boundaryFaceNodeMappings) 
   nullify(grp%localFaceTangentArray)
   grp%numberOfBoundaryCVs = grp%brp%numberOfBoundaryFaces

   ! read sides
   read(INFILE) grp%numberOfSides
   allocate(grp%sideIndexArray(grp%numberOfSides,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfSides
    read(INFILE) grp%sideIndexArray(i,1:2)
   end do

   ! read side weights
   allocate(grp%sideWeightsArray(grp%numberOfSides,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfSides
    read(INFILE) grp%sideWeightsArray(i,1:2)
   end do
  ! read side lengths into sideLengthsArray where the last component is
  ! the length of side and the two first are components of side tangent
  ! vector
  allocate(grp%sideLengthArray(grp%numberOfSides,3),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  do i=1,grp%numberOfSides
   read(INFILE) grp%sideLengthArray(i,1:2)
  end do

  ! calculate length of sides
  do i=1,grp%numberOfSides
   dx = grp%sideLengthArray(i,1)
   dy = grp%sideLengthArray(i,2)
   grp%sideLengthArray(i,3) = sqrt(dx*dx+dy*dy)
   if(grp%sideLengthArray(i,3)==0.0) then
    write(*,*) "ERROR: side ",i," has zero length"
   end if
  end do

  read(INFILE) grp%numberOfNodes
  print * , "Number of nodes :", grp%numberOfNodes
  allocate(grp%nodeVolume(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  do i=1,grp%numberOfNodes
   read(INFILE) grp%nodeVolume(i)
  end do
  ! read coordinates

  allocate(grp%coordinates(grp%numberOfNodes,2),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  do i=1,grp%numberOfNodes
   read(INFILE) grp%coordinates(i,1:2)
  end do

  allocate(grp%wallLength(grp%brp%numberOfBoundaryFaces),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  grp%wallLength = 0.0
  do i=1,grp%numberOfSides
   if(grp%sideIndexArray(i,1)<=grp%brp%numberOfBoundaryFaces.and.&
      grp%sideIndexArray(i,2)<=grp%brp%numberOfBoundaryFaces) then
    grp%wallLength(grp%sideIndexArray(i,1)) =&
     grp%wallLength(grp%sideIndexArray(i,1)) + grp%sideLengthArray(i,3)
    grp%wallLength(grp%sideIndexArray(i,2)) =&
     grp%wallLength(grp%sideIndexArray(i,2)) + grp%sideLengthArray(i,3)
   end if
  end do
  grp%wallLength = 0.5*grp%wallLength

  allocate(grp%brp%sideLengthArray(grp%brp%numberOfBoundaryFaces),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  do i=1,grp%brp%numberOfBoundaryFaces
   ind1 = grp%brp%faceIndexArray(i,1)
   ind2 = grp%brp%faceIndexArray(i,2)

   dx = grp%coordinates(ind2,1) - grp%coordinates(ind1,1)
   dy = grp%coordinates(ind2,2) - grp%coordinates(ind1,2)
 
   grp%brp%sideLengthArray(i) = sqrt(dx*dx+dy*dy)
  end do

  print * , 'find separation fields'
  call findSeparationFields(grp,ivd)
  print * , 'set IO coordinates'
  call setIOCoordinates(grp%brp,grp%coordinates)

!  if(ivd%turbulenceModel>0) then 
   print * , 'find node distance from wall'
   call findNodeDistanceFromWall(grp,ivd)
!  else
!   nullify(grp%wallDistance)
!   nullify(grp%wallDistanceFaceArray)
!  end if

  ! place subroutines that need coordinates here
 
!  allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
!  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
!  allocate(grp%wallDistanceFaceArray(grp%numberOfNodes),stat=allocateStatus)
!  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
!  do i=1,grp%numberOfNodes
!   read(INFILE) grp%wallDistanceFaceArray(i),grp%wallDistance(i)
!   write(*,*) i,grp%wallDistanceFaceArray(i),grp%wallDistance(i)
!  end do

  ! make node connectivity array
  allocate(grp%nodeConnectivityArray(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%nodeConnectivityArray = 0
  do i=1,grp%numberOfSides
   ind1 = grp%sideIndexArray(i,1)
   ind2 = grp%sideIndexArray(i,2)
   grp%nodeConnectivityArray(ind1) = grp%nodeConnectivityArray(ind1) + 1
   grp%nodeConnectivityArray(ind2) = grp%nodeConnectivityArray(ind2) + 1
  end do
  allocate(grp%AfterAndBefore(grp%numberOfNodes,2),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  grp%AfterAndBefore = 0
  
  allocate(grp%hoc(grp%numberOfSides,4),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%hoc = 0.

  allocate(grp%hoc2(grp%numberOfSides,4),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%hoc2 = 0.

!temporary memory allocation
  allocate(grp%hoc3(grp%numberOfSides,4),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%hoc3 = 0.

  print * , "make surface plot data"
  call makeSurfacePlotData(grp)

  else ! coarser meshes

   ! read boundary mapping which decides which control volumes in the
   ! coarse grid a boundary face shall contribute to (the boundary faces
   ! are not agglomerated to avoid problems with different boundary
   ! conditions etc.)
   allocate(grp%boundaryFaceNodeMappings(grp%brp%numberOfBoundaryFaces),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%brp%numberOfBoundaryFaces
    read(INFILE) ind1,grp%boundaryFaceNodeMappings(ind1)
   end do

   allocate(grp%localFaceTangentArray(grp%brp%numberOfBoundaryFaces,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
   grp%localFaceTangentArray = 0.0

   grp%brp%numberOfCoarseBoundaryFaces = 0
   do i =1,grp%brp%numberOfBoundaryFaces
    ip = grp%brp%faceIndicatorArray(i)
    tx = -1.0*grp%brp%faceTangentArray(ip,1)
    ty = -1.0*grp%brp%faceTangentArray(ip,2)
    coarseNodeNumber = grp%boundaryFaceNodeMappings(ip)
    grp%brp%numberOfCoarseBoundaryFaces = max(grp%brp%numberOfCoarseBoundaryFaces,coarseNodeNumber)
    ww = grp%brp%faceWeightNorms(ip)
    grp%localFaceTangentArray(coarseNodeNumber,1) = grp%localFaceTangentArray(coarseNodeNumber,1) + ww*tx
    grp%localFaceTangentArray(coarseNodeNumber,2) = grp%localFaceTangentArray(coarseNodeNumber,2) + ww*ty
   end do
 
   do i=1,grp%brp%numberOfCoarseBoundaryFaces
    grp%localFaceTangentArray(i,:) = grp%localFaceTangentArray(i,:)/&
         (sqrt(sum(grp%localFaceTangentArray(i,:)*grp%localFaceTangentArray(i,:))))
   end do

   nullify(grp%brp%sideLengthArray)

   read(INFILE) grp%numberOfBoundaryCVs
 
   ! read sides
   read(INFILE) grp%numberOfSides
   allocate(grp%sideIndexArray(grp%numberOfSides,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfSides
    read(INFILE) grp%sideIndexArray(i,1:2)
   end do

   ! read side weights
   allocate(grp%sideWeightsArray(grp%numberOfSides,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfSides
    read(INFILE) grp%sideWeightsArray(i,:)
   end do

   nullify(grp%tripNodeFieldIndexes)
   nullify(grp%tripNodeFieldDistances)

   ! read side lengths and inverted side lengths for coarse grid
   allocate(grp%sideLengthArray(grp%numberOfSides,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfSides
    read(INFILE) grp%sideLengthArray(i,1:3)
   end do
   allocate(grp%inverseSideLengthArray(grp%numberOfSides,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfSides
    read(INFILE) grp%inverseSideLengthArray(i,1:3)
   end do

   ! read lumped mass matrix
   read(INFILE) grp%numberOfNodes
   allocate(grp%nodeVolume(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%coordinates(grp%numberOfNodes,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   do i=1,grp%numberOfNodes
    read(INFILE) grp%coordinates(i,:),grp%nodeVolume(i)
   end do

!   allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
!   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
!   allocate(grp%wallDistanceFaceArray(grp%numberOfNodes),stat=allocateStatus)
!   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
!   do i=1,grp%numberOfNodes
!    read(INFILE) grp%wallDistanceFaceArray(i),grp%wallDistance(i)
!    write(*,*) i,grp%wallDistanceFaceArray(i),grp%wallDistance(i)
!   end do

   ! make node connectivity array
   allocate(grp%nodeConnectivityArray(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   grp%nodeConnectivityArray = 0
   do i=1,grp%numberOfSides
    ind1 = grp%sideIndexArray(i,1)
    ind2 = grp%sideIndexArray(i,2)
    grp%nodeConnectivityArray(ind1) = grp%nodeConnectivityArray(ind1) + 1
    grp%nodeConnectivityArray(ind2) = grp%nodeConnectivityArray(ind2) + 1
   end do

   ! read prolongation operator (coarse to fine)
   allocate(grp%pod,stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to create prolongation operator "
   grp%pod%numberOfNodes = grp%numberOfNodes
   call readProlongationOperatorData(grp%pod,INFILE)

   ! read restriction operator (fine to coarse)
   allocate(grp%rod,stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to create restriction operator "
   grp%rod%numberOfNodes = grp%numberOfNodes
   call readRestrictionOperatorData(grp%rod,INFILE)
   nullify(grp%wallLength)
   
!   if(ivd%turbulenceModel>0) then 
    call findWallDistanceForCoarseMesh(grp,hgrp)
!   else
!    nullify(grp%wallDistance)
!    nullify(grp%wallDistanceFaceArray)
!   end if

    nullify(grp%plotNodePairs)
    nullify(grp%plotNodeDistance)
    grp%numberOfPlotNodes = 0


  ! read coordinates

   allocate(grp%coordinates(grp%numberOfNodes,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

   call mapFineToCoarse(hgrp%coordinates,grp,grp%coordinates)


   grp%sideLengthArray = 0.0
!   allocate(grp%sideLengthArray(grp%numberOfSides,3),stat=allocateStatus)
!   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  ! calculate length of sides
   do i=1,grp%numberOfSides
    ind1 = grp%sideIndexArray(i,1)
    ind2 = grp%sideIndexArray(i,2)
    dx = grp%coordinates(ind2,1)-grp%coordinates(ind1,1)
    dy = grp%coordinates(ind2,2)-grp%coordinates(ind1,2)
    grp%sideLengthArray(i,1) = dx
    grp%sideLengthArray(i,2) = dy
    grp%sideLengthArray(i,3) = sqrt(dx*dx+dy*dy)
    if(grp%sideLengthArray(i,3)==0.0) then
     write(*,*) "ERROR: side ",i," has zero length"
    end if
   end do


  end if

! do some allocation

  ! unknowns at this and previous time levels

  allocate(grp%u(grp%numberOfNodes,5),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%uprev(grp%numberOfNodes,5),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%dissipation(grp%numberOfNodes,4),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate dissipation field"

  allocate(grp%rhs(grp%numberOfNodes,5),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

 ! pressure field
  allocate(grp%p(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate pressure field"

  allocate(grp%wakeDataArray(0:2*ivd%numberOfWakeSections+1,4),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate pressure field"

  allocate(grp%localTimeSteps(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  allocate(grp%localInviscidTimeSteps(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  allocate(grp%laminarViscosity(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate viscosity field"
  grp%laminarViscosity = 0.0

  allocate(grp%vorticity(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate vorticity field"

  allocate(grp%divergence(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate vorticity field"


  allocate(grp%tripSourceArray(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate vorticity field"

 
!  if(isInitial) then 
   allocate(grp%wallStress(grp%brp%numberOfBoundaryFaces,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
!  else
!   nullify(grp%wallStress) 
!  end if

  ! if initial grid, allocate help arrays, else point to original 
  if (isInitial) then 
   nullify(grp%sourceTerm)
   nullify(grp%turbMGSourceTerm)
   allocate(grp%nodeHelpArray(grp%numberOfNodes,16),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate help field"
   allocate(grp%isAdiabaticBoundary(grp%brp%numberOfBoundaryFaces),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate help field"
   grp%isAdiabaticBoundary = .false.
   ist = grp%brp%faceIndicatorArray(-10)
   ien = grp%brp%faceIndicatorArray(-9)
   do ib =ist,ien
    ip = grp%brp%faceIndicatorArray(ib)
    grp%isAdiabaticBoundary(ip) = .false.   
   end do

!   allocate(grp%sideHelpArray(grp%numberOfSides,4),stat=allocateStatus)
!   if(allocateStatus/=0)&
!      STOP "ERROR: Not enough memory to allocate help field"
   nullify(grp%sideHelpArray)
  else
   ! allocate FAS source term
   allocate(grp%sourceTerm(grp%numberOfNodes,5),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate source term field"
   allocate(grp%turbMGSourceTerm(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate source term field"

   ! to avoid using memory to create help variables for coarse grids
   grp%nodeHelpArray => finegrp%nodeHelpArray
   grp%sideHelpArray => finegrp%sideHelpArray
   grp%isAdiabaticBoundary => finegrp%isAdiabaticBoundary
  end if


! set Jameson timesteps
  
  allocate(grp%JamesonCoefficients(3),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%JamesonCoefficients(1) = 0.6
  grp%JamesonCoefficients(2) = 0.6 
  grp%JamesonCoefficients(3) = 1.0 
  
  grp%residualScalingFactor = 1.0/float(grp%numberOfNodes)

! if(.not.associated(grp%sourceTerm)) then 
! open(15,file="wall_distance.plt",form='formatted',status='unknown')
! do i=1,grp%numberOfNodes
!  write(15,*) i,grp%wallDistance(i),0,0,0,0,0,0
! end do 
! close(15)
! end if
! Let's do a test
!if(.false.) then 
  grp%u(:,1:2) = 0.0
  do i=1,grp%numberOfSides
   ind1 = grp%sideIndexArray(i,1)
   ind2 = grp%sideIndexArray(i,2)
   grp%u(ind1,1:2) = grp%u(ind1,1:2) + grp%sideWeightsArray(i,:)
   grp%u(ind2,1:2) = grp%u(ind2,1:2) - grp%sideWeightsArray(i,:)
  end do  

  do i=1,grp%brp%numberOfBoundaryFaces
   ind1 = grp%brp%faceIndexArray(i,1)
   ind2 = grp%brp%faceIndexArray(i,2)

   grp%u(ind1,1:2) = grp%u(ind1,1:2) + 2.0*grp%brp%faceWeightsArray(i,:)
   grp%u(ind2,1:2) = grp%u(ind2,1:2) + 2.0*grp%brp%faceWeightsArray(i,:)
  end do 

  do i=1,grp%numberOfNodes
   if(abs(grp%u(i,1))>1.0e-10.or.abs(grp%u(i,2))>1.0e-10) then 
!    write(*,*) "Side error: ",i,grp%u(i,1:2)
   end if
  end do 

 ! calculate the coefficients of the high-order 
 ! discretisations of the variables involved 
 ! in the sides of the boundary layers
! call calculateHOCoefficients(grp(1),ivd)

!end if 

 grp%brp%engineInletAreas = 0.0
 grp%brp%engineInletNormals = 0.0
 if(grp%gridNumber==1) then
  do i=1,grp%brp%numberOfEngineInletSides
   ind1 = grp%brp%engineInletSideIndexes(i,1)
   ind2 = grp%brp%engineInletSideIndexes(i,2)
   wx = grp%brp%engineInletSideCoefficients(i,:)
   grp%brp%engineInletAreas = grp%brp%engineInletAreas + sqrt(sum(wx*wx)) 
   grp%brp%engineInletNormals = grp%brp%engineInletNormals + wx 
  end do
 end if
 

 end subroutine readGridComputationData
!------------------------------------------------------------------------
 subroutine calculateRoeMatrix(i,rnx,rny,rhor,rhol,uxr,uxl,vyr,vyl,hr,hl,&
            epsr,epsl,gamma,z,grp)
 ! used in setBoundaryConditions

 IMPLICIT NONE

 integer :: i
 real :: rnx,rny,rhor,rhol,uxr,uxl,vyr,vyl,hr,hl,gamma,epsr,epsl,z
 type(GridSolverData) :: grp

 real :: gam1,di,diinv,ui,vi,hi,q2,ci2,ci,ucap,vcap,cx,cy,gcap,rlam1,rlam2,rlam3,rlam4
 real :: zz,rinv11,rinv12,rinv13,rinv14,rinv21,rinv22,rinv23,rinv24,rinv31,rinv32,rinv33
 real :: rinv34,rinv41,rinv42,rinv43,rinv44,r11,r12,r13,r14,r21,r22,r23,r24,r31,r32,r33
 real :: r34,r41,r42,r43,r44,dr,du,dv,de,w1,w2,w3,w4,roes1,roes2,roes3,roes4,fakm

 gam1 = gamma-1

 di=sqrt(rhor/rhol)
 diinv=1.0/(di+1.0)
 ui=(di*uxr+uxl)*diinv
 vi=(di*vyr+vyl)*diinv
 hi=(di*hr+hl)*diinv

 q2=ui*ui+vi*vi        
 ci2=gam1*(hi-q2/2.) ! removed ekini
 ci=sqrt(ci2)
 ucap= ui*rnx+vi*rny
 vcap=-ui*rny+vi*rnx
 cx=ci*rnx
 cy=ci*rny
 gcap=gam1/ci2

 rlam1=abs(ucap)
 rlam2=abs(ucap)
 rlam3=abs(ucap+ci)
 rlam4=abs(ucap-ci)

 zz=z

!if(rlam1.lt.zz) rlam1=0.5*(rlam1*rlam1/zz+zz)
!if(rlam2.lt.zz) rlam2=0.5*(rlam2*rlam2/zz+zz)
!if(rlam3.lt.zz) rlam3=0.5*(rlam3*rlam3/zz+zz)
!if(rlam4.lt.zz) rlam4=0.5*(rlam4*rlam4/zz+zz)

rinv11=1.0-gcap*q2/2.
rinv12=gcap*ui
rinv13=gcap*vi
rinv14=-gcap
rinv21=-vcap
rinv22=-rny
rinv23=rnx
rinv24=0.
rinv31=gcap*q2/4.-ucap/(2.*ci)
rinv32=(rnx/ci-gcap*ui)/2.
rinv33=(rny/ci-gcap*vi)/2.
rinv34=gcap/2.
rinv41=gcap*q2/4.+ucap/(2.*ci)
rinv42=(-rnx/ci-gcap*ui)/2.
rinv43=(-rny/ci-gcap*vi)/2.
rinv44=gcap/2.

r11=1.
r12=0.
r13=1.
r14=1.
r21=ui
r22=-rny
r23=ui+rnx*ci
r24=ui-rnx*ci
r31=vi
r32=rnx
r33=vi+rny*ci
r34=vi-rny*ci
r41=q2/2. ! removed ekini
r42=vcap
r43=hi+ucap*ci
r44=hi-ucap*ci

dr = rhor       - rhol
du = rhor*uxr   - rhol*uxl
dv = rhor*vyr   - rhol*vyl
de = epsr       - epsl

w1=rinv11*dr+rinv12*du+rinv13*dv+rinv14*de
w2=rinv21*dr+rinv22*du+rinv23*dv+rinv24*de
w3=rinv31*dr+rinv32*du+rinv33*dv+rinv34*de
w4=rinv41*dr+rinv42*du+rinv43*dv+rinv44*de

w1=w1*rlam1
w2=w2*rlam2
w3=w3*rlam3
w4=w4*rlam4

roes1=r11*w1+r12*w2+r13*w3+r14*w4
roes2=r21*w1+r22*w2+r23*w3+r24*w4
roes3=r31*w1+r32*w2+r33*w3+r34*w4 
roes4=r41*w1+r42*w2+r43*w3+r44*w4 

fakm = 1.0
grp%nodeHelpArray(i,1) = roes1*fakm
grp%nodeHelpArray(i,2) = roes2*fakm
grp%nodeHelpArray(i,3) = roes3*fakm
grp%nodeHelpArray(i,4) = roes4*fakm

 end subroutine calculateRoeMatrix
!-----------------------------------------------------------------------
 subroutine calculateLiftAndDrag(grp,ivd,lift,drag,momentum,pressureLift,pressureDrag,frictionDrag)
 ! as the name says
 IMPLICIT NONE
 
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real :: lift,drag,momentum,pressureDrag,pressureLift,frictionDrag

 integer :: i,ist,ien,ip
 real :: flowTangent(2),flowNormal(2),wallTangent(2),wallNormal(2),r(2) 
 real :: p,tau_t,tau_n,PI,alphaRad,sum

 lift = 0.0
 drag = 0.0
 momentum = 0.0

 PI = 4.0*atan(1.0)
 alphaRad = ivd%alpha*PI/180. 

 flowTangent(1) = cos(alphaRad)
 if (abs(flowTangent(1)) < 1.d-15) then
    flowTangent(1) = 0.D0
 endif
 flowTangent(2) = sin(alphaRad) 
 if (abs(flowTangent(2)) < 1.d-15) then
    flowTangent(2) = 0.D0
 endif
 flowNormal(1) = - flowTangent(2)
 flowNormal(2) = flowTangent(1)

 ! inviscid faces
 ist = grp%brp%faceIndicatorArray(-20)
 ien = grp%brp%faceIndicatorArray(-19)
 do i=ist,ien 
 ip = grp%brp%faceIndicatorArray(i)
  wallTangent(1) = grp%brp%faceTangentArray(ip,1)
  wallTangent(2) = grp%brp%faceTangentArray(ip,2)
  wallNormal(1) = -wallTangent(2)
  wallNormal(2) = wallTangent(1) 
  p = grp%p(ip)  
  r = grp%coordinates(ip,:)-ivd%momentumPoint
  lift = lift - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowNormal(1)+wallNormal(2)*flowNormal(2))
  drag = drag - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowTangent(1)+wallNormal(2)*flowTangent(2))
  momentum = momentum - 2.0*grp%wallLength(ip)*p*(r(1)*wallNormal(2)-r(2)*wallNormal(1))

  pressureLift = pressureLift - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowNormal(1)+wallNormal(2)*flowNormal(2))
  pressureDrag = pressureDrag - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowTangent(1)+wallNormal(2)*flowTangent(2))

  if(ivd%ReynoldsNumber>0.0) then
  lift = lift - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowNormal(1)+grp%wallStress(ip,2)*flowNormal(2)) 
  drag = drag - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowTangent(1)+grp%wallStress(ip,2)*flowTangent(2))
  frictionDrag = frictionDrag - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowTangent(1)+grp%wallStress(ip,2)*flowTangent(2))
  momentum = momentum - 2.0*grp%wallLength(ip)*(r(1)*grp%wallStress(ip,2)-r(2)*grp%wallStress(ip,1))
   end if
 end do

 ! isothermal viscous faces
 ist = grp%brp%faceIndicatorArray(-12)
 ien = grp%brp%faceIndicatorArray(-11)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  wallTangent(1) = grp%brp%faceTangentArray(ip,1)
  wallTangent(2) = grp%brp%faceTangentArray(ip,2)
  wallNormal(1) = -wallTangent(2)
  wallNormal(2) = wallTangent(1) 
  p = grp%p(ip)
  r = grp%coordinates(ip,:)-ivd%momentumPoint
  lift = lift - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowNormal(1)+wallNormal(2)*flowNormal(2))
  drag = drag - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowTangent(1)+wallNormal(2)*flowTangent(2))
  momentum = momentum - 2.0*grp%wallLength(ip)*p*(r(1)*wallNormal(2)-r(2)*wallNormal(1))
  pressureLift = pressureLift - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowNormal(1)+wallNormal(2)*flowNormal(2))
  pressureDrag = pressureDrag - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowTangent(1)+wallNormal(2)*flowTangent(2))

  if(ivd%ReynoldsNumber>0.0) then
  lift = lift - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowNormal(1)+grp%wallStress(ip,2)*flowNormal(2)) 
  drag = drag - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowTangent(1)+grp%wallStress(ip,2)*flowTangent(2)) 
  frictionDrag = frictionDrag - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowTangent(1)+grp%wallStress(ip,2)*flowTangent(2))
  momentum = momentum - 2.0*grp%wallLength(ip)*(r(1)*grp%wallStress(ip,2)-r(2)*grp%wallStress(ip,1))
  end if
 end do

 ! adiabatic viscous faces
 ist = grp%brp%faceIndicatorArray(-10)
 ien = grp%brp%faceIndicatorArray(-9)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  wallTangent(1) = grp%brp%faceTangentArray(ip,1)
  wallTangent(2) = grp%brp%faceTangentArray(ip,2)
  wallNormal(1) = -wallTangent(2)
  wallNormal(2) = wallTangent(1) 
  p = grp%p(ip)
  r = grp%coordinates(ip,:)-ivd%momentumPoint
  lift = lift - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowNormal(1)+wallNormal(2)*flowNormal(2))
  drag = drag - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowTangent(1)+wallNormal(2)*flowTangent(2))
  momentum = momentum - 2.0*grp%wallLength(ip)*p*(r(1)*wallNormal(2)-r(2)*wallNormal(1))
  pressureLift = pressureLift - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowNormal(1)+wallNormal(2)*flowNormal(2))
  pressureDrag = pressureDrag - 2.0*grp%wallLength(ip)*p*(wallNormal(1)*flowTangent(1)+wallNormal(2)*flowTangent(2))

  if(ivd%ReynoldsNumber>0.0) then
  lift = lift - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowNormal(1)+grp%wallStress(ip,2)*flowNormal(2)) 
  drag = drag - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowTangent(1)+grp%wallStress(ip,2)*flowTangent(2)) 
  frictionDrag = frictionDrag - 2.0*grp%wallLength(ip)*(grp%wallStress(ip,1)*flowTangent(1)+grp%wallStress(ip,2)*flowTangent(2))
  momentum = momentum - 2.0*grp%wallLength(ip)*(r(1)*grp%wallStress(ip,2)-r(2)*grp%wallStress(ip,1))
  end if
 end do

! open(15,file="skin_friction.plt",form='formatted',status='unknown')
! do i=1,grp%brp%numberOfBoundaryFaces
!  wallTangent(1) = grp%brp%faceTangentArray(i,1)
!  wallTangent(2) = grp%brp%faceTangentArray(i,2)
!  if(wallTangent(1) >= 0.0) then 
!   write(15,*) i,grp%coordinates(i,1),-2.0*(grp%wallStress(i,1)*wallTangent(1) + grp%wallStress(i,2)*wallTangent(2))
!  else
!   write(15,*) i,grp%coordinates(i,1),2.0*(grp%wallStress(i,1)*wallTangent(1) + grp%wallStress(i,2)*wallTangent(2))
!  end if
! end do 
! close(15)

 end subroutine calculateLiftAndDrag
!-----------------------------------------------------------------------
 subroutine findSeparationFields(grp,ivd)
 IMPLICIT NONE
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 
 logical,pointer :: isNodeInField(:)
 real :: currentCoordinate(2),tripCoordinate(2),currentDistance,maximumDistance 
 integer :: i,j,k,maximumNode,allocateStatus

 ! find number of separation points

  grp%numberOfTripNodes=0
  do i=1,10
   if(ivd%tripNodes(i)>0) then 
    grp%numberOfTripNodes = grp%numberOfTripNodes + 1
   else
    exit
   end if
  end do
  write(*,*) "Number of separation points: ",grp%numberOfTripNodes
  ! now place the ivd%sizeOfSeparationField closest nodes to each separation node in a field 

  allocate(grp%tripNodeFieldIndexes(grp%numberOfTripNodes,0:ivd%sizeOfSeparationField),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%tripNodeFieldDistances(grp%numberOfTripNodes,ivd%sizeOfSeparationField),&
                          stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  allocate(isNodeInField(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  isNodeInField = .false.
  grp%tripNodeFieldDistances = 1.0e10
  grp%tripNodeFieldIndexes = 1
  grp%tripNodeFieldIndexes(:,0) = ivd%sizeOfSeparationField
  do i=1,grp%numberOfTripNodes
   write(*,*) ivd%tripNodes(i)
   tripCoordinate = grp%coordinates(ivd%tripNodes(i),:)
   maximumDistance = 1e10
   maximumNode = 1 
   do j=1,grp%numberOfNodes
    if(.not.isNodeInField(j)) then 
     currentCoordinate = grp%coordinates(j,:) 
     currentDistance = sum((currentCoordinate-tripCoordinate)*(currentCoordinate-tripCoordinate))
     if(currentDistance<maximumDistance) then 
      isNodeInField(grp%tripNodeFieldIndexes(i,maximumNode)) = .false.
      grp%tripNodeFieldIndexes(i,maximumNode) = j
      grp%tripNodeFieldDistances(i,maximumNode) = currentDistance
      isNodeInField(j) = .true.
      ! find new furthermost node
      maximumDistance = grp%tripNodeFieldDistances(i,1)
      maximumNode = 1 
      do k=2,ivd%sizeOfSeparationField
       if(maximumDistance<grp%tripNodeFieldDistances(i,k)) then 
        maximumDistance = grp%tripNodeFieldDistances(i,k)
        maximumNode = k
       end if 
      end do
     end if
    end if
   end do
   grp%tripNodeFieldDistances(i,:) = sqrt(grp%tripNodeFieldDistances(i,:))
  end do 
 
! open(15,file="turbdist.plt",form='formatted',status='unknown')
 do k=1,grp%numberOfNodes
  currentDistance = 0.0
  do i=1,grp%numberOfTripNodes
   do j=1,ivd%sizeOfSeparationField
    if(grp%tripNodeFieldIndexes(i,j)==k) then 
     currentDistance = grp%tripNodeFieldDistances(i,j) 
    end if
   end do
  end do
!  write(15,*) k,100000.0*currentDistance,0.0,0.0,0.0,0.0     
 end do
! close(15)

 grp%sizeOfSeparationFields = ivd%sizeOfSeparationField

 allocate(grp%tripWallLength(grp%numberOfTripNodes),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFields"

 do i=1,grp%numberOfTripNodes
  grp%tripWallLength(i) = grp%wallLength(ivd%tripNodes(i))
 end do

 end subroutine findSeparationFields
!-----------------------------------------------------------------------
 subroutine findSepFieldsForCoarseMesh(finegrp,coarsegrp,ivd)
 IMPLICIT NONE

 type(GridSolverData) :: finegrp,coarsegrp
 type(InputVariablesData) :: ivd

 integer,pointer :: separationPointIndexes(:,:),numberOfCoarseNodesForSF(:),numberOfFineNodesForCoarse(:,:)
 integer :: i,j,k,cnt,allocateStatus
 logical :: alreadyCounted
 type(LinkedInteger),pointer :: currentInteger

 allocate(separationPointIndexes(finegrp%numberOfTripNodes,finegrp%sizeOfSeparationFields),&
                                                              stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFieldsForCoarseMesh"

 do i=1,coarsegrp%numberOfNodes
  currentInteger=>coarsegrp%pod%prolongationArray(i)%first
  do while(associated(currentInteger))
   do j=1,finegrp%numberOfTripNodes
    do k=1,finegrp%sizeOfSeparationFields
     if(currentInteger%int==finegrp%tripNodeFieldIndexes(j,k)) then
      ! found a coarse mesh separation node
      separationPointIndexes(j,k) = i
     end if
    end do
   end do
   currentInteger=>currentInteger%next
  end do
 end do

 ! separationPointIndexes(i,j) now contains the coarse mesh index of
 ! the fine mesh separation field node (i,j)

 ! count how many coarse mesh nodes separation field i contains


 allocate(numberOfCoarseNodesForSF(finegrp%numberOfTripNodes),&
                                                              stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFieldsForCoarseMesh"
 numberOfCoarseNodesForSF = 0

 do i=1,finegrp%numberOfTripNodes
  do j=1,finegrp%sizeOfSeparationFields
   if(separationPointIndexes(i,j)==0) exit
   alreadyCounted = .false.
   do k=1,j-1
    if(separationPointIndexes(i,j)==separationPointIndexes(i,k)) then

     alreadyCounted = .true.
     exit
    end if
   end do
   if(.not.alreadyCounted) then
    numberOfCoarseNodesForSF(i) = numberOfCoarseNodesForSF(i) + 1
   end if
  end do
 end do

 ! now find distances and create coarsegrp%tripNodeFieldIndexes

 coarsegrp%numberOfTripNodes = finegrp%numberOfTripNodes
 coarsegrp%sizeOfSeparationFields = maxval(numberOfCoarseNodesForSF)
 allocate(coarsegrp%tripNodeFieldIndexes(coarsegrp%numberOfTripNodes,coarsegrp%sizeOfSeparationFields),&
                                                              stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFieldsForCoarseMesh"
 allocate(coarsegrp%tripNodeFieldDistances(coarsegrp%numberOfTripNodes,coarsegrp%sizeOfSeparationFields),&
                                                              stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFieldsForCoarseMesh"

 allocate(numberOfFineNodesForCoarse(finegrp%numberOfTripNodes,finegrp%sizeOfSeparationFields),&
                                                              stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFieldsForCoarseMesh"

 numberOfFineNodesForCoarse = 0
 coarsegrp%tripNodeFieldIndexes = 0
 coarsegrp%tripNodeFieldDistances = 0.0

 do i=1,finegrp%numberOfTripNodes
  cnt = 0
  do j=1,finegrp%sizeOfSeparationFields
   if(separationPointIndexes(i,j)==0) exit
   alreadyCounted = .false.
   do k=1,cnt
    if(separationPointIndexes(i,j)==coarsegrp%tripNodeFieldIndexes(i,k)) then
     alreadyCounted = .true.
      coarsegrp%tripNodeFieldDistances(i,k) = coarsegrp%tripNodeFieldDistances(i,k) +&
                                                finegrp%tripNodeFieldDistances(i,j)
      numberOfFineNodesForCoarse(i,k) = numberOfFineNodesForCoarse(i,k) + 1
     exit
    end if
   end do
   if(.not.alreadyCounted) then
    cnt = cnt + 1
    coarsegrp%tripNodeFieldIndexes(i,cnt) = separationPointIndexes(i,j)
    coarsegrp%tripNodeFieldDistances(i,cnt) = finegrp%tripNodeFieldDistances(i,j)
    if(numberOfFineNodesForCoarse(i,cnt)>0) STOP&
       "ERROR: Something's wrong in findSeparationFieldsForCoarseMesh -1"
     numberOfFineNodesForCoarse(i,cnt) = 1
   end if
  end do
 end do

 ! average out distances

 do i=1,coarsegrp%numberOfTripNodes
  do j=1,coarsegrp%sizeOfSeparationFields
   if(coarsegrp%tripNodeFieldIndexes(i,j)==0) exit
   coarsegrp%tripNodeFieldDistances(i,j) = coarsegrp%tripNodeFieldDistances(i,j)/float(numberOfFineNodesForCoarse(i,j))
  end do
 end do

 deallocate(separationPointIndexes,stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Couldn't deallocate in findSeparationFieldsForCoarseMesh"
 deallocate(numberOfFineNodesForCoarse,stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Couldn't deallocate in findSeparationFieldsForCoarseMesh"
 deallocate(numberOfCoarseNodesForSF,stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Couldn't deallocate in findSeparationFieldsForCoarseMesh"


 coarsegrp%rhs = 0.0
 do i=1,finegrp%numberOfTripNodes
 do j=1,coarsegrp%sizeOfSeparationFields
  coarsegrp%rhs(coarsegrp%tripNodeFieldIndexes(i,j),1) = coarsegrp%tripNodeFieldDistances(i,j)
 end do
 end do


! call mapCoarseToFine(coarsegrp%rhs,coarsegrp,finegrp,finegrp%rhs)
! do i=1,finegrp%numberOfNodes
!  write(888,*) i,finegrp%rhs(i,1),0.0,0.0,0.0,0.0,0.0,0.0
! end do
! close(888)

 allocate(coarsegrp%tripWallLength(coarsegrp%numberOfTripNodes),&
                                                              stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in findSeparationFieldsForCoarseMesh"

 if(associated(finegrp%wallLength)) then
  do i=1,coarsegrp%numberOfTripNodes
   coarsegrp%tripWallLength(i) = finegrp%wallLength(ivd%tripNodes(i))
  end do
 else
  do i=1,coarsegrp%numberOfTripNodes
   coarsegrp%tripWallLength(i) = finegrp%tripWallLength(i)
  end do
 end if
 end subroutine findSepFieldsForCoarseMesh
!-----------------------------------------------------------------------
 subroutine writeGridResidual(grp,iterationNumber)
 IMPLICIT NONE
 
 type(GridSolverData) :: grp
 integer :: iterationNumber

 integer :: i
 real :: res(4)

 do i=1,4
  res(i) = sqrt(sum((grp%rhs(:,i)-grp%uprev(:,i))*(grp%rhs(:,i)-grp%uprev(:,i))))
 end do

 write(*,*) iterationNumber,res

 end subroutine writeGridResidual 
!-----------------------------------------------------------------------
 subroutine makeSurfacePlotData(grp)
 IMPLICIT NONE

 type(GridSolverData) :: grp

 integer :: i,j,ind1,ind2,ist,ien,ip,cnt,allocateStatus
 real :: r1(2),r2(2),innerProd

 real,pointer :: directionArray(:)
 integer,pointer :: boundaryFlag(:)

 integer :: searchNode,initialSearch

 ! first count wall nodes

 write(*,*) "Making surface plot data..."
 cnt = 0


 allocate(boundaryFlag(grp%brp%numberOfBoundaryFaces),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in makeSurfacePlotData"

 grp%numberOfPlotNodes = 0
 boundaryFlag = 0 

 ist = grp%brp%faceIndicatorArray(-20)
 ien = grp%brp%faceIndicatorArray(-19)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  if(.not.boundaryFlag(ip)) then
   boundaryFlag(ip) = 1
   grp%numberOfPlotNodes = grp%numberOfPlotNodes + 1
  end if
 end do

 ist = grp%brp%faceIndicatorArray(-12)
 ien = grp%brp%faceIndicatorArray(-11)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  if(.not.boundaryFlag(ip)) then 
   boundaryFlag(ip) = 1 
   grp%numberOfPlotNodes = grp%numberOfPlotNodes + 1
  end if
 end do

 ist = grp%brp%faceIndicatorArray(-10)
 ien = grp%brp%faceIndicatorArray(-9)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  if(.not.boundaryFlag(ip)) then
   boundaryFlag(ip) = 1 
   grp%numberOfPlotNodes = grp%numberOfPlotNodes + 1
  end if
 end do

 allocate(grp%plotNodePairs(grp%numberOfPlotNodes,2),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in makeSurfacePlotData"

 allocate(grp%plotNodeDistance(grp%numberOfPlotNodes),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in makeSurfacePlotData"

 do i=1,grp%brp%numberOfBoundaryFaces
  if(boundaryFlag(i)==1) then 
   boundaryFlag(i)=2
   cnt = cnt + 1
   grp%plotNodePairs(cnt,1) = i 
   searchNode = i
   do while(searchNode>0) 
    initialSearch = searchNode
    do j=1,grp%brp%numberOfBoundaryFaces
     if(grp%brp%faceIndexArray(j,1)==searchNode.and.boundaryFlag(grp%brp%faceIndexArray(j,2))==1) then 
      boundaryFlag(grp%brp%faceIndexArray(j,2))=2
      cnt = cnt + 1
      grp%plotNodePairs(cnt,1) = grp%brp%faceIndexArray(j,2)
      searchNode = grp%brp%faceIndexArray(j,2) 
     else if(grp%brp%faceIndexArray(j,2)==searchNode.and.boundaryFlag(grp%brp%faceIndexArray(j,1))==1) then 
      boundaryFlag(grp%brp%faceIndexArray(j,1))=2
      cnt = cnt + 1
      grp%plotNodePairs(cnt,1) = grp%brp%faceIndexArray(j,1)
      searchNode = grp%brp%faceIndexArray(j,1)
     end if
    end do 
    if(initialSearch==searchNode) then 
     ! couldn't find next node - terminate
     searchNode = 0
    end if
   end do
  end if 
 end do

 if(cnt.ne.grp%numberOfPlotNodes) then 
  write(*,*) cnt,grp%numberOfPlotNodes
  STOP "ERROR: Something's gone wrong in makeSurfacePlotData - 1"
 end if

 deallocate(boundaryFlag,stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: makeSurfacePlotData couldn't deallocate"

 allocate(directionArray(grp%numberOfPlotNodes),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in makeSurfacePlotData"

! find nodes normal to plot nodes

 directionArray = 1000.0

 do i=1,grp%numberOfSides
  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)
!  if(ind1.le.grp%brp%numberOfBoundaryFaces.and.ind2>grp%brp%numberOfBoundaryFaces) then  
  if(ind1.le.grp%brp%numberOfBoundaryFaces) then  
   ! see if ind1 is plot node
   do j=1,grp%numberOfPlotNodes
    if(ind1==grp%plotNodePairs(j,1)) then 
     r1 = grp%coordinates(ind2,:) - grp%coordinates(ind1,:)
     r1 = r1/sqrt(sum(r1*r1))
     r2(1) = grp%brp%faceTangentArray(ind1,2)
     r2(2) = -grp%brp%faceTangentArray(ind1,1)
     innerProd = sum(r1*r2)
     if(innerProd<directionArray(j)) then 
      directionArray(j) = innerProd
      grp%plotNodePairs(j,2) = ind2
     end if
    end if
   end do
  end if
 end do

 deallocate(directionArray,stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: makeSurfacePlotData couldn't deallocate"

 do i=1,grp%numberOfPlotNodes
  r1 = grp%coordinates(grp%plotNodePairs(i,2),:) - grp%coordinates(grp%plotNodePairs(i,1),:)
  grp%plotNodeDistance(i) = sqrt(sum(r1*r1))
 end do


 write(*,*) "number of plot nodes: ",grp%numberOfPlotNodes
! do i=1,grp%numberOfPlotNodes
!  write(*,*) "b: ",grp%plotNodePairs(i,1:2),grp%plotNodeDistance(i)
! end do

 end subroutine makeSurfacePlotData
!-----------------------------------------------------------------------
! subroutine writeGridResults(grp,OUTFILE,OUTFILEt,ivd)
 subroutine writeGridResults(grp,OUTFILE,ivd)
 
 IMPLICIT NONE

 type(GridSolverData) :: grp
 integer :: OUTFILE
 type(InputVariablesData) :: ivd

! integer :: i,j,is,in
! real :: x0,y0,xn,yn,um,ds

integer :: i
 
 do i=1,grp%numberOfNodes
  write(OUTFILE,'(I5,5E14.4)') i,grp%u(i,1),grp%u(i,2:4)/grp%u(i,1),grp%p(i)
 end do

!  do i=1,grp%numberOfNodes
!   is = i
!   x0 = grp%coordinates(i,1)
!   y0 = grp%coordinates(i,2)
!   in = is
!   if(grp%AfterAndBefore(i,1).ne.0.and.grp%AfterAndBefore(i,2).eq.0)then
! 30  xn = grp%coordinates(in,1)
!     yn = grp%coordinates(in,2)
!     ds = sqrt((x0-xn)**2+(y0-yn)**2)
!     um = sqrt((grp%u(in,2)/grp%u(in,1))**2+(grp%u(in,3)/grp%u(in,1))**2)
!     in = grp%AfterAndBefore(in,1)
!     if(in.ne.0) goto 30
!   end if
!  end do
!  close(888)
! open(999,file='turbulence.res',form='formatted',status='unknown') 
! if(ivd%turbulenceModel>0) then 
!  do i=1,grp%numberOfNodes
!   write(OUTFILE,'(I5,5E25.15)') i,grp%u(i,1),grp%u(i,2:4)/grp%u(i,1),grp%u(i,5)
!   write(999,'(I5,5E14.4)') i,grp%u(i,5),grp%u(i,5)/grp%u(i,1),grp%wallDistance(i),0,0
!  end do
! else
!  do i=1,grp%numberOfNodes
!   write(OUTFILE,'(I5,5E14.4)') i,grp%u(i,1),grp%u(i,2:4)/grp%u(i,1),0.0
!  end do
! end if
! close(999)

 end subroutine writeGridResults
!------------------------------------------------------------------------
subroutine writeGridResults2(grp,OUTFILE)
 
 IMPLICIT NONE

 type(GridSolverData) :: grp
 integer :: OUTFILE

 integer :: ist, ien, ib, ip

 ist = grp%brp%faceIndicatorArray(-6)
 ien = grp%brp%faceIndicatorArray(-5)

 do ib =ist,ien
  ip = grp%brp%faceIndicatorArray(ib)
  write(OUTFILE,'(I5,5E25.15)') ip,grp%u(ip,1),grp%u(ip,2:4)/grp%u(ip,1),grp%p(ip)
 end do

 end subroutine writeGridResults2
!------------------------------------------------------------------------
 subroutine writeSurface(grp,OUTFILE,ivd)
 
 IMPLICIT NONE

 type(GridSolverData) :: grp
 integer :: OUTFILE
 type(InputVariablesData) :: ivd

 integer :: i,j,ist,ien,ip
 real :: wallTangent(2),st,cp,p,pinf


 pinf = ivd%inflowField(5,1)

 ist = grp%brp%faceIndicatorArray(-20)
 ien = grp%brp%faceIndicatorArray(-19)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  p = (ivd%gamma-1.0)*(grp%u(ip,4)-0.5*sum(grp%u(ip,2:3)*grp%u(ip,2:3))/grp%u(ip,1))
  cp = 2.0*(pinf-p)
  if(wallTangent(1) >= 0.0) then
   st = 0.0
   write(OUTFILE,*) grp%coordinates(ip,1),st,cp
  else
   st = 0.0
   write(OUTFILE,*) grp%coordinates(ip,1),st,cp
  end if
 end do

 ist = grp%brp%faceIndicatorArray(-12)
 ien = grp%brp%faceIndicatorArray(-11)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  wallTangent(1) = grp%brp%faceTangentArray(ip,1)
  wallTangent(2) = grp%brp%faceTangentArray(ip,2)
  p = (ivd%gamma-1.0)*(grp%u(ip,4)-0.5*sum(grp%u(ip,2:3)*grp%u(ip,2:3))/grp%u(ip,1))
  cp = 2.0*(pinf-p) 
  if(wallTangent(1) >= 0.0) then
   st = -2.0*(grp%wallStress(ip,1)*wallTangent(1) + grp%wallStress(ip,2)*wallTangent(2))
   write(OUTFILE,*) grp%coordinates(ip,1),st,cp
  else
   st = 2.0*(grp%wallStress(ip,1)*wallTangent(1) + grp%wallStress(ip,2)*wallTangent(2))
   write(OUTFILE,*) grp%coordinates(ip,1),st,cp
  end if
 end do

 ist = grp%brp%faceIndicatorArray(-10)
 ien = grp%brp%faceIndicatorArray(-9)
 do i=ist,ien
  ip = grp%brp%faceIndicatorArray(i)
  wallTangent(1) = grp%brp%faceTangentArray(ip,1)
  wallTangent(2) = grp%brp%faceTangentArray(ip,2)
  p = (ivd%gamma-1.0)*(grp%u(ip,4)-0.5*sum(grp%u(ip,2:3)*grp%u(ip,2:3))/grp%u(ip,1))
  cp = 2.0*(pinf-p) 
  if(wallTangent(1) >= 0.0) then
   st = -2.0*(grp%wallStress(ip,1)*wallTangent(1) + grp%wallStress(ip,2)*wallTangent(2))
   write(OUTFILE,*) grp%coordinates(ip,1),st,cp
  else
   st = 2.0*(grp%wallStress(ip,1)*wallTangent(1) + grp%wallStress(ip,2)*wallTangent(2))
   write(OUTFILE,*) grp%coordinates(ip,1),st,cp
  end if
 end do
 end subroutine writeSurface
!------------------------------------------------------------------------
 subroutine writeSurface2(grp,OUTFILE,ivd)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 integer :: OUTFILE
 type(InputVariablesData) :: ivd

 integer :: i,j,ist,ien,ind1,ind2,ind3,ind4
 real :: wallTangent(2),st,cp,p,pinf,u1(2),u2(2),u3(2),u4(2),uu1,uu2,uu3,uu4,stn
 real :: compScaling


 pinf = ivd%inflowField(5,1)
 do i=1,grp%numberOfPlotNodes
  ind1 = grp%plotNodePairs(i,1)
  ind2 = grp%AfterAndBefore(ind1,1)
  if(ind2.eq.0) then
    ind2 = grp%plotNodePairs(i,2)
  else
    ind3 = grp%AfterAndBefore(ind2,1)
    ind4 = grp%AfterAndBefore(ind3,1)
  end if
  p = (ivd%gamma-1.0)*(grp%u(ind1,4)-0.5*sum(grp%u(ind1,2:3)*grp%u(ind1,2:3))/grp%u(ind1,1))
  cp = 2.0*(pinf-p)
  u1 = grp%u(ind1,2:3)/grp%u(ind1,1)
  u2 = grp%u(ind2,2:3)/grp%u(ind2,1)
  uu1 = sqrt(sum(u1*u1))
  uu2 = sqrt(sum(u2*u2))
  st = (uu2-uu1)/grp%plotNodeDistance(i)
  st = 2.0*grp%laminarViscosity(ind1)*st
  stn= st
  if(ind3.ne.0) then
    u3 = grp%u(ind3,2:3)/grp%u(ind3,1)
    u4 = grp%u(ind4,2:3)/grp%u(ind4,1)
    uu3 = sqrt(sum(u3*u3))
    uu4 = sqrt(sum(u4*u4))
    stn = (-11.0*uu1+18.0*uu2-9.0*uu3+2.0*uu4)/grp%plotNodeDistance(i)/6.0
    stn = 2.0*grp%laminarViscosity(ind1)*stn
  end if

!  compScaling = 0.5*(grp%u(ind1,1)+grp%u(ind2,1))
!  st = st*compScaling
  
  if(ind1==120) write(*,*) "A: ",grp%laminarViscosity(ind1),st,grp%u(ind1,1:4)
  write(OUTFILE,'(13E17.8)') grp%coordinates(ind1,1:2),st,stn,cp,grp%u(ind1,1:4),uu1,uu2,uu3,uu4
 end do

 end subroutine writeSurface2
!------------------------------------------------------------------------
end module GridSolver
