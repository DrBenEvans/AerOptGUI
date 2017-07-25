!**********************************************
!*aggl2d_def.f90                              *
!*                                            *
!*Module definition file for program          *
!*                                            * 
!*   aggl2d v 2.0                             * 
!*                                            *
!*                                            *
!*Made by                                     *
!*Kaare A Sorensen,                           *
!*20.08.98-                                   *
!**********************************************
!*Description:                                *
!* This file includes procedures for file     *
!* communication, construction of side based  *
!* datastructure from elements, coloring for  *
!* for vectorization and construction of      *
!* coarser grids from the given original one  *
!* by agglomeration - for use in multigrid    *
!**********************************************


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module RectangularElements
 integer,parameter,private:: NSD = 2

 type :: SREData ! data structure for a single element    
  integer :: pointIndexes(4) ! indexes to points in element 
!  real :: area 
 end type SREData

 type RectangularElementData ! data structure for element array
   type(SREData), pointer :: elements(:) 
   integer ::  numberOfElements 
 end type RectangularElementData

 contains

!-------------------------------------------------------------------------
 subroutine constructRectElArray(n,rep)
! Sets up a REData type for use
  
  IMPLICIT NONE 

  integer :: n
  type(RectangularElementData) :: rep

  integer :: allocateStatus

! allocate memory for element array

  if(n>0) then  
   allocate(rep%elements(n),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create element array"
  else
   nullify(rep%elements)
  end if
  rep%numberOfElements=n  
 
 end subroutine constructRectElArray
!-------------------------------------------------------------------------
 subroutine readRectElData(INFILE,rep)
! Reads element data 
 IMPLICIT NONE

 integer :: INFILE
 type(RectangularElementData) :: rep

 integer :: dummy,i,ip
 integer :: nel 
 nel = rep%numberOfElements

  do ip=1,nel
   read(INFILE,*)dummy, (rep%elements(ip)%pointIndexes(i),i=1,4)   
  end do 

 end subroutine readRectElData
!-------------------------------------------------------------------------
  function getRectElement(i,rep) result(c)
  
  integer :: i
  type(RectangularElementData) :: rep 
  type(SREData),pointer :: c

  c => rep%elements(i)

 end function getRectElement
!-------------------------------------------------------------------------
end module RectangularElements


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


module TriangularElements
 integer,parameter,private:: NSD = 2

 type :: STEData ! data structure for a single element
  integer :: pointIndexes(3) ! indexes to points in element
!  real :: area
 end type STEData

 type TriangularElementData ! data structure for element array
   type(STEData), pointer :: elements(:)
   integer ::  numberOfElements
 end type TriangularElementData

 contains

! GENERAL REMARK:  To increase modularity of code the same standard set
!                  of subroutines and functions should be implemented
!                  in other type of elements, like f.ex. in three dimensions
!                  As a general rule, all public routines and variables
!                  need to be defined in a new element module

!-------------------------------------------------------------------------
 subroutine constructTriElArray(n,tep)
! Sets up a TEData type for use

  IMPLICIT NONE

  integer :: n
  type(TriangularElementData) :: tep

  integer :: allocateStatus

! allocate memory for element array

  if(n>0) then
   allocate(tep%elements(n),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create element array"
  else
   nullify(tep%elements)
  end if

  tep%numberOfElements=n

 end subroutine constructTriElArray
!-------------------------------------------------------------------------
 subroutine readTriElData(INFILE,tep)
! Reads element data
 IMPLICIT NONE

 integer :: INFILE
 type(TriangularElementData) :: tep

 integer :: dummy,i,ip
 integer :: nel
 nel = tep%numberOfElements

  do ip=1,nel
   read(INFILE,*)dummy, (tep%elements(ip)%pointIndexes(i),i=1,3)
  end do

 end subroutine readTriElData
!-------------------------------------------------------------------------
  function getTriElement(i,tep) result(c)

  integer :: i
  type(TriangularElementData) :: tep
  type(STEData),pointer :: c

  c => tep%elements(i)

 end function getTriElement
!-------------------------------------------------------------------------
end module TriangularElements

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


module CoordinateRegister
! keeps track of coordinate points
! since the coarser grids use the same points, they will only contain
! pointers to this register

 type CoordinateRegisterData ! coordinate register data
  real, pointer :: points(:,:)
  integer :: numberOfPoints,NSD
 end type CoordinateRegisterData
 
 integer, parameter, private :: NSD = 2

 contains

!-------------------------------------------------------------------------
 subroutine constructCoordinateRegister(np,nsd,crp)
! sets up register

 IMPLICIT NONE

 integer :: np,nsd
 type(CoordinateRegisterData) :: crp

 integer allocateStatus

 crp%NSD = nsd
 crp%numberOfPoints = np

 allocate(crp%points(np,nsd),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create point array"

 end subroutine constructCoordinateRegister
!-------------------------------------------------------------------------
 subroutine readPointData(INFILE,crp)
! reads data from file

 IMPLICIT NONE

 integer :: INFILE
 type(CoordinateRegisterData) :: crp

 integer ip,i,dummy
 
 integer NSD,NP

 NSD = crp%NSD
 NP = crp%numberOfPoints
 crp%points = 0.0

 do ip=1,NP
  read(INFILE,*)dummy, (crp%points(ip,i),i=1,NSD)
 end do
!crp%points(:,1) = maxval(crp%points(:,1))-crp%points(:,1)


 end subroutine readPointData
!------------------------------------------------------------------------- 
 subroutine writeCoordinateData(OUTFILE,crp)
! writes coordinate data for computation file

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 integer :: OUTFILE

 integer :: i,j

 do i=1,crp%numberOfPoints
  write(OUTFILE) crp%points(i,1:NSD)
 end do

 end subroutine writeCoordinateData
!-------------------------------------------------------------------------
 function getCoor(i,crp) result (c)
! returns coordinates related to a point index
 IMPLICIT NONE

 integer :: i
 type(CoordinateRegisterData) :: crp

 real :: c(NSD)

 c = crp%points(i,:)

 end function getCoor
!-------------------------------------------------------------------------
 subroutine setCoor(i,c,crp) 
! returns coordinates related to a point index
 IMPLICIT NONE

 integer :: i
 type(CoordinateRegisterData) :: crp
 real :: c(NSD)

 crp%points(i,:) = c
 end subroutine setCoor
!-------------------------------------------------------------------------
end module CoordinateRegister

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module SideModule

 type SideData
  integer :: sideIndexes(2) ! mapping between side indexes and point indexes
  real, pointer :: sideCoefficients(:)
!  real :: sideCoefficients(2)

  integer,pointer :: neighbouringElements(:) ! pointer to neighbouring elements,
                                             ! used in visualization
  integer,pointer :: oldIndexes(:)
  real,pointer :: sideLength(:)
  real,pointer :: invertedSideLength(:)
!  real :: sideCoefficients(2)
 end type SideData

 contains

!--------------------------------------------------------------------------
  subroutine sideConstruct(n,sd)
! initializes side

  IMPLICIT NONE
 
  integer :: n
  type(SideData),pointer :: sd

  integer :: i,allocateStatus

  allocate(sd%sideCoefficients(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create side"

  do i=1,n
   sd%sideCoefficients(i) = 0.
  end do 

  end subroutine sideConstruct
!-------------------------------------------------------------------------

end module SideModule

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module VisualizationModule
! does the visualization bit
 type VisualizationSide
  integer :: sideIndexes(2)
  integer :: neighbouringElements(2)
  integer :: controlVolumes(2)
 end type VisualizationSide

 type VisualizationSidePointer
  type(VisualizationSide),pointer :: vs
  type(VisualizationSidePointer),pointer :: next
 end type VisualizationSidePointer

 type VisualizationSidePointers
  type(VisualizationSidePointer),pointer :: first,last 
 end type VisualizationSidePointers

 type VisualizationData
  type(VisualizationSidePointers),pointer :: vsp(:)
  integer :: numberOfVisualizationSides 
 end type VisualizationData

end module VisualizationModule

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module SideRegister
! register for sides
 use SideModule

 type SideRegisterData
  type(SideData), pointer :: sideData(:)
 end type SideRegisterData

end module SideRegister

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module SuperSide
 type SuperSideData
  real :: SSCoefficients(2) ! holds sum of side coefficients 
  integer :: numberOfSides
  integer :: SSindexes(2) ! links to neighbouring control volumes
  ! linked list of sides
 end type SuperSideData

end module SuperSide

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module BoundaryFace2D
! module for a side on boundary

 type BoundaryFace2DData
  integer :: faceIndexes(2)
  integer,pointer :: oldFaceIndexes(:)
  integer :: indicator
  real,pointer :: faceCoefficients(:)
 end type BoundaryFace2DData

end module BoundaryFace2D

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module BoundaryRegister
! register for boundaries
 use BoundaryFace2D
 use CoordinateRegister
 
 integer,parameter,private:: NSD = 2
 
 type LinkedIntegerPointer
  integer :: index
  type(LinkedIntegerPointer), pointer :: next
 end type  

 type LinkedFace
  type(BoundaryFace2DData),pointer :: bd
  type(LinkedFace),pointer :: next
 end type LinkedFace

 type LinkedInteger
  integer :: col   
  type(LinkedInteger),pointer :: next
 end type LinkedInteger

 type LinkedIntegerArray
  integer,pointer :: arr(:)   
  type(LinkedIntegerArray),pointer :: next
 end type LinkedIntegerArray

 type LinkedIntegerP
  type(LinkedInteger),pointer :: first
 end type LinkedIntegerP

 type LinkedIntegerArrayP
  type(LinkedIntegerArray),pointer :: first
 end type LinkedIntegerArrayP

 type BoundaryColorData
  type(LinkedFace),pointer :: first,last
  type(BoundaryColorData),pointer :: next
  integer :: numberOfSides,colorNumber
 end type BoundaryColorData 

 type BoundaryColorList
  type(BoundaryColorData),pointer :: first,last
  integer :: numberOfColors
 end type BoundaryColorList 
 
 type BoundaryRegisterData
  type(BoundaryFace2DData), pointer :: faces(:)
  integer :: numberOfBoundaryFaces,numberOfBoundaryNodes 
  integer :: numberOfColors
  real,pointer :: faceTangents(:,:)
  type(LinkedIntegerPointer),pointer :: firstTE,lastTE  ! points to first and last 
                                                        ! instance in trailing edge chain
  integer,pointer :: boundaryNodeIndicators(:)
  type(LinkedIntegerP),pointer :: colorsOccupied(:)
  integer,pointer :: internalOutflowRegister(:)  ! for internal outflow chains
  integer,pointer :: IORegisterLengths(:)
  integer :: numberOfIORegisters,numberOfIONodes
  integer :: numberOfTrailingEdges
  integer,pointer :: startOfColors(:)  ! points to beginning index of color in faces array
  
  !added with engine inlet
  integer :: numberOfEngineInletSides
  integer :: numberOfEngineInletSideLimit
  integer,pointer :: engineInletSideIndexes(:,:)
  real,pointer :: engineInletSideCoefficients(:,:)
  
 end type BoundaryRegisterData

 contains

!------------------------------------------------------------------------
 subroutine constructBoundaryRegister(n,brp)
! intializes module

  IMPLICIT NONE
 
  integer :: n ! number of boundary faces
  type(BoundaryRegisterData) :: brp

  integer :: i,allocateStatus

  allocate(brp%faces(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

  allocate(brp%faceTangents(n,2),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

  allocate(brp%boundaryNodeIndicators(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

  do i=1,n
!   allocate(brp%faces(i)%indicator,stat=allocateStatus)
!   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
   allocate(brp%faces(i)%faceCoefficients(2),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
   brp%faces(i)%faceCoefficients = 0.0 
   nullify(brp%faces(i)%oldFaceIndexes)
  end do

  brp%faceTangents = 0.0


! the number of boundary points equals the number of boundary 
! faces if only closed boundaries are used
! the boundary nodes are always the first nodes in the file, that is 
! nodes (1,numberOfBoundaryFaces) are boundary nodes

  brp%numberOfBoundaryFaces = n
  brp%numberOfBoundaryNodes = n
 end subroutine constructBoundaryRegister
!-------------------------------------------------------------------------
 subroutine setUpNewBoundaryData(brp,newbrp)
 ! sets up boundary data for agglomerated grid
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp,newbrp

 integer :: i,allocateStatus

 allocate(newbrp%faces(brp%numberOfBoundaryFaces),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(newbrp%boundaryNodeIndicators(brp%numberOfBoundaryFaces),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 do i = 1,brp%numberOfBoundaryFaces
  newbrp%faces(i)%faceCoefficients => brp%faces(i)%faceCoefficients
!  newbrp%faces(i)%indicator => brp%faces(i)%indicator
  newbrp%boundaryNodeIndicators(i) = brp%boundaryNodeIndicators(i)
  if(associated(brp%faces(i)%oldFaceIndexes)) then 
   newbrp%faces(i)%oldFaceIndexes => brp%faces(i)%oldFaceIndexes
  else 
  ! newbrp belongs to second grid
   allocate(newbrp%faces(i)%oldFaceIndexes(2),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
   newbrp%faces(i)%oldFaceIndexes = brp%faces(i)%faceIndexes
  end if
  ! the face indexes cannot be pointed over since they are different for agglomerated grids
 end do

 newbrp%faceTangents => brp%faceTangents
 newbrp%firstTE => brp%firstTE
 newbrp%lastTE => brp%lastTE

 nullify(newbrp%colorsOccupied)

 newbrp%numberOfBoundaryFaces = brp%numberOfBoundaryFaces
 newbrp%numberOfTrailingEdges = brp%numberOfTrailingEdges

 newbrp%numberOfIORegisters = 0
 newbrp%numberOfIONodes = 0
 nullify(newbrp%IORegisterLengths)
 nullify(newbrp%internalOutflowRegister)

 end subroutine setUpNewBoundaryData
!-------------------------------------------------------------------------
 subroutine readBoundData(INFILE,brp)
! Reads boundary data form file 
 IMPLICIT NONE

 integer :: INFILE
 type(BoundaryRegisterData) :: brp

 integer :: i,ip

 integer :: nbd

 nbd = brp%numberOfBoundaryFaces
  do ip=1,nbd
   read(INFILE,*)brp%faces(ip)%faceIndexes(1), brp%faces(ip)%faceIndexes(2),brp%faces(ip)%indicator
!   if(brp%faces(ip)%indicator==7) write(*,*) "eet: ",ip,brp%faces(ip)%faceIndexes
  end do 

 
 end subroutine readBoundData
!-------------------------------------------------------------------------
 subroutine writeBoundaryData(OUTFILE,brp)
! writes coordinate data for computation file
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp 
 integer :: OUTFILE

 integer :: i,j,allocateStatus,ind

 integer,pointer :: boundaryIndicatorArray(:)
 integer :: startIndex(8),currentIndex
 type(LinkedIntegerPointer),pointer :: currentInteger

 allocate(boundaryIndicatorArray(-20:brp%numberOfBoundaryFaces),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeBoundaryData out of memory"


 boundaryIndicatorArray = 0
 ! first run through all faces to count the number of each BC
 do i=1,brp%numberOfBoundaryFaces
!  ind = brp%faces(i)%indicator
!  ind = brp%boundaryNodeIndicators(i)
! ind = brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1))
  currentInteger=>brp%firstTE
  if(associated(brp%faces(i)%oldFaceIndexes)) then
   currentIndex = brp%faces(i)%oldFaceIndexes(1) 
  ind = brp%boundaryNodeIndicators(brp%faces(i)%oldfaceIndexes(1))
  if(ind>8) STOP "ERROR: Boundary indicator not defined" 
  else
   currentIndex = brp%faces(i)%faceIndexes(1)
  ind = brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1))
  if(ind>8) STOP "ERROR: Boundary indicator not defined" 
  end if
  do while(associated(currentInteger)) ! check if trailing edge
!   if(currentInteger%index==i) then 
   if(currentInteger%index==currentIndex) then 
    goto 76 
   end if
   currentInteger=>currentInteger%next
  end do
  boundaryIndicatorArray(-21+2*ind) = boundaryIndicatorArray(-21+2*ind) + 1
  boundaryIndicatorArray(0) = boundaryIndicatorArray(0)+1 
76 continue
 end do 
 ! set negative indexes

 ind = 0 
 do i=1,8
  boundaryIndicatorArray(-22+2*i) = ind+1 ! start index
  ind = ind + boundaryIndicatorArray(-21+2*i) 
  boundaryIndicatorArray(-21+2*i) = ind ! end index
 end do

 ! check the number of types of farfield BC's

 if(boundaryIndicatorArray(-14)>=boundaryIndicatorArray(-13)) then 
  boundaryIndicatorArray(-1) = 1
 else
  boundaryIndicatorArray(-1) = 2
 end if

 ! set up start indices
 
 do i=1,8
  startIndex(i) = boundaryIndicatorArray(-22+2*i)
 end do

 ! now insert the node number belonging to each type of BC
 do i=1,brp%numberOfBoundaryFaces 
!  ind = brp%faces(i)%indicator
!  ind = brp%boundaryNodeIndicators(i)
! ind = brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1))
  currentInteger=>brp%firstTE
  if(associated(brp%faces(i)%oldFaceIndexes)) then
   currentIndex = brp%faces(i)%oldFaceIndexes(1) 
  ind = brp%boundaryNodeIndicators(brp%faces(i)%oldfaceIndexes(1))
  else
   currentIndex = brp%faces(i)%faceIndexes(1)
  ind = brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1))
  end if
  do while(associated(currentInteger))
   if(currentInteger%index==currentIndex) then 
    goto 79 ! trailing edge has no entry in matrix
   end if
   currentInteger=>currentInteger%next
  end do
  boundaryIndicatorArray(startIndex(ind)) = currentIndex 
  startIndex(ind) = startIndex(ind)+1 
79 continue
 end do

 ! finally set trailing edges at the end 

 currentInteger=>brp%firstTE
 do while(associated(currentInteger)) 
  boundaryindicatorArray(startIndex(7)) = currentInteger%index
  startIndex(7) = startIndex(7)+1
  currentInteger=>currentInteger%next
 end do 

 ! write face indexes to file

 do i=1,brp%numberOfBoundaryFaces 
  write(OUTFILE) brp%faces(i)%faceIndexes,brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1))
 end do

 ! write boundary indicator array to file
 do i=-20,brp%numberOfBoundaryFaces
  write(OUTFILE) boundaryIndicatorArray(i)
!  write(291,*) boundaryIndicatorArray(i)
 end do

 ! write face coefficients to file

 do i=1,brp%numberOfBoundaryFaces 
  write(OUTFILE) brp%faces(i)%faceCoefficients
 end do

 ! write face tangents to file

 do i=1,brp%numberOfBoundaryFaces 
  write(OUTFILE) brp%faceTangents(i,:) 
 end do

 deallocate(boundaryIndicatorArray,stat=allocateStatus)
 if(allocateStatus>0) STOP "ERROR: could not deallocate in writeBoundaryData"

 write(OUTFILE) brp%numberOfIORegisters

! do i=1,brp%numberOfIORegisters
!  write(OUTFILE) brp%IORegisterLengths(i)
! end do
 write(OUTFILE) brp%numberOfIONodes
 do i=1,brp%numberOfIONodes
  write(OUTFILE) brp%internalOutflowRegister(i)
 end do

 write(OUTFILE) brp%numberOfEngineInletSides
 write(OUTFILE) brp%engineInletSideIndexes(1:brp%numberOfEngineInletSides,1:2)
 write(OUTFILE) brp%engineInletSideCoefficients(1:brp%numberOfEngineInletSides,1:2)
  
 end subroutine writeBoundaryData
!-------------------------------------------------------------------------
 subroutine makeBoundaryTangents(brp,crp)
 use CoordinateRegister
 use BoundaryFace2D
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp
 type(CoordinateRegisterData) :: crp

 real :: x1(2),x2(2),nn(2) 
 integer :: i,j,nob
 type(BoundaryFace2DData) :: bfp
 real :: norm

 nob = brp%numberOfBoundaryFaces

! can now assume that coordinates are set up correctly from file
 brp%faceTangents = 0.0
 do i=1,nob
! get coordinates of side points

  bfp = brp%faces(i)

  x1 = getCoor(bfp%faceIndexes(1),crp)
  x2 = getCoor(bfp%faceIndexes(2),crp)

  nn = x2-x1
  brp%faceTangents(bfp%faceIndexes(1),:) = brp%faceTangents(bfp%faceIndexes(1),:) + nn
  brp%faceTangents(bfp%faceIndexes(2),:) = brp%faceTangents(bfp%faceIndexes(2),:) + nn
 end do
 
 ! normalize tangents

 do i=1,nob
  nn = brp%faceTangents(i,:)
  norm = sqrt(sum(nn*nn)) 
  brp%faceTangents(i,:) = brp%faceTangents(i,:)/norm
 end do

 end subroutine makeBoundaryTangents
!-------------------------------------------------------------------------
 subroutine searchForTrailingEdges(brp,crp)
! uses a double loop to find nodes which belong to elements 
! creating a sharp angle between them
! The indexes of trailing edge nodes are put in a linked list
 
 use CoordinateRegister
 use BoundaryFace2D

 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp
 type(CoordinateRegisterData) :: crp

 integer :: i,j,nob,allocateStatus
 real :: x1(2),x2(2),x3(2),v1(2),v2(2)

 nob = brp%numberOfBoundaryFaces
 nullify(brp%firstTE)
 nullify(brp%lastTE)
 brp%numberOfTrailingEdges = 0

 do i=1,nob
  do j=1,nob
   if(brp%faces(i)%faceIndexes(1) == brp%faces(j)%faceIndexes(2)) then
    x1 = getCoor(brp%faces(i)%faceIndexes(2),crp)
    x2 = getCoor(brp%faces(i)%faceIndexes(1),crp)
    x3 = getCoor(brp%faces(j)%faceIndexes(1),crp)              
    v1 = x1 - x2
    v2 = x3 - x2
    v1 = v1/sqrt(dot_product(v1,v1))
    v2 = v2/sqrt(dot_product(v2,v2))
! the angle between v1 and v2 is now arccos(dot_product(v1,v2))
! we do not allow this angle to be smaller than, say, arccos(0.9)
! (which is about 25.8 degrees)
    if(dot_product(v1,v2) > 0.9) then 
     ! have found a trailing edge
     if(associated(brp%lastTE))  then
     ! have allready found TE's      
      allocate(brp%lastTE%next,stat=allocateStatus)
      if (allocateStatus /= 0) STOP&
                             "ERROR: Not enough memory to create boundary array"
      brp%lastTE => brp%lastTE%next
!      brp%lastTE%index = i
      brp%lastTE%index = brp%faces(i)%faceIndexes(1)
      nullify(brp%lastTE%next)
     else
      allocate(brp%firstTE,stat=allocateStatus)
      if (allocateStatus /= 0) STOP&
                             "ERROR: Not enough memory to create boundary array" 
      brp%lastTE => brp%firstTE
!      brp%lastTE%index = i 
      brp%lastTE%index = brp%faces(i)%faceIndexes(1)
      nullify(brp%lastTE%next)  
     end if
     brp%numberOfTrailingEdges = brp%numberOfTrailingEdges + 1
    end if   
   end if
  end do
 end do  

 write(*,*) "Number of trailing edges found: ",brp%numberOfTrailingEdges
 end subroutine searchForTrailingEdges
!--------------------------------------------------------------------------
 subroutine makeInternalOutflowRegister(brp,rep,tep)
 ! for internal outflow (for example for flow in a channel), a register containing
 ! node numbers of connected nodes is found. This can be used by the solver 
 ! for interpolation to boundary 
 
 use CoordinateRegister
 use RectangularElements
 use TriangularElements 

 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp
 type(RectangularElementData) :: rep
 type(TriangularElementData) :: tep

 integer :: i,j,k,l,m,ind1,ind2,ind3,numberOfConnectedElements,numberOfStartNodes
 integer :: nodeNumber,searchNodeNumber,allocateStatus,i1,i2,buff,numberOfIOChains
 integer :: numberOfConnectedNodes,elementNumber,chainLengthAddition,counter1,counter2
 logical :: elementIsAsItShould,hasChanged
 type(LinkedIntegerP),pointer :: elementChain
 type(LinkedInteger),pointer :: lastElement,currentElement,previousElement
 type(LinkedIntegerArrayP),pointer :: startNodes
 type(LinkedIntegerArray),pointer :: lastNode,currentNode,previousNode

 write(*,*) "Making internal outflow register..."


 nullify(elementChain)

 numberOfConnectedElements = 0 
 numberOfStartNodes = 0
 chainLengthAddition = 0
 do i=1,rep%numberOfElements
  do j=1,4
   nodeNumber = rep%elements(i)%pointIndexes(j) 
   if(nodeNumber.le.brp%numberOfBoundaryFaces) then ! node is on boundary 
    if(brp%boundaryNodeIndicators(nodeNumber)==7) then ! node is internal outflow

     ! check if element has two internal outflow nodes, if so 
     ! add to chain
     if(j==1) then 
      ind1 = 2
      ind2 = 3
      ind3 = 4
     else if(j==2) then 
      ind1 = 1
      ind2 = 3
      ind3 = 4
     else if(j==3) then 
      ind1 = 1
      ind2 = 2
      ind3 = 4
     else
      ind1 = 1
      ind2 = 2
      ind3 = 3
     end if 

     elementIsAsItShould = .false.  
     
     counter1 = 1
     counter2 = 1
     searchNodeNumber = rep%elements(i)%pointIndexes(ind1) 
     if(searchNodeNumber.le.brp%numberOfBoundaryFaces) then 
      counter1 = counter1 + 1
      if(brp%boundaryNodeIndicators(searchNodeNumber)==7) then 
       counter2 = counter2 + 1
       elementIsAsItShould = .true.
      end if
     end if
    
     searchNodeNumber = rep%elements(i)%pointIndexes(ind2) 
     if(searchNodeNumber.le.brp%numberOfBoundaryFaces) then 
      counter1 = counter1 + 1
      if(brp%boundaryNodeIndicators(searchNodeNumber)==7) then
       counter2 = counter2 + 1
       elementIsAsItShould = .true.
      end if
     end if
    
     searchNodeNumber = rep%elements(i)%pointIndexes(ind3) 
     if(searchNodeNumber.le.brp%numberOfBoundaryFaces) then 
      counter1 = counter1 + 1
      if(brp%boundaryNodeIndicators(searchNodeNumber)==7) then
       counter2 = counter2 + 1
       elementIsAsItShould = .true.
      end if
     end if

     if(counter1==2.and.counter2==2) then 
      elementIsAsItShould = .true.
     else if(counter2==1) then 
      elementIsAsItShould = .false.
     else if(counter1==1) then 
      chainLengthAddition = chainLengthAddition + 1
      elementIsAsItShould = .true.
     else if(counter1==3.and.counter2==2) then
      chainLengthAddition = chainLengthAddition - 1
      elementIsAsItShould = .false.
     else if(counter1==4) then 
      write(*,*) "Warning: Rectangular element ",i," is cornered"
      STOP "stopped because of cornered element"
     end if    

     if(elementIsAsItShould) then 
     ! element is on IO boundary and not the first in chain 
     ! add to element list
      numberOfConnectedElements = numberOfConnectedElements + 1
      if(associated(lastElement)) then
       allocate(lastElement%next,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       lastElement => lastElement%next
       nullify(lastElement%next)
       lastElement%col = i
      else ! first element in list
       allocate(elementChain,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       allocate(elementChain%first,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       lastElement=>elementChain%first
       nullify(lastElement%next)
       lastElement%col = i
      end if
     else
     ! element is first or last in IO chain
     ! add to start list
      numberOfStartNodes = numberOfStartNodes + 1
      if(associated(startNodes)) then
       allocate(lastNode%next,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       allocate(lastNode%next%arr(2),stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       lastNode => lastNode%next
       nullify(lastNode%next)
      else
       ! first start node
       allocate(startNodes,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       allocate(startNodes%first,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       allocate(startNodes%first%arr(2),stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory"
       lastNode => startNodes%first
       nullify(lastNode%next)
      end if
      i1 = 0
      i2 = 0
      if(counter2==1) then 
       i1 = j+1
       if(i1==5) i1 = 1
       i2 = j-1
       if(i2==0) i2 = 4
      else if(counter1==3.and.counter2==2) then 
       ! first find internal node
       do m=1,4
        if(rep%elements(i)%pointIndexes(m)>brp%numberOfBoundaryFaces) then  
         i2 = m  
         exit     
        end if
       end do
       do m=1,4
        if(rep%elements(i)%pointIndexes(m).le.brp%numberOfBoundaryFaces) then  
         if(brp%boundaryNodeIndicators(rep%elements(i)%pointIndexes(m))/=7) then 
          i1 = m
          exit 
         end if
        end if
       end do 
      end if
      ! make sure that i1 points to first node in chain - i.e. is on boundary
      if(rep%elements(i)%pointindexes(i1)>brp%numberOfBoundaryFaces) then 
       buff = i2
       i2 = i1
       i1 = buff 
      end if   
 
      lastNode%arr(1) = rep%elements(i)%pointIndexes(i1)
      lastNode%arr(2) = rep%elements(i)%pointIndexes(i2) 
     end if
     exit 
    end if
   end if 
  end do
 end do


 do i=1,tep%numberOfElements
  do j=1,3
   nodeNumber = tep%elements(i)%pointIndexes(j) 
   if(nodeNumber.le.brp%numberOfBoundaryFaces) then ! node is on boundary 
    if(brp%boundaryNodeIndicators(nodeNumber)==7) then ! node is internal outflow
     ! check if element has two non - internal outflow nodes, if so 
     ! add to chain
     if(j==1) then 
      ind1 = 2
      ind2 = 3
     else if(j==2) then 
      ind1 = 1
      ind2 = 3
     else
      ind1 = 1
      ind2 = 2
     end if 
   
     elementIsAsItShould = .true.  
      
     do k=1,3
      if(k/=j) then
       ! now find all nodes connected to internal outflow node
       searchNodeNumber = tep%elements(i)%pointIndexes(k)
       if(searchNodeNumber.le.brp%numberOfBoundaryFaces) then ! connected node is external 
        elementIsAsItShould = .false. 
        if(brp%boundaryNodeIndicators(searchNodeNumber)/=7) then 
        ! element has two nodes on boundary but only one is IO, this is therefore
        ! the start of a chain
         numberOfStartNodes = numberOfStartNodes + 1
         if(associated(startNodes)) then 
          allocate(lastNode%next,stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
          allocate(lastNode%next%arr(2),stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
          lastNode => lastNode%next
          nullify(lastNode%next)
         else 
          ! first start node
          allocate(startNodes,stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
          allocate(startNodes%first,stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
          allocate(startNodes%first%arr(2),stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
          lastNode => startNodes%first
          nullify(lastNode%next)
         end if 
         lastNode%arr(1) = searchNodeNumber 
         if(searchNodeNumber==tep%elements(i)%pointIndexes(ind1)) then 
          lastNode%arr(2) = tep%elements(i)%pointIndexes(ind2)
         else
          lastNode%arr(2) = tep%elements(i)%pointIndexes(ind1) 
         end if                 
        end if
       end if 
      end if
     end do 
  
     if(elementIsAsItShould) then 
      numberOfConnectedElements = numberOfConnectedElements + 1
      if(associated(lastElement)) then 
       allocate(lastElement%next,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
       lastElement => lastElement%next
       nullify(lastElement%next)
       lastElement%col = i + rep%numberOfElements 
      else ! first element in list
       allocate(elementChain,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
       allocate(elementChain%first,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
       lastElement=>elementChain%first
       nullify(lastElement%next)
       lastElement%col = i + rep%numberOfElements
      end if
     end if 
     exit
    end if
!    exit
   end if
  end do
 end do 


 if(mod(numberOfStartNodes,2)/=0) &
 STOP "ERROR: Something's gone terribly wrong in makeInternalOutflowRegister"
 if(numberOfStartNodes>0) then 
  numberOfIOChains = max(numberOfStartNodes/2,1) 
 else
  numberOfIOChains = 0 
 end if
 
 brp%numberOfIORegisters = numberOfIOChains
 if(numberOfIOChains>0) then 
  allocate(brp%IORegisterLengths(numberOfIOChains),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 
 else
  nullify(brp%IORegisterLengths) 
 end if

 write(*,*) "number of IO chains: ",numberOfIOChains
 
 ! we now have a list of elements that contain all the nodes connected
 ! to internal outflow boundaries, now lets create an ordered chain of these
 ! to make interpolation fast

! lastNode=>startNodes%first
! do while(associated(lastNode)) 
!  write(*,*) "sverre: ",lastNode%arr
!  lastNode => lastNode%next
! end do

!  j = 0
!  currentElement=>elementChain%first
!  do while(associated(currentElement))
!   j = j + 1
!   write(*,*) "Jonas: ",currentElement%col,j
!   currentElement => currentElement%next
!  end do
 

 if(numberOfConnectedElements+numberOfStartNodes>0) then 
!  write(*,*) "cLA: ",chainLengthAddition
!  write(*,*) startNodes%first%arr
!  write(*,*) startNodes%first%next%arr
!  numberOfConnectedNodes = numberOfConnectedElements + 3*numberOfStartNodes/2 + chainLengthAddition
  numberOfConnectedNodes = numberOfConnectedElements + 3*numberOfStartNodes/2
!  numberOfConnectedNodes = numberOfConnectedElements + numberOfStartNodes
  brp%numberOfIONodes = numberOfConnectedNodes
 write(*,*) "aal: ",numberOfConnectedElements,brp%numberOfIONodes
  allocate(brp%internalOutflowRegister(numberOfConnectedNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister out of memory" 



  j = 0
  do i=1,numberOfIOChains
   j = j+1
   ! set first chain 
   if(numberOfStartNodes>0) then 
    brp%internalOutflowRegister(j) = startNodes%first%arr(1)
    write(*,*) "q1: ",brp%internalOutflowRegister(j)
    j = j+1
    brp%internalOutflowRegister(j) = startNodes%first%arr(2)
    write(*,*) "q2: ",brp%internalOutflowRegister(j)
    lastNode=>startNodes%first
    startNodes%first => startNodes%first%next
    deallocate(lastNode%arr,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
    deallocate(lastNode,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
   else
    ! this means that the chain is connected - this deserves a comment 
    ! to user
    write(*,*) "Warning: Internal outflow boundary is connected"
    ! just start at a random element
    elementNumber = elementChain%first%col
    currentElement => elementChain%first
    elementChain%first => currentElement%next
    deallocate(currentElement,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
    
    if(elementNumber.le.rep%numberOfElements) then 
     ind1 = 0
     ind2 = 0 
     do k=1,4 
      if(rep%elements(elementNumber)%pointIndexes(k)>brp%numberOfBoundaryFaces) then 
       if(ind1==0) then 
        ind1 = rep%elements(elementNumber)%pointIndexes(k)
       else
        ind2 = rep%elements(elementNumber)%pointIndexes(k)
        exit
       end if
      end if
     end do
    else
     elementNumber = elementNumber - tep%numberOfElements
     ind1 = 0
     ind2 = 0     
     do k=1,3
      if(tep%elements(elementNumber)%pointIndexes(k)>brp%numberOfBoundaryFaces) then 
       if(ind1==0) then 
        ind1 = tep%elements(elementNumber)%pointIndexes(k)
       else
        ind2 = tep%elements(elementNumber)%pointIndexes(k)
        exit
       end if
      end if
     end do
    end if     
    brp%internalOutflowRegister(j) = ind1
    j = j+1
    brp%internalOutflowRegister(j) = ind2
   end if
 

   ! find element in which last element in chain is found


   if(associated(elementChain)) then 
    currentElement=>elementChain%first
    nullify(previousElement)
    do while(associated(currentElement))
     elementNumber = currentElement%col
     if(elementNumber.le.rep%numberOfElements) then 
      ind1 = 0
      ind2 = 0     
      do k=1,4
       if(rep%elements(elementNumber)%pointIndexes(k)>brp%numberOfBoundaryFaces) then 
        if(ind1==0) then 
         ind1 = rep%elements(elementNumber)%pointIndexes(k)
        else
         ind2 = rep%elements(elementNumber)%pointIndexes(k)
         exit
        end if
       end if
      end do
     else
      elementNumber = elementNumber - rep%numberOfElements
      ind1 = 0
      ind2 = 0     
      do k=1,3
       if(tep%elements(elementNumber)%pointIndexes(k)>brp%numberOfBoundaryFaces) then 
        if(ind1==0) then 
         ind1 = tep%elements(elementNumber)%pointIndexes(k)
        else
         ind2 = tep%elements(elementNumber)%pointIndexes(k)
         exit
        end if
       end if
      end do
     end if

     hasChanged = .false.
     if(ind1==brp%internalOutflowRegister(j)) then 
      j = j + 1
      brp%internalOutflowRegister(j) = ind2
      hasChanged = .true.
     else if(ind2==brp%internalOutflowRegister(j)) then 
      j = j + 1
      brp%internalOutflowRegister(j) = ind1
      hasChanged = .true.
     end if
   
     if(hasChanged) then  
      if(associated(previousElement)) then 
       previousElement%next => currentElement%next
      else
       elementChain%first => currentElement%next
      end if  
      deallocate(currentElement,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
      if(associated(previousElement)) then  
       currentElement => previousElement%next
       previousElement => currentElement
       currentElement => currentElement%next
      else
       currentElement => elementChain%first
      end if
     else 
      currentElement => currentElement%next
     end if
    end do
   end if

  ! terminate chain
   currentNode => startNodes%first
   nullify(previousNode)
   do while(associated(currentNode)) 
    ind1 = currentNode%arr(1)
    ind2 = currentNode%arr(2)
    write(*,*) "aae: ",ind1,ind2,brp%internalOutflowRegister(j),j
!    if(ind1==brp%internalOutflowRegister(j)) then
!!    if(ind1==brp%internalOutflowRegister(j).or.ind2==brp%internalOutflowRegister(j)) then
!     j = j+1 
!     brp%internalOutflowRegister(j) = ind2
!     if(associated(previousNode)) then 
!      previousNode%next => currentNode%next
!     else
!      startNodes%first => startNodes%first%next
!     end if
!     deallocate(currentNode%arr,stat=allocateStatus)
!     if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
!     deallocate(currentNode,stat=allocateStatus)
!     if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
!     exit
!    else if(ind2==brp%internalOutflowRegister(j)) then 
!!    else  
     if(ind2==brp%internalOutflowRegister(j)) then
     j = j+1 
     brp%internalOutflowRegister(j) = ind1



     if(associated(previousNode)) then
      previousNode%next => currentNode%next
     else
      startNodes%first => startNodes%first%next
     end if
     deallocate(currentNode%arr,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"
     deallocate(currentNode,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeInternalOutflowRegister couldn't deallocate"



     exit 
    end if
    previousNode => currentNode
    currentNode => currentNode%next
   end do
   brp%IORegisterLengths(i) = j
  end do  
 else
  nullify(brp%internalOutflowRegister)
 end if 

! write(*,*) "IO register:"
! do i=1,brp%numberOfIONodes
!  write(*,*) i,brp%internalOutflowRegister(i) 
! end do


 end subroutine makeInternalOutflowRegister
!-------------------------------------------------------------------------
subroutine makeEngineInflowData(crp,brp)

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp ! contains coordinate register for grid
 type(BoundaryRegisterData) :: brp

 integer :: i,p1,p2,allocateStatus
 double precision :: x1(2),x2(2)
 
! added for engine inlet
 brp%numberOfEngineInletSideLimit = 0
 do i=1,brp%numberOfBoundaryFaces
  if(brp%faces(i)%indicator == 8) then
   brp%numberOfEngineInletSideLimit = brp%numberOfEngineInletSideLimit + 1
  end if
 end do


 allocate(brp%engineInletSideIndexes(brp%numberOfEngineInletSideLimit,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(brp%engineInletSideCoefficients(brp%numberOfEngineInletSideLimit,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 brp%numberOfEngineInletSides = 0
 brp%engineInletSideCoefficients = 0.0
 brp%engineInletSideIndexes = 0

 write(*,*) "Making engine inlet registers..."

 do i=1,brp%numberOfBoundaryFaces
  if(brp%faces(i)%indicator == 8) then
   p1=brp%faces(i)%faceIndexes(1)
   p2=brp%faces(i)%faceIndexes(2)
   x1 = getCoor(p1,crp)
   x2 = getCoor(p2,crp)   
   call registerEngineInletSide(p1,p2,x2-x1,brp)
  end if
 end do


 write(*,*) "Number of engine inlet sides: ",brp%numberOfEngineInletSides

 end subroutine makeEngineInflowData
 !-------------------------------------------------------------------------
 subroutine registerEngineInletSide(p1,p2,Cxx,brp)
 IMPLICIT NONE

 integer :: p1,p2
 double precision :: Cxx(2)
 type(BoundaryRegisterData) :: brp

 double precision :: C12(2)
 integer :: xp1,xp2,i,sideNumber,sp1,sp2

 C12 = Cxx

 if(p1>p2) then
  xp1 = p2
  xp2 = p1
 else
  xp1 = p1
  xp2 = p2
 end if

 ! find out if side exists
 sideNumber = 1
 do i=1,brp%numberOfEngineInletSideLimit
  sp1 = brp%engineInletSideIndexes(i,1)
  sp2 = brp%engineInletSideIndexes(i,2)
  if(sp1>0) then
   if(xp1==sp1.and.xp2==sp2) then
    exit
   else
    sideNumber = sideNumber + 1
   end if
  else
   sideNumber = -sideNumber
   brp%numberOfEngineInletSides = brp%numberOfEngineInletSides + 1
   exit
  end if
 end do

 if(sideNumber>0) then
  ! side exists
   brp%engineInletSideCoefficients(sideNumber,1) = brp%engineInletSideCoefficients(sideNumber,1) + C12(2)
   brp%engineInletSideCoefficients(sideNumber,2) = brp%engineInletSideCoefficients(sideNumber,2) - C12(1)
 else
  ! create new side
  sideNumber = -sideNumber
  if(sideNumber.le.brp%numberOfEngineInletSideLimit) then
   brp%engineInletSideIndexes(sideNumber,1) = xp1
   brp%engineInletSideIndexes(sideNumber,2) = xp2
   brp%engineInletSideCoefficients(sideNumber,1) = C12(2)
   brp%engineInletSideCoefficients(sideNumber,2) = -C12(1)    
  else
   STOP "ERROR: Something's wrong in registerEngineInletSide - 1"
  end if
 end if

 end subroutine registerEngineInletSide 
!--------------------------------------------------------------------------
 subroutine colorBoundaryFaces(brp)
 ! colors boundary faces for vectorization 
 ! if boundary is correctly set up, maximum three colors are needed
 ! in the fine grid, in the agglomerated grids however number of colors
 ! are only bounded by the number og boundary faces
 ! This subroutine is heavily complicated for two dimensions, but 
 ! may be warranted for 3D applications

 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp ! contains grid data  

 integer :: optNumbOfFacesPerColor,numberOfFaces
 integer :: i,j,k,allocateStatus,slack
 logical :: notFinished

 type(LinkedFace),pointer :: currentFace, previousFace, nextFace
 type(BoundaryColorData),pointer :: currentColor,searchColor
 type(LinkedInteger),pointer :: swapInteger
 type(BoundaryColorList),pointer :: colors

 write(*,*) "Coloring boundary faces..."

 allocate(brp%colorsOccupied(brp%numberOfBoundaryFaces),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 

 allocate(colors,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 

 do i=1,brp%numberOfBoundaryFaces
  nullify(brp%colorsOccupied(i)%first)
 end do

! start coloring
! a linked list of colors is created, each color has a linked
! list of faces in it 
 
 do i=1,brp%numberOfBoundaryFaces !run through all points
!  write(*,*) 'Node no.: ', i
  allocate(currentFace,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 
  allocate(currentFace%bd,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 
   currentFace%bd = brp%faces(i) ! makes a copy  
   currentColor => colors%first
   do while(associated(currentColor).and.associated(currentFace))
    if(faceFitInColor(currentFace,currentColor%colorNumber,brp)) then
! face fits in color   
     currentColor%last%next=>currentFace
     do k=1,2
      j = currentFace%bd%faceIndexes(k)  
      swapInteger => brp%colorsOccupied(j)%first
      allocate(brp%colorsOccupied(j)%first,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 
      brp%colorsOccupied(j)%first%col = currentColor%colorNumber
      brp%colorsOccupied(j)%first%next => swapInteger 
     end do   
     currentColor%last=>currentFace
     currentColor%numberOfSides = currentColor%numberOfSides+1
     nullify(currentColor%last%next)
     nullify(currentFace)
     exit
    else   
! switch color
     currentColor => currentColor%next 
    end if 
   end do ! while associated(currentColor)
   if(associated(currentFace)) then 
! have run through all colors, but current side hasn't been placed
! this means that the current side couldn't be fitted into one of the 
! existing colors, must then create new color

    if(associated(colors%first)) then 
     allocate(colors%last%next,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 
     colors%last => colors%last%next
    else
     ! at the beginning, no colors allocated
     allocate(colors%first,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 
     colors%last => colors%first 
     colors%numberOfColors = 0
    end if
    
    nullify(colors%last%next)

     do k=1,2
      j = currentFace%bd%faceIndexes(k)  
      swapInteger => brp%colorsOccupied(j)%first
      allocate(brp%colorsOccupied(j)%first,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 
      brp%colorsOccupied(j)%first%col = colors%numberOfColors+1
      brp%colorsOccupied(j)%first%next => swapInteger 
     end do   
! this color is empty so the side must fit
    colors%last%first=>currentFace
    colors%last%last=>currentFace
    colors%last%numberOfSides = 1
    colors%numberOfColors = colors%numberOfColors + 1   
    colors%last%colorNumber = colors%numberOfColors
    nullify(colors%last%last%next)
   end if
 end do !number of nodes 

! try to get the color sizes similar 

 optNumbOfFacesPerColor = brp%numberOfBoundaryFaces/colors%numberOfColors
 brp%numberOfColors = colors%numberOfColors

 currentColor => colors%first
 do while(associated(currentColor)) ! run through all colors
  if(currentColor%colorNumber == 1) then ! try to avoid this 
   slack = brp%numberOfBoundaryFaces-optNumbOfFacesPerColor*brp%numberOfColors
  else
   slack = 0
  end if
  currentFace => currentColor%first
  nullify(previousFace)
  nextFace => currentFace%next
  do while ((currentColor%numberOfSides>(optNumbOfFacesPerColor+slack)).and.&
            (associated(currentFace)))
   searchColor => colors%first 
   do while (associated(searchColor)) 
    if (searchColor%numberOfSides<optNumbOfFacesPerColor) then
     if(faceFitInColor(currentFace,searchColor%colorNumber,brp)) then
     ! transfer side from currentColor to searchColor
      searchColor%last%next => currentFace
      searchColor%last => currentFace

      if(associated(previousFace)) then 
       previousFace%next=>currentFace%next ! jump over currentSide in 
                                           ! currentColor
      else ! must mean that currentFace is first in color
       currentColor%first=>currentFace%next
      end if
      nullify(currentFace%next) 
      currentColor%numberOfSides = currentColor%numberOfSides-1
      searchColor%numberOfSides = searchColor%numberOfSides+1   
     ! fix colorsOccupied register
      do k=1,2
       j = currentFace%bd%faceIndexes(k)  
       swapInteger=>brp%colorsOccupied(j)%first
       do 
        if(.not.associated(swapInteger)) then
         STOP "Something wrong with color register"
        end if
        if (swapInteger%col == currentColor%colorNumber) then
         ! swap sides
         swapInteger%col = searchColor%colorNumber
         exit
        end if
        swapInteger=>swapInteger%next 
       end do
      end do ! k=1,2 
      currentFace=>previousFace ! this will lead to that the previous side
                                ! of next iteration will equal the one in this  
      exit
     end if ! side fit in color
    end if ! searchColor%numberOfSides<optNumbOfSidesPerColor  
    searchColor => searchColor%next
   end do ! while associated(searchColor)
   previousFace => currentFace
   currentFace => nextFace 
   if(associated(currentFace)) then 
    nextFace => currentFace%next
   end if
  end do
  currentColor => currentColor%next
 end do ! while associated(currentColor)
 
 write(*,*) 'Number of face colors : ',colors%numberOfColors

! finally place faces, ordered by their color, into brp%faces 

 i = 1
 j = 1
 allocate(brp%startOfColors(colors%numberOfColors),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces out of memory" 

 currentColor=>colors%first
 do while(associated(currentColor))
  currentFace=>currentColor%first
  do while(associated(currentFace))
   brp%faces(i) = currentFace%bd
   previousFace => currentFace
   currentFace => currentFace%next
   deallocate(previousFace,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces couldn't deallocate"
   i = i+1 
  end do
  brp%startOfColors(j) = i  
  j = j + 1
  searchColor => currentColor
  currentColor => currentColor%next
  deallocate(searchColor,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces couldn't deallocate"
 end do 

 deallocate(colors,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: colorBoundaryFaces couldn't deallocate"

 end subroutine colorBoundaryFaces 
!-------------------------------------------------------------------------
 logical function faceFitInColor(face,colorNumber,brp)
! Find out whether side is allowable in color or not
! Side is only allowable in color if it has no nodes
! in common with the sides allready in color and
! if the color is not full. The color is full if 
! the number of sides in it exceeds 
! sidesInColors + numberOfColors
 IMPLICIT NONE

 type(LinkedFace) :: face 
 integer :: colorNumber
 type(BoundaryRegisterData) :: brp

 type(LinkedInteger),pointer :: current
 integer :: k

  faceFitInColor = .true.
  do k=1,2
   current => brp%colorsOccupied(face%bd%faceIndexes(k))%first 
   do while(associated(current))
    if (current%col == colorNumber) then ! does not fit 
     faceFitInColor = .false.
     goto 55
    end if 
    current => current%next
   end do
55  end do 
 end function faceFitInColor
!-------------------------------------------------------------------------
end module BoundaryRegister

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module Grid
! holds variables and procedures for a single grid
 use SideRegister
 use CoordinateRegister
 use TriangularElements
 use RectangularElements
 use BoundaryRegister
 use VisualizationModule

 type LinkedSide ! a side data type for use in linked lists
  type(SideData), pointer :: sd
  type(LinkedSide), pointer :: next
 end type LinkedSide
 
 type SideDataPointer
  type(SideData),pointer :: sd
 end type

 type FirstPointer
  type(LinkedSide), pointer :: first,last 
 end type FirstPointer

 type LinkedDoubleInteger
   integer :: index(2)
   type(LinkedDoubleInteger),pointer :: next
 end type LinkedDoubleInteger

 type LinkedDoubleIntegerP
  type(LinkedDoubleInteger),pointer :: first,last 
 end type LinkedDoubleIntegerP

 type ColorData
  type(LinkedSide),pointer :: first,last
  type(ColorData),pointer :: next
  integer :: numberOfSides,colorNumber
 end type ColorData 

 type ColorList
  type(ColorData),pointer :: first,last
  integer :: numberOfColors
 end type ColorList 
 
 type ControlVolume
  type(FirstPointer),pointer :: sides
 end type ControlVolume


 type GridData
  type(BoundaryRegisterData),pointer :: brp ! holds boundary data in grid
  real,pointer :: controlVolumeAreas(:)
  type(SideDataPointer),pointer :: srp(:) ! to hold side or superside register 
  type(LinkedIntegerP),pointer :: colorsOccupied(:) ! holds information on
                                                    ! which colors index
                                                    ! is used in
  type(VisualizationData),pointer :: vdp
  type(ColorList),pointer :: colors ! holds color information
  real,pointer :: controlVolumeNodeSize(:) ! contains area of control volumes
                                           ! to be used for intergrid mappings
!  real,pointer :: basisFunctionIntegration(:) ! contains three times the integral of
                                              ! the basis function over entire domain,
                                              ! which is (three times the) support area 
  real,pointer :: wallDistance(:) ! distance from any point to closest wall, used in 
                                  ! turbulence modeling
  integer,pointer :: wallDistanceFaceIndex(:) ! index of closest face


  real,pointer :: coordinates(:,:)
  type(FirstPointer),pointer :: sidePointers(:)
  type(FirstPointer),pointer :: visualizationSidePointers(:)
  type(FirstPointer),pointer :: superSidePointers(:)
  type(LinkedIntegerP),pointer :: controlVolumeNodePointers(:)
  integer,pointer :: boundaryFaceNodeMappings(:) ! figures out to which control
                                                 ! volumes a face belongs 
  integer :: numberOfSides,numberOfSuperSides,numberOfNodes
  integer :: optimalNumberOfColors,sidesInColors,numberOfCVSides
  real :: directionalityParameter,minimumAspectRatio 
  integer :: numberOfBoundaryCVs
 end type GridData

 contains

!-------------------------------------------------------------------------
 subroutine setBoundaryFirst2(grp,crp,rep,tep,brp)
 ! subroutine to make sure that the first nodes are on boundary

 IMPLICIT NONE

 type(GridData) :: grp
 type(CoordinateRegisterData) :: crp
 type(RectangularElementData) :: rep
 type(TriangularElementData) :: tep
 type(BoundaryRegisterData) :: brp

 integer :: i,jstart,kstart,j,k,ind1,nbf,switchNode,l,allocateStatus,searchNode
 real :: coor(2)
 integer,pointer :: swapList(:,:)
 integer :: numberOfSwitches,ind
 logical :: itWentAsPlanned
 integer,pointer :: nodeSwitchRegister(:)
 logical,pointer :: nodeOK(:)
 integer :: count1,count2,buff

 type(STEData),pointer :: triTemp  ! single triangular element pointer
 type(SREData),pointer :: rectTemp  ! single rectangular element pointer

 write(*,*) "Reordering..."

 allocate(nodeOK(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"

 allocate(swapList(brp%numberOfBoundaryNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"

 brp%numberOfBoundaryNodes = brp%numberOfBoundaryFaces

 nodeOK = .true.
 nodeOK(1:brp%numberOfBoundaryFaces) = .false.
 count1 = 0
 count2 = 0

 do i=1,brp%numberOfBoundaryFaces
  do j=1,2
   ind = brp%faces(i)%faceIndexes(j)
   if(ind.le.brp%numberOfBoundaryFaces) then
    if(.not.nodeOK(ind)) count2 = count2 + 1
    nodeOK(ind) = .true.
   else
    if(nodeOK(ind)) count1 = count1 + 1
    nodeOK(ind) = .false.
   end if
  end do
 end do

 ! now look for nodes that are larger

 numberOfSwitches = 0

 do i=1,brp%numberOfBoundaryFaces
  do j=1,2
   ind = brp%faces(i)%faceIndexes(j)
   if(.not.nodeOK(ind)) then
    ! find node to switch with
    itWentAsPlanned = .false.
    do k=1,brp%numberOfBoundaryFaces
     if(.not.nodeOK(k)) then
      numberOfSwitches = numberOfSwitches + 1
      swapList(numberOfSwitches,1) = ind
      swapList(numberOfSwitches,2) = k
      nodeOK(k) = .true.
      nodeOK(ind) = .true.
      itWentAsPlanned = .true.
      exit
     end if
    end do
    if(.not.itWentAsPlanned) write(*,*) "Switches: ",numberOfSwitches
    if(.not.itWentAsPlanned) STOP "ERROR: Something's gone wrong in setBoundaryFirst - 1"
   end if
  end do
 end do

 deallocate(nodeOK,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst couldn't deallocate"

! switch coordinates

 do i=1,numberOfSwitches
  coor = getCoor(swapList(i,2),crp)
  call setCoor(swapList(i,2),getCoor(swapList(i,1),crp),crp)
  call setCoor(swapList(i,1),coor,crp)
 end do

 ! to speed things up, make a full node register

 allocate(nodeSwitchRegister(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"

 do i=1,grp%numberOfNodes
  nodeSwitchRegister(i) = i
 end do

 do i=1,numberOfSwitches
  buff = nodeSwitchRegister(swapList(i,1))
  nodeSwitchRegister(swapList(i,1)) = nodeSwitchRegister(swapList(i,2))
  nodeSwitchRegister(swapList(i,2)) = buff
 end do

! switch element indexes

 do j=1,rep%numberOfElements
  rectTemp => getRectElement(j,rep)
  do k=1,4
   rectTemp%pointIndexes(k) = nodeSwitchRegister(rectTemp%pointIndexes(k))
  end do
 end do

 do j=1,tep%numberOfElements
  triTemp => getTriElement(j,tep)
  do k=1,3
   triTemp%pointIndexes(k) = nodeSwitchRegister(triTemp%pointIndexes(k))
  end do
 end do

! finally do switching in boundary register

 do j=1,brp%numberOfBoundaryFaces
  do k=1,2
   brp%faces(j)%faceIndexes(k) = nodeSwitchRegister(brp%faces(j)%faceIndexes(k))
  end do
 end do 


 ! write switch list

! write(*,*) "Writing switch file..."
! open(21,file="switch.reg",form='formatted',status='unknown')
! write(21,*) numberOfSwitches
! do i=1,numberOfSwitches
!  write(21,*) swapList(i,:)
! end do
! close(21)

! deallocate(swapList,stat=allocateStatus)
! if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst couldn't deallocate"
 
 do i=1,brp%numberOfBoundaryFaces
  brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1)) = brp%faces(i)%indicator
 end do

 ! write new base file

! write(*,*) "Writing base file..."
! open(27,file="base.plt",form='formatted',status='unknown')

! write(27,*) 1
! write(27,*) "title"
! write(27,*) "nelem   npoin   nboun nrect ntrian"
!! write(27,*) tep%numberOfElements+rep%numberOfElements,&
!! write(27,*) tep%numberOfElements+2*rep%numberOfElements,&
!!             crp%numberOfPoints,brp%numberOfBoundaryFaces,&
!!             rep%numberOfElements,tep%numberOfElements
! write(27,*) "connectivities"
! do i=1,rep%numberOfElements
!  rectTemp => getRectElement(i,rep)
!!  write(27,*) i,rectTemp%pointIndexes,0
!   write(27,*) 2*i-1,rectTemp%pointIndexes(1:3),0,0
!   write(27,*) 2*i,rectTemp%pointIndexes(1),rectTemp%pointIndexes(3),&
!                 rectTemp%pointIndexes(4),0,0
! end do
! do i=1,tep%numberOfElements
!  triTemp => getTriElement(i,tep)
!  write(27,*) i+2*rep%numberOfElements,triTemp%pointIndexes,0
! end do
! write(27,*) "coordinates"
! do i=1,crp%numberOfPoints
!  write(27,'(I8,5E20.8)') i,getCoor(i,crp),0,0
! end do
! write(27,*) "unknowns"
! do i=1,crp%numberOfPoints
!  write(27,*) "6*0.0"
! end do
! write(27,*) "boundary faces"
! do i=1,brp%numberOfBoundaryFaces
!  write(27,*) brp%faces(i)%faceIndexes,0,brp%boundaryNodeIndicators(i),brp%boundaryNodeIndicators(i)
! end do
! close(27)
 

 end subroutine setBoundaryFirst2
!-------------------------------------------------------------------------
 subroutine setBoundaryFirst(grp,crp,rep,tep,brp)
 ! subroutine to make sure that the first nodes are on boundary

 IMPLICIT NONE

 type(GridData) :: grp
 type(CoordinateRegisterData) :: crp
 type(RectangularElementData) :: rep
 type(TriangularElementData) :: tep
 type(BoundaryRegisterData) :: brp

 integer :: i,jstart,kstart,j,k,ind1,nbf,switchNode,l,allocateStatus,searchNode
 real :: coor(2)
 logical,pointer :: pointTaken(:)
 logical :: notOnBoundary,foundOne
 integer,pointer :: swapList(:,:) 
 integer :: numberOfSwitches,ind
 logical :: itWentAsPlanned

 type(STEData),pointer :: triTemp  ! single triangular element pointer
 type(SREData),pointer :: rectTemp  ! single rectangular element pointer

 write(*,*) "Reordering..."
 
 nbf = brp%numberOfBoundaryFaces

 allocate(pointTaken(nbf),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"
 allocate(swapList(brp%numberOfBoundaryFaces,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"

 numberOfSwitches = 0
 searchNode = nbf
 ! first fix boundary
 do i=1,nbf
  if(brp%faces(i)%faceIndexes(1)/=i) then
   if(brp%faces(i)%faceIndexes(1).le.nbf) then
    numberOfSwitches = numberOfSwitches + 1
    ! find node to switch with
    searchNode = searchNode + 1
    foundOne = .false.
    do while(.not.foundOne.and.searchNode.le.crp%numberOfPoints)
     foundOne = .true.
     do j=1,nbf
      if(brp%faces(j)%faceIndexes(1)==searchNode) then
       foundOne = .false.
       searchNode = searchNode + 1
       exit
      end if
     end do
    end do
    if(.not.foundOne) STOP "ERROR: Something's wrong in setBoundaryFirst - 1"
    swapList(numberOfSwitches,1) = brp%faces(i)%faceIndexes(1)
    swapList(numberOfSwitches,2) = searchNode
   end if
  end if
 end do

 do i=1,nbf
  if(brp%faces(i)%faceIndexes(1)/=i) then
   numberOfSwitches = numberOfSwitches + 1
   if(brp%faces(i)%faceIndexes(1)>nbf) then
    swapList(numberOfSwitches,1) = brp%faces(i)%faceIndexes(1)
   else
    ! find new node number
    k = 0
    do j=1,numberOfSwitches
     if(swapList(j,1)==brp%faces(i)%faceIndexes(1)) then
      k = swapList(j,2)
      exit
     end if
    end do
    if(k==0) k = brp%faces(i)%faceIndexes(1)
    swapList(numberOfSwitches,1) = k
   end if
   swapList(numberOfSwitches,2) = i
  end if
 end do

! switch coordinates

 do i=1,numberOfSwitches 
  coor = getCoor(swapList(i,2),crp) 
  call setCoor(swapList(i,2),getCoor(swapList(i,1),crp),crp)
  call setCoor(swapList(i,1),coor,crp) 
 end do
! switch element indexes

 do i=1,numberOfSwitches
  do j=1,rep%numberOfElements
   rectTemp => getRectElement(j,rep) 
   do k=1,4
    if(rectTemp%pointIndexes(k) == swapList(i,1)) then
     rectTemp%pointIndexes(k) = swapList(i,2) 
    else if(rectTemp%pointIndexes(k) == swapList(i,2)) then
     rectTemp%pointIndexes(k) = swapList(i,1) 
    end if
   end do
  end do
 end do
 
 do i=1,numberOfSwitches
  do j=1,tep%numberOfElements
   triTemp => getTriElement(j,tep) 
   do k=1,3
    if(triTemp%pointIndexes(k) == swapList(i,1)) then
     triTemp%pointIndexes(k) = swapList(i,2) 
    else if(triTemp%pointIndexes(k) == swapList(i,2)) then
     triTemp%pointIndexes(k) = swapList(i,1) 
    end if
   end do
  end do
 end do

! finally do switching in boundary register
 do i=1,numberOfSwitches
  do j=1,nbf 
   do k=1,2
    if(brp%faces(j)%faceIndexes(k)==swapList(i,1)) then 
     brp%faces(j)%faceIndexes(k) = swapList(i,2) 
    end if
   end do
  end do
 end do

! set boundary indicator array

 do i=1,nbf
  brp%boundaryNodeIndicators(brp%faces(i)%faceIndexes(1)) = brp%faces(i)%indicator
 end do

 write(*,*) "BBX: ",brp%boundaryNodeIndicators(271)

 ! write switch list

! write(*,*) "Writing switch file..."
! open(21,file="switch.reg",form='formatted',status='unknown')
! write(21,*) numberOfSwitches 
! do i=1,numberOfSwitches
!  write(21,*) swapList(i,:) 
! end do
! close(21)

! deallocate(swapList,stat=allocateStatus)
! if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst couldn't deallocate"
  

 ! write new base file
 
! write(*,*) "Writing base file..."
! open(27,file="base.plt",form='formatted',status='unknown')

! write(27,*) 1
! write(27,*) "title"
! write(27,*) "nelem   npoin   nboun nrect ntrian"
!! write(27,*) tep%numberOfElements+rep%numberOfElements,&
!! write(27,*) tep%numberOfElements+2*rep%numberOfElements,&
!!             crp%numberOfPoints,brp%numberOfBoundaryFaces,&
!!             rep%numberOfElements,tep%numberOfElements
! write(27,*) "connectivities"
! do i=1,rep%numberOfElements
!  rectTemp => getRectElement(i,rep)
!!  write(27,*) i,rectTemp%pointIndexes,0 
!   write(27,*) 2*i-1,rectTemp%pointIndexes(1:3),0,0
!   write(27,*) 2*i,rectTemp%pointIndexes(1),rectTemp%pointIndexes(3),&
!                 rectTemp%pointIndexes(4),0,0
! end do
! do i=1,tep%numberOfElements
!  triTemp => getTriElement(i,tep)
!  write(27,*) i+2*rep%numberOfElements,triTemp%pointIndexes,0 
! end do
! write(27,*) "coordinates"
! do i=1,crp%numberOfPoints
!  write(27,'(I5,5E14.4)') i,getCoor(i,crp),0,0
! end do 
! write(27,*) "unknowns"
! do i=1,crp%numberOfPoints
!  write(27,*) "6*0.0"
! end do
! write(27,*) "boundary faces"
! do i=1,brp%numberOfBoundaryFaces
!  write(27,*) brp%faces(i)%faceIndexes,0,brp%boundaryNodeIndicators(i),brp%boundaryNodeIndicators(i)
! end do
! close(27)
 end subroutine setBoundaryFirst
!-------------------------------------------------------------------------
 subroutine setUpBoundary(crp,brp,grp,rep,tep)
! calculates boundary normals for each node on boundary
! Boundary indicator convention:
!  1: inviscid wall
!  2: symmetry surface
!  3-4: far field (can have two different)
!  5: isothermal viscous wall
!  6: adiabatic viscous wall  
!  7: internal outflow (f.ex in channel flow)
!  8: engine inlet % Manon

 IMPLICIT NONE
 
 type(CoordinateRegisterData) :: crp ! contains coordinate register for grid
 type(BoundaryRegisterData) :: brp  ! contains boundary face information in 2D          
 type(GridData) :: grp ! holds grid data 
 type(RectangularElementData) :: rep
 type(TriangularElementData) :: tep
 
 integer :: i 


! first create boundary tangents 
! can now assume that coordinates are set up correctly from file

  call makeBoundaryTangents(brp,crp)

! find internal outflow chains 

  call makeInternalOutflowRegister(brp,rep,tep)

! for engine inlet (8)
  call makeEngineInflowData(crp,brp)
 
! now search for trailing edges

  call searchForTrailingEdges(brp,crp)

! color boundary 

!  call colorBoundaryFaces(brp)

  grp%numberOfBoundaryCVs = brp%numberOfBoundaryFaces

 end subroutine setUpBoundary
!-------------------------------------------------------------------------
 subroutine makeCalculationData(crp,rep,tep,grp,brp,doVisualization)

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp ! contains coordinate register for grid
 type(RectangularElementData) :: rep ! contains element register for rectangular grid
 type(TriangularElementData) :: tep ! contains element register for triangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 integer :: allocateStatus,i


  grp%numberOfNodes = crp%numberOfPoints 

! allocate

  allocate(grp%sidePointers(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 

  allocate(grp%controlVolumeNodeSize(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 

  if(doVisualization) then 
   allocate(grp%vdp,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 
   allocate(grp%vdp%vsp(grp%numberOfNodes),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 

   do i=1,grp%numberOfNodes
    nullify(grp%vdp%vsp(i)%first)
    nullify(grp%vdp%vsp(i)%last)
   end do
  end if

! initialize

  do i=1,grp%numberOfNodes
   nullify(grp%sidePointers(i)%first)
   nullify(grp%sidePointers(i)%last)
  end do
  
  grp%numberOfSides = 0
  grp%controlVolumeNodeSize = 0
  write(*,*) "Processing rectangular elements..."
  call makeCalcDataForRect(crp,rep,grp,brp,doVisualization)
  write(*,*) "Processing triangular elements..."
  call makeCalcDataForTri(crp,tep,grp,brp,doVisualization,rep%numberOfElements)

  write(*,*) "Processing engine inlet faces..."
  call makeEngineInflowData(crp,brp)


  write(*,*) "Number of sides in fine grid: ",grp%numberOfSides
  
  
 
 end subroutine makeCalculationData
!-------------------------------------------------------------------------
 subroutine makeCalcDataForRect(crp,rep,grp,brp,doVisualization)
 ! makes side coefficients from rectangular elements and splits
 ! data structure into side based

 IMPLICIT NONE
 
 type(CoordinateRegisterData) :: crp ! contains coordinate register for grid
 type(RectangularElementData) :: rep ! contains element register for rectangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 type(LinkedSide),pointer :: current
 type(VisualizationSidePointer),pointer :: currentVis

 type(SREData),pointer :: temp  ! single triangular element pointer
 integer, parameter :: NSD = 2
 integer, parameter :: NODES_IN_ELEMENT = 4 
 
 integer ::  nel ! number of elements in REP
 integer :: i,j,k,l,allocateStatus,total,p1,p2,p3,p4,sp1,sp2,ptemp,elmind1,elmind2 
 integer :: cj,ck,cl
 real :: xm(NSD) ! to contain element median
 real :: xs(NSD) ! to contain side median
 real :: xc(NSD) ! to contain side coefficient
 real :: xcs(NSD)  
 real :: xf(NSD),v1(NSD),v2(NSD)
 real :: x1(NSD),x2(NSD),x3(NSD),x4(NSD) ! contains side edge coordinates
 real :: C12(NSD),D12(NSD)
 logical :: isAntiClockwise
 integer :: in1,in2,count1,count2
 real :: elementArea,sign

 logical,pointer :: validate(:)

! Do not know how many sides there are in grid, therefore incorporate
! a linked list to dynamically create new sides as we go along. 

  nel = rep%numberOfElements
 
  count1 = 0
  count2 = 0

  do i = 1,nel ! do for each element
  temp => getRectElement(i,rep) 
! find centroid in element 
  xm = 0. 
  do j = 1,NODES_IN_ELEMENT
   xm = xm + getCoor(temp%pointIndexes(j),crp)
  end do
  xm = xm/NODES_IN_ELEMENT
  p1=temp%pointIndexes(1)
  p2=temp%pointIndexes(2)
  p3=temp%pointIndexes(3)
  p4=temp%pointIndexes(4)
  x1 = getCoor(p1,crp)
  x2 = getCoor(p2,crp)
  x3 = getCoor(p3,crp)
  x4 = getCoor(p4,crp)
  elementArea = 0.5*abs((x2(1)-x1(1))*(x3(2)-x1(2))-(x2(2)-x1(2))*(x3(1)-x1(1)))&
              + 0.5*abs((x3(1)-x1(1))*(x4(2)-x1(2))-(x3(2)-x1(2))*(x4(1)-x1(1))) 

! run through points in element to create sides 
  do j = 1,NODES_IN_ELEMENT -1
   do k = j+1,NODES_IN_ELEMENT
    if(abs(j-k)/=2) then 
    if(j==1.and.k==4) then ! get numbering to follow right hand rule
     cj = k
     ck = j 
    else
     cj = j
     ck = k
    end if
    cl = ck + 2 
    if(cl>4) cl = cl - 4 
    p1=temp%pointIndexes(cj)
    p2=temp%pointIndexes(ck)

    
   if(p1 > p2) then ! sort indices to always have smallest index first
     ptemp = p1
     p1 = p2
     p2 = ptemp
     ptemp = cj
     cj=ck
     ck=ptemp
     sign = -1.0
    else
     sign = 1.0
    end if 
    current => findSide(p1,p2,grp)
    if(doVisualization) then 
     currentVis => findVisSide(p1,p2,grp) 
    end if

    ! set controlVolumeNodesSize array as areas of CV around node
    l = 6-j-k
    x1 = getCoor(p1,crp)
    x2 = getCoor(p2,crp)
    xm = getCoor(temp%pointIndexes(1),crp)+&
         getCoor(temp%pointIndexes(2),crp)+&
         getCoor(temp%pointIndexes(3),crp)+&
         getCoor(temp%pointIndexes(4),crp)
    xm = 0.25*xm
    v1 = xm - x1
    v2 = 0.5*(x2-x1)
    grp%controlVolumeNodeSize(p1) = grp%controlVolumeNodeSize(p1)+&
           0.5*abs(v1(1)*v2(2)-v2(1)*v1(2))

    v1 = xm - x2
    v2 = 0.5*(x1-x2)
    grp%controlVolumeNodeSize(p2) = grp%controlVolumeNodeSize(p2)+&
           0.5*abs(v1(1)*v2(2)-v2(1)*v1(2))

    if(associated(current)) then
!    side already exists, do not create new side but add index in element
!    register if visualization is needed
     if(doVisualization) then 
      currentVis%vs%neighbouringElements(2) = i 
     end if
     count1 = count1 + 1
    else

!   create new side
     allocate(current,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
     allocate(current%sd,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
     call sideConstruct(NSD,current%sd)

! insert data into side

     current%sd%sideIndexes(1) = p1  
     current%sd%sideIndexes(2) = p2 
    
     current%sd%sideCoefficients = 0.0
     allocate(current%sd%sideLength(2),stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
     current%sd%sideLength = x2 - x1 
     nullify(current%sd%invertedSideLength)
     
     count2 = count2 + 1
! For visualization, the old indices of each side is needed since agglomerated
! sides don't have physical nodes and thus can't directly be plotted
! Boundary sides need their original nodes since they are connected to 
! a global register (boundaries aren't agglomerated). In addition this
! yields a very fast method to check if side is on boundary since the 
! node numbers in original grid starts with boundary nodes 
 
     if(doVisualization) then ! do visualization
      allocate(currentVis,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
      allocate(currentVis%vs,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
      currentVis%vs%sideIndexes(1) = p1 
      currentVis%vs%sideIndexes(2) = p2 
      currentVis%vs%neighbouringElements(1) = i
      currentVis%vs%neighbouringElements(2) = 0 
      currentVis%vs%controlVolumes(1) = p1
      currentVis%vs%controlVolumes(2) = p2 
     end if

     if((p1<=brp%numberOfBoundaryFaces).and.(p2<=brp%numberOfBoundaryFaces)) then  
      ! side is on boundary
      allocate(current%sd%oldIndexes(2),stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
      current%sd%oldIndexes(1) = p1
      current%sd%oldIndexes(2) = p2
     else
      nullify(current%sd%oldIndexes)
     end if

     if(associated(grp%sidePointers(p1)%first)) then
      grp%sidePointers(p1)%last%next=>current 
     else
      grp%sidePointers(p1)%first=>current
     end if
     nullify(current%next)
     grp%sidePointers(p1)%last => current
     grp%numberOfSides = grp%numberOfSides + 1
 
     if(doVisualization) then 
      if(associated(grp%vdp%vsp(p1)%first)) then
       grp%vdp%vsp(p1)%last%next => currentVis
      else
       grp%vdp%vsp(p1)%first => currentVis
      end if
      nullify(currentVis%next)
      grp%vdp%vsp(p1)%last => currentVis
     end if

    end if ! associated(current)

! now calculate the side weights
     xc = 0.5*(x1+x2)
     C12(1) = xm(2)-xc(2)
     C12(2) = xc(1)-xm(1)
     current%sd%sideCoefficients = current%sd%sideCoefficients + 0.5*sign*C12 
    end if
   end do
  end do
 end do ! for each element 

! set face weights

! do i=1,brp%numberOfBoundaryFaces
!  x1 = getCoor(brp%faces(i)%faceIndexes(1),crp)
!  x2 = getCoor(brp%faces(i)%faceIndexes(2),crp) 
!  ! insert face data
!  brp%faces(i)%faceCoefficients(1) = 0.25*(x2(2)-x1(2))
!  brp%faces(i)%faceCoefficients(2) = 0.25*(x1(1)-x2(1))
! end do

 if(doVisualization) then 
  grp%vdp%numberOfVisualizationSides = grp%numberOfSides
 end if 
 
! do i=1,grp%numberOfNodes
!  if(grp%controlVolumeNodeSize(i)>0.001) write(*,*) i,": ",grp%controlVolumeNodeSize(i)
! end do

 end subroutine makeCalcDataForRect
!-------------------------------------------------------------------------
 subroutine makeCalcDataForTri(crp,tep,grp,brp,doVisualization,nrel)
! makes side coefficients from a given triangulation
! and splits data stucture into side-based

 IMPLICIT NONE
 
 type(CoordinateRegisterData) :: crp ! contains coordinate register for grid
 type(TriangularElementData) :: tep ! contains element register for triangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization
 integer :: nrel ! number of rectangular elements allready created

 type(LinkedSide),pointer :: current,iter
 type(VisualizationSidePointer),pointer :: currentVis

 type(STEData),pointer :: temp  ! single triangular element pointer

 integer, parameter :: NSD = 2
 integer, parameter :: NODES_IN_ELEMENT = 3
 
 integer ::  nel ! number of elements in TEP
 integer :: i,j,k,l,allocateStatus,total,p1,p2,p3,sp1,sp2,ptemp,elmind1,elmind2 
 integer :: cj,ck,cl
 real :: xm(NSD) ! to contain element median
 real :: xs(NSD) ! to contain side median
 real :: xc(NSD) ! to contain side coefficient
 real :: xcs(NSD)  
 real :: xf(NSD)
 real :: x1(NSD),x2(NSD), x3(NSD) ! contains side edge coordinates
 real :: C12(NSD),D12(NSD)
 logical :: isAntiClockwise
 integer :: in1,in2,count1,count2
 real :: elementArea,sign


! Do not know how many sides there are in grid, therefore incorporate
! a linked list to dynamically create new sides as we go along. 

  nel = tep%numberOfElements

  do i = 1,nel ! do for each element
  temp => getTriElement(i,tep) 
! find centroid in element 
  xm = 0. 
  do j = 1,NODES_IN_ELEMENT
   xm = xm + getCoor(temp%pointIndexes(j),crp)
  end do 
  xm = xm/NODES_IN_ELEMENT
  p1=temp%pointIndexes(1)
  p2=temp%pointIndexes(2)
  p3=temp%pointIndexes(3)
  x1 = getCoor(p1,crp)
  x2 = getCoor(p2,crp)
  x3 = getCoor(p3,crp)
  elementArea = 0.5*abs((x2(1)-x1(1))*(x3(2)-x1(2))-(x2(2)-x1(2))*(x3(1)-x1(1))) 
! run through points in element to create sides 
  do j = 1,NODES_IN_ELEMENT -1
   do k = j+1,NODES_IN_ELEMENT
    if(j==1.and.k==3) then ! get numbering to follow right hand rule
     cj = k
     ck = j 
    else
     cj = j
     ck = k
    end if
    cl = 6 - j - k
    p1=temp%pointIndexes(cj)
    p2=temp%pointIndexes(ck)

   if(p1 > p2) then ! sort indices to always have smallest index first
     ptemp = p1
     p1 = p2
     p2 = ptemp
     ptemp = cj
     cj=ck
     ck=ptemp
     sign = -1.0
    else
     sign = 1.0
    end if 


    current => findSide(p1,p2,grp)
    if(doVisualization) then 
     currentVis => findVisSide(p1,p2,grp) 
    end if

    ! set controlVolumeNodesSize array as areas of CV around node
    l = 6-j-k
    x1 = getCoor(p1,crp)
    x2 = getCoor(p2,crp)
    x3 = getCoor(temp%pointIndexes(cl),crp)
    grp%controlVolumeNodeSize(p1) = grp%controlVolumeNodeSize(p1)+&
           (0.25/3.)*abs((x2(1)+x3(1)-2*x1(1))*(x2(2)-x1(2))-&
                    (x2(2)+x3(2)-2*x1(2))*(x2(1)-x1(1)))


    grp%controlVolumeNodeSize(p2) = grp%controlVolumeNodeSize(p2)+&
           (0.25/3.)*abs((x1(1)+x3(1)-2*x2(1))*(x1(2)-x2(2))-&
                    (x1(2)+x3(2)-2*x2(2))*(x1(1)-x2(1)))

    if(associated(current)) then
!    side already exists, do not create new side but add index in element
!    register if visualization is needed
     if(doVisualization) then 
     ! current%sd%neighbouringElements(2) = i 
      currentVis%vs%neighbouringElements(2) = i + nrel
     end if
     count1 = count1 + 1
    else

!   create new side
 
     allocate(current,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
     allocate(current%sd,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
     call sideConstruct(NSD,current%sd)

! insert data into side

     current%sd%sideIndexes(1) = p1  
     current%sd%sideIndexes(2) = p2 
    
     current%sd%sideCoefficients = 0.0
     allocate(current%sd%sideLength(2),stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
     current%sd%sideLength = x2 - x1 
     nullify(current%sd%invertedSideLength)
     count2 = count2 + 1
! For visualization, the old indices of each side is needed since agglomerated
! sides don't have physical nodes and thus can't directly be plotted
! Boundary sides need their original node since they are connected to 
! a global register (boundaries aren't agglomerated). In addition this
! yields a very fast method to check if side is on boundary since the 
! node numbers in original grid starts with boundary nodes 
 
     if(doVisualization) then ! do visualization
      allocate(currentVis,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
      allocate(currentVis%vs,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
      currentVis%vs%sideIndexes(1) = p1 
      currentVis%vs%sideIndexes(2) = p2 
      currentVis%vs%neighbouringElements(1) = i + nrel
      currentVis%vs%neighbouringElements(2) = 0 
      currentVis%vs%controlVolumes(1) = p1
      currentVis%vs%controlVolumes(2) = p2 
     end if
!     else if((p1<=brp%numberOfBoundaryFaces).and.(p2<=brp%numberOfBoundaryFaces)) then  
     if((p1<=brp%numberOfBoundaryFaces).and.(p2<=brp%numberOfBoundaryFaces)) then  
      ! side is on boundary
      allocate(current%sd%oldIndexes(2),stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
      current%sd%oldIndexes(1) = p1
      current%sd%oldIndexes(2) = p2
     else
      nullify(current%sd%oldIndexes)
     end if
 
     if(associated(grp%sidePointers(p1)%first)) then
      grp%sidePointers(p1)%last%next=>current 
     else
      grp%sidePointers(p1)%first=>current
     end if
     nullify(current%next)
     grp%sidePointers(p1)%last => current
     grp%numberOfSides = grp%numberOfSides + 1

     if(doVisualization) then 
      if(associated(grp%vdp%vsp(p1)%first)) then
       grp%vdp%vsp(p1)%last%next => currentVis
      else
       grp%vdp%vsp(p1)%first => currentVis
      end if
      nullify(currentVis%next)
      grp%vdp%vsp(p1)%last => currentVis
     end if

    end if ! associated(current)

! now calculate the side weights

    xm = (x1+x2+x3)/3.0
    xc = 0.5*(x1+x2)
    C12(1) = xm(2)-xc(2)
    C12(2) = xc(1)-xm(1)
    current%sd%sideCoefficients = current%sd%sideCoefficients + 0.5*sign*C12 
   end do
  end do
 end do ! for each element 

 k = 0
 do i=1,grp%numberOfNodes
  current => grp%sidePointers(i)%first
  do while(associated(current)) 
   k = k + 1 
   current => current%next
  end do
 end do

! set boundary face weights

 do i=1,brp%numberOfBoundaryFaces
  x1 = getCoor(brp%faces(i)%faceIndexes(1),crp)
  x2 = getCoor(brp%faces(i)%faceIndexes(2),crp) 

  ! insert face data

  brp%faces(i)%faceCoefficients(1) = 0.125*(x2(2)-x1(2))
  brp%faces(i)%faceCoefficients(2) = 0.125*(x1(1)-x2(1))
 end do

 if(doVisualization) then 
  grp%vdp%numberOfVisualizationSides = grp%numberOfSides
 end if 

 end subroutine makeCalcDataForTri 
!-------------------------------------------------------------------------
 subroutine setColors(grp,sidePointers,sizeOfSidePointers,numberOfNodes)
! routine to divide the side array into independent sets
! (ie. sets where a node does not occur twice) to 
! enable efficient vectorization

 IMPLICIT NONE

 type(GridData) ::  grp ! contains grid data  
 type(FirstPointer) :: sidePointers(:)
 integer :: sizeOfSidePointers,numberOfNodes
 

 integer,allocatable :: nSidesArray(:) 
 integer :: optNumbOfSidesPerColor,numberOfSides
 integer :: i,j,k,allocateStatus,slack
 logical :: notFinished

 type(LinkedSide),pointer :: currentSide, previousSide, nextSide
 type(ColorData),pointer :: currentColor,searchColor
 type(LinkedInteger),pointer :: swapInteger

 write(*,*) "Coloring..."

 numberOfSides = sizeOfSidePointers

! first find out how many colors are needed - the largest
! amount of sides attached to one node

 allocate(nSidesArray(numberOfSides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 

 nSidesArray = 0 !initialize
 do i=1,numberOfNodes
  currentSide => sidePointers(i)%first
  do while(associated(currentSide))
   nSidesArray(currentSide%sd%sideIndexes(1)) =&
                 nSidesArray(currentSide%sd%sideIndexes(1)) + 1
   nSidesArray(currentSide%sd%sideIndexes(2)) =&
                 nSidesArray(currentSide%sd%sideIndexes(2)) + 1
   currentSide => currentSide%next
  end do
 end do

 grp%optimalNumberOfColors = maxval(nSidesArray)

 deallocate(nSidesArray,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setColors couldn't deallocate"

 allocate(grp%colorsOccupied(numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setColors out of memory" 

 allocate(grp%colors,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setColors out of memory" 

 do i=1,numberOfNodes
  nullify(grp%colorsOccupied(i)%first)
 end do

! start coloring
! a linked list of colors is created, each color has a linked
! list of sides in it 
 do i=1,numberOfNodes !run through all points
  currentSide => sidePointers(i)%first 
  do while(associated(currentSide)) !run through all sides in point
! try to find a place for currentSide in allready existing colors
   currentColor => grp%colors%first
   do while(associated(currentColor).and.associated(currentSide))
    if(sideFitInColor(currentSide,currentColor%colorNumber,grp)) then
! side fits in color   
     currentColor%last%next=>currentSide
     do k=1,2
      j = currentSide%sd%sideIndexes(k)  
      swapInteger => grp%colorsOccupied(j)%first
      allocate(grp%colorsOccupied(j)%first,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: setColors out of memory" 
      grp%colorsOccupied(j)%first%col = currentColor%colorNumber
      grp%colorsOccupied(j)%first%next => swapInteger 
     end do   
 
     currentColor%last=>currentSide
     currentSide=>currentSide%next
     currentColor%numberOfSides = currentColor%numberOfSides+1
     nullify(currentColor%last%next)
     currentColor => grp%colors%first ! makes the loop start from first color again
    else   
! switch color
     currentColor => currentColor%next 
    end if 
   end do ! while associated(currentColor) and associated(currentSide)
   if(associated(currentSide)) then 
! have run through all colors, but current side hasn't been placed
! this means that the current side couldn't be fitted into one of the 
! existing colors, must then create new color

    if(associated(grp%colors%first)) then 
     allocate(grp%colors%last%next,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: setColors out of memory" 
     grp%colors%last => grp%colors%last%next
    else
     ! at the beginning, no colors allocated
     allocate(grp%colors%first,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: setColors out of memory" 
     grp%colors%last => grp%colors%first 
     grp%colors%numberOfColors = 0
    end if
    
    nullify(grp%colors%last%next)

     do k=1,2
      j = currentSide%sd%sideIndexes(k)  
      swapInteger => grp%colorsOccupied(j)%first
      allocate(grp%colorsOccupied(j)%first,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: setColors out of memory" 
      grp%colorsOccupied(j)%first%col = grp%colors%numberOfColors+1
      grp%colorsOccupied(j)%first%next => swapInteger 
     end do   
! this color is empty so the side must fit
    grp%colors%last%first=>currentSide
    grp%colors%last%last=>currentSide
    grp%colors%last%numberOfSides = 1
    grp%colors%numberOfColors = grp%colors%numberOfColors + 1   
    grp%colors%last%colorNumber = grp%colors%numberOfColors
    currentSide=>currentSide%next
    nullify(grp%colors%last%last%next)
   end if
  end do ! while associated(currentSide)
  nullify(sidePointers(i)%first)
 end do !number of nodes 

! try to get the color sizes similar 

 optNumbOfSidesPerColor = numberOfSides/grp%optimalNumberOfColors
 currentColor => grp%colors%first
 do while(associated(currentColor)) ! run through all colors
  currentSide => currentColor%first
  nullify(previousSide)
  nextSide => currentSide%next
  if(currentColor%colorNumber == 1) then ! try to avoid this 
   slack = numberOfSides-optNumbOfSidesPerColor*grp%optimalNumberOfColors
  else
   slack = 0
  end if
  do while ((currentColor%numberOfSides>(optNumbOfSidesPerColor+slack)).and.&
            (associated(currentSide)))
   searchColor => grp%colors%first 
   do while (associated(searchColor)) 
    if (searchColor%numberOfSides<optNumbOfSidesPerColor) then
     if(sideFitInColor(currentSide,searchColor%colorNumber,grp)) then
     ! transfer side from currentColor to searchColor
      searchColor%last%next => currentSide
      searchColor%last => currentSide

      if(associated(previousSide)) then 
       previousSide%next=>currentSide%next ! jump over currentSide in 
                                           ! currentColor
      else ! must mean that currentSide is first in color
       currentColor%first=>currentSide%next
      end if
      nullify(currentSide%next) 
      currentColor%numberOfSides = currentColor%numberOfSides-1
      searchColor%numberOfSides = searchColor%numberOfSides+1   
     ! fix colorsOccupied register
      do k=1,2
       j = currentSide%sd%sideIndexes(k)  
       swapInteger=>grp%colorsOccupied(j)%first
       do 
        if(.not.associated(swapInteger)) then
         STOP "Something wrong with color register"
        end if
        if (swapInteger%col == currentColor%colorNumber) then
         ! swap sides
         swapInteger%col = searchColor%colorNumber
         exit
        end if
        swapInteger=>swapInteger%next 
       end do
      end do ! k=1,2 
      currentSide=>previousSide ! this will lead to that the previous side
                                ! of next iteration will equal the one in this  
      exit
     end if ! side fit in color
    end if ! searchColor%numberOfSides<optNumbOfSidesPerColor  
    searchColor => searchColor%next
   end do ! while associated(searchColor)
   previousSide => currentSide
   currentSide => nextSide 
   if(associated(currentSide)) then 
    nextSide => currentSide%next
   end if
  end do
   if(.not.associated(currentSide)) then
!    STOP "Could not simplify color register"
   ! much smarter to avoid this problem by listing colors 
   ! and color sizes for solver to read. 
   end if
  currentColor => currentColor%next
 end do ! while associated(currentColor)

! finally place sides, ordered by their color, in a single array

 allocate(grp%srp(numberOfSides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 
 
 i = 1
 currentColor=>grp%colors%first
 do while(associated(currentColor))
  currentSide=>currentColor%first
  do while(associated(currentSide))
   grp%srp(i)%sd => currentSide%sd
   previousSide => currentSide
   currentSide => currentSide%next
!   deallocate(previousSide,stat=allocateStatus)
!   if (allocateStatus /= 0) STOP "ERROR: setColor couldn't deallocate"
   i = i+1 
  end do
!  searchColor => currentColor
  currentColor => currentColor%next
!  deallocate(searchColor,stat=allocateStatus)
!  if (allocateStatus /= 0) STOP "ERROR: CentroidCVgenerator couldn't deallocate"
 end do 

 write(*,*) 'Number of colors: ',grp%colors%numberOfColors

 end subroutine setColors 
!-------------------------------------------------------------------------
subroutine prepareForAgglomeration(grp)
! sets up grid correctly if sides are not colored

IMPLICIT NONE

type(GridData) :: grp

type(LinkedSide),pointer :: currentSide
type(firstPointer),pointer :: sidePointers(:)
integer :: i,j,allocateStatus

if (associated(grp%superSidePointers)) then 
 ! the fine grid is an agglomerated one
 sidePointers => grp%superSidePointers
else
 ! the fine grid is the original grid
 sidePointers => grp%sidePointers
end if

allocate(grp%srp(grp%numberOfSides),stat=allocateStatus)
if (allocateStatus /= 0) STOP "ERROR: CentroidCVGenerator out of memory" 
 
 j = 1
 do i=1,grp%numberOfNodes 
  currentSide=>sidePointers(i)%first
  do while(associated(currentSide))
   grp%srp(j)%sd => currentSide%sd
   currentSide => currentSide%next
   j = j+1 
  end do
  nullify(sidePointers(i)%first)
 end do 
end subroutine prepareForAgglomeration
!-------------------------------------------------------------------------
 logical function sideFitInColor(side,colorNumber,grp)
! Find out whether side is allowable in color or not
! Side is only allowable in color if it has no nodes
! in common with the sides allready in color and
! if the color is not full. The color is full if 
! the number of sides in it exceeds 
! sidesInColors + numberOfColors
! use SideRegister

 IMPLICIT NONE

 type(LinkedSide) :: side
 integer :: colorNumber
 type(GridData) :: grp

 type(LinkedInteger),pointer :: current
 integer :: k

  sideFitInColor = .true.
  do k=1,2
   current => grp%colorsOccupied(side%sd%sideIndexes(k))%first 
   do while(associated(current))
    if (current%col == colorNumber) then ! does not fit 
     sideFitInColor = .false.
     goto 55
    end if 
    current => current%next
   end do
55  end do 
!   end do
 end function sideFitInColor
!-------------------------------------------------------------------------
 function findSide(p1,p2,grp) result(c)
! looks for side in linked list
! use SideRegister

 IMPLICIT NONE

 type(GridData) :: grp
 integer :: p1,p2,numberOfSides
 type(LinkedSide),pointer :: first
 type(LinkedSide),pointer :: c

 integer :: i,ip
 type(LinkedSide),pointer :: current

 nullify(c)

 if(associated(grp%sidePointers(p1)%first)) then
 current => grp%sidePointers(p1)%first
 
 do while(associated(current)) 
  if(current%sd%sideIndexes(2) == p2) then
    ! side is found   
    c => current 
    exit 
  end if 
  current => current%next
 end do
 end if 
 
end function findSide
!-------------------------------------------------------------------------
  function findBoundarySide(p1,p2,grp) result(c)
! looks for side in linked list
! use SideRegister

 IMPLICIT NONE

 type(GridData) :: grp
 integer :: p1,p2,numberOfSides
 type(LinkedSide),pointer :: first
 type(LinkedSide),pointer :: c

 integer :: i,ip
 type(LinkedSide),pointer :: current

 nullify(c)

 if(associated(grp%sidePointers(p1)%first)) then
 current => grp%sidePointers(p1)%first
 
 do while(associated(current)) 
  if(current%sd%sideIndexes(2) == p2) then
    ! side is found   
    c => current 
    exit 
  end if 
  current => current%next
 end do
 end if 
 
end function findBoundarySide
!-------------------------------------------------------------

 function findVisSide(p1,p2,grp) result(c)

 IMPLICIT NONE

 type(GridData) :: grp
 integer :: p1,p2
 type(VisualizationSidePointer),pointer :: first
 type(VisualizationSidePointer),pointer :: c

 integer :: i,ip
 type(VisualizationSidePointer),pointer :: currentVis

 nullify(c)

 if(associated(grp%vdp%vsp(p1)%first)) then
 currentVis => grp%vdp%vsp(p1)%first
 
 do while(associated(currentVis)) 
  if(currentVis%vs%sideIndexes(2) == p2) then
    ! side is found   
    c => currentVis 
    exit 
  end if 
  currentVis => currentVis%next
 end do
 end if 
 
end function findVisSide
!-------------------------------------------------------------------------
 subroutine completeSideRegister(grp,sidePointers,sizeOfSidePointers)
! Prepares the fine grid for agglomeration by extending the 
! sidePointer register such that each node has a complete listing
! of the sides connected to it, not just the sides where the
! lowest node equals the current point index

 IMPLICIT NONE
 
 type(GridData) :: grp 
 type(FirstPointer) :: sidePointers(:)
 integer :: sizeOfSidePointers

 type(LinkedSide),pointer :: currentSide
 integer :: i,ind1,ind2,allocateStatus,k 
 integer :: ind(2)
 integer :: numberOfSides  

 numberOfSides = sizeOfSidePointers
  

  do i=1,numberOfSides
   ind = grp%srp(i)%sd%sideIndexes
   do k=1,2
    if(associated(sidePointers(ind(k))%first)) then 
     allocate(sidePointers(ind(k))%last%next,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: completeSideRegister out of memory" 
     sidePointers(ind(k))%last => sidePointers(ind(k))%last%next
     nullify(sidePointers(ind(k))%last%next)
     sidePointers(ind(k))%last%sd => grp%srp(i)%sd
    else
     allocate(sidePointers(ind(k))%first,stat=allocateStatus)
     if (allocateStatus /= 0) STOP "ERROR: completeSideRegister out of memory" 
     sidePointers(ind(k))%last => sidePointers(ind(k))%first
     nullify(sidePointers(ind(k))%last%next)
     sidePointers(ind(k))%last%sd => grp%srp(i)%sd
    end if 
   end do
  end do
 end subroutine completeSideRegister
!-------------------------------------------------------------------------
 subroutine agglomerateGrid(grp,newgrp,crp,doVisualization,mergeLonesomeNodes,isInitial)
! agglomerates given grid

  IMPLICIT NONE
 
  type(GridData) :: grp,newgrp 
  type(CoordinateRegisterData) :: crp
  logical :: doVisualization,mergeLonesomeNodes,isInitial

  integer :: i,j,k,l,m,p,n,s,ind,ind2,sind1,sind2,minOldInd,ind1,minInd
  integer :: allocateStatus,numberOfNewNodes,numberOfBoundaryFaces
  integer :: numberOfLonesomeNodes, sideNumber, maxInd
  logical,pointer :: nodesAgglomerated(:)
  type(LinkedSide),pointer :: currentSide,nextSide,searchSide,previousSide 
  type(LinkedInteger),pointer :: currentNode,searchNode,swapNode
  type(VisualizationSidePointer),pointer :: currentVisSide,newVisSide

  type(LinkedIntegerP),pointer :: newNodes(:) ! generalized nodes created by 
                                              ! agglomeration process

  integer,pointer :: newControlVolumeMappings(:)
  type(firstPointer),pointer :: sidePointers(:)
  type(LinkedDoubleInteger),pointer :: currentDoubleIndex
  type(LinkedDoubleIntegerP),pointer :: donutList

  real :: x1(2),x2(2),xm(2),dx(2),d12,averageCoefficient,currentSideCoefficient
  real :: localAverageCoefficient,maximumCoefficient,minimumCoefficient
  integer :: numberOfSidesAttachedToNode,numberOfDonuts
  integer :: indx1,indx2,indx,swapIndex,numberOfNewBoundaryCVs
  logical :: breakFree

  numberOfBoundaryFaces = grp%brp%numberOfBoundaryFaces
 
! The simplest scheme is the one used by Lallemand et al. where
! the seed nodes are taken randomly
  
  allocate(nodesAgglomerated(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

  allocate(newNodes(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

! initialize
  do i=1,grp%numberOfNodes
   nodesAgglomerated(i) = .false. ! this indicates that the node hasn't been agglomerated
   nullify(newNodes(i)%first)
  end do 

  nullify(donutList)  
  numberOfNewNodes = 0
  newgrp%numberOfCVSides = 0

  if(associated(grp%superSidePointers)) then 
   ! the fine grid is an agglomerated one
   sidePointers => grp%superSidePointers
  else
   ! the fine grid is the original grid
   sidePointers => grp%sidePointers
  end if

! Agglomerate the nodes. Seed nodes are selected and all neighbouring nodes
! (i.e connected to the seed node by a side) that have not yet been 
! agglomerated and satisfy the directionality condition
! are agglomerated into a new control volume. Lonesome 
! seed nodes are merged into one of its neighbours if wanted 

 
  do i=1,grp%numberOfNodes
   if(.not.nodesAgglomerated(i)) then ! checks if node is already agglomerated 
   ! agglomerate node
    nodesAgglomerated(i) = .true. 
    numberOfNewNodes = numberOfNewNodes + 1

   ! run through associated sides to mark the other nodes  
    allocate(newNodes(i)%first,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
    nullify(newNodes(i)%first%next) 
    newNodes(i)%first%col = i

    ! first run through sides attached to node to calculate average side coefficient
    maximumCoefficient = -1.0
    minimumCoefficient = 1.0E16
    numberOfSidesAttachedToNode = 0
    averageCoefficient = 0.0  
    currentSide => sidePointers(i)%first
    do while(associated(currentSide)) 
     ! use side coefficients as measure
!     averageCoefficient = averageCoefficient +& 
!              sqrt(sum(currentSide%sd%sideCoefficients*currentSide%sd%sideCoefficients))
     currentSideCoefficient = &
              sqrt(sum(currentSide%sd%sideCoefficients*currentSide%sd%sideCoefficients))
     averageCoefficient = averageCoefficient + currentSideCoefficient
     maximumCoefficient = max(maximumCoefficient,currentSideCoefficient)
     minimumCoefficient = min(minimumCoefficient,currentSideCoefficient)
     numberOfSidesAttachedToNode = numberOfSidesAttachedToNode + 1 
     currentSide => currentSide%next
    end do
    averageCoefficient = averageCoefficient/dble(numberOfSidesAttachedToNode)
!    averageCoefficient = averageCoefficient
    
    currentSide => sidePointers(i)%first
    currentNode => newNodes(i)%first
    do while(associated(currentSide))
     currentSideCoefficient = sqrt(sum(currentSide%sd%sideCoefficients*&
                                       currentSide%sd%sideCoefficients))

     localAverageCoefficient = (averageCoefficient*numberOfSidesAttachedToNode-&
                               currentSideCoefficient)/(numberOfSidesAttachedToNode-1.0)
!     if(minimumCoefficient/maximumCoefficient<grp%minimumAspectRatio) then 
!      localAverageCoefficient = maximumCoefficient
!     else
!      ! control volume is not streched enough to use directional agglomeration
!      localAverageCoefficient = 0.0
!     end if

     ind = currentSide%sd%sideIndexes(1)
     if (ind == i) then
      ind = currentSide%sd%sideIndexes(2) 
     end if
     ! ind now points to other node
     if ((.not.nodesAgglomerated(ind)).and.&
!         currentSideCoefficient > grp%directionalityParameter*localAverageCoefficient) then  
         currentSideCoefficient > grp%directionalityParameter*averageCoefficient) then  
      allocate(currentNode%next,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
      currentNode => currentNode%next
      nullify(currentNode%next)
      currentNode%col = ind 
      nodesAgglomerated(ind) = .true. 
     end if
     currentSide => currentSide%next
    end do ! while associated(currentSide)

    if((.not.associated(newNodes(i)%first%next)).and.&
!      ((mergeLonesomeNodes).or.(localAverageCoefficient==0.0))) then 
                                                (mergeLonesomeNodes)) then 
    ! This means that new control volume only has one node, 
    ! the node is lonesome. 
    ! This control volume is therefore merged with one of its
    ! neighbours (to make it feel better)
     breakFree = .false.
     searchSide => sidePointers(i)%first
     do while(associated(searchSide).and.(.not.breakFree))
      if(searchSide%sd%sideIndexes(1)==i) then 
       ind = searchSide%sd%sideIndexes(2)
      else if(searchSide%sd%sideIndexes(2)==i) then 
       ind = searchSide%sd%sideIndexes(1)
      else
       STOP "ERROR: Something's wrong in grid agglomeration - 1"
      end if
      ! Search through sides connected to ind to find seed 
      currentSide => sidePointers(ind)%first 
      do while(associated(currentSide).and.(.not.breakFree))
       ind2 = currentSide%sd%sideIndexes(1)
       if (ind2 == ind) then
        ind2 = currentSide%sd%sideIndexes(2) 
       end if   
       if(i/=ind2.and.associated(newNodes(ind2)%first)) then 
        ! yes, we have found a seed node connected to ind
        allocate(currentNode,stat=allocateStatus)
        if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
        ! place new node after seed node (i.e. as number two in list)
        currentNode%next=>newNodes(ind2)%first%next
        newNodes(ind2)%first%next => currentNode
        currentNode%col = i
        numberOfNewNodes = numberOfNewNodes - 1
        nullify(newNodes(i)%first)

!         newNodes(i)%first%next => newNodes(ind2)%first%next
!         newNodes(ind2)%first%next => newNodes(i)%first
!         nullify(newNodes(i)%first)
        !exit
        breakFree = .true.
       end if
       currentSide => currentSide%next
      end do
      searchSide => searchSide%next
     end do
    end if ! .not.associated(newNodes...  
   end if
  end do


! make newControlVolumeMappings array
 allocate(newControlVolumeMappings(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
 k = 0 
 do i=1,grp%numberOfNodes
  currentNode => newNodes(i)%first
  if(associated(currentNode)) then 
   k = k+1 
  end if 
  do while(associated(currentNode)) 
   newControlVolumeMappings(currentNode%col) = k 
   currentNode => currentNode%next
  end do 
 end do 
! search for surrounded CV's (donuts)
 k = 0
 numberOfDonuts = 0
 do i=1,numberOfNewNodes
  k = k+1
  currentNode => newNodes(k)%first
  do while(.not.associated(currentNode))
   k = k+1
   currentNode => newNodes(k)%first
  end do  
  ! am I a donut?
  breakFree = .false.
  ind1 = newControlVolumeMappings(currentNode%col)
  ind2 = -1  
  searchNode => currentNode
  do while(associated(searchNode).and.(.not.breakFree))
   currentSide => grp%sidePointers(searchNode%col)%first
   do while(associated(currentSide).and.(.not.breakFree)) 
    indx1 = newControlVolumeMappings(currentSide%sd%sideIndexes(1))
    indx2 = newControlVolumeMappings(currentSide%sd%sideIndexes(2))
    if(indx1 /= indx2) then 
    ! side is not internal
     if(indx1==ind1) then 
      indx = indx2
     else if(indx2==ind1) then 
      indx = indx1
     else
      STOP "ERROR: Something's wrong in agglomerateGrid - 6"
     end if
     ! indx now points to outer side  
     if((ind2 > 0).and.(ind2/=indx)) then 
      breakFree = .true.
     else  
      ind2 = indx
     end if 
    else if(indx1/=ind1) then 
     STOP "ERROR: Something's wrong in agglomerateGrid - 7"
    end if      

    currentSide => currentSide%next
   end do
   searchNode => searchNode%next
  end do 

  if(ind2<0) then 
   write(*,*) "WARNING: It appears that the coarse grid has only one node"
  end if

  if(.not.breakFree) then 
  ! yes, I am a donut - put me in the donut chain
   if(.not.associated(donutList)) then 
    allocate(donutList,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
   end if

   allocate(currentDoubleIndex,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
   currentDoubleIndex%index(1) = k 
   currentDoubleIndex%index(2) = ind2  
   nullify(currentDoubleIndex%next)
   if(associated(donutList%first)) then 
    donutList%last%next => currentDoubleIndex
    donutList%last => currentDoubleIndex
   else
    donutList%first => currentDoubleIndex 
    donutList%last => currentDoubleIndex
   end if  
   numberOfDonuts = numberOfDonuts + 1
  end if
 end do

! now lets remove the donuts
 if(numberOfDonuts>0) then 
  if(numberOfNewNodes-numberOfDonuts<=1) then
   STOP "To many grids to be generated, coarse grid has only one node!"
  end if

  k = 0
  l = 0
   do 
   k = k+1
   if(k>grp%numberOfNodes) STOP "ERROR: Something's wrong in agglomerateGrid - 9"
   currentNode => newNodes(k)%first
   do while(.not.associated(currentNode))
    k = k+1
    if(k>grp%numberOfNodes) STOP "ERROR: Something's wrong in agglomerateGrid - 10"
    currentNode => newNodes(k)%first
   end do  
   ! is node k connected to a donut? 
   currentDoubleIndex => donutList%first
   do while(associated(currentDoubleIndex)) 
    if(currentDoubleIndex%index(2)==newControlVolumeMappings(k)) then     
    ! yes it is - engulf it 
    ! keep seed node as first in chain
     swapNode => newNodes(currentDoubleIndex%index(1))%first
     do while(associated(swapNode))
      searchNode => swapNode%next
      swapNode%next => newNodes(k)%first%next
      newNodes(k)%first%next => swapNode 
      swapNode => searchNode
     end do 
     nullify(newNodes(currentDoubleIndex%index(1))%first)
     l = l+1
    end if
    currentDoubleIndex => currentDoubleIndex%next
   end do
   if(l==numberOfDonuts) exit
  end do 

  write(*,*) numberOfDonuts,"surrounded control volumes were removed"
  numberOfNewNodes = numberOfNewNodes - numberOfDonuts

 ! make newControlVolumeMappings array again
  deallocate(newControlVolumeMappings,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid couldn't deallocate"
  allocate(newControlVolumeMappings(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
  k = 0 
  do i=1,grp%numberOfNodes
   currentNode => newNodes(i)%first
   if(associated(currentNode)) then 
    k = k+1 
   end if 
   do while(associated(currentNode)) 
    newControlVolumeMappings(currentNode%col) = k 
    currentNode => currentNode%next
   end do 
  end do 
 end if ! (numberOfDonuts>0) ...

! We now have an array, newNodes, where only the indexes of the seed nodes
! point to something. These indexes contain a chain of the nodes agglomerated
! to it. The task now is to find common sides of these control volumes,
! create them and place them into a list of control volumes of size
! numberOfNewNodes  

 if(numberOfNewNodes <= 1) then 
  STOP "To many grids to be generated, coarse grid has only one node!"
 end if

! map boundary control volumes (i.e. volumes that include boundary nodes) to
! beginning of array

nodesAgglomerated = .false. ! now used to indicate which volumes are on boundary 
numberOfNewBoundaryCVs = 0
do i=1,grp%numberOfNodes
 currentNode => newNodes(i)%first
 do while(associated(currentNode)) 
  if(currentNode%col<=grp%numberOfBoundaryCVs) then 
   nodesAgglomerated = .true.
   numberOfNewBoundaryCVs = numberOfNewBoundaryCVs + 1
   exit
  end if 
  currentNode => currentNode%next
 end do 
end do

! now swap chains to get boundary first

if(numberOfNewNodes>numberOfNewBoundaryCVs) then 
 k = 0
 i = 0 
 do while(k<numberOfNewBoundaryCVs) 
  i = i + 1
  currentNode => newNodes(i)%first
  if(associated(currentNode)) then  
   k = k + 1
   if(.not.nodesAgglomerated(i)) then ! swapping candidate 
    ! now find someone to swap with
    swapIndex = 0
    do j = i+1,grp%numberOfNodes
     if(nodesAgglomerated(j)) then 
      swapIndex = j 
      write(*,*) "swapped"
      nodesAgglomerated(j) = .false.
      exit
     end if
    end do
    if(swapIndex>0) then 
     swapNode => newNodes(j)%first
     newNodes(i)%first => swapNode 
     newNodes(j)%first => currentNode
    else
     STOP "ERROR: Somethings wrong in agglomerateGrid - 10"
    end if
   end if
  end if
 end do
else
 write(*,*) "WARNING: New grid consists only of boundary control volumes"
end if

newgrp%numberOfBoundaryCVs = numberOfNewBoundaryCVs

deallocate(nodesAgglomerated,stat=allocateStatus)
if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid couldn't deallocate"

write(*,*) "Number of nodes in coarse grid: ",numberOfNewNodes
write(*,*) "Number of boundary nodes in coarse grid: ",numberOfNewBoundaryCVs

newgrp%numberOfNodes = numberOfNewNodes
newgrp%numberOfSides = 0

! do visualization

if(doVisualization) then 
 allocate(newgrp%vdp,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
 allocate(newgrp%vdp%vsp(numberOfNewNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

 do i=1,numberOfNewNodes
  nullify(newgrp%vdp%vsp(i)%first)
  nullify(newgrp%vdp%vsp(i)%last)
 end do 
 
 k = 0
 do i=1,grp%numberOfNodes
  currentVisSide => grp%vdp%vsp(i)%first 
  do while(associated(currentVisSide)) 
   ind1 = currentVisSide%vs%controlVolumes(1)
   ind2 = currentVisSide%vs%controlVolumes(2)  
   ind1 = newControlVolumeMappings(ind1)
   ind2 = newControlVolumeMappings(ind2)
   if(ind1/=ind2) then
   ! keep side
    k = k+1
    allocate(newVisSide,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
    allocate(newVisSide%vs,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
    newVisSide%vs%sideIndexes = currentVisSide%vs%sideIndexes
    newVisSide%vs%neighbouringElements = currentVisSide%vs%neighbouringElements
    newVisSide%vs%controlVolumes(1) = ind1
    newVisSide%vs%controlVolumes(2) = ind2
    minInd = min(ind1,ind2)
    nullify(newVisSide%next)
    if(associated(newgrp%vdp%vsp(minInd)%first)) then 
     newgrp%vdp%vsp(minInd)%last%next => newVisSide
    else
     newgrp%vdp%vsp(minInd)%first => newVisSide
    end if
    newgrp%vdp%vsp(minInd)%last => newVisSide
   end if
   currentVisSide => currentVisSide%next
  end do
 end do

 newgrp%vdp%numberOfVisualizationSides = k 
end if
deallocate(newControlVolumeMappings,stat=allocateStatus)
if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid couldn't deallocate"

! For each node in control volume, run through every side
! connected to node to find out which sides are internal 
! (internal means side spans between to nodes of the control volume)
! The internal sides are deleted, while the external sides are kept
! in the sidePointers register of the new grid. 
! The superSidePointers
! register is created by summing up the side coefficients of sides
! bordering to the same pair of control volumes. The indexes
! of a superside is the neighbouring control volume numbers  

 allocate(newgrp%controlVolumeNodePointers(numberOfNewNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

 allocate(newgrp%sidePointers(numberOfNewNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

 allocate(newgrp%boundaryFaceNodeMappings(numberOfBoundaryFaces),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

 newgrp%boundaryFaceNodeMappings = 0 

 do i=1,numberOfNewNodes
  nullify(newgrp%controlVolumeNodePointers(i)%first)
  nullify(newgrp%sidePointers(i)%first)
  nullify(newgrp%sidePointers(i)%last)
 end do

 j = 1 
 k = 0
 l = 0 
 m = 0 
 do i=1,grp%numberOfNodes
  if(associated(newNodes(i)%first)) then
   currentNode => newNodes(i)%first 
   ! first link the grid control volume register to node chain
   newgrp%controlVolumeNodePointers(j)%first => currentNode 
   ! now, lets find out which sides to transfer to new grid, and 
   ! which to deallocate
   do while(associated(currentNode))
    ind = currentNode%col 

    currentSide => sidePointers(ind)%first ! sidePointers(ind) points to all
                                           ! sides connected to node ind 
    do while(associated(currentSide))
     nextSide => currentSide%next 
     ind2 = currentSide%sd%sideIndexes(1)
     if (ind2 == ind) then
      ind2 = currentSide%sd%sideIndexes(2) 
     end if  
     ! we have now made sure that ind and ind2 point to the two
     ! different nodes in currentSide
     if(ind2>0) then 
      ! if ind2 < 0 then the side has allready been found to be external
      ! in another control volume
      searchNode => newNodes(i)%first 
      do while(associated(searchNode))
       if(searchNode%col == ind2) then 
        ! find and delete side pointer in other node
        !searchSide => grp%sidePointers(ind2)%first
        searchSide => sidePointers(ind2)%first
        nullify(previousSide)
        do while(associated(searchSide))
         if(searchSide%sd%sideIndexes(1)==ind.or.searchSide%sd%sideIndexes(2)==ind)then 
        !  deallocate(searchSide%sd,stat=allocateStatus)
        !  if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid couldn't deallocate"
          if(associated(previousSide)) then 
           previousSide%next => searchSide%next
          else
           !grp%sidePointers(ind2)%first => searchSide%next
           sidePointers(ind2)%first => searchSide%next
          end if 
          deallocate(searchSide%sd,stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR:  agglomerateGrid couldn't deallocate"
          deallocate(searchSide,stat=allocateStatus)
          if (allocateStatus /= 0) STOP "ERROR:  agglomerateGrid couldn't deallocate"
          exit
         end if
         previousSide => searchSide
         searchSide => searchSide%next
        end do 
        nullify(currentSide)  
        k = k + 1   
        exit
       end if 
       searchNode => searchNode%next
      end do 
      end if ! ind2 > 0

     ! if currentSide still exists, it is external  
     if(associated(currentSide)) then 
      if(associated(newgrp%sidePointers(j)%first)) then   
       allocate(newgrp%sidePointers(j)%last%next,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR:  agglomerateGrid out of memory"
       newgrp%sidePointers(j)%last => newgrp%sidePointers(j)%last%next
       newgrp%sidePointers(j)%last%sd => currentSide%sd
      else
      ! first side in newgrp%sidePointers(j)
       allocate(newgrp%sidePointers(j)%first,stat=allocateStatus)
       if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory"
       newgrp%sidePointers(j)%last => newgrp%sidePointers(j)%first
       newgrp%sidePointers(j)%last%sd => currentSide%sd
      end if
      nullify(newgrp%sidePointers(j)%last%next)
!      deallocate(currentSide,stat=allocateStatus)
!      if (allocateStatus /= 0) STOP "ERROR: 1 agglomerateGrid couldn't deallocate"
      if(ind2>0) then 
       ! this means that side is untouched
       ! now check whether to invert side weight
       if (ind /= currentSide%sd%sideIndexes(1)) then 
        ! this means that even though side is untouched and
        ! the side weights should be directed from this control volume,
        ! the current side weight is directed the wrong way.  
        ! that is, the node number associated with this control volume
        ! is larger than the other node on the side, i.e. side index
        ! number 2 is attached to this control volume
        newgrp%sidePointers(j)%last%sd%sideCoefficients =&
           -1.0*newgrp%sidePointers(j)%last%sd%sideCoefficients
        newgrp%sidePointers(j)%last%sd%sideLength(1:2) =&
           -1.0*newgrp%sidePointers(j)%last%sd%sideLength(1:2) 
        if(associated(newgrp%sidePointers(j)%last%sd%invertedSideLength)) then 
         newgrp%sidePointers(j)%last%sd%invertedSideLength(1:2) =&
          -1.0*newgrp%sidePointers(j)%last%sd%invertedSideLength(1:2) 
        end if
       end if 
       ! -j is placed in first side index
       ! since sides are always to have smallest index first
       newgrp%sidePointers(j)%last%sd%sideIndexes(1) = -j
       newgrp%numberOfSides = newgrp%numberOfSides + 1
      else
       ! side has been touched before, set index number two
       newgrp%sidePointers(j)%last%sd%sideIndexes(2) = -j
!       do p=1,2 
!        if(currentSide%sd%neighbouringElements(p) > 0) then 
!         newgrp%numberOfCVSides = newgrp%numberOfCVSides + 1
!        end if
!       end do
      end if 
      l = l + 1
     else
     end if ! associated(currentSide...
    currentSide => nextSide 
    end do
    currentNode => currentNode%next 
   end do ! while(associated(currentNode)) 
   j = j + 1
  end if
 end do

!write(*,*) "Sides killed: ",k
!write(*,*) "Sides created: ",l

! set controlVolumeNodeSize and basisFunctionIntegration for
! agglomerated grid

 allocate(newgrp%controlVolumeNodeSize(numberOfNewNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
 allocate(newgrp%coordinates(numberOfNewNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

! allocate(newgrp%basisFunctionIntegration(numberOfNewNodes),stat=allocateStatus)
! if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

 newgrp%controlVolumeNodeSize = 0.0
! newgrp%basisFunctionIntegration = 0.0

 do i=1,numberOfNewNodes
  currentNode => newgrp%controlVolumeNodePointers(i)%first
  newgrp%coordinates(i,:) = 0.0
!  if(isInitial) then 
!   newgrp%coordinates(i,:) = getCoor(currentNode%col,crp)
!  else
!   newgrp%coordinates(i,:) = grp%coordinates(currentNode%col,:)
!  end if
  do while(associated(currentNode)) 
   newgrp%controlVolumeNodeSize(i) = newgrp%controlVolumeNodeSize(i) +&
             grp%controlVolumeNodeSize(currentNode%col)

  if(isInitial) then
   newgrp%coordinates(i,:) = newgrp%coordinates(i,:) + grp%controlVolumeNodeSize(currentNode%col)*getCoor(currentNode%col,crp)
  else
   newgrp%coordinates(i,:) = newgrp%coordinates(i,:) + grp%controlVolumeNodeSize(currentNode%col)*grp%coordinates(currentNode%col,:)
  end if


!   newgrp%basisFunctionIntegration(i) = newgrp%basisFunctionIntegration(i) +&
!             grp%basisFunctionIntegration(currentNode%col)
   currentNode => currentNode%next
  end do
  newgrp%coordinates(i,:) = newgrp%coordinates(i,:)/newgrp%controlVolumeNodeSize(i)
 end do

k = 0
numberOfLonesomeNodes = 0
do i=1,newgrp%numberOfNodes
 currentNode => newgrp%controlVolumeNodePointers(i)%first
 k = 0 
 do while(associated(currentNode))
  k = k + 1
  currentNode => currentNode%next
 end do
! if(k<2) write(*,*) "lonesome node"
 if(k<2) numberOfLonesomeNodes = numberOfLonesomeNodes + 1
end do

if(numberOfLonesomeNodes>0) then 
 write(*,*) "There are",numberOfLonesomeNodes, "lonesome nodes in new grid"
end if

k = 0

! set side indexes positive
do i=1,newgrp%numberOfNodes
 currentSide => newgrp%sidePointers(i)%first
 do while(associated(currentSide))
  k = k + 1
  currentSide%sd%sideIndexes = abs(currentSide%sd%sideIndexes)
  currentSide => currentSide%next
 end do
end do

! fix boundary face node mappings
if(.not.associated(grp%boundaryFaceNodeMappings)) then 
! grp is finest grid
 k = 0 
 do i=1,grp%numberOfNodes
  currentNode => newNodes(i)%first
  if(associated(currentNode)) then 
   k = k + 1
  end if 
  do while(associated(currentNode))
   if(currentNode%col<=grp%brp%numberOfBoundaryFaces) then 
    ! means that current node is on boundary
    newgrp%boundaryFaceNodeMappings(currentNode%col) = k 
   end if
   currentNode => currentNode%next
  end do 
 end do
else
 do i=1,grp%brp%numberOfBoundaryFaces
  k = 0 
  do j=1,grp%numberOfNodes
   currentNode => newNodes(j)%first
   if(associated(currentNode)) then 
    k = k + 1
   end if
   do while(associated(currentNode))
    if(currentNode%col==grp%boundaryFaceNodeMappings(i)) then 
     newgrp%boundaryFaceNodeMappings(i) = k 
     exit
    end if
    currentNode=>currentNode%next
   end do
  end do 
 end do
end if

! fix boundary face indexes
do i=1,newgrp%brp%numberOfBoundaryFaces
! minOldInd = newgrp%brp%faces(i)%oldFaceIndexes(1)
 if(associated(grp%brp%faces(i)%oldFaceIndexes)) then 
  minOldInd = grp%brp%faces(i)%oldFaceIndexes(1)
 else
  minOldInd = grp%brp%faces(i)%faceIndexes(1)
 end if
 newgrp%brp%faces(i)%oldFaceIndexes(1) = minOldInd
 newgrp%brp%faces(i)%faceIndexes(1) = &
    newgrp%boundaryFaceNodeMappings(minOldInd) 
! minOldInd = newgrp%brp%faces(i)%oldFaceIndexes(2)
 if(associated(grp%brp%faces(i)%oldFaceIndexes)) then 
  minOldInd = grp%brp%faces(i)%oldFaceIndexes(2)
 else
  minOldInd = grp%brp%faces(i)%faceIndexes(2)
 end if
 newgrp%brp%faces(i)%oldFaceIndexes(2) = minOldInd
 newgrp%brp%faces(i)%faceIndexes(2) = &
    newgrp%boundaryFaceNodeMappings(minOldInd) 
end do

write(*,*) "Sides in coarse grid: ",newgrp%numberOfSides

! now create super sides

! The superSidePointers register is created by summing up the side coefficients of sides
! bordering to the same pair of control volumes. The indexes
! of a superside is the neighbouring control volume numbers  
 
 newgrp%numberOfSuperSides = 0

 allocate(newgrp%superSidePointers(newgrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 

 do i=1,newgrp%numberOfNodes
  nullify(newgrp%superSidePointers(i)%first)
  nullify(newgrp%superSidePointers(i)%last)
 end do
 

 do i=1,newgrp%numberOfNodes
  currentSide => newgrp%sidePointers(i)%first
  do while(associated(currentSide))
   ind = currentSide%sd%sideIndexes(1)
   if(ind == i) then
    ind = currentSide%sd%sideIndexes(2)
   end if
   if (i < ind) then ! create only sides once (supersides are to be colored)
    nextSide => isSideRegistered(newgrp%superSidePointers(i)%first,ind)
    if(associated(nextSide)) then 
    ! this means that a side between the two control volumes ind,i have been 
    ! created earlier, we just have to add the coefficients to the existing one
     nextSide%sd%sideCoefficients = nextSide%sd%sideCoefficients +&
                                     currentSide%sd%sideCoefficients
     if(associated(currentSide%sd%invertedSideLength)) then  
     ! third or higher grid
      nextSide%sd%sideLength = nextSide%sd%sideLength +&
                                     currentSide%sd%sideLength
      nextSide%sd%invertedSideLength = nextSide%sd%invertedSideLength +&
                                     currentSide%sd%invertedSideLength
     else
     ! second grid
      dx = currentSide%sd%sideLength
      nextSide%sd%sideLength(1:2) = nextSide%sd%sideLength(1:2)+dx
      nextSide%sd%sideLength(3) = nextSide%sd%sideLength(3) +&
                            sqrt(sum(dx*dx)) 
!      nextSide%sd%invertedSideLength(1:2) = nextSide%sd%invertedSideLength(1:2) +&
!           1.0/dx
      nextSide%sd%invertedSideLength(1:2) = nextSide%sd%invertedSideLength(1:2) +&
           dx
      nextSide%sd%invertedSideLength(3) = nextSide%sd%invertedSideLength(3) +&
           1.0/sqrt(sum(dx*dx))
     end if
    else
    ! side between ind and i does not exist yet in superSidePointer register,
    ! must therefore create a new one
     if(associated(newgrp%superSidePointers(i)%first)) then 
      allocate(newgrp%superSidePointers(i)%last%next,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
      newgrp%superSidePointers(i)%last => newgrp%superSidePointers(i)%last%next
      nullify(newgrp%superSidePointers(i)%last%next)
     else
      allocate(newgrp%superSidePointers(i)%first,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
      newgrp%superSidePointers(i)%last => newgrp%superSidePointers(i)%first
      nullify(newgrp%superSidePointers(i)%last%next) 
     end if
     newgrp%superSidePointers(i)%last%sd => currentSide%sd
     if(.not.associated(currentSide%sd%invertedSideLength)) then 
      ! this means that the current grid is the second grid. 
      ! to preserve memory the side length arrays aren't stored for 
      ! the finest grid since it is easy to calculate in this case
      dx = currentSide%sd%sideLength
      ! need to reallocate sideLength since it is only of size 2 in fine grid
      deallocate(currentSide%sd%sideLength,stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid couldn't deallocate" 
      allocate(currentSide%sd%sideLength(3),stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
      allocate(currentSide%sd%invertedSideLength(3),stat=allocateStatus)
      if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
      currentSide%sd%sideLength(1:2) = dx
      currentSide%sd%sideLength(3) = sqrt(sum(dx*dx)) 
      currentSide%sd%invertedSideLength(1:2) = dx
      currentSide%sd%invertedSideLength(3) = 1.0/currentSide%sd%sideLength(3)
     end if
     newgrp%numberOfSuperSides = newgrp%numberOfSuperSides + 1
    end if ! associated(nextSide) 
   end if ! i<ind
   currentSide => currentSide%next
  end do    
 end do

 write(*,*) "Number of supersides: ",newgrp%numberOfSuperSides 

 allocate(newControlVolumeMappings(newgrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateGrid out of memory" 
 newControlVolumeMappings = 0
 do i=1,newgrp%numberOfNodes
  currentSide => newgrp%superSidePointers(i)%first
  do while(associated(currentSide)) 
   newControlVolumeMappings(currentSide%sd%sideIndexes(1)) = &
     newControlVolumeMappings(currentSide%sd%sideIndexes(1)) + 1 
   newControlVolumeMappings(currentSide%sd%sideIndexes(2)) = &
     newControlVolumeMappings(currentSide%sd%sideIndexes(2)) + 1 
   currentSide => currentSide%next
  end do 
 end do
 do i=1,newgrp%numberOfNodes
  if(newControlVolumeMappings(i) <= 1) write(*,*) "CV nr. ",i," has ",newControlVolumeMappings(i), "sides" 
 end do
 
 
 
  ! set up engine inflow boundaries

 allocate(newgrp%brp%engineInletSideIndexes(grp%brp%numberOfEngineInletSides,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateBoundary out of memory"
 allocate(newgrp%brp%engineInletSideCoefficients(grp%brp%numberOfEngineInletSides,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateBoundary out of memory"

 newgrp%brp%numberOfEngineInletSides = 0

 do i=1,grp%brp%numberOfEngineInletSides
  ind1 = newgrp%boundaryFaceNodeMappings(grp%brp%engineInletSideIndexes(i,1))
  ind2 = newgrp%boundaryFaceNodeMappings(grp%brp%engineInletSideIndexes(i,2))
  minInd = min(ind1,ind2)
  maxInd = max(ind1,ind2)

  ! check whether side already exists

  sideNumber = 0
  do j=1,newgrp%brp%numberOfEngineInletSides
   if(newgrp%brp%engineInletSideIndexes(j,1)==minInd.and.newgrp%brp%engineInletSideIndexes(j,2)==maxInd) then
    sideNumber=j
    exit
   end if
  end do

  if(sideNumber>0) then
   newgrp%brp%engineInletSideCoefficients(sideNumber,:) = newgrp%brp%engineInletSideCoefficients(sideNumber,:) + &
    grp%brp%engineInletSideCoefficients(i,:)
  else
   newgrp%brp%numberOfEngineInletSides = newgrp%brp%numberOfEngineInletSides + 1
   sideNumber = newgrp%brp%numberOfEngineInletSides
   newgrp%brp%engineInletSideCoefficients(sideNumber,:) = grp%brp%engineInletSideCoefficients(i,:)
   newgrp%brp%engineInletSideIndexes(sideNumber,1) = minInd
   newgrp%brp%engineInletSideIndexes(sideNumber,2) = maxInd
   newgrp%brp%engineInletSideIndexes(sideNumber,3) = grp%brp%engineInletSideIndexes(i,3)
  end if
 end do

 write(*,*) "Number of engine inlet sides: ",newgrp%brp%numberOfEngineInletSides
 
 
 
 end subroutine agglomerateGrid
!-------------------------------------------------------------------------
 function isSideRegistered(firstSide,ind) result(c)
 ! finds out if side with index ind is found in chain with starting 
 ! side firstSide 
! use SideRegister

 type(LinkedSide),pointer :: c

 type(LinkedSide),pointer :: firstSide  
 integer :: ind
 
 type(LinkedSide),pointer :: currentSide

 nullify(c)
 currentSide => firstSide 
 do while(associated(currentSide))
  if((currentSide%sd%sideIndexes(1)==ind).or.(currentSide%sd%sideIndexes(2)==ind)) then 
   c => currentSide 
   exit
  end if 
  currentSide => currentSide%next 
 end do

 end function isSideRegistered
!-------------------------------------------------------------------------
 subroutine getNodeVolume(crp,rep,tep,grp)
 ! integrates the test function over entire domain (to get lumped mass matrix coefficient)
 ! the result is multiplied with three
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 type(RectangularElementData) :: rep
 type(TriangularElementData) :: tep
 type(GridData) :: grp

 integer :: allocateStatus

! allocate(grp%basisFunctionIntegration(grp%numberOfNodes),stat=allocateStatus)
! if (allocateStatus /= 0) STOP "ERROR: getNodeVolume out of memory"
! grp%basisFunctionIntegration = 0.0 

! call getNodeVolumeForRect(crp,rep,grp)
! call getNodeVolumeForTri(crp,tep,grp)
 end subroutine getNodeVolume
!-------------------------------------------------------------------------
 subroutine getNodeVolumeForRect(crp,rep,grp)
 ! finds area of control volume surrounding node
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 type(RectangularElementData) :: rep
 type(GridData) :: grp

 integer :: allocateStatus,i,j,i1,i2,i3,i4
 real :: x1(2),x2(2),x3(2),x4(2),xm(2),x12(2),x23(2),x34(2),x41(2),r1(2),r2(2),area

 ! run through elements

 do i=1,rep%numberOfElements
  ! calculate area of element
  i1 = rep%elements(i)%pointIndexes(1)
  i2 = rep%elements(i)%pointIndexes(2)
  i3 = rep%elements(i)%pointIndexes(3)
  i4 = rep%elements(i)%pointIndexes(4)
  x1 = crp%points(i1,:) 
  x2 = crp%points(i2,:) 
  x3 = crp%points(i3,:) 
  x4 = crp%points(i4,:) 
  xm = 0.25*(x1+x2+x3+x4)
  x12 = 0.5*(x1+x2)
  x23 = 0.5*(x2+x3)
  x34 = 0.5*(x3+x4)
  x41 = 0.5*(x4+x1)

  ! control volume area

  r1 = xm - x1
  r2 = x41 - x1
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x1
  r2 = x12 - x1
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
! grp%basisFunctionIntegration(i1) =& 
!                       grp%basisFunctionIntegration(i1) + area

  r1 = xm - x2
  r2 = x23 - x2
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x2
  r2 = x12 - x2
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
!  grp%basisFunctionIntegration(i2) =& 
!                       grp%basisFunctionIntegration(i2) + area

  r1 = xm - x3
  r2 = x23 - x3
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x3
  r2 = x34 - x3
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
!  grp%basisFunctionIntegration(i3) =& 
!                       grp%basisFunctionIntegration(i3) + area

  r1 = xm - x4
  r2 = x41 - x4
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x4
  r2 = x34 - x4
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
!  grp%basisFunctionIntegration(i4) =& 
!                       grp%basisFunctionIntegration(i4) + area
 end do 

 end subroutine getNodeVolumeForRect
!-------------------------------------------------------------------------
 subroutine getNodeVolumeForTri(crp,tep,grp)
 ! finds area of control volume surrounding node
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 type(TriangularElementData) :: tep
 type(GridData) :: grp

 integer :: allocateStatus,i,j,i1,i2,i3
 real :: x1(2),x2(2),x3(2),xm(2),x12(2),x23(2),x31(2),r1(2),r2(2),area

 ! run through elements

 do i=1,tep%numberOfElements
  ! calculate area of element
  i1 = tep%elements(i)%pointIndexes(1)
  i2 = tep%elements(i)%pointIndexes(2)
  i3 = tep%elements(i)%pointIndexes(3)
  x1 = crp%points(i1,:) 
  x2 = crp%points(i2,:) 
  x3 = crp%points(i3,:) 
  xm = (x1+x2+x3)/3.
  x12 = 0.5*(x1+x2)
  x23 = 0.5*(x2+x3)
  x31 = 0.5*(x3+x1)

  ! control volume area 

  r1 = xm - x1
  r2 = x31 - x1
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x1
  r2 = x12 - x1
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))

!  grp%basisFunctionIntegration(i1) =& 
!                       grp%basisFunctionIntegration(i1) + area

  r1 = xm - x2
  r2 = x23 - x2
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x2
  r2 = x12 - x2
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
!  grp%basisFunctionIntegration(i2) =& 
!                       grp%basisFunctionIntegration(i2) + area

  r1 = xm - x3
  r2 = x31 - x3
  area = 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
  r1 = xm - x3
  r2 = x23 - x3
  area = area + 0.5*abs(r2(1)*r1(2)-r2(2)*r1(1))
!  grp%basisFunctionIntegration(i3) =& 
!                       grp%basisFunctionIntegration(i3) + area
 end do 

 end subroutine getNodeVolumeForTri
!-------------------------------------------------------------------------
 subroutine findNodeDistanceFromWall(crp,grp)
 ! finds the shortest distance from the nodes to a wall boundary
 ! this is for the moment only done on the finest grid but can 
 ! easily be implemented for agglomerated grid using weihghted averages
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 type(GridData) :: grp

 integer :: i,j,indicator,faceIndex,allocateStatus,nodeOnFace
 real :: minimumDistance,dist1,dist2,dist3,xp(2),x1(2),x2(2),tmin
 logical :: gridIsViscid

 write(*,*) "Finding node distance from wall..."

 ! first check if there are viscid walls in grid
 gridIsViscid = .false.
 do i=1,grp%brp%numberOfBoundaryFaces
  indicator = grp%brp%boundaryNodeIndicators(i)
  if(indicator==1.or.indicator==5.or.indicator==6) then ! boundary face is viscous wall
!  if(indicator==5.or.indicator==6) then ! boundary face is viscous wall
   gridIsViscid = .true. 
   exit
  end if
 end do


 allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
 allocate(grp%wallDistanceFaceIndex(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"

 grp%wallDistance = 0.0
 grp%wallDistanceFaceIndex = 0

 if(gridIsViscid) then 
 do i=1,grp%numberOfNodes
  xp = getCoor(i,crp)
  faceIndex = 0
  if(xp(1)>1.0) then  ! C
   minimumDistance = xp(2)*xp(2) ! distance from wake ! C 
  else ! C
   minimumDistance = 1.0e10
  end if ! C
  
  do j=1,grp%brp%numberOfBoundaryFaces 
   indicator = grp%brp%boundaryNodeIndicators(j)
   if(indicator==1.or.indicator==5.or.indicator==6) then ! boundary face is viscous wall
!   if(indicator==5.or.indicator==6) then ! boundary face is viscous wall
    x1 = getCoor(grp%brp%faces(j)%faceIndexes(1),crp)
    x2 = getCoor(grp%brp%faces(j)%faceIndexes(2),crp)
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
!     faceIndex = j 
     faceIndex = grp%brp%faces(j)%faceIndexes(nodeOnFace)
    end if
   end if
  end do
!  if(faceIndex==0) STOP "ERROR: Something's wrong in findNodeDistanceFromWall"
  grp%wallDistance(i) = sqrt(minimumDistance)
  grp%wallDistanceFaceIndex(i) = faceIndex
 end do
 else
  write(*,*) " There are no viscid walls present" 
 end if 
 end subroutine findNodeDistanceFromWall
!-------------------------------------------------------------------------
subroutine writeVisFile1(crp,grp,tep,brp,OUTFILE)
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 type(GridData) :: grp
 type(TriangularElementData) :: tep
 type(BoundaryRegisterData) :: brp
 integer :: OUTFILE

 integer :: i,nel,ncoor,e1,e2,c1,c2,c3,k,allocateStatus
 logical,pointer :: isElProcessed(:)
 type(VisualizationSidePointer),pointer :: currentVisSide
 write(*,*) "Writing visualization file..."

 ! write header
 write(OUTFILE,*) "1"
 write(OUTFILE,*) "visualization file from Aggl2d v2.1" 
 write(OUTFILE,*) "ne   np   nb" 

 ! count elements 

 allocate(isElProcessed(tep%numberOfElements),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: getNodeVolume out of memory"

 isElProcessed = .false.
 nel = 0
 ncoor = 0
 do i=1,grp%numberOfNodes
  if(associated(grp%vdp%vsp(i)%first)) then 
  currentVisSide => grp%vdp%vsp(i)%first
  do while(associated(currentVisSide)) 
   e1 = currentVisSide%vs%neighbouringElements(1)
   e2 = currentVisSide%vs%neighbouringElements(2)
   if(e1>0) then 
    if(.not.isElProcessed(e1)) then 
     nel = nel + 1 
     isElProcessed(e1) = .true. 
    end if
   end if
   if(e2>0) then 
    if(.not.isElProcessed(e2)) then 
     nel = nel + 1 
     isElProcessed(e2) = .true. 
    end if
   end if
   currentVisSide => currentVisSide%next
  end do 
  end if
 end do
 write(OUTFILE,*) nel,crp%numberOfPoints,grp%brp%numberOfBoundaryFaces

 write(OUTFILE,*) "connectivities"
 if(isElProcessed(2093)) write(*,*) "det var da som ..."
 k = 0
 do i=1,tep%numberOfElements
  if(isElProcessed(i)) then 
   k = k+1
   write(OUTFILE,*) k,tep%elements(i)%pointIndexes
  end if
 end do

 write(OUTFILE,*) "coordinates"
 
 k = 0
 do i=1,crp%numberOfPoints
!  if(isCoorProcessed(i)) then 
!   k = k+1
!   write(OUTFILE,*) k,getCoor(i,crp)
!  end if
  write(OUTFILE,*) i,getCoor(i,crp)
 end do

 write(OUTFILE,*) "unknowns"

 do i=1,crp%numberOfPoints
  write(OUTFILE,'(I7,6E17.7)') i,0.0,0.0,0.0,0.0,0.0,0.0
 end do

 write(OUTFILE,*) "boundaries"

 do i=1,brp%numberOfBoundaryFaces
  write(OUTFILE,*) brp%faces(i)%faceIndexes,0,0,0 
 end do 

 deallocate(isElProcessed,stat=allocateStatus)
 if(allocateStatus>0) STOP "ERROR: could not deallocate in writeVisFile"

end subroutine writeVisFile1
!-------------------------------------------------------------------------
subroutine writeVisFile2(crp,grp,rep,tep,brp,OUTFILE)
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 type(GridData) :: grp
 type(RectangularElementData) :: rep
 type(TriangularElementData) :: tep
 type(BoundaryRegisterData) :: brp
 integer :: OUTFILE

 integer :: i,j,nel,ncoor,nconn,nbs,nbf,e1,e2,c1,c2,c3,k,l,allocateStatus
 integer,pointer :: isElProcessed(:)
 real,pointer :: elementCentroid(:,:)
 real :: x(2)
 integer :: nrep ! number of rectangular elements

 type(VisualizationSidePointer),pointer :: currentVisSide
 write(*,*) "Writing visualization file..."

 ! write header
 write(OUTFILE,*) "1"
 write(OUTFILE,*) "visualization file from Aggl2d v2.1" 
 write(OUTFILE,*) "ne   np   nb" 

 ! count elements 
 allocate(isElProcessed(tep%numberOfElements+rep%numberOfElements),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeVisFile out of memory"

 nrep = rep%numberOfElements
 isElProcessed = 0 
 nel = 0
 nbs = 0
 do i=1,grp%numberOfNodes
  if(associated(grp%vdp%vsp(i)%first)) then 
  currentVisSide => grp%vdp%vsp(i)%first
  do while(associated(currentVisSide)) 
   e1 = currentVisSide%vs%neighbouringElements(1)
   e2 = currentVisSide%vs%neighbouringElements(2)
   if(e1>0) then 
    if(isElProcessed(e1)==0) then 
     nel = nel + 1 
     isElProcessed(e1) = nel 
    end if
   else
    nbs = nbs+1
   end if
   if(e2>0) then 
    if(isElProcessed(e2)==0) then 
     nel = nel + 1 
     isElProcessed(e2) = nel 
    end if
   else
    nbs = nbs+1
   end if
   currentVisSide => currentVisSide%next
  end do 
  end if
 end do

 nbf = brp%numberOfBoundaryFaces
 ncoor = nbf+nel+grp%vdp%numberOfVisualizationSides  
 nconn = 2*grp%vdp%numberOfVisualizationSides - nbs
 write(OUTFILE,*) nconn,ncoor,grp%brp%numberOfBoundaryFaces

 allocate(elementCentroid(nel,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeVisFile out of memory"

 write(OUTFILE,*) "connectivities"

 k = 0
 l = 0
 do i=1,grp%numberOfNodes
  currentVisSide => grp%vdp%vsp(i)%first
  do while(associated(currentVisSide)) 
   l = l+1 
   if(currentVisSide%vs%neighbouringElements(1)>0) then
    k = k+1
    write(OUTFILE,*) k,nbf+nel+l,nbf+isElProcessed(currentVisSide%vs%neighbouringElements(1)),nbf+nel+l
   end if
   if(currentVisSide%vs%neighbouringElements(2)>0) then
    k = k+1
    write(OUTFILE,*) k,nbf+nel+l,nbf+isElProcessed(currentVisSide%vs%neighbouringElements(2)),nbf+nel+l
   end if
   currentVisSide => currentVisSide%next
  end do
 end do
 write(OUTFILE,*) "coordinates"

 ! boundaryPoints

 do i=1,nbf
  write(OUTFILE,*) i,getCoor(i,crp)
 end do

 ! element centroids

 do i=1,rep%numberOfElements
  if(isElProcessed(i)>0) then 
   x = 0.0
   do j=1,4
    x = x + getCoor(rep%elements(i)%pointIndexes(j),crp) 
   end do 
   elementCentroid(isElProcessed(i),:) = x/4.0
  end if
 end do

 do i=1,tep%numberOfElements
  if(isElProcessed(i+nrep)>0) then 
   x = 0.0
   do j=1,3
    x = x + getCoor(tep%elements(i)%pointIndexes(j),crp) 
   end do 
   elementCentroid(isElProcessed(i+nrep),:) = x/3.0
  end if
 end do

 do i=1,nel
  write(OUTFILE,*) nbf+i,elementCentroid(i,:)
 end do

 ! then visualization side midpoints
 l = 0
 do i=1,grp%numberOfNodes
  currentVisSide => grp%vdp%vsp(i)%first
  do while(associated(currentVisSide)) 
   l = l+1 
   x = getCoor(currentVisSide%vs%sideIndexes(1),crp)
   x = x + getCoor(currentVisSide%vs%sideIndexes(2),crp)
   x = 0.5*x
   write(OUTFILE,*) nbf+nel+l,x
   currentVisSide => currentVisSide%next
  end do
 end do

 write(OUTFILE,*) "unknowns"

 do i=1,ncoor
!  write(OUTFILE,*) i,0.0,0.0,0.0,0.0,0.0,0.0
  write(OUTFILE,'(I7,6E17.7)') i,0.0,0.0,0.0,0.0,0.0,0.0
 end do

 write(OUTFILE,*) "boundaries"

 do i=1,brp%numberOfBoundaryFaces
  write(OUTFILE,*) brp%faces(i)%faceIndexes,0,0,0 
 end do 

 deallocate(isElProcessed,stat=allocateStatus)
 if(allocateStatus>0) STOP "ERROR: could not deallocate in writeVisFile"
 deallocate(elementCentroid,stat=allocateStatus)
 if(allocateStatus>0) STOP "ERROR: could not deallocate in writeVisFile"
end subroutine writeVisFile2
!-------------------------------------------------------------------------
 subroutine writeComputationFile(isInitial,crp,finegrp,grp,brp,OUTFILE)
! write file for solver input
!use CoordinateRegister
!use BoundaryRegister

 IMPLICIT NONE

 logical :: isInitial 
 type(GridData) :: finegrp,grp
 type(CoordinateRegisterData) :: crp
 type(BoundaryRegisterData) :: brp
 integer :: OUTFILE 
 integer,pointer :: helpArray(:)

 integer :: i,j,ind1,allocateStatus
 real :: x1(2),x2(2)
 real :: help
 type(LinkedInteger),pointer :: currentNode

 if(isInitial) then
! write face data
  write(OUTFILE) brp%numberOfBoundaryFaces
  call writeBoundaryData(OUTFILE,brp) 

! write sides
  write(OUTFILE) grp%numberOfSides
  do i=1,grp%numberOfSides
   write(OUTFILE) grp%srp(i)%sd%sideIndexes(1:2) 
  end do
  do i=1,grp%numberOfSides
   write(OUTFILE) grp%srp(i)%sd%sideCoefficients(1:2)
  end do
 
  do i=1,grp%numberOfSides
   write(OUTFILE) grp%srp(i)%sd%sideLength(1:2)
  end do

  write(OUTFILE) grp%numberOfNodes   
! write control volume area 
  do i=1,grp%numberOfNodes
!   write(OUTFILE) grp%basisFunctionIntegration(i) 
   write(OUTFILE) grp%controlVolumeNodeSize(i) 
!   write(*,*) grp%basisFunctionIntegration(i)/grp%controlVolumeNodeSize(i) 
!   if(abs(grp%basisFunctionIntegration(i)/grp%controlVolumeNodeSize(i)-1.0)>0.01) then
!    write(*,*) i, grp%basisFunctionIntegration(i),grp%controlVolumeNodeSize(i)
!   else
!    write(*,*) 
!   end if
  end do

  call writeCoordinateData(OUTFILE,crp) 

! write node distance from wall
  do i=1,grp%numberOfNodes
!  write(OUTFILE) grp%wallDistanceFaceIndex(i),grp%wallDistance(i)
  end do

 else
  write(OUTFILE) brp%numberOfBoundaryFaces
  call writeBoundaryData(OUTFILE,brp) 

  do i=1,brp%numberOfBoundaryFaces
   ind1 = brp%faces(i)%oldFaceIndexes(1)
   write(OUTFILE) ind1,grp%boundaryFaceNodeMappings(ind1)
  end do

  write(OUTFILE) grp%numberOfBoundaryCVs
  write(OUTFILE) grp%numberOfSuperSides   
  do i=1,grp%numberOfSuperSides
   write(OUTFILE) grp%srp(i)%sd%sideIndexes(1:2) 
  end do

  do i=1,grp%numberOfSuperSides
   write(OUTFILE) grp%srp(i)%sd%sideCoefficients(:)
  end do

  do i=1,grp%numberOfSuperSides
   write(OUTFILE) grp%srp(i)%sd%sideLength
  end do  
  do i=1,grp%numberOfSuperSides
   write(OUTFILE) grp%srp(i)%sd%invertedSideLength
  end do  

  write(OUTFILE) grp%numberOfNodes   

! write lumped mass matrix 
  do i=1,grp%numberOfNodes
!   write(OUTFILE) grp%basisFunctionIntegration(i) 
   write(OUTFILE) grp%coordinates(i,:),grp%controlVolumeNodeSize(i) 
  end do

! write node distance from wall

 allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
 allocate(grp%wallDistanceFaceIndex(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"

 do i=1,grp%numberOfNodes 
  ! first loop through all nodes in control volume to find total area
  currentNode=>grp%controlVolumeNodePointers(i)%first
  grp%wallDistanceFaceIndex(i) =&
                                  finegrp%wallDistanceFaceIndex(currentNode%col)
  do while(associated(currentNode))
   grp%wallDistance(i) = finegrp%wallDistance(currentNode%col)&
       *finegrp%controlVolumeNodeSize(currentNode%col)/grp%controlVolumeNodeSize(i)
   currentNode=>currentNode%next 
  end do 
 end do

 do i=1,grp%numberOfNodes
  if(grp%wallDistanceFaceIndex(i)>0) then 
!   write(OUTFILE) brp%faces(grp%wallDistanceFaceIndex(i))%faceIndexes(1),&
!                  grp%wallDistance(i)
  else
!   write(OUTFILE) 0,grp%wallDistance(i) 
  end if 
 end do


! write mappings 
  call writeMappingToFile(OUTFILE,crp,finegrp,grp,brp,1) 

 end if

 end subroutine writeComputationFile
!-------------------------------------------------------------------------
 subroutine writeMappingToFile(OUTFILE,crp,finegrp,grp,brp,scheme)
! writes inter-grid transformations to computation file

 IMPLICIT NONE

 integer :: OUTFILE
 type(CoordinateRegisterData) :: crp
 type(GridData) :: finegrp,grp
 type(BoundaryRegisterData) :: brp
 integer :: scheme ! decides which mapping scheme to use

 integer :: i,j
 type(LinkedInteger),pointer :: currentNode

! CONVENTION: 
! scheme = 1: injection for coarse to fine (smoothing can be applied at solver level)
!             linear interpolation for fine to coarse

 if(scheme ==1) then 
! first write association between grids, that is, which nodes in the 
! finer grid are associated with the different supernodes in current grid 
! i.e. write out the newNodes chain array
! this data can be used as is for insertion, or form the basis for some weighted
! mapping 
 j = 0
 do i=1,grp%numberOfNodes
 currentNode=>grp%controlVolumeNodePointers(i)%first 
 if(.not.associated(currentNode)) write(*,*) i
  do while(associated(currentNode)) 
   write(OUTFILE) currentNode%col
   currentNode => currentNode%next
  end do 
 write(OUTFILE) 0
 j = j + 1
 end do  
! fine to coarse

! use linear interpolation

 do i=1,grp%numberOfNodes 
  ! first loop through all nodes in control volume to find total area
  currentNode=>grp%controlVolumeNodePointers(i)%first
  do while(associated(currentNode))
   write(OUTFILE)&
       finegrp%controlVolumeNodeSize(currentNode%col)/grp%controlVolumeNodeSize(i)
   currentNode=>currentNode%next 
  end do 
  write(OUTFILE) -1.0 
 end do
 else
  STOP "ERROR: Unknown mapping scheme: " 
 end if 

 end subroutine writeMappingToFile
!-------------------------------------------------------------------------

end module Grid 

