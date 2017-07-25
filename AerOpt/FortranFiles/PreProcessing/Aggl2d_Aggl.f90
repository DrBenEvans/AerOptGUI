!**********************************************
!*aggl2d_aggl.f90                             *
!*                                            *
!*Main program module                         *
!*                                            * 
!*   aggl2d v 2.0                             * 
!*                                            *
!*                                            *
!*Made by                                     *
!*Kaare A Sorensen,                           *
!*20.08.98-                                   *
!**********************************************
!*Description:                                *
!* This file controls the datastructure       *
!* creation and agglomeration for the pre-    *
!* processor.                                 *
!**********************************************
 
module Agglomerator2d
! main program module, handles everything

! Include the definitions of other modules. 
! These could be changed if other types are 
! wanted
 use BoundaryRegister ! handles boundary faces
 use TriangularElements ! handles elements
 use RectangularElements
 use CoordinateRegister ! provides coordinates for nodes in grid
 use Grid ! handles a grid   
 
 IMPLICIT NONE ! to disallow implicit variable declarations


 integer,parameter:: NSD = 2 ! number of space dimensions

 type(CoordinateRegisterData),pointer :: crp
 type(TriangularElementData),pointer :: tep
 type(RectangularElementData),pointer :: rep
 type(GridData),pointer :: grp(:)

! integer, parameter :: MAX_NAME_LENGTH = 80
 character :: fileName*80 ! Input file. Contains a fine triangular
                          ! grid to be agglomerated and processed 
                          ! for side based computations

 character :: compFile*80 ! Output file to be used for computations
                          ! does not contain any element or control 
                          ! volume data, only side coefficients and
                          ! boundary information

 character :: visFile*80  ! Output file for visualization purposes,
                          ! contains element or control volume 
                          ! information

 integer :: npoint,nelem,nbound  ! number of points elements and boundary sides in 
                                 ! fine grid

 integer :: ngrids ! number of grids to produce (1: just fine grid)
 real :: directionalityParameter,minimumAspectRatio 
 logical :: doVisualization,mergeLonesomeNodes,isHybrid
 integer :: visualizationMode,numberOfDirectionalAggl
 integer :: numberOfRectEl,numberOfTriEl
 contains

!--------------------------------------------------------------------------
  subroutine main()
  ! runs the show, calls other functions
 
  call communicate() ! communicates with user - calls procedures to get and
                     ! set up data from file
  call readFineGrid() ! reads fine grid from file and prepares data structure
  call setUpInitialGrid() ! prepares fine grid data structure for agglomeration
  call agglomerate() ! conducts agglomoration
  write(*,*) "Ha en god dag!" ! "Have a nice day!" in Norwegian
  end subroutine main
!--------------------------------------------------------------------------
  subroutine communicate()

  IMPLICIT NONE
  
  integer :: allocateStatus
  
      
  write(*,*) "" 
  write(*,'(A)') "************************************************"
  write(*,'(A)') "    WELCOME TO AGGL2D - 2D GRID AGGLOMERATOR    "  
  write(*,'(A)') "************************************************"
  write(*,*) "" 
  write(*,'(A)',advance="no") "Enter input filename: "  
  read(*,'(A)') fileName
  write(*,'(A)',advance="no") "Is grid hybrid (.true./.false.)/(T/F): "  
  read(*,*) isHybrid
  write(*,'(A)',advance="no") "Number of grids to generate: "  
  read(*,*) ngrids 
  write(*,'(A)',advance="no") "Directionality parameter: "  
  read(*,*) directionalityParameter 
  if(directionalityParameter>0.0) then 
   write(*,'(A)',advance="no") "Number of directional agglomerations: "  
   read(*,*) numberOfDirectionalAggl 
  end if
!  write(*,'(A)',advance="no") "Minimum aspect ratio: "  
!  read(*,*) minimumAspectRatio 
!  write(*,'(A)',advance="no") "Merge lonesome nodes (.true./.false.): "  
!  read(*,*) mergeLonesomeNodes 
  mergeLonesomeNodes = .true.
  write(*,'(A)',advance="no") "Visualization mode (0: none): "  
  read(*,*) visualizationMode
  doVisualization = visualizationMode==1.or.visualizationMode==2 
  
  allocate(grp(ngrids),stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary register"

  end subroutine communicate  
!--------------------------------------------------------------------------
  subroutine agglomerate()

! Only need to store the current and previous
! grids, the rest can be written to file and deleted 
! from memory. The two first alone will then decide 
! the size limitations for the agglomerator

  integer :: i,visNameLength,compNameLength,allocateStatus
  integer :: VISOUTFILE,COMPOUTFILE  ! visualization and calculation output
                                     ! files 

! set directionality parameter
  do i=1,ngrids
   if(i-1<numberOfDirectionalAggl) then 
    grp(i)%directionalityParameter = directionalityParameter
   else
    grp(i)%directionalityParameter = 0.0 
!    mergeLonesomeNodes = .true.
   end if 
   grp(i)%minimumAspectRatio = minimumAspectRatio
  end do

  COMPOUTFILE = -1 
  write(*,'(A)',advance="no") "Computation output file: "  
  read(*,'(A)') compFile

  compNameLength = nameLen(compFile)
  if(compNameLength > 0) then
   write(*,*) "Writing grid data #",1," to file..."
   COMPOUTFILE = 26

   open(COMPOUTFILE,file=compFile(1:compNameLength),form='unformatted',status='unknown')
   write(COMPOUTFILE) ngrids

   call writeComputationFile(.true.,crp,grp(1),grp(1),grp(1)%brp,COMPOUTFILE)
! this must be done before agglomeration as agglomerator deletes
! sides not included in coarse grid 
  end if

   if(doVisualization) then
    write(*,'(A)',advance="no") "Visualization output file for grid: "
    read(*,'(A)') visFile
    visNameLength = nameLen(visFile)
   else
    visNameLength = 0
   end if
   if(visNameLength > 0) then
    VISOUTFILE = 27
    open(VISOUTFILE,file=visFile(1:visNameLength),form='formatted',status='unknown')
    if(visualizationMode==1) then
     call writeVisFile1(crp,grp(1),tep,grp(1)%brp,VISOUTFILE)
    else if(visualizationMode==2) then
     call writeVisFile2(crp,grp(1),rep,tep,grp(1)%brp,VISOUTFILE)
    end if
    close(VISOUTFILE)
   end if

  VISOUTFILE = -1 
  do i=1,ngrids-1
   write(*,*) ""
   write(*,*) "Grid #",i+1,":" 
   if(doVisualization) then 
    write(*,'(A)',advance="no") "Visualization output file for grid: "  
    read(*,'(A)') visFile
    visNameLength = nameLen(visFile)
   else
    visNameLength = 0
   end if
   if(visNameLength > 0) then
    VISOUTFILE = 27
    open(VISOUTFILE,file=visFile(1:visNameLength),form='formatted',status='unknown')
   end if
   write(*,*) "Agglomerating grid #",i
   allocate(grp(i+1)%brp,stat=allocateStatus)
   if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary register"
   call setUpNewBoundaryData(grp(i)%brp,grp(i+1)%brp)
   call agglomerateGrid(grp(i),grp(i+1),crp,doVisualization,mergeLonesomeNodes,i==1)
!   call colorBoundaryFaces(grp(i+1)%brp)
!   call setColors(grp(i+1),grp(i+1)%superSidePointers,&
!                  grp(i+1)%numberOfSuperSides,grp(i+1)%numberOfNodes)
   call prepareForAgglomeration(grp(i+1))
   call completeSideRegister(grp(i+1),grp(i+1)%superSidePointers,&
                                  grp(i+1)%numberOfSuperSides) 

   if(visNameLength>0) then 
!    call writeVisualizationFile(crp,grp(i+1),tep,grp(i+1)%brp,VISOUTFILE)
    if(visualizationMode==1) then 
     call writeVisFile1(crp,grp(i+1),tep,grp(1)%brp,VISOUTFILE)
    else if(visualizationMode==2) then 
     call writeVisFile2(crp,grp(i+1),rep,tep,grp(1)%brp,VISOUTFILE)
    end if
    close(VISOUTFILE)
   end if

   if(compNameLength>0) then 
    call writeComputationFile(.false.,crp,grp(i),grp(i+1),grp(i+1)%brp,COMPOUTFILE)
   end if
  end do

  if(compNameLength>0) then 
   close(COMPOUTFILE)
  end if

  end subroutine agglomerate
!--------------------------------------------------------------------------
  subroutine setUpInitialGrid()
! This procedure sets up the initial grid for use. 
! The initial grid is given from the grid generator
! as a set of triangular (tetrahedral) elements, 
! this is first transformed into a side-based structure.
! The data structure is modified to accomodate the agglomoration
! techniques applied later. 
!  call setBoundaryFirst(grp(1),crp,rep,tep,grp(1)%brp)
  call setBoundaryFirst2(grp(1),crp,rep,tep,grp(1)%brp)
  call makeCalculationData(crp,rep,tep,grp(1),grp(1)%brp,doVisualization)
  call getNodeVolume(crp,rep,tep,grp(1))
  call setUpBoundary(crp,grp(1)%brp,grp(1),rep,tep)
!  call setColors(grp(1),grp(1)%sidePointers,grp(1)%numberOfSides,grp(1)%numberOfNodes)
  nullify(grp(1)%superSidePointers)
  call prepareForAgglomeration(grp(1))
  call completeSideRegister(grp(1),grp(1)%sidePointers,grp(1)%numberOfSides) 
  call findNodeDistanceFromWall(crp,grp(1))
  end subroutine setUpInitialGrid 
!-------------------------------------------------------------------------- 
  subroutine readFineGrid()
! communicates with user and reads fine grid

  IMPLICIT NONE

  integer, parameter :: INFILE = 25 
  integer :: nameLength,i,ip,allocateStatus,dummy,num

  nameLength = nameLen(fileName)
  open(INFILE,file=fileName(1:nameLength),form='formatted',status='old')
  write(*,*) 'File '//fileName(1:nameLength)//' opened...' 
  write(*,*) "" 

! Start reading file

  rewind(INFILE)
  read(INFILE,*) num
  do i=1,num
   read(INFILE,*) 
  end do 
  read(INFILE,*) 
  if(isHybrid) then 
   read(INFILE,*) numberOfRectEl, numberOfTriEl,npoint, nbound
   nelem = numberOfRectEl + numberOfTriEl 
   write(*,*) ""
   write(*,'(''Quadrilaterals:'',i8)') numberOfRectEl
   write(*,'(''Triangles     :'',i8)') numberOfTriEl
   write(*,'(''Points        :'',i8)') npoint
   write(*,'(''Boundary faces:'',i8)') nbound
   write(*,*) ""
  else
   read(INFILE,*) numberOfTriEl,npoint, nbound
   numberOfRectEl = 0
   nelem = numberOfRectEl + numberOfTriEl
   write(*,*) ""
   write(*,'(''Elements      :'',i8)') nelem
   write(*,'(''Points        :'',i8)') npoint
   write(*,'(''Boundary faces:'',i8)') nbound
   write(*,*) ""
  end if
  read(INFILE,*) 

  
! initiates data structure for input

  grp%numberOfNodes = npoint

  allocate(rep,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create element register"
  allocate(tep,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create element register"


  call constructRectElArray(numberOfRectEl,rep) ! calls the constructor of the current element 
                                                ! module 
  call constructTriElArray(numberOfTriEl,tep)  

  allocate(crp,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create coordinate register"

  call constructCoordinateRegister(npoint,NSD,crp) ! calls the constructor of the point
                                                ! register

  allocate(grp(1)%brp,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary register"
    
  call constructBoundaryRegister(nbound,grp(1)%brp) ! calls the constructor of the current 
                                      ! boundary module

! read element indices 
  write(*,*) "Reading elements..." 
  call readRectElData(INFILE,rep) ! calls element module to read it's own data 
  call readTriElData(INFILE,tep) 
! read coordinates for points 
  read(INFILE,*)
  write(*,*) "Reading points..."
  call readPointData(INFILE,crp) ! calls point register module to read data 

!read unknowns (useless here)
!  read(INFILE,*) 
!  do i=1,npoint
!   read(INFILE,*) 
!  end do


! read boundary face indices and indicators
  write(*,*) "Reading boundary faces..."
  read(INFILE,*)
  call readBoundData(INFILE,grp(1)%brp) ! calls point register module to read data 
  close(INFILE)
  write(*,*) ""

  end subroutine readFineGrid
!--------------------------------------------------------------------------
  integer function nameLen(fn)
  IMPLICIT NONE
 
  character*80 :: fn
  integer :: i

  do i = 80,1,-1
   nameLen = i
   if(fn(i:i)/=' ') GOTO 77 ! EXIT
   end do 
   nameLen = 0
77  end function nameLen
!--------------------------------------------------------------------------
end module Agglomerator2d

