    module Smoothing
    
    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    use FDGD
    
    contains
    
    subroutine SmoothingLinear(CNDisp)
    
        ! Variables
        implicit none
        integer :: ismooth, NoSweeps, NoPmove, NoIterFDGD, isweep, i
        double precision, dimension(:), allocatable :: ddn_before, ddn_after, nx, ny, x, y, xnew, ynew, CNxfinal, CNyfinal, CNxsmooth, CNysmooth, beta
        double precision, dimension(maxDoF) :: CNDisp
        double precision :: conv, initialResidual, res, convergence, smoothfactor, smoothfactormax
        
        ! Body of SubSmoothing
        convergence = -3
        initialResidual = 0.0
        conv = 0.0
        NoPmove = size(MovingGeomIndex, dim = 1)
        NoSweeps = NoPmove
        smoothfactor = 0.005
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xnew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ynew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxfinal(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNyfinal(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxsmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNysmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(nx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ny(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        x = RD%coord(orderedBoundaryIndex,1)
        y = RD%coord(orderedBoundaryIndex,2)
        
        ! CN of smoothed surface
        CNxsmooth = x(CN_ind)
        CNysmooth = y(CN_ind)
        
        ! Final CN
        CNxfinal = x(CN_ind) + CNDisp(1:IV%NoCN)
        CNyfinal = y(CN_ind) + CNDisp((IV%NoCN+1):(2*IV%NoCN))
        
        do while (conv > convergence)            
                        
            ! 0. Evaluate maximum smoothing and adjust parameters
            beta(1) = sqrt((x(2)-x(NoPmove))**2 + (y(2)-y(NoPmove))**2)/2
            do i = 2, NoPmove-1
                beta(i) = sqrt((x(i+1)-x(i-1))**2 + (y(i+1)-y(i-1))**2)/2
            end do
            beta(NoPmove) = sqrt((x(1)-x(NoPmove-1))**2 + (y(1)-y(NoPmove-1))**2)/2
            smoothfactormax = minval(beta)/2
            if (smoothfactor > smoothfactormax) then
                print *, 'The smoothfactor is too big and has been reduced from', smoothfactor, 'to', smoothfactormax
                Nosweeps = floor(NoPmove*(smoothfactor/smoothfactormax))
                print *, 'The number of sweeps changed from', NoPmove, 'to', Nosweeps
                smoothfactor = smoothfactormax
            end if
            
            ! 1. Calculate Second Derivative before       
            allocate(ddn_before(NoPmove),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
            call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny)
            
            !2. Move Boundary Nodes (linear)
            call LinearMotion(x, y, CNxfinal, CNyfinal)
            
            do isweep = 1, NoSweeps
                !3. Calculate Second Derivative after
                allocate(ddn_after(NoPmove),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
                call CalcSecondDerivative(ddn_after, NoPmove, x, y, nx, ny)
                
                !4. Apply smoothing
                do ismooth = 1, NoPmove
                    xnew(ismooth) = x(ismooth) + smoothfactor*(ddn_after(ismooth) - ddn_before(ismooth))*nx(ismooth)
                    ynew(ismooth) = y(ismooth) + smoothfactor*(ddn_after(ismooth) - ddn_before(ismooth))*ny(ismooth)  
                end do
                x = xnew
                y = ynew
            end do
            
            !5. Check for convergence 
            res = sqrt(sum((x(CN_ind) - CNxsmooth)**2 + (y(CN_ind) - CNysmooth)**2))
            if (initialResidual == 0) then
                initialResidual = res
            end if
            conv = log(res/initialResidual)/log(10.0)
            CNxsmooth = x(CN_ind) 
            CNysmooth = y(CN_ind)
        end do
        
        !6. if converged: Move CN locations and bound back onto desired position
        call LinearMotion(x, y, CNxfinal, CNyfinal)
        
        ! Hand over Coordinates to temporary coord
        RD%coord_temp(orderedBoundaryIndex,1) = x
        RD%coord_temp(orderedBoundaryIndex,2) = y
    
    end subroutine SmoothingLinear
    
    recursive subroutine SmoothingFDGD(CNDisp, NoIterFDGD)
    
        ! Variables
        implicit none
        integer :: ismooth, NoSweeps, NoPmove, NoIterFDGD, isweep, intersect, i
        double precision, dimension(:), allocatable :: ddn_before, ddn_after, nx, ny, x, y, xnew, ynew, CNxsmooth, CNysmooth, beta
        double precision, dimension(maxDoF) :: CNDisp
        double precision :: conv, initialResidual, res, convergence, smoothfactor, smoothfactormax
        
        ! Body of SubSmoothing
        smoothfactor = 0.000001
        convergence = -3
        initialResidual = 0.0
        conv = 0.0
        NoPmove = size(MovingGeomIndex, dim = 1)
        NoSweeps = NoPmove
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xnew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ynew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxsmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNysmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(nx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ny(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ddn_before(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ddn_after(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "

        ! Extract actual boundary order from boundary faces
        x = RD%coord(orderedBoundaryIndex,1)
        y = RD%coord(orderedBoundaryIndex,2)
        
        ! CN of smoothed surface
        CNxsmooth = x(CN_indordered)
        CNysmooth = y(CN_indordered)
                
        do while (conv > convergence)
            
            ! 0. Evaluate maximum smoothing and adjust parameters
            beta(1) = sqrt((x(2)-x(NoPmove))**2 + (y(2)-y(NoPmove))**2)/2
            do i = 2, NoPmove-1
                beta(i) = sqrt((x(i+1)-x(i-1))**2 + (y(i+1)-y(i-1))**2)/2
            end do
            beta(NoPmove) = sqrt((x(1)-x(NoPmove-1))**2 + (y(1)-y(NoPmove-1))**2)/2
            smoothfactormax = minval(beta)/2
            if (smoothfactor > smoothfactormax) then
                print *, 'The smoothfactor is too big and has been reduced from', smoothfactor, 'to', smoothfactormax
                Nosweeps = floor(Nosweeps*(smoothfactor/smoothfactormax))
                print *, 'The number of sweeps changed from', NoPmove, 'to', Nosweeps
                smoothfactor = smoothfactormax
            end if

            ! 1. Calculate Second Derivative before       
            call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny)
            
            !2. Pre-Processing of starting Geometry
            call PreMeshingBoundary()
            
            !3. Set new DelaunayCoord
            call RelocateCN(CNDisp, NoIterFDGD)
            
            !4. Move Boundary Mesh
            call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1)) 
            x = RD%coord_temp(orderedBoundaryIndex,1)
            y = RD%coord_temp(orderedBoundaryIndex,2)
            
            do isweep = 1, NoSweeps
                !5. Calculate Second Derivative after
                call CalcSecondDerivative(ddn_after, NoPmove, x, y, nx, ny)
                
                !6. Apply smoothing
                do ismooth = 1, NoPmove
                    xnew(ismooth) = x(ismooth) + beta(ismooth)*smoothfactor*(ddn_after(ismooth) - ddn_before(ismooth))*nx(ismooth)
                    ynew(ismooth) = y(ismooth) + beta(ismooth)*smoothfactor*(ddn_after(ismooth) - ddn_before(ismooth))*ny(ismooth)  
                end do
                x = xnew
                y = ynew
            end do
            RD%coord_temp(orderedBoundaryIndex,1) = x
            RD%coord_temp(orderedBoundaryIndex,2) = y
   
            !7. Check for convergence 
            res = sqrt(sum((x(CN_indordered) - CNxsmooth)**2 + (y(CN_indordered) - CNysmooth)**2))
            if (abs(initialResidual - res) < 10e-10) then
                conv = -4
            end if 
            if (initialResidual == 0) then
                initialResidual = res
            end if
            if (conv /= -4) then
                conv = log(res/initialResidual)/log(10.0)
            end if
                
            CNxsmooth = x(CN_indordered) 
            CNysmooth = y(CN_indordered)
        end do
        
        !8. if converged: Move CN locations and bound back onto desired position
        call PreMeshingBoundary()
        call RelocateCN(CNDisp, NoIterFDGD)
        call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))
                
        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        !call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect)
        intersect = 1
        ! Move Domain Nodes
        if (intersect == 1) then
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                CNDisp = CNDisp*(1.0/NoIterFDGD)
        else
            NoIterFDGD = NoIterFDGD*2
            call SmoothingFDGD(CNDisp, NoIterFDGD)
        end if
    
    end subroutine SmoothingFDGD
    
    subroutine CalcSecondDerivative(ddn, NoPmove, x, y, nx, ny)
    
        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, sx, sy, ddn, nx, ny
    
        ! Body of CalcSecondDerivative        
        allocate(sx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sy(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        
        ! Body of ComputeNormal
        nx(1) = -(y(2)-y(NoPmove))
        ny(1) = x(2)-x(NoPmove)
        do ip = 2, NoPmove-1
            nx(ip) = -(y(ip+1)-y(ip-1))
            ny(ip) = x(ip+1)-x(ip-1)
        end do
        nx(NoPmove) = -(y(1)-y(NoPmove-1))
        ny(NoPmove) = x(1)-x(NoPmove-1)

        ! Normalize
        do i = 1, NoPmove
            magnitude = sqrt(nx(i)**2 + ny(i)**2)
            nx(i) = nx(i)/magnitude
            ny(i) = ny(i)/magnitude
        end do
            
        ! Compute the s vectors
        sx = ny
        sy = -nx
             
        ! Compute Gradient
        dx1 = (x(1) - x(NoPmove))*sx(1) + (y(1) - y(NoPmove))*sy(1)
        dx2 = (x(2) - x(1))*sx(1) + (y(2) - y(1))*sy(1)
        ddn(1) = 2*(dx1*(x(2)*nx(1) + y(2)*ny(1)) - &
                (dx1+dx2)*(x(1)*nx(1) + y(1)*ny(1)) + &
                dx2*(x(NoPmove)*nx(1) + y(NoPmove)*ny(1)))/(dx1*dx2*(dx1 + dx2))
        do i = 2, NoPmove-1
            dx1 = (x(i) - x(i-1))*sx(i) + (y(i) - y(i-1))*sy(i)
            dx2 = (x(i+1) - x(i))*sx(i) + (y(i+1) - y(i))*sy(i)
            ddn(i) = 2*(dx1*(x(i+1)*nx(i) + y(i+1)*ny(i)) - &
                    (dx1+dx2)*(x(i)*nx(i) + y(i)*ny(i)) + &
                    dx2*(x(i-1)*nx(i) + y(i-1)*ny(i)))/(dx1*dx2*(dx1 + dx2))                     
        end do
        dx1 = (x(NoPmove) - x(NoPmove-1))*sx(NoPmove) + (y(NoPmove) - y(NoPmove-1))*sy(NoPmove)
        dx2 = (x(1) - x(NoPmove))*sx(NoPmove) + (y(1) - y(NoPmove))*sy(NoPmove)
        ddn(NoPmove) = 2*(dx1*(x(1)*nx(NoPmove) + y(1)*ny(NoPmove)) - &
                (dx1+dx2)*(x(NoPmove)*nx(NoPmove) + y(NoPmove)*ny(NoPmove)) + &
                dx2*(x(NoPmove-1)*nx(NoPmove) + y(NoPmove-1)*ny(NoPmove)))/(dx1*dx2*(dx1 + dx2))
        
    
    end subroutine CalcSecondDerivative
    
    subroutine LinearMotion(x, y, CNxfinal, CNyfinal)
    
        ! Variables
        implicit none
        integer :: i, j
        double precision :: CNdist
        double precision, dimension(maxDoF) :: CNDisp
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal
    
        ! Body of LinearMotion
        do i = 1, IV%NoCN
            CNDisp(i) = CNxfinal(i) - x(CN_ind(i))
            CNDisp(i+IV%NoCN) = CNyfinal(i) - y(CN_ind(i))
            x(CN_ind(i)) = x(CN_ind(i)) + CNDisp(i)
            y(CN_ind(i)) = y(CN_ind(i)) + CNDisp(i+IV%NoCN)
        end do
        
        ! Move boundary points linearly
        do i = 1, IV%NoCN          
            if (i /= IV%NoCN) then
                CNdist = abs(x(CN_ind(i+1)) - x(CN_ind(i)))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    y(j) = y(j) + CNDisp(i+IV%NoCN) - abs(x(CN_ind(i)) - x(j))*(1/CNdist)*CNDisp(i+IV%NoCN)
                    y(j) = y(j) + CNDisp(i+IV%NoCN+1) - abs(x(CN_ind(i+1)) - x(j))*(1/CNdist)*CNDisp(i+IV%NoCN+1)
                end do
                CNdist = abs(y(CN_ind(i+1)) - y(CN_ind(i)))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    x(j) = x(j) + CNDisp(i) - abs(y(CN_ind(i)) - y(j))*(1/CNdist)*CNDisp(i)
                    x(j) = x(j) + CNDisp(i+1) - abs(y(CN_ind(i+1)) - y(j))*(1/CNdist)*CNDisp(i+1)
                end do
            else
                CNdist = x(CN_ind(1)) - x(CN_ind(i))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    y(j) = y(j) + CNDisp(i+IV%NoCN) - abs(x(CN_ind(i)) - x(j))*(1/CNdist)*CNDisp(i+IV%NoCN)
                    y(j) = y(j) + CNDisp(1+IV%NoCN) - abs(x(CN_ind(1)) - x(j))*(1/CNdist)*CNDisp(1+IV%NoCN)
                end do                  
                CNdist = y(CN_ind(1)) - y(CN_ind(i))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    x(j) = x(j) + CNDisp(i) - abs(y(CN_ind(i)) - y(j))*(1/CNdist)*CNDisp(i)
                    x(j) = x(j) + CNDisp(1) - abs(y(CN_ind(1)) - y(j))*(1/CNdist)*CNDisp(1)
                end do
            end if
        end do    


    end subroutine LinearMotion
    
    end module Smoothing