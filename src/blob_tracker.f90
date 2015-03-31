  module blob_tracker_mod
   implicit none
!=====================================================================
! --- NAMELIST
!=====================================================================

   logical :: track_blobs = .true.
! contour finding parameters:
   real    :: crit_val     = -1.0 !units of data field in the code
   real    :: crit_area    = 1e-2 !units of total domain area, \in [0,1]
   real    :: diff_dist    = 1e-1 !units of total domain size: distinguish between contours
   real    :: max_dist     = 1.0  !units of total domain size: distinguish open contours
! block finding parameters:
   real    :: crit_dist    = 0.5  !units of radius 
   real    :: crit_dur     = 5.0  !units of time
   integer :: maxdur       = 100  !max number of time steps for block life
! generally used parameters:
   integer :: maxnum       = 100  ! max number of contours and blobs per time step
   integer :: maxpts       = 5000  !max number of points per contour and blob
!=====================================================================
   integer :: numBlobs, numConts, prev_numBlobs, numCands, numBlocks
   integer :: totBlobs
   integer,dimension(:),allocatable :: numPts,blockCandDuration
   real                ,allocatable :: blobArea(:),blobCentr(:,:),contours(:,:,:)
   real                ,allocatable :: prev_Area(:),prev_Centrs(:,:)
   real                ,allocatable :: blockCandCentrs(:,:,:),blockCandAreas(:,:)&
        &,blockCandTime(:,:)
   real                ,allocatable :: blockAreas(:,:),blockCentrs(:,:,:),blockTimes(:,:)
   real                             :: prev_time
   real     :: pi
   parameter ( pi = 3.14159 )
   logical :: vecout_called,blob_called
   real :: xlast,ylast,dx,dy
   integer :: nx,ny

   namelist /block_tracker_nml/ track_blobs, &
                                crit_val, crit_area,  crit_dist, crit_dur, &
                                maxnum, maxpts, &
                                diff_dist, max_dist

!=====================================================================
    contains

!######################################################################
      subroutine init_blobs
        implicit none
        
! read namelist
        open(8, FILE='blocks.nml')
        read(8, nml=block_tracker_nml)
        close(8)

        totBlobs = 0
        numBlobs = 0
        numConts = 0
        prev_numBlobs = 0
        numCands = 0
        if( track_blobs )then
           allocate(numPts(maxnum),blobArea(maxnum),blobCentr(2,maxnum))
           allocate(contours(2,maxpts,maxnum)) !(x,y),points per contour,contours
           allocate(prev_Area(maxnum),prev_Centrs(2,maxnum))
           allocate(blockCandCentrs(2,maxdur,maxnum),blockCandDuration(maxnum))
           allocate(blockCandAreas(maxdur,maxnum),blockCandTime(maxdur,maxnum))
           allocate(blockAreas(maxdur,maxnum),blockCentrs(2,maxdur,maxnum)&
                &,blockTimes(maxdur,maxnum))
        endif
        numPts   = 0
        blockCandDuration = 0

        vecout_called = .false.
        blob_called = .false.
        
      end subroutine init_blobs

!######################################################################

      subroutine find_blocks(data,time)
        implicit none
        real,dimension(:,:),intent(in) :: data
        real               ,intent(in) :: time
        integer :: i,j

        call find_blob_candidates(data,time)

        if ( blob_called .and. numBlobs > 0 ) then ! not first time step
           call find_block_candidates(time)
        else
           blob_called = .true.
        endif

        call check_block_candidates
              
      end subroutine find_blocks

!######################################################################
      subroutine find_blob_candidates(data,time)
        implicit none
        real,    dimension(:,:),intent(in) :: data  !2D gridded data
        real                   ,intent(in) :: time  !purely for diagnostics
!
        integer, dimension(:,:),allocatable :: candidates, picked
        integer :: minind(2),i,j,k,l
        real,dimension(:  ),allocatable :: x,y
        real    :: normalize(2),area,centroid(2)
        logical :: last = .true.

        nx = size(data,1)
        ny = size(data,2)
        allocate(x(nx),y(ny),candidates(nx,ny),picked(nx,ny))
        !global x and y will be within [0,1]
        dx = 1./(nx-1)
        dy = 1./(ny-1)
        do i=1,nx
           x(i) = (i-1)*dx
        enddo
        do j=1,ny
           y(j) = (j-1)*dy
        enddo
        ! find contours of crit_val
        if(minval(data) < crit_val)then
           call conrec(data,1,nx,1,ny,x,y,1,crit_val)
        else
           return
        endif
        ! find characteristics of contours
        do i=1,numConts
           call compute_blob(data,contours(1,1:numPts(i),i),contours(2,1:numPts(i),i),area,centroid)
           if( area > crit_area .and. numBlobs < maxnum)then !we've got a block candidate
              numBlobs = numBlobs+1
              blobArea(numBlobs) = area
              blobCentr(:,numBlobs) = centroid
           endif
        enddo
        
        ! diagnostics
        totBlobs = totBlobs + numBlobs
        do i=1,numBlobs
           do j=1,numPts(i)
              write(16,"(es11.3,',',i3,2(',',es11.3))") time,i,contours(:,j,i)
           enddo
           write(17,"(es11.3,',',i3,3(',',es11.3))") time,i,blobCentr(:,i),blobArea(i)
        enddo

        ! re-initialize after this
        numConts = 0
        numPts   = 0
        vecout_called = .false.
        
      end subroutine find_blob_candidates

!######################################################################

      subroutine find_block_candidates(time)
        implicit none
        real,intent(in) :: time
        real    :: radius,prev_locs(2,1,maxnum),dists(maxnum)
        integer :: i,j,jnd(1)
        logical :: foundMatch
        

        do i=1,numBlobs
           radius = sqrt(blobArea(i)/pi)
           ! check current blobs against last position of block candidates
           if ( numCands > 0 ) then
              do j=1,numCands
                 dists(j) = sqrt( ( blobCentr(1,i) - blockCandCentrs(1,blockCandDuration(j),j) )**2&
                      & + ( blobCentr(2,i) - blockCandCentrs(2,blockCandDuration(j),j) )**2 )
              enddo
              ! now find the closest candidate
              jnd = minloc(dists(1:numCands))
              j = jnd(1)
              foundMatch = .false.
              if ( dists(j) < radius*crit_dist ) then ! overlap with a block candidate
                 blockCandDuration(j) = blockCandDuration(j) + 1
                 blockCandCentrs(:,blockCandDuration(j),j) = blobCentr(:,i)
                 blockCandAreas(blockCandDuration(j),j) = blobArea(i)
                 blockCandTime(blockCandDuration(j),j) = time
                 foundMatch = .true.
              endif
           endif
           if ( .not. foundMatch .and. prev_numBlobs > 0 ) then !check among previous blobs 
              prev_locs(:,1,:)  = prev_Centrs
              call compute_dists(blobCentr(1,i),blobCentr(2,i),prev_locs(:,:,1:prev_numBlobs),dists(1:prev_numBlobs))
              jnd = minloc(dists(1:prev_numBlobs))
              j = jnd(1)
              if ( dists(j) < radius*crit_dist ) then ! we have a new block candidate
                 numCands = numCands + 1
                 blockCandDuration(numCands) = 2
                 blockCandCentrs(:,1,numCands) = prev_Centrs(:,j)
                 blockCandAreas(1,numCands) = prev_Area(j)
                 blockCandTime(1,numCands) = prev_time
                 blockCandCentrs(:,2,numCands) = blobCentr(:,i)
                 blockCandAreas(2,numCands) = blobArea(i)
                 blockCandTime(2,numCands) = time
              endif
           endif
        enddo

        ! initialize next time step
        prev_Centrs   = blobCentr
        prev_Area     = blobArea
        prev_time     = time
        prev_numBlobs = numBlobs
        numBlobs = 0

      end subroutine find_block_candidates

!######################################################################

      subroutine check_block_candidates
        implicit none
        integer i,j,d
        
        blockAreas = 0.
        blockCentrs= 0.
        blockTimes = 0.
        numBlocks  = 0
        do i=1,numCands
           d = blockCandDuration(i)
           if ( d >= crit_dur ) then !we have a block
              numBlocks = numBlocks + 1
              blockAreas(1:d,numBlocks) = blockCandAreas(1:d,i)
              blockCentrs(:,1:d,numBlocks) = blockCandCentrs(:,1:d,i)
              blockTimes(1:d,numBlocks) = blockCandTime(1:d,i)
           endif
        enddo

      end subroutine check_block_candidates
!######################################################################

      subroutine compute_blob(data,xcont,ycont,area,centroid)
        implicit none

        real,dimension(:,:)   ,intent(in) :: data
        real,dimension(:)     ,intent(in) :: xcont,ycont
        real                  ,intent(out):: area
        real,dimension(2)     ,intent(out):: centroid

        real,allocatable,dimension(:,:) :: f
        real :: dA,N,rsq,xc,yc,maxdist,max_x,min_x,max_y,min_y
        integer :: is,js,iss,jss,lx,ly,i,j
        
        max_x = maxval(xcont)
        min_x = minval(xcont)
        max_y = maxval(ycont)
        min_y = minval(ycont)

        maxdist = max(max_x - min_x, max_y - min_y)
        if ( maxdist > max_dist ) then !this is an open contour
           centroid = 0.
           area     = 0.
           return
        endif

        is = max(1 ,floor  (min_x*nx))
        iss= min(nx,ceiling(max_x*nx))
        js = max(1 ,floor  (min_y*ny))
        jss= min(ny,ceiling(max_y*ny))

        lx=iss-is+1
        ly=jss-js+1
        allocate(f(lx,ly))
        !compute normalization constant
        N=0.
        do j=js,jss
           do i=is,iss
              f(i-is+1,j-js+1) = min(0., data(i,j) - crit_val)
              N = N + f(i-is+1,j-js+1)
           enddo
        enddo
        !compute centroid (in global coordinates)
        xc=0.
        yc=0.
        do j=js,jss
           do i=is,iss
             xc = xc + (i-1)*dx*f(i-is+1,j-js+1)
             yc = yc + (j-1)*dy*f(i-is+1,j-js+1)
          enddo
       enddo
       xc = xc/N
       yc = yc/N
       centroid(1) = xc
       centroid(2) = yc
       !compute equivalent radius square
       rsq=0.
       do j=js,jss
          do i=is,iss
             rsq = rsq + ( ((i-1)*dx - xc)**2 + ((j-1)*dy - yc)**2 )*f(i-is+1,j-js+1)
          enddo
       enddo
       rsq = rsq/N
       !finally, get equivalent area
       area = pi*rsq

      end subroutine compute_blob

!######################################################################

!
!======================================================================
!
!     This is a sample vector output routine. For a local environment
!     either replace the VECOUT call in the main line, or better
!     replace the contents of this subroutine between the *'s shown.
!
!     There is often the requirement to distinguish each contour
!     line with a different colour or a different line style. This
!     can be done in many ways using the contour values z for a
!     particular line segment.
!
      subroutine vecout(data,x1,y1,x2,y2,z,last)
        implicit none
        real,dimension(nx,ny),intent(in) :: data
        real                 ,intent(in) :: x1,y1,x2,y2,z
        logical              ,intent(in) :: last
        real :: dist,thisarea,area,centroid(2)
        real :: dists(maxnum)
!!$     real,parameter :: eps=1.e-16
        integer :: minxy(1),ind,dummy
        !
        if(.not. vecout_called)then !fist time this is called
           numConts = 1
           numPts   = 0
           numPts(1)= 1
           contours = 2.
           contours(1,1,1) = 0.5*(x1+x2)
           contours(2,1,1) = 0.5*(y1+y2)
           dists = 2.
           vecout_called = .true.
           return
        endif
        
        call compute_dists(x1,y1,contours(:,:,1:numConts),dists(1:numConts))
       
        if ( minval(dists(1:numConts)) > diff_dist ) then ! we have a new contour
           if ( numConts < maxnum ) then
              numConts = numConts + 1
              numPts(numConts) = 1
              contours(1,1,numConts) = 0.5*(x1+x2)
              contours(2,1,numConts) = 0.5*(y1+y2)
           else
              print*,'WARNING: MAXIMUM NUMBER OF CONTOURS REACHED'
           endif
        else ! this is part of an existing contour
           minxy = minloc(dists(1:numConts)) ! the contour with minimum distance
           ind = minxy(1)
           if ( numPts(ind) < maxpts ) then
              numPts(ind) = numPts(ind) + 1
              contours(1,numPts(ind),ind) = 0.5*(x1+x2)
              contours(2,numPts(ind),ind) = 0.5*(y1+y2)
           else
              print*,'WARNING: MAXIMUMM POINTS PER CONTOUR REACHED'
           endif
        endif
        
      end subroutine vecout

!######################################################################
      

      subroutine compute_dists(x,y,pos_vec,dists)
        implicit none
        real                 ,intent(in) :: x,y
        real,dimension(:,:,:),intent(in) :: pos_vec !(x,y),points,contours
        real,dimension(:    ),intent(out):: dists
        real,dimension(size(pos_vec,2)) :: locdists

        integer i,j,s

        s = size(pos_vec,2)
        locdists = 2.

        ! go through all existing contours and determine minimum distance of x,y
        do i=1,size(dists)
           where ( pos_vec(1,:,i) <= 1. .and. pos_vec(2,:,i) <= 1. )
              locdists(:) = sqrt( (x - pos_vec(1,:,i))**2 + (y - pos_vec(2,:,i))**2 )
           endwhere
           dists(i) = minval(locdists)
        enddo
      end subroutine compute_dists
!######################################################################
        
      subroutine end_blobs
        implicit none
        integer n,d
        
        if( track_blobs )then
           open(unit=8,file='blocks.csv')
           do n=1,numBlocks
              do d=1,maxnum
                 if( blockTimes(d,n) .ne. 0. )then
                    write(8,101) n,blockTimes(d,n),blockCentrs(1,d,n), &
                         blockCentrs(2,d,n), blockAreas(d,n)
                 endif
              enddo
           enddo
           close(8)

101        format(1x,i3,4(',',es14.6))

           deallocate(numPts,blobArea,blobCentr,contours,prev_Area,prev_Centrs)
           deallocate(blockCandCentrs,blockCandDuration)
           deallocate(blockCandAreas,blockCandTime,blockAreas,blockCentrs,blockTimes)
        endif

      end subroutine end_blobs

!#######################################################################

    end module blob_tracker_mod
