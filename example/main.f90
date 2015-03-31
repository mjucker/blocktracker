program main
  use blob_tracker_mod, only: init_blobs,find_blocks,totBlobs,numBlocks
  implicit none
  
  integer :: t                     ! time integration index
  integer :: nx=100, ny = 150      ! number of grid points in x and y
  integer :: nt=30                 ! number of integration steps
  integer :: start_t=10, stop_t=20 ! start and stop integration step where blocks appear
  real    :: sigma=0.1             ! Gaussian block widths
  real    :: noise_amp = 0.        ! background noise amplitude
  integer :: nsecx,nsecy
  real,dimension(:,:),allocatable :: x,y,data
  integer :: i,j,margBlobs,margBlocks
  real    :: xloc,yloc

  namelist /blocks_nml/ nx, ny, nt, start_t, stop_t, sigma, noise_amp

  
  !read namelist
  open(11, FILE='input.nml')
  read(11, nml=blocks_nml)

  allocate(x(nx,ny),y(nx,ny),data(nx,ny))

  !create (x,y) grid
  do j=1,ny
     yloc = (j-1.)/(ny-1.)
     do i=1,nx
        x(i,j) = (i-1.)/(nx-1.)
        y(i,j) = yloc
     enddo
  enddo
  
  call init_blobs

  call random_seed

  margBlobs = 0
  margBlocks= 0
  do t=1,nt
     ! start with some random noise
     call random_number(data)
     data = noise_amp*2*(0.5 - data)
     if ( t .ge. start_t .and. t .le. stop_t ) then
        ! add one Gausian
        xloc = 1.*(t - start_t)/(stop_t - start_t)
        yloc = min(0.5,xloc)
        data = data - exp(-((x-xloc)/sigma)**2 - ((y-yloc)/sigma)**2)
        ! add another Gaussian
        xloc = 1. - 1.*(t - start_t)/(stop_t - start_t)
        yloc = 0.75
        data = data - exp(-((x-xloc)/sigma)**2 - ((y-yloc)/sigma)**2)
        ! add another Gaussian
        xloc = 1.*(t - start_t)/(stop_t - start_t)
        yloc = 1. - xloc
        data = data - exp(-((x-xloc)/sigma)**2 - ((y-yloc)/sigma)**2)
        do j=1,ny
           do i=1,nx
              write(16,'(i3,3es14.6)')t,x(i,j),y(i,j),data(i,j)
           enddo
        enddo
     endif
     call find_blocks(data,real(t)/nt)
     margBlobs = totBlobs - margBlobs
     margBlocks= numBlocks- margBlocks
     write(*,'(a,i3,a,i3,a,i4)')'Found ',margBlobs,' blobs, ',margBlocks,' blocks at time ',t
     margBlobs = totBlobs
     margBlocks= numBlocks
  enddo
  write(*,'(a,i3,a)')'Found a total of ',numBlocks,' blocks.'


end program main
