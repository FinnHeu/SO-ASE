
program main
  ! To smooth Topography over the model grid by taking the weighted
  ! mean of neighbour values. 
  !
  ! The smoothing starts from the shallowest grid point to the deepest
  ! grid point. The smoothed value of a center node will be the weighted
  ! mean value from the neighbour nodes including the center node itself.
  ! One can specify the weight for the center node. The basic weight is
  ! the number of surrounding neighbor nodes; the effective weight will
  ! be the basis weight times the specified relative_center_weight value.
  ! One can also specify how many iterations to do such
  ! smoothing.
  !
  ! Input files:
  ! mesh and depth files (nod2d.out, elem2d.out, nodhn.out)
  !
  ! Qiang, 19.05.2011
  !------------------------------------------------------------------------------

  implicit none

  integer                              :: i, j, m, n, k, a, tr(3),tet(4)
  integer                              :: counter, el, ml(1), cnt
  integer                              :: num_iteration
  integer                              :: nod2d, elem2d
  integer, allocatable, dimension(:,:) :: elem2d_nodes
  integer, allocatable, dimension(:)   :: ind, ind2d
  integer, dimension(100)              :: aux=0
  real(kind=8)                         :: dsum, h_min, h_max, relative_center_weight,rad
  real(kind=8)                         :: min_water_depth, max_water_depth, coeff1
  real(kind=8), allocatable            :: dep(:), dep_old(:), xynod(:,:),lon(:),lat(:)
  logical                              :: set_min_depth=.false.
  logical                              :: set_max_depth=.false. 
  logical			       :: rotated_grid=.true. 
  logical			       :: special_region=.true.
  character(100)                        :: meshdir='./mesh/'

  type addresstype
     integer                                   :: nmb
     integer(KIND=4), dimension(:), pointer    :: addresses
  end type addresstype
  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D 
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D

  !------------------------------------------------------------
  ! user specification starts here

  meshdir='./mesh_RTopo2.0.4_30sec_SOCAV_v4/01_raw/'

  rotated_grid=.false.  !true or false depends on how you specify the special region: tricky thing

  special_region=.false.  !go to the code to do modification!

  relative_center_weight=2.0    ! 2! >=1
  num_iteration= 1 !3  !2

  set_min_depth=.true.
  min_water_depth=30.0

  set_max_depth=.true.
  max_water_depth=6500.

  !user specification ends here
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! read mesh
  write(*,*) 'Reading the mesh...'
  write(*,*) 'The mesh is: ', meshdir

  open(11, file=trim(meshdir)//'nod2d.out')
  read(11,*) nod2d
  allocate(xynod(nod2d,2),ind2d(nod2d))
  do n=1,nod2d
     read(11,*) i, xynod(n,1),xynod(n,2),ind2d(n)
  enddo
  close(11)

  allocate(lon(nod2d),lat(nod2d))
  rad=3.141592653589793/180.0
  xynod=xynod*rad
  if(rotated_grid) then
  do n=1,nod2d
     call r2g(lon(n), lat(n), xynod(n,1),xynod(n,2))
  end do
  else
  	lon=xynod(:,1)
	lat=xynod(:,2)
  endif
  lon=lon/rad
  lat=lat/rad

  open(12, file=trim(meshdir)//'elem2d.out')
  read(12,*) elem2d
  allocate(elem2d_nodes(3,elem2d))
  do n=1,elem2d
     read(12,*) elem2d_nodes(1:3,n)
  end do
  close(12)

  allocate(dep(nod2d), dep_old(nod2d))
  open(13,file=trim(meshdir)//'nodhn.out')
  read(13,*) dep
  close(13)

  dep=-dep

  if(set_min_depth) then
     min_water_depth=-abs(min_water_depth)
     do k=1,nod2d
        if(dep(k)>min_water_depth) dep(k)=min_water_depth
     end do
  end if
  h_min=abs(maxval(dep))

  if(set_max_depth) then
     max_water_depth=-abs(max_water_depth)
     do k=1,nod2d
        if(dep(k)<max_water_depth) dep(k)=max_water_depth
     end do
  end if
  h_max=abs(minval(dep))

  write(*,*) 'Mesh and depth files are read'
  if(set_min_depth .or. set_max_depth) then
     write(*,*) 'After applying the min/max depth limitation,'
  end if
  write(*,*) 'The minimum depth is ', h_min
  write(*,*) 'The maximum depth is ', h_max

  !----------------------------------------------------------
  ! Builds nod_in_elem2D

  allocate(ind(nod2D))

  ind=0
  do j=1,elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(nod2D))
  nod_in_elem2D%nmb=ind    
  do j=1,nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,elem2D   
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
     do k=1,3
        nod_in_elem2D(tr(k))%addresses(ind(tr(k)))=j
     end do
  end do

  ! Builds nghbr_nod2D
  allocate(nghbr_nod2D(nod2D))
  ind=0
  do j=1, nod2D
     counter=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if(a==j) cycle    ! the neighbour array of a node does not contain itself!!
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=counter
     allocate(nghbr_nod2D(j)%addresses(counter))

     ! we need to sort array aux(1:counter)
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod2D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do

  !-------------------------------------------------------------
  ! smooth topography

  write(*,*) 'Smoothing topography  ...'

  do i=1, num_iteration
     dep_old=dep

     do j=1,nod2D

        if (mod(j, 10000) == 0) then
           write(*,'(I10,A,I10)') j, '/', nod2D
        end if

        ml=maxloc(dep_old)
        m=ml(1)
        dep_old(m)=-99999.0

  	coeff1=1.0
        !if(lat(m)>82.0 .and. dep(m)>=maxval(dep(nghbr_nod2d(m)%addresses))) then
        !if(lat(m)>83.0) then
	if(special_region .and. lat(m)>63.0 .and. lat(m)<80.0 .and. (lon(m)<-158.0 .or. lon(m)>178.0)) then
		coeff1=3.0
	endif

        dsum=0.0
        cnt=nghbr_nod2d(m)%nmb
        do k=1, cnt
           n=nghbr_nod2d(m)%addresses(k)
           dsum=dsum+dep(n)
        end do

        dsum=dsum + real(coeff1*relative_center_weight*(cnt-1.0)-1.0)*dep(m)

        dep(m)=dsum/real(cnt-1+(cnt-1)*coeff1*relative_center_weight)
     end do

     write(*,*) 'The ',i, ' iteration finished.'

  end do

  !-------------------------------------------------------------
  !save modified/smoothed depth

  open(13, file=trim(meshdir)//'depth_smooth.out')
  do m=1,nod2d
     write(13,'(1f8.0)') dep(m)
  end do
  close(13)
  write(*,*) 'The smoothed depth is saved.'
  write(*,*) 'It is suggested to compare the new depth with the original depth.'

end program main

  !
  !----------------------------------------------------------------
  !

subroutine r2g(lon, lat, rlon, rlat)

  implicit none

  real(kind=8)        :: rotate_matrix(3,3)
  real(kind=8)      :: al, be, ga, rad
  real(kind=8)      :: xr, yr, zr, xg, yg, zg
  real(kind=8), intent(out)      :: lon, lat
  real(kind=8), intent(in)       :: rlon, rlat

  real(kind=8)         	:: alphaEuler=50. 		![degree] Euler angles, convention:
  real(kind=8)         	:: betaEuler=15.  		![degree] first around z, then around new x,
  real(kind=8)		:: gammaEuler=-90.		![degree] then around new z.
  
  rad=3.141592653589793/180.0

  al=alphaEuler*rad
  be=betaEuler*rad
  ga=gammaEuler*rad

  ! rotation matrix
  rotate_matrix(1,1)=cos(ga)*cos(al)-sin(ga)*cos(be)*sin(al)
  rotate_matrix(1,2)=cos(ga)*sin(al)+sin(ga)*cos(be)*cos(al)
  rotate_matrix(1,3)=sin(ga)*sin(be)
  rotate_matrix(2,1)=-sin(ga)*cos(al)-cos(ga)*cos(be)*sin(al)
  rotate_matrix(2,2)=-sin(ga)*sin(al)+cos(ga)*cos(be)*cos(al)
  rotate_matrix(2,3)=cos(ga)*sin(be)
  rotate_matrix(3,1)=sin(be)*sin(al) 
  rotate_matrix(3,2)=-sin(be)*cos(al)  
  rotate_matrix(3,3)=cos(be)


  ! Rotated Cartesian coordinates:
  xr=cos(rlat)*cos(rlon)
  yr=cos(rlat)*sin(rlon)
  zr=sin(rlat)

  ! Geographical Cartesian coordinates:
  xg=rotate_matrix(1,1)*xr + rotate_matrix(2,1)*yr + rotate_matrix(3,1)*zr
  yg=rotate_matrix(1,2)*xr + rotate_matrix(2,2)*yr + rotate_matrix(3,2)*zr  
  zg=rotate_matrix(1,3)*xr + rotate_matrix(2,3)*yr + rotate_matrix(3,3)*zr  

  ! Geographical coordinates:
  lat=asin(zg)
  if(yg==0. .and. xg==0.) then
     lon=0.0     ! exactly at the poles
  else
     lon=atan2(yg,xg)
  end if
end subroutine r2g

