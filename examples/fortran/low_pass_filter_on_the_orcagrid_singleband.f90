program low_pass_filter_on_the_orcagrid
  USE netcdf
  !USE ensdam_spharea
  USE ensdam_sphylm
  IMPLICIT none

  CHARACTER(len=256) :: domcfg, infile, outfile, varname, arg
  INTEGER :: lmin, lmax
  INTEGER :: nlon, nlat
  INTEGER :: array_size, array_shape(2)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spectrum
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_1d, y_1d, ranfield_1d

  INTEGER :: i, j, l, rank

  INTEGER :: is, ncid, idx, idy, idv, idlon, idlat
  REAL(KIND=8) :: latmin, latmax, dlatmax


  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x, y, ranfield, area, lon, lat, mask


  ! Initialize parallel computation
#if defined MPI
  integer, save :: mpi_code
  include "mpif.h"
  call mpi_init(mpi_code)
  !call mpi_comm_size(mpi_comm_world,nproc,mpi_code)
  call mpi_comm_rank(mpi_comm_world,rank,mpi_code)
#endif

  ! Get command line arguments
  if (command_argument_count() < 4) then
      print *, "Usage: ./low_pass_filter_on_the_orcagrid domcfg.nc input.nc output.nc varname"
      stop 1
  endif

  call get_command_argument(1, domcfg)
  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)
  call get_command_argument(4, varname)
  call get_command_argument(5, arg); read(arg, *) lmin
  call get_command_argument(6, arg); read(arg, *) lmax

  IF (rank==0) THEN
     print *, "Input file:   ", trim(infile)
     print *, "Output file:  ", trim(outfile)
     print *, "Variable:     ", trim(varname)
     print *, "lmin/lmax:     ", lmin, lmax
  ENDIF

  IF (rank==0) THEN 
     CALL read_var_2d(TRIM(domcfg), "e1t", x)
     CALL read_var_2d(TRIM(domcfg), "e2t", y)
     array_shape = SHAPE(x)
  ENDIF

  ! Bcast the size of x,y
  CALL MPI_Bcast(array_shape, 2, MPI_INT, 0, mpi_comm_world, mpi_code)
  
  IF (rank/=0) THEN
     ALLOCATE(x(array_shape(1), array_shape(2)))
     ALLOCATE(y(array_shape(1), array_shape(2)))
  ENDIF

  CALL MPI_Bcast(x, SIZE(x), MPI_DOUBLE, 0, mpi_comm_world, mpi_code)
  CALL MPI_Bcast(y, SIZE(y), MPI_DOUBLE, 0, mpi_comm_world, mpi_code)

  ! Caluclate the area of each grid cell
  ALLOCATE(mask, mold=x)
  ALLOCATE(area, mold=x)
  area = x * y
  area = area/SUM(area)

  ! DEALLOCATE scale factors
  DEALLOCATE(x, y)

  ! Read field and coords
  IF (rank==0) THEN 
     call read_var_2d(TRIM(infile), TRIM(varname), ranfield)
     call read_var_2d(TRIM(infile), "nav_lon", lon)
     call read_var_2d(TRIM(infile), "nav_lat", lat)
  ENDIF
  IF (rank/=0) THEN
     ALLOCATE(ranfield(array_shape(1), array_shape(2)))
     ALLOCATE(lon(array_shape(1), array_shape(2)))
     ALLOCATE(lat(array_shape(1), array_shape(2)))
  ENDIF

  CALL MPI_Bcast(ranfield, SIZE(x), MPI_DOUBLE, 0, mpi_comm_world, mpi_code)
  CALL MPI_Bcast(lon, SIZE(y), MPI_DOUBLE, 0, mpi_comm_world, mpi_code)
  CALL MPI_Bcast(lat, SIZE(y), MPI_DOUBLE, 0, mpi_comm_world, mpi_code)

  IF (rank==0) THEN 
     print *, "printing shapes of used arrays"
     print *, "e1t size: ", shape(x)
     print *, "e2t size: ", shape(y)
     print *, "sla size: ", shape(ranfield)
     print *, "area size: ", shape(area)
  ENDIF
  
  ! Weight input field with area of each grid cell
  mask(:,:) = 1.0
  WHERE(ABS(ranfield)>100) mask = 0.
  ranfield = ranfield * area * mask
  
  ! Initialize computation of spherical harmonics (up to degree lmax)
  latmin = -90. ; latmax = 90. ; dlatmax = 0.25
  CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )

  ! Compute spectrum of the input field up to degree lmax
  array_size = SIZE(area)
  array_shape = SHAPE(area)
  ALLOCATE(x_1d(array_size), y_1d(array_size), ranfield_1d(array_size))
  ALLOCATE(spectrum(0:lmax,-lmax:lmax))
  x_1d=RESHAPE(lon,(/array_size/))
  y_1d=RESHAPE(lat,(/array_size/))
  ranfield_1d=RESHAPE(ranfield,(/array_size/))
  call proj_ylm(spectrum,ranfield_1d,x_1d,y_1d)

  ! Transform back into original space (up to degree lmax)
  !call back_ylm(spectrum,ranfield_1d,x_1d,y_1d)
  call back_ylm_loc(spectrum,ranfield_1d,x_1d,y_1d, 80, lmax)
  ranfield=RESHAPE(ranfield_1d,(/SHAPE(area)/))
  ranfield = ranfield * mask  !/ area

  IF (rank==0) THEN  
     ! Write output file (in NetCDF)
     is = NF90_CREATE(TRIM(outfile),NF90_CLOBBER,ncid)
     is = NF90_DEF_DIM(ncid,'x',array_shape(1),idx)
     is = NF90_DEF_DIM(ncid,'y',array_shape(2),idy)
     is = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx,idy/),idlon)
     is = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idx,idy/),idlat)
     is = NF90_DEF_VAR(ncid,'ranfield',NF90_FLOAT,(/idx,idy/),idv)
     is = NF90_ENDDEF(ncid)
     is = NF90_PUT_VAR(ncid,idlon,lon)
     is = NF90_PUT_VAR(ncid,idlat,lat)
     is = NF90_PUT_VAR(ncid,idv,ranfield)
     is = NF90_CLOSE(ncid)
  ENDIF
contains

  subroutine read_var_2d(filename, varname, field)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: varname
    real(kind=8), allocatable, intent(out) :: field(:,:)

    integer :: ncid, varid, status, xtype
    integer :: dimids(NF90_MAX_VAR_DIMS), ndims
    integer :: n1, n2

    status = NF90_OPEN(filename, NF90_NOWRITE, ncid)
    if (status /= NF90_NOERR) stop "Error opening file: "//trim(filename)

    status = NF90_INQ_VARID(ncid, varname, varid)
    if (status /= NF90_NOERR) stop "Variable not found: "//trim(varname)

    status = NF90_INQUIRE_VARIABLE(ncid, varid, ndims = ndims, dimids=dimids)
    if (ndims < 2) stop "Variable "//trim(varname)//" is not 2D!"

    status = NF90_INQUIRE_DIMENSION(ncid, dimids(1), len=n1)
    status = NF90_INQUIRE_DIMENSION(ncid, dimids(2), len=n2)

    allocate(field(n1,n2))

    status = NF90_GET_VAR(ncid, varid, field)
    if (status /= NF90_NOERR) stop "Error reading variable: "//trim(varname)

    status = NF90_CLOSE(ncid)
  end subroutine read_var_2d

end
