! This based on pres_temp_4D_wr.f90 but changes
! to have the right use of time for ParaView
! and refactored for use in the LES.
! This file copyright 2014 Wim Vanderbauwhede
! pres_temp_4D_wr.f90  Copyright 2006 University Corporation for Atmospheric Research/Unidata.

module module_LES_write_netcdf

  use netcdf
  use common_sn
  implicit none
  integer :: ncid, ncid_p, ncid_u,ncid_v,ncid_w,ncid_uvwsum
  integer :: pres_varid, pres_varid_p
  integer :: vel_x_varid_u, vel_y_varid_v, vel_z_varid_w
  integer :: velsum_x_varid_uvwsum, velsum_y_varid_uvwsum, velsum_z_varid_uvwsum
  ! These settings tell netcdf to write one timestep of data. (The
  ! setting of start(4) inside the loop below tells netCDF which
  ! timestep to write.)
  integer, parameter :: NDIMS = 4
  integer, parameter :: NTIMESTEPS=20 ! FIXME: should be taken from global macro
  ! For p only: 0:ip+2,0:jp+2,0:kp+1
  integer, parameter :: NLVLS_P_UV = kp+2, NLATS_P_UVW = jp+3, NLONS_P = ip+3
  ! For u,v
  integer, parameter :: NLONS_UVW = ip+2
  ! For w
  integer, parameter :: NLVLS_W = kp+3
  ! For usum/vsum/wsum
  integer, parameter :: NLVLS_UVWSUM = kp+1, NLATS_UVWSUM = jp+1, NLONS_UVWSUM = ip+1

  integer :: start(NDIMS), count(NDIMS)
  integer :: count_p(NDIMS), count_u(NDIMS),count_v(NDIMS),count_w(NDIMS), count_uvwsum(NDIMS)

  integer :: lvl_dimid, lon_dimid, lat_dimid, time_series_dimid
  integer :: lvl_dimid_p, lon_dimid_p, lat_dimid_p, time_series_dimid_p
  integer :: lvl_dimid_u, lon_dimid_u, lat_dimid_u, time_series_dimid_u
  integer :: lvl_dimid_v, lon_dimid_v, lat_dimid_v, time_series_dimid_v
  integer :: lvl_dimid_w, lon_dimid_w, lat_dimid_w, time_series_dimid_w
  integer :: lvl_dimid_uvwsum, lon_dimid_uvwsum, lat_dimid_uvwsum, time_series_dimid_uvwsum

  integer :: dimids(NDIMS)
  integer :: dimids_p(NDIMS)
  integer :: dimids_u(NDIMS)
  integer :: dimids_v(NDIMS)
  integer :: dimids_w(NDIMS)
  integer :: dimids_uvwsum(NDIMS)

 ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "./LES_output.nc"
  character (len = *), parameter :: FILE_NAME_P = "./LES_output_p.nc"
  character (len = *), parameter :: FILE_NAME_U = "./LES_output_u.nc"
  character (len = *), parameter :: FILE_NAME_V = "./LES_output_v.nc"
  character (len = *), parameter :: FILE_NAME_W = "./LES_output_w.nc"
  character (len = *), parameter :: FILE_NAME_UVWSUM = "./LES_output_uvwsum.nc"
  
  character (len = *), parameter :: LVL_NAME = "level"
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: TIME_SERIES_NAME = "time"

  ! These program variables hold the latitudes and longitudes.
  real :: times(NTIMESTEPS)
  real :: lats_p_uvw(NLATS_P_UVW), lons_p(NLONS_P),levs_p_uv(NLVLS_P_UV)
  real :: lons_uvw(NLONS_UVW),levs_w(NLVLS_W)
  real :: lats_uvwsum(NLATS_UVWSUM), lons_uvwsum(NLONS_UVWSUM),levs_uvwsum(NLVLS_UVWSUM)

  integer :: lon_varid, lat_varid, time_series_varid
  integer :: lon_varid_p, lat_varid_p, time_series_varid_p
  integer :: lon_varid_u, lat_varid_u, time_series_varid_u
  integer :: lon_varid_v, lat_varid_v, time_series_varid_v
  integer :: lon_varid_w, lat_varid_w, time_series_varid_w
  integer :: lon_varid_uvwsum, lat_varid_uvwsum, time_series_varid_uvwsum

  ! We will create netCDF variables for p/uvw/uvwsum
  character (len = *), parameter :: PRES_NAME="pressure"
  character (len = *), parameter :: VEL_X_NAME="velocity_x"
  character (len = *), parameter :: VEL_Y_NAME="velocity_y"
  character (len = *), parameter :: VEL_Z_NAME="velocity_z"
  character (len = *), parameter :: VELSUM_X_NAME="velocity_sum_x"
  character (len = *), parameter :: VELSUM_Y_NAME="velocity_sum_y"
  character (len = *), parameter :: VELSUM_Z_NAME="velocity_sum_z"

  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PRES_UNITS = "hPa"
  character (len = *), parameter :: VEL_UNITS = "m/s"
  character (len = *), parameter :: LAT_UNITS = "degrees_north" ! 35° 0' 0" N / 135° 45' 0" E 
  ! 35.0, 135.7 - 34.9, 135.8 => .1 degree is 150, so we need to do 34.9+0.1*lat/NLATS, 135.7+0.1*lon/NLONS
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: TIME_SERIES_UNITS = "hours since 2014-08-05 16:00:00"
  
  
  save

  contains

  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  subroutine init_netcdf_file()       
      ! Loop indices
      integer :: lat, lon,t
  
     ! Create the file.
print *,'nf90_create()'
      call check( nf90_create( FILE_NAME, NF90_CLOBBER,  ncid) )

      call check( nf90_create( FILE_NAME_P, NF90_CLOBBER,  ncid_p) )
      call check( nf90_create( FILE_NAME_U, NF90_CLOBBER,  ncid_u) )
      call check( nf90_create( FILE_NAME_V, NF90_CLOBBER,  ncid_v) )
      call check( nf90_create( FILE_NAME_W, NF90_CLOBBER,  ncid_w) )
      call check( nf90_create( FILE_NAME_UVWSUM, NF90_CLOBBER,  ncid_uvwsum) )

      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
print *,'nf90_def_dim()'
      call check( nf90_def_dim(ncid, LVL_NAME, NLVLS_P_UV, lvl_dimid) )
      call check( nf90_def_dim(ncid, LAT_NAME, NLATS_P_UVW, lat_dimid) )
      call check( nf90_def_dim(ncid, LON_NAME, NLONS_P, lon_dimid) )
      call check( nf90_def_dim(ncid, TIME_SERIES_NAME, NTIMESTEPS, time_series_dimid) ) !  NF90_UNLIMITED -> WV: can't get that to work

      call check( nf90_def_dim(ncid_p, LVL_NAME, NLVLS_P_UV, lvl_dimid_p) )
      call check( nf90_def_dim(ncid_p, LAT_NAME, NLATS_P_UVW, lat_dimid_p) )
      call check( nf90_def_dim(ncid_p, LON_NAME, NLONS_P, lon_dimid_p) )
      call check( nf90_def_dim(ncid_p, TIME_SERIES_NAME, NTIMESTEPS, time_series_dimid_p) ) !  NF90_UNLIMITED

      call check( nf90_def_dim(ncid_u, LVL_NAME, NLVLS_P_UV, lvl_dimid_u) )
      call check( nf90_def_dim(ncid_u, LAT_NAME, NLATS_P_UVW, lat_dimid_u) )
      call check( nf90_def_dim(ncid_u, LON_NAME, NLONS_UVW, lon_dimid_u) )
      call check( nf90_def_dim(ncid_u, TIME_SERIES_NAME, NTIMESTEPS, time_series_dimid_u) ) !  NF90_UNLIMITED

      call check( nf90_def_dim(ncid_v, LVL_NAME, NLVLS_P_UV, lvl_dimid_v) )
      call check( nf90_def_dim(ncid_v, LAT_NAME, NLATS_P_UVW, lat_dimid_v) )
      call check( nf90_def_dim(ncid_v, LON_NAME, NLONS_UVW, lon_dimid_v) )
      call check( nf90_def_dim(ncid_v, TIME_SERIES_NAME, NTIMESTEPS, time_series_dimid_v) ) !  NF90_UNLIMITED

      call check( nf90_def_dim(ncid_w, LVL_NAME, NLVLS_W, lvl_dimid_w) )
      call check( nf90_def_dim(ncid_w, LAT_NAME, NLATS_P_UVW, lat_dimid_w) )
      call check( nf90_def_dim(ncid_w, LON_NAME, NLONS_UVW, lon_dimid_w) )
      call check( nf90_def_dim(ncid_w, TIME_SERIES_NAME, NTIMESTEPS, time_series_dimid_w) ) !  NF90_UNLIMITED

      call check( nf90_def_dim(ncid_uvwsum, LVL_NAME, NLVLS_UVWSUM, lvl_dimid_uvwsum) )
      call check( nf90_def_dim(ncid_uvwsum, LAT_NAME, NLATS_UVWSUM, lat_dimid_uvwsum) )
      call check( nf90_def_dim(ncid_uvwsum, LON_NAME, NLONS_UVWSUM, lon_dimid_uvwsum) )
      call check( nf90_def_dim(ncid_uvwsum, TIME_SERIES_NAME, NTIMESTEPS, time_series_dimid_uvwsum) ) !  NF90_UNLIMITED

      ! Define the coordinate variables. We will only define coordinate
      ! variables for lat and lon.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (lat_dimid) and
      ! similarly for (lon_dimid).
print*,'nf90_def_var()'
      call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
      call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
      call check( nf90_def_var(ncid, TIME_SERIES_NAME, NF90_DOUBLE, time_series_dimid,  time_series_varid) )

      call check( nf90_def_var(ncid_p, LAT_NAME, NF90_REAL, lat_dimid_p, lat_varid_p) )
      call check( nf90_def_var(ncid_p, LON_NAME, NF90_REAL, lon_dimid_p, lon_varid_p) )
      call check( nf90_def_var(ncid_p, TIME_SERIES_NAME, NF90_DOUBLE, time_series_dimid_p,  time_series_varid_p) )

      call check( nf90_def_var(ncid_u, LAT_NAME, NF90_REAL, lat_dimid_u, lat_varid_u) )
      call check( nf90_def_var(ncid_u, LON_NAME, NF90_REAL, lon_dimid_u, lon_varid_u) )
      call check( nf90_def_var(ncid_u, TIME_SERIES_NAME, NF90_DOUBLE, time_series_dimid_u,  time_series_varid_u) )

      call check( nf90_def_var(ncid_v, LAT_NAME, NF90_REAL, lat_dimid_v, lat_varid_v) )
      call check( nf90_def_var(ncid_v, LON_NAME, NF90_REAL, lon_dimid_v, lon_varid_v) )
      call check( nf90_def_var(ncid_v, TIME_SERIES_NAME, NF90_DOUBLE, time_series_dimid_v,  time_series_varid_v) )

      call check( nf90_def_var(ncid_w, LAT_NAME, NF90_REAL, lat_dimid_w, lat_varid_w) )
      call check( nf90_def_var(ncid_w, LON_NAME, NF90_REAL, lon_dimid_w, lon_varid_w) )
      call check( nf90_def_var(ncid_w, TIME_SERIES_NAME, NF90_DOUBLE, time_series_dimid_w,  time_series_varid_w) )

      call check( nf90_def_var(ncid_uvwsum, LAT_NAME, NF90_REAL, lat_dimid_uvwsum, lat_varid_uvwsum) )
      call check( nf90_def_var(ncid_uvwsum, LON_NAME, NF90_REAL, lon_dimid_uvwsum, lon_varid_uvwsum) )
      call check( nf90_def_var(ncid_uvwsum, TIME_SERIES_NAME, NF90_DOUBLE, time_series_dimid_uvwsum,  time_series_varid_uvwsum) )

print*,'nf90_put_att()'

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
      call check( nf90_put_att(ncid, time_series_varid, UNITS, TIME_SERIES_UNITS) )

      call check( nf90_put_att(ncid_p, lat_varid_p, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid_p, lon_varid_p, UNITS, LON_UNITS) )
      call check( nf90_put_att(ncid_p, time_series_varid_p, UNITS, TIME_SERIES_UNITS) )

      call check( nf90_put_att(ncid_u, lat_varid_u, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid_u, lon_varid_u, UNITS, LON_UNITS) )
      call check( nf90_put_att(ncid_u, time_series_varid_u, UNITS, TIME_SERIES_UNITS) )

      call check( nf90_put_att(ncid_w, lat_varid_w, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid_w, lon_varid_w, UNITS, LON_UNITS) )
      call check( nf90_put_att(ncid_w, time_series_varid_w, UNITS, TIME_SERIES_UNITS) )

      call check( nf90_put_att(ncid_v, lat_varid_v, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid_v, lon_varid_v, UNITS, LON_UNITS) )
      call check( nf90_put_att(ncid_v, time_series_varid_v, UNITS, TIME_SERIES_UNITS) )

      call check( nf90_put_att(ncid_uvwsum, lat_varid_uvwsum, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid_uvwsum, lon_varid_uvwsum, UNITS, LON_UNITS) )
      call check( nf90_put_att(ncid_uvwsum, time_series_varid_uvwsum, UNITS, TIME_SERIES_UNITS) )


      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ lon_dimid, lat_dimid, lvl_dimid, time_series_dimid /)

      dimids_p = (/ lon_dimid_p, lat_dimid_p, lvl_dimid_p, time_series_dimid_p /)
      dimids_u = (/ lon_dimid_u, lat_dimid_u, lvl_dimid_u, time_series_dimid_u /)
      dimids_v = (/ lon_dimid_v, lat_dimid_v, lvl_dimid_v, time_series_dimid_v /)
      dimids_w = (/ lon_dimid_w, lat_dimid_w, lvl_dimid_w, time_series_dimid_w /)
      dimids_uvwsum = (/ lon_dimid_uvwsum, lat_dimid_uvwsum, lvl_dimid_uvwsum, time_series_dimid_uvwsum /)
print*,'nf90_def_var(p,u,v,w,uvwsum)'
      ! Define the netCDF variables for the pressure and velocity data.
      call check( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )

      call check( nf90_def_var(ncid_p, PRES_NAME, NF90_REAL, dimids_p, pres_varid_p) )
      call check( nf90_def_var(ncid_u, VEL_X_NAME, NF90_REAL, dimids_u, vel_x_varid_u) )
      call check( nf90_def_var(ncid_v, VEL_Y_NAME, NF90_REAL, dimids_v, vel_y_varid_v) )
      call check( nf90_def_var(ncid_w, VEL_Z_NAME, NF90_REAL, dimids_w, vel_z_varid_w) )
      call check( nf90_def_var(ncid_uvwsum, VELSUM_X_NAME, NF90_REAL, dimids, velsum_x_varid_uvwsum) )
      call check( nf90_def_var(ncid_uvwsum, VELSUM_Y_NAME, NF90_REAL, dimids, velsum_y_varid_uvwsum) )
      call check( nf90_def_var(ncid_uvwsum, VELSUM_Z_NAME, NF90_REAL, dimids, velsum_z_varid_uvwsum) )
print*,'nf90_put_att(UNITS)'
      ! Assign units attributes to the netCDF variables.
      call check( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )

      call check( nf90_put_att(ncid_p, pres_varid_p, UNITS, PRES_UNITS) )
      call check( nf90_put_att(ncid_u, vel_x_varid_u, UNITS, VEL_UNITS) )
      call check( nf90_put_att(ncid_v, vel_y_varid_v, UNITS, VEL_UNITS) )
      call check( nf90_put_att(ncid_w, vel_z_varid_w, UNITS, VEL_UNITS) )
      call check( nf90_put_att(ncid_uvwsum, velsum_x_varid_uvwsum, UNITS, VEL_UNITS) )
      call check( nf90_put_att(ncid_uvwsum, velsum_y_varid_uvwsum, UNITS, VEL_UNITS) )
      call check( nf90_put_att(ncid_uvwsum, velsum_z_varid_uvwsum, UNITS, VEL_UNITS) )

print*,'nf90_enddef()'
      ! End define mode.
      call check( nf90_enddef(ncid) )

      call check( nf90_enddef(ncid_p) )
      call check( nf90_enddef(ncid_u) )
      call check( nf90_enddef(ncid_v) )
      call check( nf90_enddef(ncid_w) )
      call check( nf90_enddef(ncid_uvwsum) )

print*,'Create lat/lon/times arrays'

       ! Create arrays for lat and long
       ! Kyoto/Uji from northwest to southeast 35.0, 135.7 - 34.9, 135.8 => .1 degree is 150, so we need to do 34.9+0.1*lat/NLATS, 135.7+0.1*lon/NLONS
       do lat = 1, NLATS_P_UVW
          lats_p_uvw(lat) = lat ! 34.9+0.1*lat/NLATS
       end do
       do lon = 1, NLONS_P
          lons_p(lon) = lon ! 135.7+0.1*lon/NLONS
       end do
       do lon = 1, NLONS_UVW
          lons_uvw(lon) = lon ! 135.7+0.1*lon/NLONS
       end do
       do lat = 1, NLATS_UVWSUM
          lats_uvwsum(lat) = lat ! 34.9+0.1*lat/NLATS
       end do
       do lon = 1, NLONS_UVWSUM
          lons_uvwsum(lon) = lon ! 135.7+0.1*lon/NLONS
       end do

       do t = 1, NTIMESTEPS
          times(t) = t
       end do

print*,'nf90_put_var()'
      ! Write the coordinate variable data. This will put the latitudes
      ! and longitudes of our data grid into the netCDF file.
      call check( nf90_put_var(ncid, lat_varid, lats_p_uvw) )
      call check( nf90_put_var(ncid, lon_varid, lons_p) )
      call check( nf90_put_var(ncid, time_series_varid, times) )
print*,'nf90_put_var(p)'
      call check( nf90_put_var(ncid_p, lat_varid_p, lats_p_uvw) )
      call check( nf90_put_var(ncid_p, lon_varid_p, lons_p) )
      call check( nf90_put_var(ncid_p, time_series_varid_p, times) )
print*,'nf90_put_var(u)'
      call check( nf90_put_var(ncid_u, lat_varid_u, lats_p_uvw) )
      call check( nf90_put_var(ncid_u, lon_varid_u, lons_uvw) )
      call check( nf90_put_var(ncid_u, time_series_varid_u, times) )
print*,'nf90_put_var(v)'
      call check( nf90_put_var(ncid_v, lat_varid_v, lats_p_uvw) )
      call check( nf90_put_var(ncid_v, lon_varid_v, lons_uvw) )
      call check( nf90_put_var(ncid_v, time_series_varid_v, times) )
print*,'nf90_put_var(w)'
      call check( nf90_put_var(ncid_w, lat_varid_w, lats_p_uvw) )
      call check( nf90_put_var(ncid_w, lon_varid_w, lons_uvw) )
      call check( nf90_put_var(ncid_w, time_series_varid_w, times) )
print*,'nf90_put_var(uvwsum)'
      call check( nf90_put_var(ncid_uvwsum, lat_varid_uvwsum, lats_uvwsum) )
      call check( nf90_put_var(ncid_uvwsum, lon_varid_uvwsum, lons_uvwsum) )
      call check( nf90_put_var(ncid_uvwsum, time_series_varid_uvwsum, times) )

print*,'Create count*/start'

      count = (/ NLONS_P, NLATS_P_UVW, NLVLS_P_UV, 1 /)

      count_p = (/ NLONS_P, NLATS_P_UVW, NLVLS_P_UV, 1 /)
      count_u = (/ NLONS_UVW, NLATS_P_UVW, NLVLS_P_UV, 1 /)
      count_v = (/ NLONS_UVW, NLATS_P_UVW, NLVLS_P_UV, 1 /)
      count_w = (/ NLONS_UVW, NLATS_P_UVW, NLVLS_W, 1 /)
      count_uvwsum = (/ NLONS_UVWSUM, NLATS_UVWSUM, NLVLS_UVWSUM, 1 /)

      start = (/ 1, 1, 1, 1 /)


!      print *, 'ncid init: ',ncid,'pres_varid: ',pres_varid
    
  end subroutine init_netcdf_file

  subroutine write_to_netcdf_file(p,u,v,w,usum,vsum,wsum,n)
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1), intent(InOut)  :: p
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: usum
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: vsum
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: wsum
        integer, intent(In) :: n

 ! The start and count arrays will tell the netCDF library where to
  ! write our data.

  ! Write the pretend data. This will write our surface pressure and
  ! surface velerature data. The arrays only hold one timestep worth
  ! of data. We will just rewrite the same data for each timestep. In
  ! a real :: application, the data would change between timesteps.
     start(4) = n
     print *, 'ncid: ',ncid,'pres_varid: ',pres_varid

     call check( nf90_put_var(ncid, pres_varid, p,  start, count) )

     call check( nf90_put_var(ncid_p, pres_varid, p,  start, count_p) )

     print *, 'ncid_u: ',ncid_u,'vel_x_varid_u: ',vel_x_varid_u
     call check( nf90_put_var(ncid_u, vel_x_varid_u, u, start, count_u) )
     call check( nf90_put_var(ncid_v, vel_y_varid_v, v, start, count_v) )
     call check( nf90_put_var(ncid_w, vel_z_varid_w, w, start, count_w) )

     call check( nf90_put_var(ncid_uvwsum, velsum_x_varid_uvwsum, usum, start, count_uvwsum) )
     call check( nf90_put_var(ncid_uvwsum, velsum_y_varid_uvwsum, vsum, start, count_uvwsum) )
     call check( nf90_put_var(ncid_uvwsum, velsum_z_varid_uvwsum, wsum, start, count_uvwsum) )

    end subroutine write_to_netcdf_file

    subroutine close_netcdf_file()
  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are really written to disk.
  call check( nf90_close(ncid) )
#ifdef VERBOSE   
  print *,"*** SUCCESS writing file ", FILE_NAME
#endif
  call check( nf90_close(ncid_p) )
#ifdef VERBOSE
  print *,"*** SUCCESS writing file ", FILE_NAME_p
#endif
  call check( nf90_close(ncid_u) )
#ifdef VERBOSE
  print *,"*** SUCCESS writing file ", FILE_NAME_U
#endif
  call check( nf90_close(ncid_v) )
#ifdef VERBOSE
  print *,"*** SUCCESS writing file ", FILE_NAME_V
#endif
  call check( nf90_close(ncid_W) )
#ifdef VERBOSE
  print *,"*** SUCCESS writing file ", FILE_NAME_W
#endif
  call check( nf90_close(ncid_uvwsum) )
#ifdef VERBOSE
  print *,"*** SUCCESS writing file ", FILE_NAME_UVWSUM
#endif    
    end subroutine close_netcdf_file

end module module_LES_write_netcdf
