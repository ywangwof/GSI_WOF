subroutine read_radar_wind_netcdf(nread,ndata,nodata,infile,lunout,obstype,twind,sis,hgtl_full,nobs)
!$$$   subprogram documentation block
!                .      .    .                                       .
!   subprogram: read_radar_wind_netcdf        read radar radial velocity in DART-like netcdf format
!   
! abstract: Read and process radar radial velocity observations in DART-like netcdf format.
!
! program history log:
!   2016-02-14  Johnson, Y. Wang, X. Wang - modify read_radar.f90 to read Vr in DART-like netcdf format
!                                           in collaboration with Carley, POC: xuguang.wang@ou.edu
!   2019-04-18  Thomas Jones	= Major changes to netcdf read portion for 2019 VR format                               
!
! program history log:
!           
!   input argument list:
!     infile   - file from which to read data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!
!   output argument list:
!     nread    - number of radar reflectivity observations read
!     ndata    - number of radar reflectivity observations retained for further processing
!     nodata   - number of radar reflectivity observations retained for further processing
!     sis      - satellite/instrument/sensor indicator
!
! Variable Definitions:
!
!  cdata_all - real - dim(maxdat,maxobs) - array holding all data for assimilation
!  dlat - real - grid relative latitude of observation (grid units)
!  dlon - real - grid relative longitude of observation (grid units)
!  maxobs - int - max number of obs converted to no precip observations
!  num_m2nopcp -int - number of missing obs 
!  num_missing - int - number of missing observations
!  num_noise - int - number of rejected noise observations
!  num_nopcp - int - number of noise obs converted to no precip observations
!  numbadtime - int - number of elevations outside time window
!  outside - logical - if observations are outside the domain -> true
!  rmins_an - real - analysis time from reference date (hour)
!  rmins_ob - real -  observation time from reference date (hour)
!  rstation_id - real - radar station id
!  thisazimuthr - real - 90deg minues the actual azimuth and converted to radians
!  thiserr - real - observation error
!  thislat - real - latitude of observation, point
!  thislon - real - longitude of observation, point
!  thisrange - real - range of observation from radar
!  thishgt - real - observation height, point
!  this_stahgt - real - radar station height (meters about sea level)
!  this_staid - char - radar station id
!  thistiltr - real- radar tilt angle (radians)
!  timeb - real - obs time (analyis relative minutes)
!  vrQC - real - radial velocity observation 
!  vr_err - real - observation error of radial velocity
!  height - real - height of observation
!  lon    - real - longitude of observation
!  lat    - real - latitude of observation
!  utime  - real - time for each observation point
!  sta_info - real - radar site information, including lon, lat, height of radar sation, 
!                    and Nyquist etc.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block

  use kinds, only: r_kind,r_double,i_kind
  use constants, only: zero,half,one,two,deg2rad,rearth,rad2deg, &
                       one_tenth,r1000,r60,r60inv,r100,r400,grav_equator, &
                       eccentricity,somigliana,grav_ratio,grav,semi_major_axis,flattening 
  use gridmod, only: regional,tll2xy,nsig,nlat,nlon
  use obsmod, only: iadatemn,mintiltvr,maxtiltvr,minobrangevr,maxobrangevr,rmesh_vr,zmesh_vr
  use obsmod,only: radar_no_thinning,missing_to_nopcp
  use convinfo, only: nconvtype,ctwind,cgross,icuse,ioctype
  use convthin, only: make3grids,map3grids,del3grids,use_all
  use deter_sfc_mod, only: deter_sfc2,deter_zsfc_model
  use jfunc, only: miter
  use mpimod, only: npe
  use netcdf

  implicit none

  include 'netcdf.inc'

! Declare passed variables
  character(len=*),intent(in   ) :: obstype,infile
  character(len=*),intent(in   ) :: sis
  real(r_kind)    ,intent(in   ) :: twind
  integer(i_kind) ,intent(in   ) :: lunout
  integer(i_kind) ,intent(inout) :: nread,ndata,nodata
  real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full
  integer(i_kind) ,dimension(npe),intent(inout) :: nobs

! Declare local parameters
  real(r_kind),parameter:: r6 = 6.0_r_kind
  real(r_kind),parameter:: r360=360.0_r_kind
  integer(i_kind),parameter:: maxdat=24_i_kind, ione = 1_i_kind, izero = 0_i_kind         ! Used in generating cdata array
  character(len=4), parameter :: radid = 'XXXX'
  
!-----------------------
!  NETCDF File Varaibles
!-----------------------
  real(r_kind), allocatable, dimension(:)       :: vrQC, vr_err, height,lon, lat, utime, &
                                                   dir1, dir2, dir3, radar_lon, radar_lat, radar_hgt, &
                                                   nyquist

!------------------
!  NETCDF-RELATED
!------------------
  INTEGER(i_kind)   :: ncdfID, status
  INTEGER(i_kind)   :: varID, dimid
  logical           :: if_input_exist
  integer(i_kind)   :: ivar, sec70, length, nn

!--Counters for diagnostics
 integer(i_kind) :: num_missing=izero,num_nopcp=izero, &      !counts 
                    numbadtime=izero, &
                    numoutside=izero, &    
                    num_m2nopcp=izero, &
                    num_noise=izero,num_limmax=izero     
  

  integer(i_kind) :: ithin,zflag,nlevz,icntpnt,klon1,klat1,kk,klatp1,klonp1
  real(r_kind) :: rmesh,xmesh,zmesh,dx,dy,dx1,dy1,w00,w01,w10,w11
  real(r_kind), allocatable, dimension(:) :: zl_thin
  real(r_kind),dimension(nsig):: hges,zges
  real(r_kind) sin2,termg,termr,termrg,zobs,hgt
  integer(i_kind) ntmp,iout,iiout,ntdrvr_thin2
  real(r_kind) crit1,timedif
  real(r_kind),parameter:: r16000 = 16000.0_r_kind
  logical :: luse
  integer(i_kind) maxout,maxdata
  integer(i_kind),allocatable,dimension(:):: isort
       
  !--General declarations
  integer(i_kind) :: ierror,i,j,k,v,na,nb,nelv,nvol, &
                     ikx,mins_an,mins_ob
  integer(i_kind) :: maxobs,nchanl,ilat,ilon,idomsfc
  
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt,&
                  ff10,sfcr,skint,zsges, &
                  dlon_radar,dlat_radar,errmax,errmin,error
  real(r_kind) :: celev,selev,thisazimuthr,t4dv, &
                  dlat,dlon,thiserr,thislon,thislat, &
                  timeb, twindm
  real(r_kind) :: rmins_an,rmins_ob                                                     
  real(r_kind),allocatable,dimension(:,:):: cdata_all
  real(r_double) rstation_id
  
  character(8) cstaid
  character(4) this_staid
  equivalence (this_staid,cstaid)
  equivalence (cstaid,rstation_id)

  logical      :: outside
    
  real(r_kind) :: minobrange,maxobrange,mintilt,maxtilt


  minobrange=minobrangevr
  maxobrange=maxobrangevr
  mintilt=mintiltvr
  maxtilt=maxtiltvr


!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS                                   !
!--------------------------------------------------------------------------------------!
   

 ithin=-1 !number of obs to keep per grid box DEFAULT: NO THINNING (-1)
 if(radar_no_thinning) then
  ithin=-1
 endif

  ikx=izero
   do i=ione,nconvtype
      if(trim(obstype) == trim(ioctype(i)) .and. abs(icuse(i))== ione) then
         ikx=i
!         radartwindow=ctwind(ikx)*r60         !Time window units converted to minutes 
                                              !  (default setting for cwp within convinfo is 0.05 hours)
         exit                                 !Exit loop when finished with initial convinfo fields     
      else if ( i==nconvtype ) then
         write(6,*) 'READ_VR: ERROR - OBSERVATION TYPE IS NOT PRESENT IN CONVINFO OR USE FLAG IS ZERO'
         write(6,*) 'READ_VR: ABORTTING read_vr.f90 - NO VR OBS READ!'
         return
      endif
  end do

  if (minobrange >= maxobrange) then
    write(6,*) 'MININMUM OB RANGE >= MAXIMUM OB RANGE FOR READING VR - PROGRAM STOPPING FROM READ_RADAR_WIND_NETCDF.F90'
    call stop2(400)
  end if
        
  !-next three values are dummy values for now
  nchanl=izero
  ilon=2_i_kind
  ilat=3_i_kind
  
  maxobs=50000000_i_kind    !value taken from read_radar.f90 

  !--Allocate cdata_all array
   allocate(cdata_all(maxdat,maxobs),isort(maxobs))
   rmesh=rmesh_vr
   zmesh=zmesh_vr


   maxout=0
   maxdata=0
   isort=0
   ntdrvr_thin2=0
   icntpnt=0
   zflag=0

   use_all=.true.
  if (ithin > 0) then
     write(6,*)'READ_RADAR_VR: ithin,rmesh :',ithin,rmesh
     use_all=.false.
     if(zflag == 0)then
        nlevz=nsig
     else
        nlevz=r16000/zmesh
     endif
     xmesh=rmesh
     call make3grids(xmesh,nlevz)
!     call make3grids2(xmesh,nlevz)

     allocate(zl_thin(nlevz))
     if (zflag == 1) then
        do k=1,nlevz
           zl_thin(k)=k*zmesh
        enddo
     endif
     write(6,*)'READ_RADAR_VR: xmesh, zflag, nlevz =', xmesh, zflag, nlevz
  endif
!!end modified for thinning

  ! CHECK IF DATA FILE EXISTS
  length     = len_trim(infile)
  inquire(file=infile(1:length), exist=if_input_exist)
  fileopen: if (if_input_exist) then

  ! OPEN NETCDF FILE
  status = nf90_open(TRIM(infile), NF90_NOWRITE, ncdfID)
  print*, '*** OPENING 88D Radial Velocity  OBS NETCDF FILE: ', infile, status

  ! Get dimension information
  status = nf90_inq_dimid(ncdfID, "index", dimid)
  status = nf90_inquire_dimension(ncdfID, dimid, len = nn)

  print*, 'NUM OBS: ', nn

  ALLOCATE( lat( nn ) )
  ALLOCATE( lon( nn ) )
  ALLOCATE( height( nn ) )   
  ALLOCATE( vrQC( nn ) )
  ALLOCATE( vr_err( nn ) )
  ALLOCATE( utime( nn ) )
  ALLOCATE( radar_lat( nn ) )
  ALLOCATE( radar_lon( nn ) )
  ALLOCATE( radar_hgt( nn ) )
  ALLOCATE( dir1( nn ) )
  ALLOCATE( dir2( nn ) )
  ALLOCATE( dir3( nn ) )
  ALLOCATE( nyquist( nn ) )

   !------------------------
   ! Get useful data arrays
   !-------------------------
   ! LAT
   status = nf90_inq_varid( ncdfID, 'lat', varID )
   status = nf90_get_var( ncdfID, varID, lat )

   ! LON
   status = nf90_inq_varid( ncdfID, 'lon', varID )
   status = nf90_get_var( ncdfID, varID, lon )

   ! HEIGHT (m) 
   status = nf90_inq_varid( ncdfID, 'height', varID )
   status = nf90_get_var( ncdfID, varID, height )

   ! VR VALUE (m / s)
   status = nf90_inq_varid( ncdfID, 'value', varID )
   status = nf90_get_var( ncdfID, varID, vrQC )

   ! VR OBSERVATION ERROR
   status = nf90_inq_varid( ncdfID, 'error_var', varID )
   status = nf90_get_var( ncdfID, varID, vr_err ) 

   ! RADAR LAT
   status = nf90_inq_varid( ncdfID, 'platform_lat', varID )
   status = nf90_get_var( ncdfID, varID, radar_lat )

   ! RADAR LON
   status = nf90_inq_varid( ncdfID, 'platform_lon', varID )
   status = nf90_get_var( ncdfID, varID, radar_lon )

   ! RADAR HEIGHT
   status = nf90_inq_varid( ncdfID, 'platform_hgt', varID )
   status = nf90_get_var( ncdfID, varID, radar_hgt )

   ! RADAR DIR 1
   status = nf90_inq_varid( ncdfID, 'platform_dir1', varID )
   status = nf90_get_var( ncdfID, varID, dir1 )

   ! RADAR DIR 2
   status = nf90_inq_varid( ncdfID, 'platform_dir2', varID )
   status = nf90_get_var( ncdfID, varID, dir2 )

   ! RADAR DIR 3
   status = nf90_inq_varid( ncdfID, 'platform_dir3', varID )
   status = nf90_get_var( ncdfID, varID, dir3 )

   ! NYQUIST VELOCITY
   status = nf90_inq_varid( ncdfID, 'platform_nyquist', varID )
   status = nf90_get_var( ncdfID, varID, nyquist )

   ! TIME
   status = nf90_inq_varid( ncdfID, 'utime', varID )
   status = nf90_get_var( ncdfID, varID, utime )

   ! CLOSE NETCDF FILE
   status = nf90_close( ncdfID )



  !-Obtain analysis time in minutes since reference date

  sec70 = 252460800  ! seconds since from 01/01/1970


  call w3fs21(iadatemn,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  rmins_an=mins_an             !convert to real number
 
  twindm = twind*60.0 !Convert namelist timewindow to minutes from hours
  
  ivar = 2
  
  do i = 1, nn

       rmins_ob = ( utime(i) - sec70 )/60 
       timeb = rmins_ob-rmins_an

       !print*, timeb

       if(abs(timeb) > abs(twindm)) then
         numbadtime=numbadtime+ione
         cycle
       end if

       thishgt = height(i) ! unit : meter
       hgt     = thishgt

       thislon = lon(i)
       thislat = lat(i)
  
       !-Check format of longitude and correct if necessary
                 
       if(thislon>=r360) thislon=thislon-r360
       if(thislon<zero ) thislon=thislon+r360
                 
       !-Convert back to radians                 
         
       thislat = thislat*deg2rad
       thislon = thislon*deg2rad
                 
       !find grid relative lat lon locations of earth lat lon
                 
       call tll2xy(thislon,thislat,dlon,dlat,outside)
       if (outside) then
         numoutside = numoutside+ione
         cycle
       endif
                                           !If observation is outside the domain
                                           ! then cycle, but don't increase range right away.
                                           ! Domain could be rectangular, so ob may be out of
                                           ! range at one end, but not the other.		     					                   		   		   
       
       ! Check if observation error > 0 (TAJ)
       if ( vr_err(i) < 0.0 ) then
         print '(A25, 2F8.2)', 'WARNING: Obs Var < 0 ', vrQC(i), vr_err(i)
         cycle
       endif              
       thiserr = sqrt(vr_err(i))

          
       ! Check if dirs are good
       if ( (dir1(i) /= dir1(i)) .or. (dir2(i) /= dir2(i)) .or. (dir3(i) /= dir3(i)) ) then
         print '(A25, 4F8.2)', 'WARNING: NaN angles, skip ', vrQC(i), dir1(i), dir2(i), dir3(i)
         cycle
       endif

       ! Get model terrain at radar station location
       ! If radar station is outside of grid, does not mean the
       !    radar obs are outside the grid - therefore no need to
       !    cycle azms.         
       radar_lon(i)=deg2rad*radar_lon(i)
       radar_lat(i)=deg2rad*radar_lat(i)
       call tll2xy(radar_lon(i),radar_lat(i),dlon_radar,dlat_radar,outside)
       call deter_zsfc_model(dlat_radar,dlon_radar,zsges)

       !  Determines land surface type based on surrounding land
       !    surface types
       t4dv=timeb*r60inv
       call deter_sfc2(thislat,thislon,t4dv,idomsfc,skint,ff10,sfcr)

       nread = nread + ione

!####################       Data thinning       ###################
       icntpnt=icntpnt+1
!reachded ok
           if(ithin > 0)then

              if(zflag == 0)then
                 klon1= int(dlon);  klat1= int(dlat)
                 dx   = dlon-klon1; dy   = dlat-klat1
                 dx1  = one-dx;     dy1  = one-dy
                 w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy
 
                 klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
                 if (klon1==0) klon1=nlon
                 klatp1=min(nlat,klat1+1); klonp1=klon1+1
                 if (klonp1==nlon+1) klonp1=1
                 do kk=1,nsig
                    hges(kk)=w00*hgtl_full(klat1 ,klon1 ,kk) +  &
                             w10*hgtl_full(klatp1,klon1 ,kk) + &
                             w01*hgtl_full(klat1 ,klonp1,kk) + &
                             w11*hgtl_full(klatp1,klonp1,kk)
                 end do
                 sin2  = sin(thislat)*sin(thislat)
                 termg = grav_equator * &
                    ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
                 termr = semi_major_axis /(one + flattening + grav_ratio -  &
                    two*flattening*sin2)
                 termrg = (termg/grav)*termr
                 do kk=1,nsig
                    zges(kk) = (termr*hges(kk)) / (termrg-hges(kk))
                    zl_thin(kk)=zges(kk)
                 end do
              endif

              zobs = hgt


              ntmp=ndata  ! counting moved to map3gridS
              timedif=abs(t4dv) !don't know about this
              crit1 = timedif/r6+half
 
              call map3grids(1,zflag,zl_thin,nlevz,thislat,thislon,&
                 zobs,crit1,ndata,iout,icntpnt,iiout,luse,.false.,.false.)


              maxout=max(maxout,iout)
              maxdata=max(maxdata,ndata)

              if (.not. luse) then
                 ntdrvr_thin2=ntdrvr_thin2+1
                 cycle
              endif
              if(iiout > 0) isort(iiout)=0
              if (ndata > ntmp) then
                 nodata=nodata+1
              endif
              isort(icntpnt)=iout
           else  ! NO THINNING
              ndata =ndata+1
              nodata=nodata+1
              iout=ndata
              isort(icntpnt)=iout
           endif

           !!end modified for thinning

            cdata_all(1,iout) = thiserr                       ! reflectivity obs error (dB) - inflated/adjusted
            cdata_all(2,iout) = dlon                          ! grid relative longitude
            cdata_all(3,iout) = dlat                          ! grid relative latitude
            cdata_all(4,iout) = thishgt                       ! obs absolute height (m)
            cdata_all(5,iout) = vrQC(i)                       ! radar radial velocity
            cdata_all(6,iout) = dir1(i)                       ! sin(azimuth)*cos(elevation)
            cdata_all(7,iout) = timeb*r60inv                  ! obs time (analyis relative hour)
            cdata_all(8,iout) = ikx                           ! type		   
            cdata_all(9,iout) = dir2(i)                       ! cos(azimuth)*cos(elevation)
            cdata_all(10,iout)= radar_hgt(i)                  ! station elevation (m)
            cdata_all(11,iout)= rstation_id                   ! station id
            cdata_all(12,iout)= icuse(ikx)                    ! usage parameter
            cdata_all(13,iout)= idomsfc                       ! dominate surface type
            cdata_all(14,iout)= skint                         ! skin temperature
            cdata_all(15,iout)= ff10                          ! 10 meter wind factor
            cdata_all(16,iout)= sfcr                          ! surface roughness
            cdata_all(17,iout)= thislon*rad2deg               ! earth relative longitude (degrees)
            cdata_all(18,iout)= thislat*rad2deg               ! earth relative latitude (degrees)
            cdata_all(19,iout)=thisrange/1000_r_kind          ! range from radar in km (used to estimate beam spread)
            cdata_all(20,iout)=zsges                          !  model elevation at radar site
            cdata_all(21,iout)=thiserr
            cdata_all(22,iout)=two                            ! Level 2 data
            cdata_all(23,iout) = nyquist(i)                   ! nyq vel (m/s)
            cdata_all(24,iout) = dir3(i)                      ! sin(elevation)

     end do    ! k

  987 continue      
  if (.not. use_all) then 
     deallocate(zl_thin) 
     call del3grids
  endif
!---all looping done now print diagnostic output

!  write(6,*)'READ_RADAR_WIND: Reached eof on radar radial velocity file'
  write(6,*)'READ_RADAR_WIND: # volumes in input file             =',nvol
  write(6,*)'READ_RADAR_WIND: # read in obs. number               =',nn
  write(6,*)'READ_RADAR_WIND: # observations outside domain       =',numoutside
  write(6,*)'READ_RADAR_WIND: # observations outside time window  =',numbadtime
  write(6,*)'LUN NUMBER: =', lunout

!---Write observation to scratch file---!
  call count_obs(ndata,maxdat,ilat,ilon,cdata_all,nobs) 
  write(lunout) obstype,sis,maxdat,nchanl,ilat,ilon
  write(lunout) ((cdata_all(k,i),k=ione,maxdat),i=ione,ndata)
 
  print*, 'FINISH READ RADIAL VELOCITY'
  
  !---------------DEALLOCATE ARRAYS-------------! 
 deallocate(cdata_all)
 deallocate(lat, lon, height, vrQC, vr_err,utime, dir1, dir2, dir3)
 deallocate(nyquist, radar_lat, radar_lon, radar_hgt)

 else  !fileopen
  write(6,*) 'READ_RADAR_WIND_NETCDF: ERROR OPENING RADAR RADIAL VELOCITY FILE: ',trim(infile),' IOSTAT ERROR: ',ierror, ' SKIPPING...'
 end if fileopen

314 continue

end subroutine read_radar_wind_netcdf
