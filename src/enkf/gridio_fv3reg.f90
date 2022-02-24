module gridio

   !========================================================================

   !$$$ Module documentation block
   !
   ! This module contains various routines to ingest and update
   ! variables from Weather Research and Forecasting (WRF) model Advanced
   ! Research WRF (ARW) and Non-hydrostatic Mesoscale Model (NMM) dynamical
   ! cores which are required by the Ensemble Kalman Filter (ENKF) currently
   ! designed for operations within the National Centers for Environmental
   ! Prediction (NCEP) Global Forecasting System (GFS)
   !
   ! prgmmr: Winterbottom        org: ESRL/PSD1       date: 2011-11-30
   !
   ! program history log:
   !
   !   2011-11-30 Winterbottom - Initial version.
   !
   !   2019-01- Ting  --  modified for fv3sar
   !   2021-04     Y. Wang & X. Wang - changes for radar DA
   ! attributes:
   !   language:  f95
   !
   !$$$

   !=========================================================================
   ! Define associated modules
   use gridinfo, only:  npts
   use kinds,    only: r_double, r_kind, r_single, i_kind
   use mpisetup, only: nproc
   use netcdf_io
   use params,   only: nlevs, cliptracers, datapath, arw, nmm, datestring
   use params,   only: nx_res,ny_res,nlevs,ntiles
   use params,   only:  pseudo_rh
   use mpeu_util, only: getindex
   use read_fv3regional_restarts,only:read_fv3_restart_data1d,read_fv3_restart_data2d
   use read_fv3regional_restarts,only:read_fv3_restart_data3d,read_fv3_restart_data4d
   use netcdf_mod,only: nc_check

   implicit none

   !-------------------------------------------------------------------------
   ! Define all public subroutines within this module
   private
   public :: readgriddata,readgriddata_pnc
   public :: writegriddata,writegriddata_pnc
   public :: writeincrement, writeincrement_pnc

   !-------------------------------------------------------------------------

contains
   ! Generic WRF read routine, calls ARW-WRF or NMM-WRF
   subroutine readgriddata(nanal1,nanal2,vars3d,vars2d,n3d,n2d,levels,ndim,ntimes,fileprefixes,filesfcprefixes,reducedgrid,vargrid,qsat)
      use constants, only:zero,one,half,fv, max_varname_length
      use gridinfo,only: eta1_ll
      use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
      use netcdf, only: nf90_inq_dimid,nf90_inq_varid
      use netcdf, only: nf90_nowrite,nf90_write,nf90_inquire,nf90_inquire_dimension
      implicit none
      integer, intent(in) :: nanal1,nanal2, n2d, n3d, ndim, ntimes
      character(len=max_varname_length), dimension(n2d), intent(in) :: vars2d
      character(len=max_varname_length), dimension(n3d), intent(in) :: vars3d
      integer, dimension(0:n3d), intent(in)        :: levels
      character(len=120), dimension(7), intent(in) :: fileprefixes
      character(len=120), dimension(7), intent(in) :: filesfcprefixes
      logical, intent(in) :: reducedgrid

      real(r_single), dimension(npts,ndim,ntimes,nanal2-nanal1+1),  intent(out) :: vargrid
      real(r_double), dimension(npts,nlevs,ntimes,nanal2-nanal1+1), intent(out) :: qsat



      ! Define local variables
      character(len=500) :: filename
      character(len=:),allocatable :: fv3filename,fv3filename_tracer,fv3filename_phyvar
      character(len=7)   :: charnanal
      integer(i_kind) file_id,file_id_t,file_id_p
      real(r_single), dimension(:,:,:), allocatable ::workvar3d,uworkvar3d,&
         vworkvar3d,tvworkvar3d,tsenworkvar3d,&
         workprsi,qworkvar3d,workvar3d_tmp
      real(r_double),dimension(:,:,:),allocatable:: qsatworkvar3d
      real(r_single), dimension(:,:),   allocatable ::pswork

      ! Define variables required for netcdf variable I/O
      character(len=12) :: varstrname


      character(len=1) char_tile
      character(len=24),parameter :: myname_ = 'fv3: getgriddata'

      ! Define counting variables
      integer :: nlevsp1
      integer :: i,j, k,nn,ntile,nn_tile0, nb,nanal,ne
      integer :: u_ind, v_ind, tv_ind,tsen_ind, q_ind, oz_ind
      integer :: w_ind, ql_ind,qr_ind,qi_ind,qs_ind,qg_ind,qh_ind,dbz_ind
      integer :: ps_ind, sst_ind
      integer :: tmp_ind
      logical :: ice

      !======================================================================
      write (6,*)"The input fileprefix, reducedgrid are not used in the current implementation"
      nlevsp1=nlevs+1
      u_ind   = getindex(vars3d, 'u')   !< indices in the state var arrays
      v_ind   = getindex(vars3d, 'v')   ! U and V (3D)
      tv_ind  = getindex(vars3d, 't')  ! Tv (3D)
      q_ind   = getindex(vars3d, 'q')   ! Q (3D)
      oz_ind  = getindex(vars3d, 'oz')  ! Oz (3D)
      tsen_ind = getindex(vars3d, 'tsen') !sensible T (3D)
!    prse_ind = getindex(vars3d, 'prse') ! pressure
      w_ind   = getindex(vars3d, 'w')   ! W (3D)
      ql_ind   = getindex(vars3d, 'ql')   ! QL (3D)
      qr_ind   = getindex(vars3d, 'qr')   ! QR (3D)
      qi_ind   = getindex(vars3d, 'qi')   ! QI (3D)
      qs_ind   = getindex(vars3d, 'qs')   ! QS (3D)
      qg_ind   = getindex(vars3d, 'qg')   ! QG (3D)
      qh_ind   = getindex(vars3d, 'qh')   ! QH (3D)
      dbz_ind   = getindex(vars3d, 'dbz')   ! DBZ (3D)

      ps_ind  = getindex(vars2d, 'ps')  ! Ps (2D)
      sst_ind = getindex(vars2d, 'sst') ! SST (2D)

      ! Initialize all constants required by routine
      allocate(workvar3d(nx_res,ny_res,nlevs))
      allocate(qworkvar3d(nx_res,ny_res,nlevs))
      allocate(qsatworkvar3d(nx_res,ny_res,nlevs))
      allocate(tvworkvar3d(nx_res,ny_res,nlevs))

      if (ntimes > 1) then
         write(6,*)'gridio/readgriddata: reading multiple backgrounds not yet supported'
         call stop2(23)
      endif
      ne = 0
      ensmemloop: do nanal=nanal1,nanal2
         ne = ne + 1

         backgroundloop: do nb=1,ntimes

            ! Define character string for ensemble member file
            if (nanal > 0) then
               write(charnanal,'(a3, i3.3)') 'mem', nanal
            else
               charnanal = 'ensmean'
            endif

            do ntile=1,ntiles
               nn_tile0=(ntile-1)*nx_res*ny_res
               write(char_tile, '(i1)') ntile

               filename = "fv3sar_tile"//char_tile//"_"//trim(charnanal)
               fv3filename=trim(adjustl(filename))//"_dynvar"
               fv3filename_tracer=trim(adjustl(filename))//"_tracer"
               fv3filename_phyvar=trim(adjustl(filename))//"_phyvar"

               !----------------------------------------------------------------------
               ! read u-component
               call nc_check( nf90_open(trim(adjustl(fv3filename)),nf90_nowrite,file_id),&
                  myname_,'open: '//trim(adjustl(fv3filename)) )

               call nc_check( nf90_open(trim(adjustl(fv3filename_tracer)),nf90_nowrite,file_id_t),&
                  myname_,'open: '//trim(adjustl(fv3filename_tracer)) )

               call nc_check( nf90_open(trim(adjustl(fv3filename_phyvar)),nf90_nowrite,file_id_p),&
                  myname_,'open: '//trim(adjustl(fv3filename_phyvar)) )



               !----------------------------------------------------------------------
               ! Update u and v variables (same for NMM and ARW)

               if (u_ind > 0) then
                  allocate(uworkvar3d(nx_res,ny_res+1,nlevs))
                  varstrname = 'u'

                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,uworkvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(u_ind-1)+k,nb,ne)=uworkvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(u_ind-1)+1, levels(u_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : u ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

                  deallocate(uworkvar3d)
               endif
               if (v_ind > 0) then
                  allocate(vworkvar3d(nx_res+1,ny_res,nlevs))
                  varstrname = 'v'
                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,vworkvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(v_ind-1)+k,nb,ne)=vworkvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(v_ind-1)+1, levels(v_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : v ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo
                  deallocate(vworkvar3d)

               endif

               if (w_ind > 0) then
                  varstrname = 'W'
                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(w_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(w_ind-1)+1, levels(w_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : w ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (q_ind > 0) then
                  varstrname = 'sphum'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,qworkvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(q_ind-1)+k,nb,ne)=qworkvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(q_ind-1)+1, levels(q_ind)
                     if (nproc .eq. 0) &
                        write(6,*) 'READFVregional : q ', &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif


               if (tv_ind > 0.or.tsen_ind>0) then
                  allocate(tsenworkvar3d(nx_res,ny_res,nlevs))
                  varstrname = 'T'
                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,tsenworkvar3d)
                  if(.not.  (q_ind > 0)) then
                     varstrname = 'sphum'
                     call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,qworkvar3d)
                  endif
                  do k=1,nlevs
                     do j=1,ny_res
                        do i=1,nx_res
                           workvar3d(i,j,k)=tsenworkvar3d(i,j,k)*(one+fv*qworkvar3d(i,j,k))
                        enddo
                     enddo
                  enddo
                  tvworkvar3d=workvar3d
                  if(tsen_ind > 0) then
                     workvar3d=tsenworkvar3d
                  endif
                  tmp_ind=max(tv_ind,tsen_ind) !then can't be both >0
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(tmp_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(tmp_ind-1)+1, levels(tmp_ind)
                     if (nproc .eq. 0)   then
                        write(6,*) 'READFVregional : t ',                           &
                        & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                     endif
                  enddo

                  if(allocated(tsenworkvar3d)) deallocate(tsenworkvar3d)
               endif



               if (oz_ind > 0) then
                  varstrname = 'o3mr'
                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(oz_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(oz_ind-1)+1, levels(oz_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : oz ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (ql_ind > 0) then
                  varstrname = 'liq_wat'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(ql_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(ql_ind-1)+1, levels(ql_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : ql ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (qr_ind > 0) then
                  varstrname = 'rainwat'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(qr_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(qr_ind-1)+1, levels(qr_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : qr ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (qi_ind > 0) then
                  varstrname = 'ice_wat'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(qi_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(qi_ind-1)+1, levels(qi_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : qi ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (qs_ind > 0) then
                  varstrname = 'snowwat'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(qs_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(qs_ind-1)+1, levels(qs_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : qs ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (qg_ind > 0) then
                  varstrname = 'graupel'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(qg_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(qg_ind-1)+1, levels(qg_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : qg ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (qh_ind > 0) then
                  varstrname = 'hailwat'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(qh_ind-1)+k,nb,ne)=workvar3d(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k = levels(qh_ind-1)+1, levels(qh_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : qh ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               if (dbz_ind > 0) then
                  allocate(workvar3d_tmp(nx_res-6,ny_res-6,nlevs))
                  varstrname = 'ref_f3d'
                  call read_fv3_restart_data3d(varstrname,fv3filename_phyvar,file_id_p,workvar3d_tmp)
                  workvar3d = 0.0_r_kind
                  workvar3d(4:nx_res-3,4:ny_res-3,1:nlevs)=workvar3d_tmp
                  deallocate(workvar3d_tmp)
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           vargrid(nn,levels(dbz_ind-1)+k,nb,ne)=workvar3d(i,j,nlevs+1-k)
                        enddo
                     enddo
                  enddo
                  do k = levels(dbz_ind-1)+1, levels(dbz_ind)
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'READFVregional : dbz ',                           &
                     & k, minval(vargrid(:,k,nb,ne)), maxval(vargrid(:,k,nb,ne))
                  enddo

               endif

               call nc_check( nf90_close(file_id),&
                  myname_,'close '//trim(fv3filename) )

               call nc_check( nf90_close(file_id_t),&
                  myname_,'close '//trim(fv3filename_tracer) )

               call nc_check( nf90_close(file_id_p),&
                  myname_,'close '//trim(fv3filename_phyvar) )

               ! set SST to zero for now
               if (sst_ind > 0) then
                  vargrid(:,levels(n3d)+sst_ind,nb,ne) = zero
               endif


               !----------------------------------------------------------------------
               ! Allocate memory for variables computed within routine

               if (ps_ind > 0) then
                  allocate(workprsi(nx_res,ny_res,nlevsp1))
                  allocate(pswork(nx_res,ny_res))
                  fv3filename=trim(adjustl(filename))//"_dynvar"
                  call nc_check( nf90_open(trim(adjustl(fv3filename)),nf90_nowrite,file_id),&
                     myname_,'open: '//trim(adjustl(fv3filename)) )
                  call read_fv3_restart_data3d('delp',fv3filename,file_id,workvar3d)
                  !print *,'min/max delp',ntile,minval(delp),maxval(delp)
                  call nc_check( nf90_close(file_id),&
                     myname_,'close '//trim(fv3filename) )
                  workprsi(:,:,nlevsp1)=eta1_ll(nlevsp1) !etal_ll is needed
                  do i=nlevs,1,-1
                     workprsi(:,:,i)=workvar3d(:,:,i)*0.01_r_kind+workprsi(:,:,i+1)
                  enddo

                  pswork(:,:)=workprsi(:,:,1)



                  nn = nn_tile0
                  do j=1,ny_res
                     do i=1,nx_res
                        nn=nn+1
                        vargrid(nn,levels(n3d)+ps_ind, nb,ne) =pswork(i,j)
                     enddo
                  enddo





                  do k=1,nlevs
                     do j=1,ny_res
                        do i=1,nx_res
                           workvar3d(i,j,k)=(workprsi(i,j,k)+workprsi(i,j,k+1))*half
                        enddo
                     enddo
                  enddo
                  ice=.true.  !tothink
                  if (pseudo_rh) then
                     call genqsat1(qworkvar3d,qsatworkvar3d,workvar3d,tvworkvar3d,ice,  &
                        nx_res*ny_res,nlevs)
                  else
                     qsatworkvar3d(:,:,:) = 1._r_double
                  endif
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           qsat(nn,k,nb,ne)=qsatworkvar3d(i,j,k)
                        enddo
                     enddo

                  enddo




                  if(allocated(workprsi))     deallocate(workprsi)
                  if(allocated(pswork))     deallocate(pswork)
                  if(allocated(tvworkvar3d)) deallocate(tvworkvar3d)
                  if(allocated(qworkvar3d)) deallocate(qworkvar3d)
                  if(allocated(qsatworkvar3d)) deallocate(qsatworkvar3d)
               endif
               !======================================================================
               ! Deallocate memory
               if(allocated(workvar3d))             deallocate(workvar3d)
            end do ! ntile loop

         end do backgroundloop ! loop over backgrounds to read in
      end do ensmemloop ! loop over ens members to read in

      return

   end subroutine readgriddata

   !========================================================================
   ! readgriddata_nmm.f90: read WRF-NMM state or control vector
   !-------------------------------------------------------------------------


   !========================================================================
   ! writegriddata.f90: write WRF-ARW or WRF-NMM analysis
   !-------------------------------------------------------------------------

   subroutine writegriddata(nanal1,nanal2,vars3d,vars2d,n3d,n2d,levels,ndim,vargrid,no_inflate_flag)
      use constants, only: zero, one,fv,half
      use gridinfo,only: eta1_ll,eta2_ll
      use params, only: nbackgrounds, anlfileprefixes, fgfileprefixes
      use params,   only: nx_res,ny_res,nlevs,ntiles,l_pres_add_saved
      use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
      use netcdf, only: nf90_inq_dimid,nf90_inq_varid
      use netcdf, only: nf90_write,nf90_write,nf90_inquire,nf90_inquire_dimension
      use write_fv3regional_restarts,only:write_fv3_restart_data1d,write_fv3_restart_data2d
      use write_fv3regional_restarts,only:write_fv3_restart_data3d,write_fv3_restart_data4d
      include 'netcdf.inc'

      !----------------------------------------------------------------------
      ! Define variables passed to subroutine
      integer, intent(in)  :: nanal1,nanal2, n2d, n3d, ndim
      character(len=*), dimension(n2d), intent(in) :: vars2d
      character(len=*), dimension(n3d), intent(in) :: vars3d
      integer, dimension(0:n3d), intent(in) :: levels
      real(r_single), dimension(npts,ndim,nbackgrounds,nanal2-nanal1+1), intent(in) :: vargrid
      logical, intent(in) :: no_inflate_flag

      !----------------------------------------------------------------------
      ! Define variables computed within subroutine
      character(len=500)  :: filename
      character(len=:),allocatable :: fv3filename,fv3filename_tracer,fv3filename_phyvar
      character(len=7)    :: charnanal

      !----------------------------------------------------------------------
      integer(i_kind) :: u_ind, v_ind, tv_ind, tsen_ind,q_ind, ps_ind,oz_ind
      integer(i_kind) :: w_ind, ql_ind,qr_ind,qi_ind,qs_ind,qg_ind,qh_ind,dbz_ind

      integer(i_kind) file_id,file_id_t,file_id_p
      real(r_single), dimension(:,:), allocatable ::pswork
      real(r_single), dimension(:,:,:), allocatable ::workvar3d,workinc3d,workinc3d2,uworkvar3d,&
         vworkvar3d,tvworkvar3d,tsenworkvar3d,&
         workprsi,qworkvar3d,qbgworkvar3d,workvar3d_tmp

      !----------------------------------------------------------------------
      ! Define variables required by for extracting netcdf variable
      ! fields
      integer :: nlevsp1
      ! Define variables required for netcdf variable I/O
      character(len=12) :: varstrname
      character(len=1) char_tile
      character(len=24),parameter :: myname_ = 'fv3: writegriddata'

      !----------------------------------------------------------------------
      ! Define counting variables
      integer :: i,j,k,nn,ntile,nn_tile0, nb,ne,nanal



      write(6,*)"anlfileprefixes, fgfileprefixes are not used in the current implementation"
      write(6,*)"the no_inflate_flag is not used in the currrent implementation ",no_inflate_flag
      !----------------------------------------------------------------------
      nlevsp1=nlevs+1

      u_ind   = getindex(vars3d, 'u')   !< indices in the state var arrays
      v_ind   = getindex(vars3d, 'v')   ! U and V (3D)
      tv_ind  = getindex(vars3d, 't')  ! Tv (3D)
      tsen_ind  = getindex(vars3d, 'tsen')  ! Tv (3D)
      q_ind   = getindex(vars3d, 'q')   ! Q (3D)

      w_ind   = getindex(vars3d, 'w')   ! W (3D)
      ql_ind   = getindex(vars3d, 'ql')   ! QL (3D)
      qr_ind   = getindex(vars3d, 'qr')   ! QR (3D)
      qi_ind   = getindex(vars3d, 'qi')   ! QI (3D)
      qs_ind   = getindex(vars3d, 'qs')   ! QS (3D)
      qg_ind   = getindex(vars3d, 'qg')   ! QG (3D)
      qh_ind   = getindex(vars3d, 'qh')   ! QH (3D)
      dbz_ind   = getindex(vars3d, 'dbz')   ! DBZ (3D)

      ps_ind  = getindex(vars2d, 'ps')  ! Ps (2D)


      !----------------------------------------------------------------------
      if (nbackgrounds > 1) then
         write(6,*)'gridio/writegriddata: writing multiple backgrounds not yet supported'
         call stop2(23)
      endif
      ne = 0
      ensmemloop: do nanal=nanal1,nanal2
         ne = ne + 1

         backgroundloop: do nb=1,nbackgrounds
            allocate(workinc3d(nx_res,ny_res,nlevs),workinc3d2(nx_res,ny_res,nlevsp1))
            allocate(workvar3d(nx_res,ny_res,nlevs))
            allocate(qworkvar3d(nx_res,ny_res,nlevs))
            allocate(qbgworkvar3d(nx_res,ny_res,nlevs))
            allocate(tvworkvar3d(nx_res,ny_res,nlevs))



            !----------------------------------------------------------------------
            ! First guess file should be copied to analysis file at scripting
            ! level; only variables updated by EnKF are changed
            write(charnanal,'(a3, i3.3)') 'mem', nanal

            !----------------------------------------------------------------------
            ! Update u and v variables (same for NMM and ARW)
            do ntile=1,ntiles
               nn_tile0=(ntile-1)*nx_res*ny_res
               write(char_tile, '(i1)') ntile
               filename = "fv3sar_tile"//char_tile//"_"//trim(charnanal)
               fv3filename=trim(adjustl(filename))//"_dynvar"

               fv3filename_tracer=trim(adjustl(filename))//"_tracer"
               fv3filename_phyvar=trim(adjustl(filename))//"_phyvar"

               !----------------------------------------------------------------------
               ! read u-component
               call nc_check( nf90_open(trim(adjustl(fv3filename)),NF_WRITE,file_id),&
                  myname_,'open: '//trim(adjustl(fv3filename)) )

               call nc_check( nf90_open(trim(adjustl(fv3filename_tracer)),nf_write,file_id_t),&
                  myname_,'open: '//trim(adjustl(fv3filename_tracer)) )

               call nc_check( nf90_open(trim(adjustl(fv3filename_phyvar)),nf_write,file_id_p),&
                  myname_,'open: '//trim(adjustl(fv3filename_phyvar)) )


               ! update CWM for WRF-NMM
               if (u_ind > 0) then
                  varstrname = 'u'
                  allocate(uworkvar3d(nx_res,ny_res+1,nlevs))

                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,uworkvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(u_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  uworkvar3d(:,1:ny_res,:)=uworkvar3d(:,1:ny_res,:)+workinc3d
                  uworkvar3d(:,ny_res+1,:)=uworkvar3d(:,ny_res,:)
                  call write_fv3_restart_data3d(varstrname,fv3filename,file_id,uworkvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,&
                           minval(uworkvar3d(:,:,k)), maxval(uworkvar3d(:,:,k))
                     END IF
                  ENDDO

                  deallocate(uworkvar3d)

               endif

               if (v_ind > 0) then
                  varstrname = 'v'
                  allocate(vworkvar3d(nx_res+1,ny_res,nlevs))

                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,vworkvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(v_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  vworkvar3d(1:nx_res,:,:)=vworkvar3d(1:nx_res,:,:)+workinc3d
                  vworkvar3d(nx_res+1,:,:)=vworkvar3d(nx_res,:,:)
                  call write_fv3_restart_data3d(varstrname,fv3filename,file_id,vworkvar3d)

                  DO k=1,nlevs
                    IF (nproc .eq. 0) THEN
                       WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,&
                          minval(vworkvar3d(:,:,k)), maxval(vworkvar3d(:,:,k))
                    END IF
                 ENDDO

                  deallocate(vworkvar3d)
               endif

               if (w_ind > 0) then
                  varstrname = 'W'

                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(w_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  workvar3d(nx_res+1,:,:)=workvar3d(nx_res,:,:)
                  call write_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if(q_ind>0) then

                  varstrname='sphum'
                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,qbgworkvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(q_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  qworkvar3d=qbgworkvar3d +workinc3d
                  where( qworkvar3d < 0.0_r_kind )
                     qworkvar3d = 0.0_r_kind
                  end where

                  workvar3d = qworkvar3d
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     if (nproc .eq. 0)                                               &
                        write(6,*) 'WRITEregional : sphum ',                           &
                     & k, minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                  enddo
               end if

               if (tv_ind > 0.or.tsen_ind>0 ) then

                  varstrname = 'T'
                  if(tsen_ind>0) then
                     do k=1,nlevs
                        nn = nn_tile0
                        do j=1,ny_res
                           do i=1,nx_res
                              nn=nn+1
                              workinc3d(i,j,k)=vargrid(nn,levels(tsen_ind-1)+k,nb,ne)
                           enddo
                        enddo
                     enddo
                     call read_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)
                     workvar3d=workvar3d+workinc3d
                     call write_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)
                  else  ! tv_ind >0
                     do k=1,nlevs
                        nn = nn_tile0
                        do j=1,ny_res
                           do i=1,nx_res
                              nn=nn+1
                              workinc3d(i,j,k)=vargrid(nn,levels(tv_ind-1)+k,nb,ne)
                           enddo
                        enddo
                     enddo

                     varstrname = 'T'
                     allocate(tsenworkvar3d(nx_res,ny_res,nlevs))
                     call read_fv3_restart_data3d(varstrname,fv3filename,file_id,tsenworkvar3d)
                     if( .not. (q_ind > 0) )then
                        varstrname = 'sphum'
                        call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,qworkvar3d)
                     end if
                     tvworkvar3d=tsenworkvar3d*(one+fv*qworkvar3d)
                     tvworkvar3d=tvworkvar3d+workinc3d
                     if(.not. ( q_ind > 0)) then
                        tsenworkvar3d=tvworkvar3d/(one+fv*qbgworkvar3d)
                     else
                        tsenworkvar3d=tvworkvar3d/(one+fv*qworkvar3d)
                     endif
                     tsenworkvar3d=tvworkvar3d/(one+fv*qworkvar3d)
                     varstrname = 'T'
                     call write_fv3_restart_data3d(varstrname,fv3filename,file_id,tsenworkvar3d)
                     do k=1,nlevs
                        if (nproc .eq. 0)                                               &
                           write(6,*) 'WRITEregional : T ',                           &
                        & k, minval(tsenworkvar3d(:,:,k)), maxval(tsenworkvar3d(:,:,k))
                     enddo

                     deallocate(tsenworkvar3d)
                  endif  !if tsens else tv

               endif


               if (oz_ind > 0) then
                  varstrname = 'o3mr'

                  call read_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(oz_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d=workvar3d+workinc3d
                  call write_fv3_restart_data3d(varstrname,fv3filename,file_id,workvar3d)

               endif

               if (ql_ind > 0) then
                  varstrname = 'liq_wat'

                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(ql_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if (qr_ind > 0) then
                  varstrname = 'rainwat'

                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(qr_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if (qi_ind > 0) then
                  varstrname = 'ice_wat'

                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(qi_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if (qs_ind > 0) then
                  varstrname = 'snowwat'

                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(qs_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if (qg_ind > 0) then
                  varstrname = 'graupel'

                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(qg_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if (qh_ind > 0) then
                  varstrname = 'hailwat'

                  call read_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,k)=vargrid(nn,levels(qh_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_tracer,file_id_t,workvar3d)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d(:,:,k)), maxval(workvar3d(:,:,k))
                     END IF
                  END DO

               endif

               if (dbz_ind > 0) then
                  varstrname = 'ref_f3d'
                  allocate(workvar3d_tmp(nx_res-6,ny_res-6,nlevs))

                  call read_fv3_restart_data3d(varstrname,fv3filename_phyvar,file_id_p,workvar3d_tmp)
                  workvar3d = 0.0_r_kind
                  workvar3d(4:nx_res-3,4:ny_res-3,1:nlevs)=workvar3d_tmp
                  where( workvar3d < 0.0_r_kind )
                     workvar3d = 0.0_r_kind
                  end where
                  do k=1,nlevs
                     nn = nn_tile0
                     do j=1,ny_res
                        do i=1,nx_res
                           nn=nn+1
                           workinc3d(i,j,nlevs+1-k)=vargrid(nn,levels(dbz_ind-1)+k,nb,ne)
                        enddo
                     enddo
                  enddo
                  workvar3d(1:nx_res,:,:)=workvar3d(1:nx_res,:,:)+workinc3d
                  workvar3d_tmp = workvar3d(4:nx_res-3,4:ny_res-3,1:nlevs)
                  where( workvar3d_tmp < 0.0_r_kind )
                     workvar3d_tmp = 0.0_r_kind
                  end where
                  call write_fv3_restart_data3d(varstrname,fv3filename_phyvar,file_id_p,workvar3d_tmp)

                  DO k=1,nlevs
                     IF (nproc .eq. 0) THEN
                        WRITE(6,*) 'WRITEregional : ',TRIM(varstrname),' k = ',k,    &
                           minval(workvar3d_tmp(:,:,k)), maxval(workvar3d_tmp(:,:,k))
                     END IF
                  END DO

                  deallocate(workvar3d_tmp)
               endif

               if (ps_ind > 0) then
                  allocate(workprsi(nx_res,ny_res,nlevsp1))
                  allocate(pswork(nx_res,ny_res))
                  varstrname = 'delp'
                  call read_fv3_restart_data3d(varstrname,filename,file_id,workvar3d)
                  !print *,'min/max delp',ntile,minval(delp),maxval(delp)
                  workprsi(:,:,nlevsp1)=eta1_ll(nlevsp1) !etal_ll is needed
                  do i=nlevs,1,-1
                     workprsi(:,:,i)=workvar3d(:,:,i)*0.01_r_kind+workprsi(:,:,i+1)
                  enddo



                  nn = nn_tile0
                  do j=1,ny_res
                     do i=1,nx_res
                        nn=nn+1
                        pswork(i,j)=vargrid(nn,levels(n3d)+ps_ind,nb,ne)
                     enddo
                  enddo
                  if(l_pres_add_saved) then
                     do k=1,nlevs+1
                        do j=1,ny_res
                           do i=1,nx_res
                              workinc3d2(i,j,k)=eta2_ll(k)*pswork(i,j)
                           enddo
                        enddo
                     enddo
                     workprsi=workprsi+workinc3d2
                  else
                     workprsi(:,:,1)=workprsi(:,:,1)+pswork
                     do k=2,nlevsp1
                        workprsi(:,:,k)=eta1_ll(k)+eta2_ll(k)*workprsi(:,:,1)
                     enddo
                  endif
                  do k=1,nlevs
                     workvar3d(:,:,k)=(workprsi(:,:,k)-workprsi(:,:,k+1))*100.0
                  enddo


                  call write_fv3_restart_data3d(varstrname,filename,file_id,workvar3d)
               endif

               call nc_check( nf90_close(file_id),&
                  myname_,'close '//trim(fv3filename) )

               call nc_check( nf90_close(file_id_t),&
                  myname_,'close '//trim(fv3filename_tracer) )

               call nc_check( nf90_close(file_id_p),&
                  myname_,'close '//trim(fv3filename_phyvar) )

               !----------------------------------------------------------------------
               ! update time stamp is to be considered NSTART_HOUR in NMM (HWRF) restart file.
               !======================================================================
            end do ! tiles
            if(allocated(workinc3d))     deallocate(workinc3d)
            if(allocated(workinc3d2))     deallocate(workinc3d2)
            if(allocated(workprsi))     deallocate(workprsi)
            if(allocated(pswork))     deallocate(pswork)
            if(allocated(tvworkvar3d)) deallocate(tvworkvar3d)
            if(allocated(qworkvar3d)) deallocate(qworkvar3d)

         end do backgroundloop ! loop over backgrounds to read in
      end do ensmemloop ! loop over ens members to read in

      if(allocated(workvar3d)) deallocate(workvar3d)

      ! Return calculated values
      return

      !======================================================================

   end subroutine writegriddata

   subroutine readgriddata_pnc(vars3d,vars2d,n3d,n2d,levels,ndim,ntimes, &
      fileprefixes,filesfcprefixes,reducedgrid,grdin,qsat)
      use constants, only: max_varname_length
      implicit none
      character(len=max_varname_length), dimension(n2d), intent(in) :: vars2d
      character(len=max_varname_length), dimension(n3d), intent(in) :: vars3d
      integer, intent(in) :: n2d, n3d
      integer, dimension(0:n3d), intent(in) :: levels
      integer, intent(in) :: ndim, ntimes
      character(len=120), dimension(7), intent(in)  :: fileprefixes
      character(len=120), dimension(7), intent(in)  :: filesfcprefixes
      logical, intent(in) :: reducedgrid
      real(r_single), dimension(npts,ndim,ntimes,1), intent(out) :: grdin
      real(r_double), dimension(npts,nlevs,ntimes,1), intent(out) :: qsat
   end subroutine readgriddata_pnc

   subroutine writegriddata_pnc(vars3d,vars2d,n3d,n2d,levels,ndim,grdin,no_inflate_flag)
      use constants, only: max_varname_length
      use params, only: nbackgrounds
      implicit none
      character(len=max_varname_length), dimension(n2d), intent(in) :: vars2d
      character(len=max_varname_length), dimension(n3d), intent(in) :: vars3d
      integer, intent(in) :: n2d,n3d,ndim
      integer, dimension(0:n3d), intent(in) :: levels
      real(r_single), dimension(npts,ndim,nbackgrounds,1), intent(inout) :: grdin
      logical, intent(in) :: no_inflate_flag
   end subroutine writegriddata_pnc

   subroutine writeincrement(nanal1,nanal2,vars3d,vars2d,n3d,n2d,levels,ndim,grdin,no_inflate_flag)
      use constants, only: max_varname_length
      use params, only: nbackgrounds
      implicit none
      integer, intent(in) :: nanal1,nanal2
      character(len=max_varname_length), dimension(n2d), intent(in) :: vars2d
      character(len=max_varname_length), dimension(n3d), intent(in) :: vars3d
      integer, intent(in) :: n2d,n3d,ndim
      integer, dimension(0:n3d), intent(in) :: levels
      real(r_single), dimension(npts,ndim,nbackgrounds,1), intent(inout) :: grdin
      logical, intent(in) :: no_inflate_flag
   end subroutine writeincrement

   subroutine writeincrement_pnc(vars3d,vars2d,n3d,n2d,levels,ndim,grdin,no_inflate_flag)
      use constants, only: max_varname_length
      use params, only: nbackgrounds
      implicit none
      character(len=max_varname_length), dimension(n2d), intent(in) :: vars2d
      character(len=max_varname_length), dimension(n3d), intent(in) :: vars3d
      integer, intent(in) :: n2d,n3d,ndim
      integer, dimension(0:n3d), intent(in) :: levels
      real(r_single), dimension(npts,ndim,nbackgrounds,1), intent(inout) :: grdin
      logical, intent(in) :: no_inflate_flag
   end subroutine writeincrement_pnc

end module gridio
