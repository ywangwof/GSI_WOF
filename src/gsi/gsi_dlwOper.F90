module gsi_dlwOper
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module gsi_dlwOper
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2018-08-10
!
! abstract: an obOper extension for dlwNode type
!
! program history log:
!   2018-08-10  j guo   - added this document block
!   2019-06-25  Hu      - add diag_radardlw for controling radar reflectivity
!                         diag file (in module obsmod).
!   2019-07-22  j Guo   - moved diag_radardlw and its description here from
!                         obsmod.
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! module interface:

  use gsi_obOper, only: obOper
  use m_dlwNode , only: dlwNode
  implicit none
  public:: dlwOper      ! data stracture
  public:: diag_radardlw

  type,extends(obOper):: dlwOper
  contains
    procedure,nopass:: mytype
    procedure,nopass:: nodeMold
    procedure:: setup_
    procedure:: intjo1_
    procedure:: stpjo1_
  end type dlwOper

!   def diag_radardlw- namelist logical to compute/write (=true) radar
!                                          reflectiivty diag files
  logical,save:: diag_radardlw=.false.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='gsi_dlwOper'
  type(dlwNode),save,target:: myNodeMold_

contains
  function mytype(nodetype)
    implicit none
    character(len=:),allocatable:: mytype
    logical,optional, intent(in):: nodetype
    mytype="[dlwOper]"
    if(present(nodetype)) then
      if(nodetype) mytype=myNodeMold_%mytype()
    endif
  end function mytype

  function nodeMold()
  !> %nodeMold() returns a mold of its corresponding obsNode
    use m_obsNode, only: obsNode
    implicit none
    class(obsNode),pointer:: nodeMold
    nodeMold => myNodeMold_
  end function nodeMold

  subroutine setup_(self, lunin, mype, is, nobs, init_pass,last_pass)
    use dlw_setup, only: setup
    use kinds, only: i_kind
    use gsi_obOper, only: len_obstype
    use gsi_obOper, only: len_isis

    use m_rhs  , only: awork => rhs_awork
    use m_rhs  , only: bwork => rhs_bwork
    use m_rhs  , only: iwork => i_dlw

    use obsmod  , only: write_diag
    use jfunc   , only: jiter
    use convinfo, only: diag_conv

    use mpeu_util, only: die
    implicit none
    class(dlwOper ), intent(inout):: self
    integer(i_kind), intent(in):: lunin
    integer(i_kind), intent(in):: mype
    integer(i_kind), intent(in):: is
    integer(i_kind), intent(in):: nobs
    logical        , intent(in):: init_pass     ! supporting multi-pass setup()
    logical        , intent(in):: last_pass     ! with incremental backgrounds.

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::setup_"

    character(len=len_obstype):: obstype
    character(len=len_isis   ):: isis
    integer(i_kind):: nreal,nchanl,ier,nele
    logical:: diagsave

    if(nobs == 0) return

    read(lunin,iostat=ier) obstype,isis,nreal,nchanl
    if(ier/=0) call die(myname_,'read(obstype,...), iostat =',ier)
    nele = nreal+nchanl

    diagsave  = write_diag(jiter) .and. diag_conv

    call setup(self%obsLL(:), self%odiagLL(:), &
        lunin,mype,bwork,awork(:,iwork),nele,nobs,is,diagsave,init_pass)

  end subroutine setup_

  subroutine intjo1_(self, ibin, rval,sval, qpred,sbias)
    use intdlwmod, only: intjo => intdlw
    use gsi_bundlemod  , only: gsi_bundle
    use bias_predictors, only: predictors
    use m_obsNode , only: obsNode
    use m_obsLList, only: obsLList_headNode
    use kinds     , only: i_kind, r_quad
    implicit none
    class(dlwOper  ),intent(in   ):: self
    integer(i_kind ),intent(in   ):: ibin
    type(gsi_bundle),intent(inout):: rval   ! (ibin)
    type(gsi_bundle),intent(in   ):: sval   ! (ibin)
    real(r_quad    ),target,dimension(:),intent(inout):: qpred  ! (ibin)
    type(predictors),target,             intent(in   ):: sbias

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::intjo1_"
    class(obsNode),pointer:: headNode

    headNode => obsLList_headNode(self%obsLL(ibin))
    call intjo(headNode, rval,sval)
    headNode => null()

  end subroutine intjo1_

  subroutine stpjo1_(self, ibin, dval,xval,pbcjo,sges,nstep,dbias,xbias)
    use stpdlwmod, only: stpjo => stpdlw
    use gsi_bundlemod, only: gsi_bundle
    use bias_predictors, only: predictors
    use m_obsNode , only: obsNode
    use m_obsLList, only: obsLList_headNode
    use kinds, only: r_quad,r_kind,i_kind
    implicit none
    class(dlwOper  ),intent(in):: self
    integer(i_kind ),intent(in):: ibin
    type(gsi_bundle),intent(in):: dval
    type(gsi_bundle),intent(in):: xval
    real(r_quad    ),dimension(:),intent(inout):: pbcjo ! (1:4)
    real(r_kind    ),dimension(:),intent(in   ):: sges
    integer(i_kind),intent(in):: nstep

    type(predictors),target, intent(in):: dbias
    type(predictors),target, intent(in):: xbias

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::stpjo1_"
    class(obsNode),pointer:: headNode

    headNode => obsLList_headNode(self%obsLL(ibin))
    call stpjo(headNode,dval,xval,pbcjo(:),sges,nstep)
    headNode => null()
  end subroutine stpjo1_

end module gsi_dlwOper
