module UrbanDynAlbMod
!----------------------------------------------------------------------- 
!DESCRIPTION:
!Dynamic Urban Albedo Input Stream Data
!The time-varing urban albedo is read in using this module (a stream)
!USES:
 use shr_strdata_mod , only : shr_strdata_type
 use shr_kind_mod    , only : r8 => shr_kind_r8, CL => shr_kind_CL
 use shr_log_mod     , only : errMsg => shr_log_errMsg
 use abortutils      , only : endrun
 use decompMod       , only : bounds_type
 use clm_varctl      , only : iulog
 use clm_varcon      , only : spval
 use LandunitType    , only : lun
 use GridcellType    , only : grc
 use mct_mod
 use landunit_varcon , only : isturb_MIN, isturb_MAX         ! isturb_MIN = 7, isturb_MAX = 9
 use clm_varctl      , only : Dynamic_UrbanAlbedoRoof, Dynamic_UrbanAlbedoImproad, Dynamic_UrbanAlbedoWall
 use clm_varpar      , only : numrad
 ! 
 implicit none
 save
 private
 
 public dynAlbinit  ! initialize transient urban albedo
 
 ! !PUBLIC TYPE
 type, public :: urbanalbtv_type
    ! urban roof albedo inputs
    real(r8), public, pointer :: dyn_alb_roof_dir        (:,:) ! dynamic lun direct  roof albedo
    real(r8), public, pointer :: dyn_alb_roof_dif        (:,:) ! dynamic lun diffuse roof albedo
    real(r8), public, pointer :: dyn_alb_improad_dir     (:,:) ! dynamic lun direct roof albedo
    real(r8), public, pointer :: dyn_alb_improad_dif     (:,:) ! dynamic lun diffuse roof albedo
    real(r8), public, pointer :: dyn_alb_wall_dir        (:,:) ! dynamic lun direct wall albedo
    real(r8), public, pointer :: dyn_alb_wall_dif        (:,:) ! dynamic lun diffuse wall albedo
    ! 
    type(shr_strdata_type)    :: sdat_urbanalbtvroof         ! urban time varying roof albedo data stream
    type(shr_strdata_type)    :: sdat_urbanalbtvimproad      ! urban time varying improad albedo data stream
    type(shr_strdata_type)    :: sdat_urbanalbtvwall         ! urban time varying wall albedo data stream
    
   contains
     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: dynAlbinit                         ! Allocate and initialize urbanalbtv
     procedure, public :: urbanalbtvroof_init                ! Initialize urban roof albedo time varying stream
     procedure, public :: urbanalbtvroof_interp              ! Interpolate urban roof alebdo time varying stream
     procedure, public :: urbanalbtvimproad_init             ! Initialize urban improad albedo time varying stream
     procedure, public :: urbanalbtvimproad_interp           ! Interpolate urban improad alebdo time varying stream
     procedure, public :: urbanalbtvwall_init                ! Initialize urban wall albedo time varying stream
     procedure, public :: urbanalbtvwall_interp              ! Interpolate urban wall alebdo time varying stream
 end type urbanalbtv_type
  
  character(30), private :: stream_var_name_roof(isturb_MIN:isturb_MAX)
  character(30), private :: stream_var_name_improad(isturb_MIN:isturb_MAX)
  character(30), private :: stream_var_name_wall(isturb_MIN:isturb_MAX)
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !----------------------------------------------------------------------- 
 contains
  !-----------------------------------------------------------------------
  subroutine dynAlbinit(this, bounds)
  ! !DESCRIPTION:
  ! Initialize data stream information for dynamic urban albedo
  ! !USES:
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use histFileMod     , only : hist_addfld2d
  ! !ARGUMENTS:
  class(urbanalbtv_type)                 :: this
  type(bounds_type)      , intent(in)    :: bounds
  ! !LOCAL VARIABLES:  
  integer             :: begl, endl
  !---------------------------------------------------------------------
  begl = bounds%begl; endl = bounds%endl                        ! beginning and ending landunit index
  ! 
  ! Allocate urbanalbtv data structures
  ! 
  allocate(this%dyn_alb_roof_dir        (begl:endl,numrad))   ; this%dyn_alb_roof_dir        (:,:) = nan
  allocate(this%dyn_alb_roof_dif        (begl:endl,numrad))   ; this%dyn_alb_roof_dif        (:,:) = nan 
  allocate(this%dyn_alb_improad_dir     (begl:endl,numrad))   ; this%dyn_alb_improad_dir     (:,:) = nan   
  allocate(this%dyn_alb_improad_dif     (begl:endl,numrad))   ; this%dyn_alb_improad_dif     (:,:) = nan
  allocate(this%dyn_alb_wall_dir        (begl:endl,numrad))   ; this%dyn_alb_wall_dir        (:,:) = nan   
  allocate(this%dyn_alb_wall_dif        (begl:endl,numrad))   ; this%dyn_alb_wall_dif        (:,:) = nan
  
  if (Dynamic_UrbanAlbedoRoof) then 
     call this%urbanalbtvroof_init(bounds)
     call this%urbanalbtvroof_interp(bounds)
     call hist_addfld2d (fname='DYNALB_ROOF_DIR', units='',      &
            avgflag='A', long_name='time varing urban roof albedo dir',  type2d='numrad', &
            ptr_lunit=this%dyn_alb_roof_dir, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
     call hist_addfld2d (fname='DYNALB_ROOF_DIF', units='',      &
            avgflag='A', long_name='time varing urban roof albedo dif',  type2d='numrad',  &
            ptr_lunit=this%dyn_alb_roof_dif, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
  end if
  
  if (Dynamic_UrbanAlbedoImproad) then    
     call this%urbanalbtvimproad_init(bounds)
     call this%urbanalbtvimproad_interp(bounds)
     call hist_addfld2d (fname='DYNALB_IMPROAD_DIR', units='',      &
            avgflag='A', long_name='time varing urban improad albedo dir',  type2d='numrad', &
            ptr_lunit=this%dyn_alb_improad_dir, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
     call hist_addfld2d (fname='DYNALB_IMPROAD_DIF', units='',      &
            avgflag='A', long_name='time varing urban improad albedo dif',  type2d='numrad',  &
            ptr_lunit=this%dyn_alb_improad_dif, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
  end if
  
  if (Dynamic_UrbanAlbedoWall) then
     call this%urbanalbtvwall_init(bounds)
     call this%urbanalbtvwall_interp(bounds)
     call hist_addfld2d (fname='DYNALB_WALL_DIR', units='',      &
            avgflag='A', long_name='time varing urban wall albedo dir',  type2d='numrad', &
            ptr_lunit=this%dyn_alb_wall_dir, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
     call hist_addfld2d (fname='DYNALB_WALL_DIF', units='',      &
            avgflag='A', long_name='time varing urban wall albedo dif',  type2d='numrad',  &
            ptr_lunit=this%dyn_alb_wall_dif, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
  end if   
  
 end subroutine dynAlbinit

 !---------------------------------------------------------------------
  
 subroutine urbanalbtvroof_init(this, bounds)
   ! !DESCRIPTION:
   ! !USES:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use spmdMod          , only : masterproc, mpicom, comp_id
   use fileutils        , only : getavu, relavu
   use shr_mpi_mod      , only : shr_mpi_bcast
   use shr_string_mod   , only : shr_string_listAppend
   use shr_strdata_mod  , only : shr_strdata_create, shr_strdata_print
   use decompMod        , only : gsmap_lnd_gdc2glo
   use domainMod        , only : ldomain
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use landunit_varcon  , only : isturb_TBD, isturb_HD, isturb_MD
   use clm_varctl       , only : NLFilename_in
   
   ! !ARGUMENTS:
   implicit none
   class(urbanalbtv_type)         :: this
   type(bounds_type), intent(in)  :: bounds
   ! 
   ! !LOCAL VARIABLES:
   integer            :: begl, endl                                  ! landunits
   integer            :: ifield                                      ! field index
   integer            :: stream_year_first_urbanalbtvroof            ! first year in urban roof albedo tv stream to use
   integer            :: stream_year_last_urbanalbtvroof             ! last year in urban roof albedo tv stream to use
   integer            :: model_year_align_urbanalbtvroof             ! align stream_year_first_urbantvroof with this model year
   integer            :: nu_nml                                      ! unit for namelist file 
   integer            :: nml_error                                   ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                                     ! domain information
   character(len=CL)  :: stream_fldFileName_urbanalbtvroof           ! urban roof albedo time-varying streams filename
   character(len=CL)  :: urbanalbtvroofmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbanalbtvroof_tintalgo = 'linear'          ! time interpolation alogrithm 
   character(*), parameter :: urbanalbtvroofString = "dyn_alb_roof_" ! base string for field string
   character(SHR_KIND_CL)  :: fldList                                ! field string
   character(*), parameter :: subName = "('urbanalbtvroof_init')"
   character(*), parameter :: F00 = "('(urbanalbtvroof_init) ',4a)"
   !-----------------------------------------------------------------------
   namelist /urbanalbtvroof_streams/       &
        stream_year_first_urbanalbtvroof,  &  
        stream_year_last_urbanalbtvroof,   &  
        model_year_align_urbanalbtvroof,   &  
        urbanalbtvroofmapalgo,             &  
        stream_fldFileName_urbanalbtvroof, &      
        urbanalbtvroof_tintalgo  
   !-----------------------------------------------------------------------                      
   begl = bounds%begl; endl = bounds%endl
   ! Default values for namelist
   stream_year_first_urbanalbtvroof  = 1      ! first year in stream to use
   stream_year_last_urbanalbtvroof   = 1      ! last  year in stream to use
   model_year_align_urbanalbtvroof   = 1      ! align stream_year_first_urbanalbtvroof with this model year
   stream_fldFileName_urbanalbtvroof = ' '
   
   ! Read urbanalbtv_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( newunit=nu_nml, file=trim(NLFilename_in), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbanalbtvroof_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbanalbtvroof_streams,iostat=nml_error) 
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbanalbtvroof_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
          call endrun(subname // ':: ERROR finding urbanalbtvroof_streams namelist')   
      end if
      close(nu_nml)
      call relavu( nu_nml )
    endif

    call shr_mpi_bcast(stream_year_first_urbanalbtvroof  , mpicom)
    call shr_mpi_bcast(stream_year_last_urbanalbtvroof   , mpicom)
    call shr_mpi_bcast(model_year_align_urbanalbtvroof   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_urbanalbtvroof , mpicom)
    call shr_mpi_bcast(urbanalbtvroof_tintalgo           , mpicom)
   

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'Attemping to read time varying urban roof albedo parameters......'
       write(iulog,'(a)') 'urbanalbtvroof_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbanalbtvroof  = ',stream_year_first_urbanalbtvroof
       write(iulog,'(a,i8)') '  stream_year_last_urbanalbtvroof   = ',stream_year_last_urbanalbtvroof
       write(iulog,'(a,i8)') '  model_year_align_urbanalbtvroof   = ',model_year_align_urbanalbtvroof
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbanalbtvroof = ',stream_fldFileName_urbanalbtvroof
       write(iulog,'(a,a)' ) '  urbanalbtvroof_tintalgo           = ',urbanalbtvroof_tintalgo
       write(iulog,*) 'Read in urbanalbtvroof namelist from:',trim(NLFilename_in)
    endif
    
    ! Initialize the cdeps data type this%sdat_urbanalbtv
    call clm_domain_mct (bounds, dom_clm)
    
    ! create the field list for urban albedo fields
    stream_var_name_roof(:)          = "NOT_SET"
    stream_var_name_roof(isturb_TBD) = urbanalbtvroofString//"TBD"
    stream_var_name_roof(isturb_HD)  = urbanalbtvroofString//"HD"
    stream_var_name_roof(isturb_MD)  = urbanalbtvroofString//"MD"   
    fldList = ""
    do ifield = isturb_MIN, isturb_MAX
       call shr_string_listAppend( fldList, stream_var_name_roof(ifield) )
    end do
    
   call shr_strdata_create(this%sdat_urbanalbtvroof,name="clmurbanalbtvroof",     &
        pio_subsystem=pio_subsystem,                   &
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,  &
        yearFirst=stream_year_first_urbanalbtvroof,        &
        yearLast=stream_year_last_urbanalbtvroof,          &
        yearAlign=model_year_align_urbanalbtvroof,         &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_urbanalbtvroof),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &
        domAreaName='area',                            &
        domMaskName='LANDMASK',                        &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_urbanalbtvroof)/) , &
        fldListFile=fldList,                           &
        fldListModel=fldList,                          &
        fillalgo='none',                               &
        mapalgo=urbanalbtvroofmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo=urbanalbtvroof_tintalgo,                     &
        taxmode='extend'                               )

   if (masterproc) then
      call shr_strdata_print(this%sdat_urbanalbtvroof,'urban time varying roof albedo')
   endif
   
        
  end subroutine urbanalbtvroof_init

  !==============================================================================
 
  subroutine urbanalbtvroof_interp(this, bounds)
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying albedo.
    ! 
    ! !USES:
    use clm_time_manager , only : get_curr_date
    use spmdMod          , only : mpicom
    use clm_instur       , only : urban_valid
    use shr_strdata_mod  , only : shr_strdata_advance
    ! 
    ! !ARGUMENTS:
    ! 
    class(urbanalbtv_type)           :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !
    logical :: found
    integer :: l       ! landunit
    integer :: glun
    integer :: ib
    integer :: ig
    integer :: g
    integer :: ip      
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    ! 
    !-----------------------------------------------------------------------
    ! 
    ! Advance sdat stream
    !
    call get_curr_date(year, mon, day, sec)
    !
    ! packing the date into an integer
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_urbanalbtvroof, mcdate, sec, mpicom,'urbanalbtvroof')

    ! Determine this%dyn_alb_roof for all landunits
    do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
             glun = lun%gridcell(l)
             ip = mct_aVect_indexRA(this%sdat_urbanalbtvroof%avs(1),trim(stream_var_name_roof(lun%itype(l))))
             ig = 0
             do g = bounds%begg, bounds%endg
                ig = ig+1
                if (g==glun) exit
             end do   
             do ib = 1,numrad         
                this%dyn_alb_roof_dir(l,ib) = this%sdat_urbanalbtvroof%avs(1)%rAttr(ip,ig)
                this%dyn_alb_roof_dif(l,ib) = this%sdat_urbanalbtvroof%avs(1)%rAttr(ip,ig)     
             end do                 
       else
             do ib = 1,numrad
                this%dyn_alb_roof_dir(l,ib) = spval
                this%dyn_alb_roof_dif(l,ib) = spval  
             end do      
       end if
    end do
    
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
             glun  = lun%gridcell(l)
          ! 
          ! Determine vector index corresponding to glun
          ! 
             ig = 0
             do g = bounds%begg,bounds%endg
                ig = ig+1
                if (g == glun) exit
             end do
             
             do ib = 1,numrad
                if ( .not. urban_valid(g) .or. (this%dyn_alb_roof_dir(l,ib) <= 0._r8) .or. (this%dyn_alb_roof_dif(l,ib) <= 0._r8)) then
                   found = .true.
                   gindx = g
                   lindx = l
                   exit
                end if
             end do     
       end if
    end do
    
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban roof data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(l)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'dyn_alb_roof_dir:  ',this%dyn_alb_roof_dir(lindx,:)
       write(iulog,*)'dyn_alb_roof_dif:  ',this%dyn_alb_roof_dif(lindx,:)
       call endrun(msg=errmsg(sourcefile, __LINE__))
    end if
  end subroutine urbanalbtvroof_interp
  !-----------------------------------------------------------------------
 
  subroutine urbanalbtvimproad_init(this, bounds)
   ! !DESCRIPTION:
   ! !USES:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use spmdMod          , only : masterproc, mpicom, comp_id
   use fileutils        , only : getavu, relavu
   use shr_mpi_mod      , only : shr_mpi_bcast
   use shr_string_mod   , only : shr_string_listAppend
   use shr_strdata_mod  , only : shr_strdata_create, shr_strdata_print
   use decompMod        , only : gsmap_lnd_gdc2glo
   use domainMod        , only : ldomain
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use landunit_varcon  , only : isturb_TBD, isturb_HD, isturb_MD
   use clm_varctl       , only : NLFilename_in
   ! !ARGUMENTS:
   implicit none
   class(urbanalbtv_type)         :: this
   type(bounds_type), intent(in)  :: bounds
   ! !LOCAL VARIABLES:
   integer            :: begl, endl                           ! landunits
   integer            :: ifield                               ! field index
   integer            :: stream_year_first_urbanalbtvimproad            ! first year in urban improad albedo tv stream to use
   integer            :: stream_year_last_urbanalbtvimproad             ! last year in urban improad albedo tv stream to use
   integer            :: model_year_align_urbanalbtvimproad             ! align stream_year_first_urbantvimproad with this model year
   integer            :: nu_nml                               ! unit for namelist file 
   integer            :: nml_error                            ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                              ! domain information
   character(len=CL)  :: stream_fldFileName_urbanalbtvimproad        ! urban improad albedo time-varying streams filename
   character(len=CL)  :: urbanalbtvimproadmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbanalbtvimproad_tintalgo = 'linear'          ! time interpolation alogrithm 
   character(*), parameter :: urbanalbtvimproadString = "dyn_alb_improad_" ! base string for field string
   character(SHR_KIND_CL)  :: fldList                         ! field string
   character(*), parameter :: subName = "('urbanalbtvimproad_init')"
   character(*), parameter :: F00 = "('(urbanalbtvimproad_init) ',4a)"
   !-----------------------------------------------------------------------
   namelist /urbanalbtvimproad_streams/       &
        stream_year_first_urbanalbtvimproad,  &  
        stream_year_last_urbanalbtvimproad,   &  
        model_year_align_urbanalbtvimproad,   &  
        urbanalbtvimproadmapalgo,             &  
        stream_fldFileName_urbanalbtvimproad, &      
        urbanalbtvimproad_tintalgo  
   !-----------------------------------------------------------------------    
   !                   
   begl = bounds%begl; endl = bounds%endl
   ! 
   ! Default values for namelist
   stream_year_first_urbanalbtvimproad  = 1      ! first year in stream to use
   stream_year_last_urbanalbtvimproad   = 1      ! last  year in stream to use
   model_year_align_urbanalbtvimproad   = 1      ! align stream_year_first_urbanalbtvimproad with this model year
   stream_fldFileName_urbanalbtvimproad = ' '
   
   ! Read urbanalbtv_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( newunit=nu_nml, file=trim(NLFilename_in), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbanalbtvimproad_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbanalbtvimproad_streams,iostat=nml_error) 
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbanalbtvimproad_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
          call endrun(subname // ':: ERROR finding urbanalbtvimproad_streams namelist')   
      end if
      close(nu_nml)
      call relavu( nu_nml )
    endif

    call shr_mpi_bcast(stream_year_first_urbanalbtvimproad  , mpicom)
    call shr_mpi_bcast(stream_year_last_urbanalbtvimproad   , mpicom)
    call shr_mpi_bcast(model_year_align_urbanalbtvimproad   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_urbanalbtvimproad , mpicom)
    call shr_mpi_bcast(urbanalbtvimproad_tintalgo           , mpicom)
   

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'Attemping to read time varying urban improad albedo parameters......'
       write(iulog,'(a)') 'urbanalbtvimproad_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbanalbtvimproad  = ',stream_year_first_urbanalbtvimproad
       write(iulog,'(a,i8)') '  stream_year_last_urbanalbtvimproad   = ',stream_year_last_urbanalbtvimproad
       write(iulog,'(a,i8)') '  model_year_align_urbanalbtvimproad   = ',model_year_align_urbanalbtvimproad
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbanalbtvimproad = ',stream_fldFileName_urbanalbtvimproad
       write(iulog,'(a,a)' ) '  urbanalbtvimproad_tintalgo           = ',urbanalbtvimproad_tintalgo
       write(iulog,*) 'Read in urbanalbtvimproad namelist from:',trim(NLFilename_in)
    endif
    
    ! Initialize the cdeps data type this%sdat_urbanalbtv
    call clm_domain_mct (bounds, dom_clm)
    
    ! create the field list for urban improad albedo fields
    stream_var_name_improad(:)          = "NOT_SET"
    stream_var_name_improad(isturb_TBD) = urbanalbtvimproadString//"TBD"
    stream_var_name_improad(isturb_HD)  = urbanalbtvimproadString//"HD"
    stream_var_name_improad(isturb_MD)  = urbanalbtvimproadString//"MD"   
    fldList = ""
    do ifield = isturb_MIN, isturb_MAX
       call shr_string_listAppend( fldList, stream_var_name_improad(ifield) )
    end do
    
   call shr_strdata_create(this%sdat_urbanalbtvimproad,name="clmurbanalbtvimproad",     &
        pio_subsystem=pio_subsystem,                   &
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,                &
        yearFirst=stream_year_first_urbanalbtvimproad,        &
        yearLast=stream_year_last_urbanalbtvimproad,          &
        yearAlign=model_year_align_urbanalbtvimproad,         &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_urbanalbtvimproad),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &
        domAreaName='area',                            &
        domMaskName='LANDMASK',                        &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_urbanalbtvimproad)/) , &
        fldListFile=fldList,                           &
        fldListModel=fldList,                          &
        fillalgo='none',                               &
        mapalgo=urbanalbtvimproadmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo=urbanalbtvimproad_tintalgo,                     &
        taxmode='extend'                               )

   if (masterproc) then
      call shr_strdata_print(this%sdat_urbanalbtvimproad,'urban time varying improad albedo')
   endif
   
  end subroutine urbanalbtvimproad_init

  !==============================================================================
 
  subroutine urbanalbtvimproad_interp(this, bounds)
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying improad albedo.
    ! 
    ! !USES:
    use clm_time_manager , only : get_curr_date
    use spmdMod          , only : mpicom
    use clm_instur       , only : urban_valid
    use shr_strdata_mod  , only : shr_strdata_advance
    ! 
    ! !ARGUMENTS:
    ! 
    class(urbanalbtv_type)           :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !
    logical :: found
    integer :: l       ! landunit
    integer :: glun
    integer :: ib
    integer :: ig
    integer :: g
    integer :: ip      
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    ! 
    !-----------------------------------------------------------------------
    ! 
    ! Advance sdat stream
    !
    call get_curr_date(year, mon, day, sec)
    !
    ! packing the date into an integer
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_urbanalbtvimproad, mcdate, sec, mpicom,'urbanalbtvimproad')

    ! Determine this%dyn_alb_improad for all landunits
    do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
             glun = lun%gridcell(l)
             ip = mct_aVect_indexRA(this%sdat_urbanalbtvimproad%avs(1),trim(stream_var_name_improad(lun%itype(l))))
             ig = 0
             do g = bounds%begg, bounds%endg
                ig = ig+1
                if (g==glun) exit
             end do    
             do ib = 1,numrad        
                this%dyn_alb_improad_dir(l,ib) = this%sdat_urbanalbtvimproad%avs(1)%rAttr(ip,ig)
                this%dyn_alb_improad_dif(l,ib) = this%sdat_urbanalbtvimproad%avs(1)%rAttr(ip,ig) 
             end do                     
       else 
             do ib = 1,numrad
                this%dyn_alb_improad_dir(l,ib) = spval
                this%dyn_alb_improad_dif(l,ib) = spval
             end do        
       end if
    end do
    
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
             glun  = lun%gridcell(l)
          ! 
          ! Determine vector index corresponding to glun
          ! 
             ig = 0
             do g = bounds%begg,bounds%endg
                ig = ig+1
                if (g == glun) exit
             end do
             
             do ib = 1,numrad
                if ( .not. urban_valid(g) .or. (this%dyn_alb_improad_dir(l,ib) <= 0._r8) .or. (this%dyn_alb_improad_dif(l,ib) <= 0._r8)) then
                   found = .true.
                   gindx = g
                   lindx = l
                   exit
                end if 
            end do     
       end if
    end do
    
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban improad data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(l)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'dyn_alb_improad_dir:  ',this%dyn_alb_improad_dir(lindx,:)
       write(iulog,*)'dyn_alb_improad_dif:  ',this%dyn_alb_improad_dif(lindx,:)
       call endrun(msg=errmsg(sourcefile, __LINE__))
    end if
  end subroutine urbanalbtvimproad_interp
  !-----------------------------------------------------------------------
  
  subroutine urbanalbtvwall_init(this, bounds)
   ! !DESCRIPTION:
   ! !USES:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use spmdMod          , only : masterproc, mpicom, comp_id
   use fileutils        , only : getavu, relavu
   use shr_mpi_mod      , only : shr_mpi_bcast
   use shr_string_mod   , only : shr_string_listAppend
   use shr_strdata_mod  , only : shr_strdata_create, shr_strdata_print
   use decompMod        , only : gsmap_lnd_gdc2glo
   use domainMod        , only : ldomain
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use landunit_varcon  , only : isturb_TBD, isturb_HD, isturb_MD
   use clm_varctl       , only : NLFilename_in
   ! !ARGUMENTS:
   implicit none
   class(urbanalbtv_type)         :: this
   type(bounds_type), intent(in)  :: bounds
   ! 
   ! !LOCAL VARIABLES:
   integer            :: begl, endl                                  ! landunits
   integer            :: ifield                                      ! field index
   integer            :: stream_year_first_urbanalbtvwall            ! first year in urban wall albedo tv stream to use
   integer            :: stream_year_last_urbanalbtvwall             ! last year in urban wall albedo tv stream to use
   integer            :: model_year_align_urbanalbtvwall             ! align stream_year_first_urbantvwall with this model year
   integer            :: nu_nml                                      ! unit for namelist file 
   integer            :: nml_error                                   ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                                     ! domain information
   character(len=CL)  :: stream_fldFileName_urbanalbtvwall           ! urban wall albedo time-varying streams filename
   character(len=CL)  :: urbanalbtvwallmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbanalbtvwall_tintalgo = 'linear'          ! time interpolation alogrithm 
   character(*), parameter :: urbanalbtvwallString = "dyn_alb_wall_" ! base string for field string
   character(SHR_KIND_CL)  :: fldList                                ! field string
   character(*), parameter :: subName = "('urbanalbtvwall_init')"
   character(*), parameter :: F00 = "('(urbanalbtvwall_init) ',4a)"
   !-----------------------------------------------------------------------
   namelist /urbanalbtvwall_streams/       &
        stream_year_first_urbanalbtvwall,  &  
        stream_year_last_urbanalbtvwall,   &  
        model_year_align_urbanalbtvwall,   &  
        urbanalbtvwallmapalgo,             &  
        stream_fldFileName_urbanalbtvwall, &      
        urbanalbtvwall_tintalgo  
   !-----------------------------------------------------------------------                      
   begl = bounds%begl; endl = bounds%endl
   ! Default values for namelist
   stream_year_first_urbanalbtvwall  = 1      ! first year in stream to use
   stream_year_last_urbanalbtvwall   = 1      ! last  year in stream to use
   model_year_align_urbanalbtvwall   = 1      ! align stream_year_first_urbanalbtvwall with this model year
   stream_fldFileName_urbanalbtvwall = ' '
   
   ! Read urbanalbtv_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( newunit=nu_nml, file=trim(NLFilename_in), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbanalbtvwall_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbanalbtvwall_streams,iostat=nml_error) 
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbanalbtvwall_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
          call endrun(subname // ':: ERROR finding urbanalbtvwall_streams namelist')   
      end if
      close(nu_nml)
      call relavu( nu_nml )
    endif

    call shr_mpi_bcast(stream_year_first_urbanalbtvwall  , mpicom)
    call shr_mpi_bcast(stream_year_last_urbanalbtvwall   , mpicom)
    call shr_mpi_bcast(model_year_align_urbanalbtvwall   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_urbanalbtvwall , mpicom)
    call shr_mpi_bcast(urbanalbtvwall_tintalgo           , mpicom)
   

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'Attemping to read time varying urban wall albedo parameters......'
       write(iulog,'(a)') 'urbanalbtvwall_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbanalbtvwall  = ',stream_year_first_urbanalbtvwall
       write(iulog,'(a,i8)') '  stream_year_last_urbanalbtvwall   = ',stream_year_last_urbanalbtvwall
       write(iulog,'(a,i8)') '  model_year_align_urbanalbtvwall   = ',model_year_align_urbanalbtvwall
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbanalbtvwall = ',stream_fldFileName_urbanalbtvwall
       write(iulog,'(a,a)' ) '  urbanalbtvwall_tintalgo           = ',urbanalbtvwall_tintalgo
       write(iulog,*) 'Read in urbanalbtvwall namelist from:',trim(NLFilename_in)
    endif
    
    ! Initialize the cdeps data type this%sdat_urbanalbtv
    call clm_domain_mct (bounds, dom_clm)
    
    ! create the field list for urban albedo fields
    stream_var_name_wall(:)          = "NOT_SET"
    stream_var_name_wall(isturb_TBD) = urbanalbtvwallString//"TBD"
    stream_var_name_wall(isturb_HD)  = urbanalbtvwallString//"HD"
    stream_var_name_wall(isturb_MD)  = urbanalbtvwallString//"MD"   
    fldList = ""
    do ifield = isturb_MIN, isturb_MAX
       call shr_string_listAppend( fldList, stream_var_name_wall(ifield) )
    end do
    
   call shr_strdata_create(this%sdat_urbanalbtvwall,name="clmurbanalbtvwall",     &
        pio_subsystem=pio_subsystem,                   &
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,         &
        yearFirst=stream_year_first_urbanalbtvwall,        &
        yearLast=stream_year_last_urbanalbtvwall,          &
        yearAlign=model_year_align_urbanalbtvwall,         &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_urbanalbtvwall),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &
        domAreaName='area',                            &
        domMaskName='LANDMASK',                        &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_urbanalbtvwall)/) , &
        fldListFile=fldList,                           &
        fldListModel=fldList,                          &
        fillalgo='none',                               &
        mapalgo=urbanalbtvwallmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo=urbanalbtvwall_tintalgo,                     &
        taxmode='extend'                               )

   if (masterproc) then
      call shr_strdata_print(this%sdat_urbanalbtvwall,'urban time varying wall albedo')
   endif
 
  end subroutine urbanalbtvwall_init

  !==============================================================================
 
  subroutine urbanalbtvwall_interp(this, bounds)
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying albedo.
    ! 
    ! !USES:
    use clm_time_manager , only : get_curr_date
    use spmdMod          , only : mpicom
    use clm_instur       , only : urban_valid
    use shr_strdata_mod  , only : shr_strdata_advance
    ! 
    ! !ARGUMENTS:
    ! 
    class(urbanalbtv_type)           :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !
    logical :: found
    integer :: l       ! landunit
    integer :: glun
    integer :: ib
    integer :: ig
    integer :: g
    integer :: ip      
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    ! 
    !-----------------------------------------------------------------------
    ! 
    ! Advance sdat stream
    !
    call get_curr_date(year, mon, day, sec)
    !
    ! packing the date into an integer
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_urbanalbtvwall, mcdate, sec, mpicom,'urbanalbtvwall')

    ! Determine this%dyn_alb_wall for all landunits
    do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
             glun = lun%gridcell(l)
             ip = mct_aVect_indexRA(this%sdat_urbanalbtvwall%avs(1),trim(stream_var_name_wall(lun%itype(l))))
             ig = 0
             do g = bounds%begg, bounds%endg
                ig = ig+1
                if (g==glun) exit
             end do     
             do ib = 1, numrad       
                this%dyn_alb_wall_dir(l,ib) = this%sdat_urbanalbtvwall%avs(1)%rAttr(ip,ig)
                this%dyn_alb_wall_dif(l,ib) = this%sdat_urbanalbtvwall%avs(1)%rAttr(ip,ig)           
             end do           
       else
             do ib = 1,numrad
                this%dyn_alb_wall_dir(l,ib) = spval
                this%dyn_alb_wall_dif(l,ib) = spval    
             end do    
       end if
    end do
    
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
             glun  = lun%gridcell(l)
          ! 
          ! Determine vector index corresponding to glun
          ! 
             ig = 0
             do g = bounds%begg,bounds%endg
                ig = ig+1
                if (g == glun) exit
             end do
             
             do ib = 1,numrad
                if ( .not. urban_valid(g) .or. (this%dyn_alb_wall_dir(l,ib) <= 0._r8) .or. (this%dyn_alb_wall_dif(l,ib) <= 0._r8)) then
                   found = .true.
                   gindx = g
                   lindx = l
                   exit
                end if  
             end do   
       end if
    end do
    
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban wall data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(l)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'dyn_alb_wall_dir:  ',this%dyn_alb_wall_dir(lindx,:)
       write(iulog,*)'dyn_alb_wall_dif:  ',this%dyn_alb_wall_dif(lindx,:)
       call endrun(msg=errmsg(sourcefile, __LINE__))
    end if
  end subroutine urbanalbtvwall_interp
end module UrbanDynAlbMod
