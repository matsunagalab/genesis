!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> C language interface for trj_analysis (zerocopy)
!! @brief   C language interface of trj_analysis
!! @authors Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_c_mod
  use, intrinsic :: iso_c_binding
  use s_trajectories_c_mod
  use trj_impl_mod
  use trj_source_mod

  use trajectory_str_mod
  use molecules_str_mod
  use string_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: trj_analysis_c
  public :: trj_analysis_lazy_c
  public :: trj_analysis_com_c
  public :: trj_analysis_com_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Type          t_trj_lazy_ctx
  !> @brief        State shared between trj_analysis_lazy_c and its guarded
  !!               body trj_analysis_lazy_body (see run_guarded in error_mod)
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  type :: t_trj_lazy_ctx
    ! inputs (copied from the bind(C) arguments)
    character(MaxFilename) :: filename_f = ''
    integer     :: trj_type = 0
    integer     :: dcd_natom_expected = 0
    type(c_ptr) :: source_selection_ptr = c_null_ptr
    integer     :: n_source_selection = 0
    type(c_ptr) :: dist_list_ptr = c_null_ptr
    integer     :: n_dist = 0
    type(c_ptr) :: angl_list_ptr = c_null_ptr
    integer     :: n_angl = 0
    type(c_ptr) :: tors_list_ptr = c_null_ptr
    integer     :: n_tors = 0
    integer     :: n_atoms = 0
    integer     :: ana_period = 1
    integer     :: n_frame = 0
    type(c_ptr) :: dist_ptr = c_null_ptr
    integer     :: dist_size = 0
    type(c_ptr) :: angl_ptr = c_null_ptr
    integer     :: angl_size = 0
    type(c_ptr) :: tors_ptr = c_null_ptr
    integer     :: tors_size = 0
    ! outputs (copied back to the bind(C) arguments by the wrapper)
    integer     :: nstru_out = 0
    integer     :: dcd_nframe_out = 0
    integer     :: dcd_natom_out = 0
    ! error state and resources released by the wrapper
    type(s_error)        :: err
    type(s_trj_source)   :: source
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)
  end type t_trj_lazy_ctx

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Type          t_trj_com_lazy_ctx
  !> @brief        State shared between trj_analysis_com_lazy_c and its
  !!               guarded body trj_analysis_com_lazy_body
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  type :: t_trj_com_lazy_ctx
    ! inputs (copied from the bind(C) arguments)
    character(MaxFilename) :: filename_f = ''
    integer     :: trj_type = 0
    integer     :: dcd_natom_expected = 0
    type(c_ptr) :: source_selection_ptr = c_null_ptr
    integer     :: n_source_selection = 0
    type(c_ptr) :: mass_ptr = c_null_ptr
    integer     :: n_atoms = 0
    integer     :: ana_period = 1
    integer     :: n_frame = 0
    type(c_ptr) :: dist_list_ptr = c_null_ptr
    integer     :: n_dist = 0
    type(c_ptr) :: angl_list_ptr = c_null_ptr
    integer     :: n_angl = 0
    type(c_ptr) :: tors_list_ptr = c_null_ptr
    integer     :: n_tors = 0
    type(c_ptr) :: cdis_atoms_ptr = c_null_ptr
    integer     :: n_cdis_atoms = 0
    type(c_ptr) :: cdis_offsets_ptr = c_null_ptr
    integer     :: n_cdis_offsets = 0
    type(c_ptr) :: cdis_pairs_ptr = c_null_ptr
    integer     :: n_cdis = 0
    type(c_ptr) :: cang_atoms_ptr = c_null_ptr
    integer     :: n_cang_atoms = 0
    type(c_ptr) :: cang_offsets_ptr = c_null_ptr
    integer     :: n_cang_offsets = 0
    type(c_ptr) :: cang_triplets_ptr = c_null_ptr
    integer     :: n_cang = 0
    type(c_ptr) :: ctor_atoms_ptr = c_null_ptr
    integer     :: n_ctor_atoms = 0
    type(c_ptr) :: ctor_offsets_ptr = c_null_ptr
    integer     :: n_ctor_offsets = 0
    type(c_ptr) :: ctor_quads_ptr = c_null_ptr
    integer     :: n_ctor = 0
    type(c_ptr) :: dist_ptr = c_null_ptr
    integer     :: dist_size = 0
    type(c_ptr) :: angl_ptr = c_null_ptr
    integer     :: angl_size = 0
    type(c_ptr) :: tors_ptr = c_null_ptr
    integer     :: tors_size = 0
    type(c_ptr) :: cdis_result_ptr = c_null_ptr
    integer     :: cdis_size = 0
    type(c_ptr) :: cang_result_ptr = c_null_ptr
    integer     :: cang_size = 0
    type(c_ptr) :: ctor_result_ptr = c_null_ptr
    integer     :: ctor_size = 0
    ! outputs (copied back to the bind(C) arguments by the wrapper)
    integer     :: nstru_out = 0
    integer     :: dcd_nframe_out = 0
    integer     :: dcd_natom_out = 0
    ! error state and resources released by the wrapper
    type(s_error)        :: err
    type(s_trj_source)   :: source
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)
    integer, allocatable :: cdis_atoms_copy(:)
    integer, allocatable :: cdis_offsets_copy(:)
    integer, allocatable :: cdis_pairs_copy(:)
    integer, allocatable :: cang_atoms_copy(:)
    integer, allocatable :: cang_offsets_copy(:)
    integer, allocatable :: cang_triplets_copy(:)
    integer, allocatable :: ctor_atoms_copy(:)
    integer, allocatable :: ctor_offsets_copy(:)
    integer, allocatable :: ctor_quads_copy(:)
  end type t_trj_com_lazy_ctx

  private :: t_trj_lazy_ctx, t_trj_com_lazy_ctx
  private :: trj_analysis_lazy_body, trj_analysis_com_lazy_body

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_c
  !> @brief        Trajectory analysis (zerocopy, pre-allocated results)
  !! @authors      Claude Code
  !! @param[in]    s_trajes_c   : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list_ptr: pointer to distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list_ptr: pointer to angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list_ptr: pointer to torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[in]    dist_ptr     : pre-allocated distance results (n_dist, n_frames)
  !! @param[in]    dist_size    : total size of dist array (n_dist * n_frames)
  !! @param[in]    angl_ptr     : pre-allocated angle results (n_angl, n_frames)
  !! @param[in]    angl_size    : total size of angl array
  !! @param[in]    tors_ptr     : pre-allocated torsion results (n_tors, n_frames)
  !! @param[in]    tors_size    : total size of tors array
  !! @param[out]   nstru_out    : actual number of structures analyzed
  !! @param[out]   status       : error status
  !! @param[out]   msg          : error message
  !! @param[in]    msglen       : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_c(s_trajes_c, ana_period, &
                            dist_list_ptr, n_dist, &
                            angl_list_ptr, n_angl, &
                            tors_list_ptr, n_tors, &
                            dist_ptr, dist_size, &
                            angl_ptr, angl_size, &
                            tors_ptr, tors_size, &
                            nstru_out, status, msg, msglen) &
        bind(C, name="trj_analysis_c")
    implicit none

    ! Arguments
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors
    type(c_ptr), value :: dist_ptr
    integer(c_int), value :: dist_size
    type(c_ptr), value :: angl_ptr
    integer(c_int), value :: angl_size
    type(c_ptr), value :: tors_ptr
    integer(c_int), value :: tors_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_trj_source) :: source
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    real(wp), pointer :: dist_f(:,:)
    real(wp), pointer :: angl_f(:,:)
    real(wp), pointer :: tors_f(:,:)
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)
    integer :: n_frames, nstru_local

    ! Initialize
    status = 0
    nstru_out = 0
    n_frames = s_trajes_c%nframe / ana_period

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    write(MsgOut,'(A)') '[STEP1] Trajectory Analysis'
    write(MsgOut,'(A)') ' '

    ! Convert C pointers to Fortran arrays (input lists)
    if (n_dist > 0 .and. c_associated(dist_list_ptr)) then
      call C_F_POINTER(dist_list_ptr, dist_list_f, [2, n_dist])
      allocate(dist_list_copy(2, n_dist))
      dist_list_copy = dist_list_f
    else
      allocate(dist_list_copy(2, 0))
    end if

    if (n_angl > 0 .and. c_associated(angl_list_ptr)) then
      call C_F_POINTER(angl_list_ptr, angl_list_f, [3, n_angl])
      allocate(angl_list_copy(3, n_angl))
      angl_list_copy = angl_list_f
    else
      allocate(angl_list_copy(3, 0))
    end if

    if (n_tors > 0 .and. c_associated(tors_list_ptr)) then
      call C_F_POINTER(tors_list_ptr, tors_list_f, [4, n_tors])
      allocate(tors_list_copy(4, n_tors))
      tors_list_copy = tors_list_f
    else
      allocate(tors_list_copy(4, 0))
    end if

    ! Create views of pre-allocated result arrays
    if (n_dist > 0 .and. c_associated(dist_ptr) .and. dist_size > 0) then
      call C_F_POINTER(dist_ptr, dist_f, [n_dist, n_frames])
    else
      nullify(dist_f)
    end if

    if (n_angl > 0 .and. c_associated(angl_ptr) .and. angl_size > 0) then
      call C_F_POINTER(angl_ptr, angl_f, [n_angl, n_frames])
    else
      nullify(angl_f)
    end if

    if (n_tors > 0 .and. c_associated(tors_ptr) .and. tors_size > 0) then
      call C_F_POINTER(tors_ptr, tors_f, [n_tors, n_frames])
    else
      nullify(tors_f)
    end if

    ! Build a memory-backed source and run the shared analysis loop.
    call init_source_memory(source, s_trajes_c%coords, s_trajes_c%pbc_boxes, &
                            s_trajes_c%natom, s_trajes_c%nframe, ana_period)

    call analyze(source, &
                 dist_list_copy, n_dist, &
                 angl_list_copy, n_angl, &
                 tors_list_copy, n_tors, &
                 dist_f, angl_f, tors_f, nstru_local)

    call finalize_source(source)

    nstru_out = nstru_local

    ! Cleanup local arrays
    deallocate(dist_list_copy)
    deallocate(angl_list_copy)
    deallocate(tors_list_copy)

  end subroutine trj_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_lazy_c
  !> @brief        Trajectory analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !! @param[in]    dcd_filename       : DCD file path (C string)
  !! @param[in]    filename_len       : length of filename
  !! @param[in]    trj_type           : trajectory type (1=COOR, 2=COOR+BOX)
  !! @param[in]    dcd_natom_expected : physical atom count in the DCD file
  !! @param[in]    source_selection_ptr : selected atom indices (1-indexed DCD)
  !! @param[in]    n_source_selection : number of selected atoms
  !! @param[in]    dist_list_ptr      : distance atom pairs (2, n_dist)
  !! @param[in]    n_dist             : number of distance measurements
  !! @param[in]    angl_list_ptr      : angle atom triplets (3, n_angl)
  !! @param[in]    n_angl             : number of angle measurements
  !! @param[in]    tors_list_ptr      : torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors             : number of torsion measurements
  !! @param[in]    n_atoms            : number of selected (logical) atoms
  !! @param[in]    ana_period         : analysis period
  !! @param[in]    n_frame            : number of result frames (columns)
  !! @param[in]    dist_ptr           : pre-allocated distance results
  !! @param[in]    dist_size          : total size of dist array
  !! @param[in]    angl_ptr           : pre-allocated angle results
  !! @param[in]    angl_size          : total size of angl array
  !! @param[in]    tors_ptr           : pre-allocated torsion results
  !! @param[in]    tors_size          : total size of tors array
  !! @param[out]   nstru_out          : actual number of structures analyzed
  !! @param[out]   dcd_nframe_out     : total frames in DCD
  !! @param[out]   dcd_natom_out      : atoms per frame in DCD
  !! @param[out]   status             : error status
  !! @param[out]   msg                : error message
  !! @param[in]    msglen             : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                 dcd_natom_expected, source_selection_ptr, &
                                 n_source_selection, &
                                 dist_list_ptr, n_dist, &
                                 angl_list_ptr, n_angl, &
                                 tors_list_ptr, n_tors, &
                                 n_atoms, ana_period, n_frame, &
                                 dist_ptr, dist_size, &
                                 angl_ptr, angl_size, &
                                 tors_ptr, tors_size, &
                                 nstru_out, dcd_nframe_out, dcd_natom_out, &
                                 status, msg, msglen) &
        bind(C, name="trj_analysis_lazy_c")
    implicit none

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    integer(c_int), value :: dcd_natom_expected
    type(c_ptr), value :: source_selection_ptr
    integer(c_int), value :: n_source_selection
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    integer(c_int), value :: n_frame
    type(c_ptr), value :: dist_ptr
    integer(c_int), value :: dist_size
    type(c_ptr), value :: angl_ptr
    integer(c_int), value :: angl_size
    type(c_ptr), value :: tors_ptr
    integer(c_int), value :: tors_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_trj_lazy_ctx), target :: c

    ! Hand every input to the guarded body through the context
    call c_filename_to_fortran(dcd_filename, filename_len, c%filename_f)
    c%trj_type             = trj_type
    c%dcd_natom_expected   = dcd_natom_expected
    c%source_selection_ptr = source_selection_ptr
    c%n_source_selection   = n_source_selection
    c%dist_list_ptr        = dist_list_ptr
    c%n_dist               = n_dist
    c%angl_list_ptr        = angl_list_ptr
    c%n_angl               = n_angl
    c%tors_list_ptr        = tors_list_ptr
    c%n_tors               = n_tors
    c%n_atoms              = n_atoms
    c%ana_period           = ana_period
    c%n_frame              = n_frame
    c%dist_ptr             = dist_ptr
    c%dist_size            = dist_size
    c%angl_ptr             = angl_ptr
    c%angl_size            = angl_size
    c%tors_ptr             = tors_ptr
    c%tors_size            = tors_size

    ! Run the body under the library-mode error guard: reading the DCD may
    ! call error_msg -> exit(1) in CLI mode, which would kill the host Python
    ! process (see run_guarded in error_mod).
    call run_guarded(trj_analysis_lazy_body, c, c%err, status, msg, msglen)

    nstru_out      = c%nstru_out
    dcd_nframe_out = c%dcd_nframe_out
    dcd_natom_out  = c%dcd_natom_out

    ! Release what the body acquired; this also runs after a longjmp. The
    ! allocatable components of c are freed automatically on return.
    call finalize_source(c%source)

  end subroutine trj_analysis_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_lazy_body
  !> @brief        Guarded body of trj_analysis_lazy_c
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_trj_lazy_ctx
  !
  !  Reports failures only through c%err; the wrapper converts them to the C
  !  status/message pair and releases the resources recorded in the context.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_lazy_body(ctx) bind(C, name="trj_analysis_lazy_body")
    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_trj_lazy_ctx), pointer :: c
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    integer, pointer :: source_selection_f(:)
    real(wp), pointer :: dist_f(:,:)
    real(wp), pointer :: angl_f(:,:)
    real(wp), pointer :: tors_f(:,:)
    integer :: nstru_local, init_status

    call c_f_pointer(ctx, c)

    ! Validate inputs
    if (c%n_source_selection /= c%n_atoms .or. &
        .not. c_associated(c%source_selection_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_lazy_c: invalid source selection")
      return
    end if

    if (c%n_frame <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_lazy_c: n_frame must be positive")
      return
    end if

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Create zero-copy views of the selection and result arrays
    call C_F_POINTER(c%source_selection_ptr, source_selection_f, &
                     [c%n_source_selection])
    if (any(source_selection_f < 1) .or. &
        any(source_selection_f > c%dcd_natom_expected)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_lazy_c: source selection out of range")
      return
    end if

    ! Convert measurement lists (already 1-indexed in selected space)
    if (c%n_dist > 0 .and. c_associated(c%dist_list_ptr)) then
      call C_F_POINTER(c%dist_list_ptr, dist_list_f, [2, c%n_dist])
      allocate(c%dist_list_copy(2, c%n_dist))
      c%dist_list_copy = dist_list_f
    else
      allocate(c%dist_list_copy(2, 0))
    end if

    if (c%n_angl > 0 .and. c_associated(c%angl_list_ptr)) then
      call C_F_POINTER(c%angl_list_ptr, angl_list_f, [3, c%n_angl])
      allocate(c%angl_list_copy(3, c%n_angl))
      c%angl_list_copy = angl_list_f
    else
      allocate(c%angl_list_copy(3, 0))
    end if

    if (c%n_tors > 0 .and. c_associated(c%tors_list_ptr)) then
      call C_F_POINTER(c%tors_list_ptr, tors_list_f, [4, c%n_tors])
      allocate(c%tors_list_copy(4, c%n_tors))
      c%tors_list_copy = tors_list_f
    else
      allocate(c%tors_list_copy(4, 0))
    end if

    ! Measurement indices must stay within the selected atom range
    if ((c%n_dist > 0 .and. (any(c%dist_list_copy < 1) .or. &
                             any(c%dist_list_copy > c%n_atoms))) .or. &
        (c%n_angl > 0 .and. (any(c%angl_list_copy < 1) .or. &
                             any(c%angl_list_copy > c%n_atoms))) .or. &
        (c%n_tors > 0 .and. (any(c%tors_list_copy < 1) .or. &
                             any(c%tors_list_copy > c%n_atoms)))) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_lazy_c: measurement index out of range")
      return
    end if

    ! Create views of pre-allocated result arrays
    if (c%n_dist > 0 .and. c_associated(c%dist_ptr) .and. c%dist_size > 0) then
      call C_F_POINTER(c%dist_ptr, dist_f, [c%n_dist, c%n_frame])
    else
      nullify(dist_f)
    end if

    if (c%n_angl > 0 .and. c_associated(c%angl_ptr) .and. c%angl_size > 0) then
      call C_F_POINTER(c%angl_ptr, angl_f, [c%n_angl, c%n_frame])
    else
      nullify(angl_f)
    end if

    if (c%n_tors > 0 .and. c_associated(c%tors_ptr) .and. c%tors_size > 0) then
      call C_F_POINTER(c%tors_ptr, tors_f, [c%n_tors, c%n_frame])
    else
      nullify(tors_f)
    end if

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for Trj Analysis'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(c%source, trim(c%filename_f), c%trj_type, &
                              c%ana_period, source_selection_f, &
                              c%n_source_selection, init_status)
    if (init_status /= 0) then
      call error_set(c%err, init_status, &
                     "trj_analysis_lazy_c: unable to initialize DCD source")
      return
    end if

    ! Return DCD info (also on the atom-count error below)
    c%dcd_nframe_out = c%source%dcd_nframe
    c%dcd_natom_out  = c%source%dcd_natom

    ! Check atom count
    if (c%source%dcd_natom /= c%dcd_natom_expected) then
      call error_set(c%err, ERROR_ATOM_COUNT, &
                     "trj_analysis_lazy_c: atom count mismatch")
      return
    end if

    ! Run the shared analysis loop (lazy loading via source abstraction)
    write(MsgOut,'(A)') '[STEP2] Trajectory Analysis (lazy loading)'
    write(MsgOut,'(A)') ' '

    call analyze(c%source, &
                 c%dist_list_copy, c%n_dist, &
                 c%angl_list_copy, c%n_angl, &
                 c%tors_list_copy, c%n_tors, &
                 dist_f, angl_f, tors_f, nstru_local)

    c%nstru_out = nstru_local

  end subroutine trj_analysis_lazy_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_com_c
  !> @brief        Trajectory analysis with COM (zerocopy, pre-allocated)
  !! @authors      Claude Code
  !! @note         All result arrays are pre-allocated by Python
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_com_c(mass_ptr, n_atoms, &
                                s_trajes_c, ana_period, &
                                dist_list_ptr, n_dist, &
                                angl_list_ptr, n_angl, &
                                tors_list_ptr, n_tors, &
                                cdis_atoms_ptr, n_cdis_atoms, &
                                cdis_offsets_ptr, n_cdis_offsets, &
                                cdis_pairs_ptr, n_cdis, &
                                cang_atoms_ptr, n_cang_atoms, &
                                cang_offsets_ptr, n_cang_offsets, &
                                cang_triplets_ptr, n_cang, &
                                ctor_atoms_ptr, n_ctor_atoms, &
                                ctor_offsets_ptr, n_ctor_offsets, &
                                ctor_quads_ptr, n_ctor, &
                                dist_ptr, dist_size, &
                                angl_ptr, angl_size, &
                                tors_ptr, tors_size, &
                                cdis_result_ptr, cdis_size, &
                                cang_result_ptr, cang_size, &
                                ctor_result_ptr, ctor_size, &
                                nstru_out, status, msg, msglen) &
        bind(C, name="trj_analysis_com_c")
    implicit none

    ! Arguments - mass array
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms

    ! Arguments - trajectories
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period

    ! Arguments - atom-based measurements
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors

    ! Arguments - COM distance
    type(c_ptr), value :: cdis_atoms_ptr
    integer(c_int), value :: n_cdis_atoms
    type(c_ptr), value :: cdis_offsets_ptr
    integer(c_int), value :: n_cdis_offsets
    type(c_ptr), value :: cdis_pairs_ptr
    integer(c_int), value :: n_cdis

    ! Arguments - COM angle
    type(c_ptr), value :: cang_atoms_ptr
    integer(c_int), value :: n_cang_atoms
    type(c_ptr), value :: cang_offsets_ptr
    integer(c_int), value :: n_cang_offsets
    type(c_ptr), value :: cang_triplets_ptr
    integer(c_int), value :: n_cang

    ! Arguments - COM torsion
    type(c_ptr), value :: ctor_atoms_ptr
    integer(c_int), value :: n_ctor_atoms
    type(c_ptr), value :: ctor_offsets_ptr
    integer(c_int), value :: n_ctor_offsets
    type(c_ptr), value :: ctor_quads_ptr
    integer(c_int), value :: n_ctor

    ! Pre-allocated output arrays
    type(c_ptr), value :: dist_ptr
    integer(c_int), value :: dist_size
    type(c_ptr), value :: angl_ptr
    integer(c_int), value :: angl_size
    type(c_ptr), value :: tors_ptr
    integer(c_int), value :: tors_size
    type(c_ptr), value :: cdis_result_ptr
    integer(c_int), value :: cdis_size
    type(c_ptr), value :: cang_result_ptr
    integer(c_int), value :: cang_size
    type(c_ptr), value :: ctor_result_ptr
    integer(c_int), value :: ctor_size

    ! Output
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables - Fortran pointers from C
    real(wp), pointer :: mass_f(:)
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    integer, pointer :: cdis_atoms_f(:)
    integer, pointer :: cdis_offsets_f(:)
    integer, pointer :: cdis_pairs_f(:)
    integer, pointer :: cang_atoms_f(:)
    integer, pointer :: cang_offsets_f(:)
    integer, pointer :: cang_triplets_f(:)
    integer, pointer :: ctor_atoms_f(:)
    integer, pointer :: ctor_offsets_f(:)
    integer, pointer :: ctor_quads_f(:)

    ! Pre-allocated result arrays (views)
    real(wp), pointer :: distance_f(:,:)
    real(wp), pointer :: angle_f(:,:)
    real(wp), pointer :: torsion_f(:,:)
    real(wp), pointer :: cdis_f(:,:)
    real(wp), pointer :: cang_f(:,:)
    real(wp), pointer :: ctor_f(:,:)

    ! Local copies for 1-indexed conversion
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)
    integer, allocatable :: cdis_atoms_copy(:)
    integer, allocatable :: cdis_offsets_copy(:)
    integer, allocatable :: cdis_pairs_copy(:)
    integer, allocatable :: cang_atoms_copy(:)
    integer, allocatable :: cang_offsets_copy(:)
    integer, allocatable :: cang_triplets_copy(:)
    integer, allocatable :: ctor_atoms_copy(:)
    integer, allocatable :: ctor_offsets_copy(:)
    integer, allocatable :: ctor_quads_copy(:)

    type(s_trj_source) :: source
    integer :: nstru_local
    integer :: num_out_frame
    integer :: n_cdis_groups, n_cang_groups, n_ctor_groups

    ! Initialize
    status = 0
    nstru_out = 0
    if (msglen > 0) msg(1) = c_null_char

    num_out_frame = s_trajes_c%nframe / ana_period

    ! Calculate number of groups
    n_cdis_groups = n_cdis_offsets - 1
    n_cang_groups = n_cang_offsets - 1
    n_ctor_groups = n_ctor_offsets - 1

    ! Get mass array
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])

    ! Create views of pre-allocated result arrays
    if (n_dist > 0 .and. c_associated(dist_ptr)) then
      call C_F_POINTER(dist_ptr, distance_f, [n_dist, num_out_frame])
    else
      allocate(distance_f(1,1))  ! Dummy
    end if

    if (n_angl > 0 .and. c_associated(angl_ptr)) then
      call C_F_POINTER(angl_ptr, angle_f, [n_angl, num_out_frame])
    else
      allocate(angle_f(1,1))  ! Dummy
    end if

    if (n_tors > 0 .and. c_associated(tors_ptr)) then
      call C_F_POINTER(tors_ptr, torsion_f, [n_tors, num_out_frame])
    else
      allocate(torsion_f(1,1))  ! Dummy
    end if

    if (n_cdis > 0 .and. c_associated(cdis_result_ptr)) then
      call C_F_POINTER(cdis_result_ptr, cdis_f, [n_cdis, num_out_frame])
    else
      allocate(cdis_f(1,1))  ! Dummy
    end if

    if (n_cang > 0 .and. c_associated(cang_result_ptr)) then
      call C_F_POINTER(cang_result_ptr, cang_f, [n_cang, num_out_frame])
    else
      allocate(cang_f(1,1))  ! Dummy
    end if

    if (n_ctor > 0 .and. c_associated(ctor_result_ptr)) then
      call C_F_POINTER(ctor_result_ptr, ctor_f, [n_ctor, num_out_frame])
    else
      allocate(ctor_f(1,1))  ! Dummy
    end if

    ! Convert atom-based lists (already 1-indexed from Python)
    if (n_dist > 0) then
      call C_F_POINTER(dist_list_ptr, dist_list_f, [2, n_dist])
      allocate(dist_list_copy(2, n_dist))
      dist_list_copy = dist_list_f
    else
      allocate(dist_list_copy(2, 1))
      dist_list_copy = 0
    end if

    if (n_angl > 0) then
      call C_F_POINTER(angl_list_ptr, angl_list_f, [3, n_angl])
      allocate(angl_list_copy(3, n_angl))
      angl_list_copy = angl_list_f
    else
      allocate(angl_list_copy(3, 1))
      angl_list_copy = 0
    end if

    if (n_tors > 0) then
      call C_F_POINTER(tors_list_ptr, tors_list_f, [4, n_tors])
      allocate(tors_list_copy(4, n_tors))
      tors_list_copy = tors_list_f
    else
      allocate(tors_list_copy(4, 1))
      tors_list_copy = 0
    end if

    ! Convert COM distance arrays (atoms already 1-indexed from Python)
    if (n_cdis_atoms > 0 .and. n_cdis > 0) then
      call C_F_POINTER(cdis_atoms_ptr, cdis_atoms_f, [n_cdis_atoms])
      call C_F_POINTER(cdis_offsets_ptr, cdis_offsets_f, [n_cdis_offsets])
      call C_F_POINTER(cdis_pairs_ptr, cdis_pairs_f, [2 * n_cdis])
      allocate(cdis_atoms_copy(n_cdis_atoms))
      allocate(cdis_offsets_copy(n_cdis_offsets))
      allocate(cdis_pairs_copy(2 * n_cdis))
      cdis_atoms_copy = cdis_atoms_f       ! Atoms already 1-indexed
      cdis_offsets_copy = cdis_offsets_f   ! Offsets stay 0-based
      cdis_pairs_copy = cdis_pairs_f       ! Group indices stay 0-based
    else
      allocate(cdis_atoms_copy(1))
      allocate(cdis_offsets_copy(2))
      allocate(cdis_pairs_copy(2))
      cdis_atoms_copy = 0
      cdis_offsets_copy = 0
      cdis_pairs_copy = 0
    end if

    ! Convert COM angle arrays (atoms already 1-indexed from Python)
    if (n_cang_atoms > 0 .and. n_cang > 0) then
      call C_F_POINTER(cang_atoms_ptr, cang_atoms_f, [n_cang_atoms])
      call C_F_POINTER(cang_offsets_ptr, cang_offsets_f, [n_cang_offsets])
      call C_F_POINTER(cang_triplets_ptr, cang_triplets_f, [3 * n_cang])
      allocate(cang_atoms_copy(n_cang_atoms))
      allocate(cang_offsets_copy(n_cang_offsets))
      allocate(cang_triplets_copy(3 * n_cang))
      cang_atoms_copy = cang_atoms_f       ! Atoms already 1-indexed
      cang_offsets_copy = cang_offsets_f
      cang_triplets_copy = cang_triplets_f
    else
      allocate(cang_atoms_copy(1))
      allocate(cang_offsets_copy(2))
      allocate(cang_triplets_copy(3))
      cang_atoms_copy = 0
      cang_offsets_copy = 0
      cang_triplets_copy = 0
    end if

    ! Convert COM torsion arrays (atoms already 1-indexed from Python)
    if (n_ctor_atoms > 0 .and. n_ctor > 0) then
      call C_F_POINTER(ctor_atoms_ptr, ctor_atoms_f, [n_ctor_atoms])
      call C_F_POINTER(ctor_offsets_ptr, ctor_offsets_f, [n_ctor_offsets])
      call C_F_POINTER(ctor_quads_ptr, ctor_quads_f, [4 * n_ctor])
      allocate(ctor_atoms_copy(n_ctor_atoms))
      allocate(ctor_offsets_copy(n_ctor_offsets))
      allocate(ctor_quads_copy(4 * n_ctor))
      ctor_atoms_copy = ctor_atoms_f       ! Atoms already 1-indexed
      ctor_offsets_copy = ctor_offsets_f
      ctor_quads_copy = ctor_quads_f
    else
      allocate(ctor_atoms_copy(1))
      allocate(ctor_offsets_copy(2))
      allocate(ctor_quads_copy(4))
      ctor_atoms_copy = 0
      ctor_offsets_copy = 0
      ctor_quads_copy = 0
    end if

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] Trajectory Analysis with COM'
    write(MsgOut,'(A)') ' '

    ! Build a memory-backed source and run the shared COM analysis loop.
    call init_source_memory(source, s_trajes_c%coords, s_trajes_c%pbc_boxes, &
                            s_trajes_c%natom, s_trajes_c%nframe, ana_period)

    call analyze_com(source, mass_f, &
                     dist_list_copy, n_dist, &
                     angl_list_copy, n_angl, &
                     tors_list_copy, n_tors, &
                     cdis_atoms_copy, cdis_offsets_copy, cdis_pairs_copy, &
                     n_cdis, n_cdis_groups, &
                     cang_atoms_copy, cang_offsets_copy, cang_triplets_copy, &
                     n_cang, n_cang_groups, &
                     ctor_atoms_copy, ctor_offsets_copy, ctor_quads_copy, &
                     n_ctor, n_ctor_groups, &
                     distance_f, angle_f, torsion_f, &
                     cdis_f, cang_f, ctor_f, nstru_local)

    call finalize_source(source)

    nstru_out = nstru_local

    ! Clean up dummy allocations and copies
    if (.not. (n_dist > 0 .and. c_associated(dist_ptr))) deallocate(distance_f)
    if (.not. (n_angl > 0 .and. c_associated(angl_ptr))) deallocate(angle_f)
    if (.not. (n_tors > 0 .and. c_associated(tors_ptr))) deallocate(torsion_f)
    if (.not. (n_cdis > 0 .and. c_associated(cdis_result_ptr))) deallocate(cdis_f)
    if (.not. (n_cang > 0 .and. c_associated(cang_result_ptr))) deallocate(cang_f)
    if (.not. (n_ctor > 0 .and. c_associated(ctor_result_ptr))) deallocate(ctor_f)

    deallocate(dist_list_copy)
    deallocate(angl_list_copy)
    deallocate(tors_list_copy)
    deallocate(cdis_atoms_copy)
    deallocate(cdis_offsets_copy)
    deallocate(cdis_pairs_copy)
    deallocate(cang_atoms_copy)
    deallocate(cang_offsets_copy)
    deallocate(cang_triplets_copy)
    deallocate(ctor_atoms_copy)
    deallocate(ctor_offsets_copy)
    deallocate(ctor_quads_copy)

  end subroutine trj_analysis_com_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_com_lazy_c
  !> @brief        COM trajectory analysis with lazy DCD loading
  !! @authors      Claude Code
  !! @note         Mirrors trj_analysis_com_c but reads frames on demand from a
  !!               DCD file via the lazy trajectory source. All measurement and
  !!               COM atom indices are expressed in the selected atom space.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_com_lazy_c(dcd_filename, filename_len, trj_type, &
                                     dcd_natom_expected, source_selection_ptr, &
                                     n_source_selection, &
                                     mass_ptr, n_atoms, ana_period, n_frame, &
                                     dist_list_ptr, n_dist, &
                                     angl_list_ptr, n_angl, &
                                     tors_list_ptr, n_tors, &
                                     cdis_atoms_ptr, n_cdis_atoms, &
                                     cdis_offsets_ptr, n_cdis_offsets, &
                                     cdis_pairs_ptr, n_cdis, &
                                     cang_atoms_ptr, n_cang_atoms, &
                                     cang_offsets_ptr, n_cang_offsets, &
                                     cang_triplets_ptr, n_cang, &
                                     ctor_atoms_ptr, n_ctor_atoms, &
                                     ctor_offsets_ptr, n_ctor_offsets, &
                                     ctor_quads_ptr, n_ctor, &
                                     dist_ptr, dist_size, &
                                     angl_ptr, angl_size, &
                                     tors_ptr, tors_size, &
                                     cdis_result_ptr, cdis_size, &
                                     cang_result_ptr, cang_size, &
                                     ctor_result_ptr, ctor_size, &
                                     nstru_out, dcd_nframe_out, dcd_natom_out, &
                                     status, msg, msglen) &
        bind(C, name="trj_analysis_com_lazy_c")
    implicit none

    ! Arguments - lazy DCD source
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    integer(c_int), value :: dcd_natom_expected
    type(c_ptr), value :: source_selection_ptr
    integer(c_int), value :: n_source_selection

    ! Arguments - mass array (selected atom space)
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    integer(c_int), value :: n_frame

    ! Arguments - atom-based measurements
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors

    ! Arguments - COM distance
    type(c_ptr), value :: cdis_atoms_ptr
    integer(c_int), value :: n_cdis_atoms
    type(c_ptr), value :: cdis_offsets_ptr
    integer(c_int), value :: n_cdis_offsets
    type(c_ptr), value :: cdis_pairs_ptr
    integer(c_int), value :: n_cdis

    ! Arguments - COM angle
    type(c_ptr), value :: cang_atoms_ptr
    integer(c_int), value :: n_cang_atoms
    type(c_ptr), value :: cang_offsets_ptr
    integer(c_int), value :: n_cang_offsets
    type(c_ptr), value :: cang_triplets_ptr
    integer(c_int), value :: n_cang

    ! Arguments - COM torsion
    type(c_ptr), value :: ctor_atoms_ptr
    integer(c_int), value :: n_ctor_atoms
    type(c_ptr), value :: ctor_offsets_ptr
    integer(c_int), value :: n_ctor_offsets
    type(c_ptr), value :: ctor_quads_ptr
    integer(c_int), value :: n_ctor

    ! Pre-allocated output arrays
    type(c_ptr), value :: dist_ptr
    integer(c_int), value :: dist_size
    type(c_ptr), value :: angl_ptr
    integer(c_int), value :: angl_size
    type(c_ptr), value :: tors_ptr
    integer(c_int), value :: tors_size
    type(c_ptr), value :: cdis_result_ptr
    integer(c_int), value :: cdis_size
    type(c_ptr), value :: cang_result_ptr
    integer(c_int), value :: cang_size
    type(c_ptr), value :: ctor_result_ptr
    integer(c_int), value :: ctor_size

    ! Output
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_trj_com_lazy_ctx), target :: c

    ! Hand every input to the guarded body through the context
    call c_filename_to_fortran(dcd_filename, filename_len, c%filename_f)
    c%trj_type             = trj_type
    c%dcd_natom_expected   = dcd_natom_expected
    c%source_selection_ptr = source_selection_ptr
    c%n_source_selection   = n_source_selection
    c%mass_ptr             = mass_ptr
    c%n_atoms              = n_atoms
    c%ana_period           = ana_period
    c%n_frame              = n_frame
    c%dist_list_ptr        = dist_list_ptr
    c%n_dist               = n_dist
    c%angl_list_ptr        = angl_list_ptr
    c%n_angl               = n_angl
    c%tors_list_ptr        = tors_list_ptr
    c%n_tors               = n_tors
    c%cdis_atoms_ptr       = cdis_atoms_ptr
    c%n_cdis_atoms         = n_cdis_atoms
    c%cdis_offsets_ptr     = cdis_offsets_ptr
    c%n_cdis_offsets       = n_cdis_offsets
    c%cdis_pairs_ptr       = cdis_pairs_ptr
    c%n_cdis               = n_cdis
    c%cang_atoms_ptr       = cang_atoms_ptr
    c%n_cang_atoms         = n_cang_atoms
    c%cang_offsets_ptr     = cang_offsets_ptr
    c%n_cang_offsets       = n_cang_offsets
    c%cang_triplets_ptr    = cang_triplets_ptr
    c%n_cang               = n_cang
    c%ctor_atoms_ptr       = ctor_atoms_ptr
    c%n_ctor_atoms         = n_ctor_atoms
    c%ctor_offsets_ptr     = ctor_offsets_ptr
    c%n_ctor_offsets       = n_ctor_offsets
    c%ctor_quads_ptr       = ctor_quads_ptr
    c%n_ctor               = n_ctor
    c%dist_ptr             = dist_ptr
    c%dist_size            = dist_size
    c%angl_ptr             = angl_ptr
    c%angl_size            = angl_size
    c%tors_ptr             = tors_ptr
    c%tors_size            = tors_size
    c%cdis_result_ptr      = cdis_result_ptr
    c%cdis_size            = cdis_size
    c%cang_result_ptr      = cang_result_ptr
    c%cang_size            = cang_size
    c%ctor_result_ptr      = ctor_result_ptr
    c%ctor_size            = ctor_size

    ! Run the body under the library-mode error guard: reading the DCD may
    ! call error_msg -> exit(1) in CLI mode, which would kill the host Python
    ! process (see run_guarded in error_mod).
    call run_guarded(trj_analysis_com_lazy_body, c, c%err, status, msg, msglen)

    nstru_out      = c%nstru_out
    dcd_nframe_out = c%dcd_nframe_out
    dcd_natom_out  = c%dcd_natom_out

    ! Release what the body acquired; this also runs after a longjmp. The
    ! allocatable components of c are freed automatically on return.
    call finalize_source(c%source)

  end subroutine trj_analysis_com_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_com_lazy_body
  !> @brief        Guarded body of trj_analysis_com_lazy_c
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_trj_com_lazy_ctx
  !
  !  Reports failures only through c%err; the wrapper converts them to the C
  !  status/message pair and releases the resources recorded in the context.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_com_lazy_body(ctx) &
        bind(C, name="trj_analysis_com_lazy_body")
    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_trj_com_lazy_ctx), pointer :: c
    integer, pointer :: source_selection_f(:)
    real(wp), pointer :: mass_f(:)
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    integer, pointer :: cdis_atoms_f(:)
    integer, pointer :: cdis_offsets_f(:)
    integer, pointer :: cdis_pairs_f(:)
    integer, pointer :: cang_atoms_f(:)
    integer, pointer :: cang_offsets_f(:)
    integer, pointer :: cang_triplets_f(:)
    integer, pointer :: ctor_atoms_f(:)
    integer, pointer :: ctor_offsets_f(:)
    integer, pointer :: ctor_quads_f(:)
    real(wp), pointer :: distance_f(:,:)
    real(wp), pointer :: angle_f(:,:)
    real(wp), pointer :: torsion_f(:,:)
    real(wp), pointer :: cdis_f(:,:)
    real(wp), pointer :: cang_f(:,:)
    real(wp), pointer :: ctor_f(:,:)
    integer :: nstru_local, init_status
    integer :: n_cdis_groups, n_cang_groups, n_ctor_groups

    call c_f_pointer(ctx, c)

    n_cdis_groups = c%n_cdis_offsets - 1
    n_cang_groups = c%n_cang_offsets - 1
    n_ctor_groups = c%n_ctor_offsets - 1

    ! Validate the source selection
    if (c%n_source_selection /= c%n_atoms .or. &
        .not. c_associated(c%source_selection_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_com_lazy_c: invalid source selection")
      return
    end if

    if (c%n_frame <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_com_lazy_c: n_frame must be positive")
      return
    end if

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    call C_F_POINTER(c%source_selection_ptr, source_selection_f, &
                     [c%n_source_selection])
    if (any(source_selection_f < 1) .or. &
        any(source_selection_f > c%dcd_natom_expected)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_com_lazy_c: source selection out of range")
      return
    end if

    ! Mass array (selected atom space)
    call C_F_POINTER(c%mass_ptr, mass_f, [c%n_atoms])

    ! Convert atom-based measurement lists (1-indexed in selected space)
    if (c%n_dist > 0 .and. c_associated(c%dist_list_ptr)) then
      call C_F_POINTER(c%dist_list_ptr, dist_list_f, [2, c%n_dist])
      allocate(c%dist_list_copy(2, c%n_dist))
      c%dist_list_copy = dist_list_f
    else
      allocate(c%dist_list_copy(2, 1))
      c%dist_list_copy = 0
    end if

    if (c%n_angl > 0 .and. c_associated(c%angl_list_ptr)) then
      call C_F_POINTER(c%angl_list_ptr, angl_list_f, [3, c%n_angl])
      allocate(c%angl_list_copy(3, c%n_angl))
      c%angl_list_copy = angl_list_f
    else
      allocate(c%angl_list_copy(3, 1))
      c%angl_list_copy = 0
    end if

    if (c%n_tors > 0 .and. c_associated(c%tors_list_ptr)) then
      call C_F_POINTER(c%tors_list_ptr, tors_list_f, [4, c%n_tors])
      allocate(c%tors_list_copy(4, c%n_tors))
      c%tors_list_copy = tors_list_f
    else
      allocate(c%tors_list_copy(4, 1))
      c%tors_list_copy = 0
    end if

    ! Convert COM distance arrays (atoms already 1-indexed in selected space)
    if (c%n_cdis_atoms > 0 .and. c%n_cdis > 0) then
      call C_F_POINTER(c%cdis_atoms_ptr, cdis_atoms_f, [c%n_cdis_atoms])
      call C_F_POINTER(c%cdis_offsets_ptr, cdis_offsets_f, [c%n_cdis_offsets])
      call C_F_POINTER(c%cdis_pairs_ptr, cdis_pairs_f, [2 * c%n_cdis])
      allocate(c%cdis_atoms_copy(c%n_cdis_atoms))
      allocate(c%cdis_offsets_copy(c%n_cdis_offsets))
      allocate(c%cdis_pairs_copy(2 * c%n_cdis))
      c%cdis_atoms_copy = cdis_atoms_f
      c%cdis_offsets_copy = cdis_offsets_f
      c%cdis_pairs_copy = cdis_pairs_f
    else
      allocate(c%cdis_atoms_copy(1))
      allocate(c%cdis_offsets_copy(2))
      allocate(c%cdis_pairs_copy(2))
      c%cdis_atoms_copy = 0
      c%cdis_offsets_copy = 0
      c%cdis_pairs_copy = 0
    end if

    ! Convert COM angle arrays
    if (c%n_cang_atoms > 0 .and. c%n_cang > 0) then
      call C_F_POINTER(c%cang_atoms_ptr, cang_atoms_f, [c%n_cang_atoms])
      call C_F_POINTER(c%cang_offsets_ptr, cang_offsets_f, [c%n_cang_offsets])
      call C_F_POINTER(c%cang_triplets_ptr, cang_triplets_f, [3 * c%n_cang])
      allocate(c%cang_atoms_copy(c%n_cang_atoms))
      allocate(c%cang_offsets_copy(c%n_cang_offsets))
      allocate(c%cang_triplets_copy(3 * c%n_cang))
      c%cang_atoms_copy = cang_atoms_f
      c%cang_offsets_copy = cang_offsets_f
      c%cang_triplets_copy = cang_triplets_f
    else
      allocate(c%cang_atoms_copy(1))
      allocate(c%cang_offsets_copy(2))
      allocate(c%cang_triplets_copy(3))
      c%cang_atoms_copy = 0
      c%cang_offsets_copy = 0
      c%cang_triplets_copy = 0
    end if

    ! Convert COM torsion arrays
    if (c%n_ctor_atoms > 0 .and. c%n_ctor > 0) then
      call C_F_POINTER(c%ctor_atoms_ptr, ctor_atoms_f, [c%n_ctor_atoms])
      call C_F_POINTER(c%ctor_offsets_ptr, ctor_offsets_f, [c%n_ctor_offsets])
      call C_F_POINTER(c%ctor_quads_ptr, ctor_quads_f, [4 * c%n_ctor])
      allocate(c%ctor_atoms_copy(c%n_ctor_atoms))
      allocate(c%ctor_offsets_copy(c%n_ctor_offsets))
      allocate(c%ctor_quads_copy(4 * c%n_ctor))
      c%ctor_atoms_copy = ctor_atoms_f
      c%ctor_offsets_copy = ctor_offsets_f
      c%ctor_quads_copy = ctor_quads_f
    else
      allocate(c%ctor_atoms_copy(1))
      allocate(c%ctor_offsets_copy(2))
      allocate(c%ctor_quads_copy(4))
      c%ctor_atoms_copy = 0
      c%ctor_offsets_copy = 0
      c%ctor_quads_copy = 0
    end if

    ! Atom-based measurement indices must stay within the selected atom range
    if ((c%n_dist > 0 .and. (any(c%dist_list_copy < 1) .or. &
                             any(c%dist_list_copy > c%n_atoms))) .or. &
        (c%n_angl > 0 .and. (any(c%angl_list_copy < 1) .or. &
                             any(c%angl_list_copy > c%n_atoms))) .or. &
        (c%n_tors > 0 .and. (any(c%tors_list_copy < 1) .or. &
                             any(c%tors_list_copy > c%n_atoms)))) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_com_lazy_c: measurement index out of range")
      return
    end if

    ! COM group atom indices must also stay within the selected atom range
    if ((c%n_cdis > 0 .and. (any(c%cdis_atoms_copy < 1) .or. &
                             any(c%cdis_atoms_copy > c%n_atoms))) .or. &
        (c%n_cang > 0 .and. (any(c%cang_atoms_copy < 1) .or. &
                             any(c%cang_atoms_copy > c%n_atoms))) .or. &
        (c%n_ctor > 0 .and. (any(c%ctor_atoms_copy < 1) .or. &
                             any(c%ctor_atoms_copy > c%n_atoms)))) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis_com_lazy_c: COM atom index out of range")
      return
    end if

    ! Create views of the pre-allocated result arrays (nullify when unused)
    if (c%n_dist > 0 .and. c_associated(c%dist_ptr) .and. c%dist_size > 0) then
      call C_F_POINTER(c%dist_ptr, distance_f, [c%n_dist, c%n_frame])
    else
      nullify(distance_f)
    end if

    if (c%n_angl > 0 .and. c_associated(c%angl_ptr) .and. c%angl_size > 0) then
      call C_F_POINTER(c%angl_ptr, angle_f, [c%n_angl, c%n_frame])
    else
      nullify(angle_f)
    end if

    if (c%n_tors > 0 .and. c_associated(c%tors_ptr) .and. c%tors_size > 0) then
      call C_F_POINTER(c%tors_ptr, torsion_f, [c%n_tors, c%n_frame])
    else
      nullify(torsion_f)
    end if

    if (c%n_cdis > 0 .and. c_associated(c%cdis_result_ptr) .and. &
        c%cdis_size > 0) then
      call C_F_POINTER(c%cdis_result_ptr, cdis_f, [c%n_cdis, c%n_frame])
    else
      nullify(cdis_f)
    end if

    if (c%n_cang > 0 .and. c_associated(c%cang_result_ptr) .and. &
        c%cang_size > 0) then
      call C_F_POINTER(c%cang_result_ptr, cang_f, [c%n_cang, c%n_frame])
    else
      nullify(cang_f)
    end if

    if (c%n_ctor > 0 .and. c_associated(c%ctor_result_ptr) .and. &
        c%ctor_size > 0) then
      call C_F_POINTER(c%ctor_result_ptr, ctor_f, [c%n_ctor, c%n_frame])
    else
      nullify(ctor_f)
    end if

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for COM Trj Analysis'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(c%source, trim(c%filename_f), c%trj_type, &
                              c%ana_period, source_selection_f, &
                              c%n_source_selection, init_status)
    if (init_status /= 0) then
      call error_set(c%err, init_status, &
                     "trj_analysis_com_lazy_c: unable to initialize DCD source")
      return
    end if

    ! Return DCD info (also on the atom-count error below)
    c%dcd_nframe_out = c%source%dcd_nframe
    c%dcd_natom_out  = c%source%dcd_natom

    if (c%source%dcd_natom /= c%dcd_natom_expected) then
      call error_set(c%err, ERROR_ATOM_COUNT, &
                     "trj_analysis_com_lazy_c: atom count mismatch")
      return
    end if

    ! Run the shared COM analysis loop (lazy loading via source abstraction)
    write(MsgOut,'(A)') '[STEP2] Trajectory Analysis with COM (lazy loading)'
    write(MsgOut,'(A)') ' '

    call analyze_com(c%source, mass_f, &
                     c%dist_list_copy, c%n_dist, &
                     c%angl_list_copy, c%n_angl, &
                     c%tors_list_copy, c%n_tors, &
                     c%cdis_atoms_copy, c%cdis_offsets_copy, c%cdis_pairs_copy, &
                     c%n_cdis, n_cdis_groups, &
                     c%cang_atoms_copy, c%cang_offsets_copy, c%cang_triplets_copy, &
                     c%n_cang, n_cang_groups, &
                     c%ctor_atoms_copy, c%ctor_offsets_copy, c%ctor_quads_copy, &
                     c%n_ctor, n_ctor_groups, &
                     distance_f, angle_f, torsion_f, &
                     cdis_f, cang_f, ctor_f, nstru_local)

    c%nstru_out = nstru_local

  end subroutine trj_analysis_com_lazy_body

end module trj_c_mod
