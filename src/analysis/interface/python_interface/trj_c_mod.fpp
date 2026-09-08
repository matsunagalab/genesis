!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   trj_c_mod
!> @brief   Python (bind(C)) entry points of trj_analysis
!! @authors Norio Takase (NT), Takaharu Mori (TM), Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!  The science lives in ta_analyze_mod (shared with the CLI program). The
!  entry points below only marshal the NumPy buffers into s_trj_measure /
!  result sinks, choose the trajectory source (memory or lazy DCD) and run
!  analyze_trj_unified under the library-mode error guard.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_c_mod
  use, intrinsic :: iso_c_binding
  use s_trajectories_c_mod
  use ta_analyze_mod
  use trj_source_mod
  use result_sink_mod

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
  !  Type          t_trj_ctx
  !> @brief        State shared between the four entry points and their common
  !!               guarded body trj_analysis_body (see run_guarded in error_mod)
  !! @authors      Claude Code
  !
  !  Measurement lists follow the Python conventions: atom indices are
  !  1-based, COM group ids and offsets are 0-based (converted in the body).
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  type :: t_trj_ctx
    ! trajectory source: in-memory (trajes_ptr) or lazy DCD (filename_f)
    type(c_ptr) :: trajes_ptr = c_null_ptr       ! caller's s_trajectories_c
    character(MaxFilename) :: filename_f = ''
    integer     :: trj_type = 0
    integer     :: dcd_natom_expected = 0
    type(c_ptr) :: source_selection_ptr = c_null_ptr
    integer     :: n_source_selection = 0
    integer     :: ana_period = 1
    integer     :: n_frame = 0                   ! result columns (lazy only)
    type(c_ptr) :: mass_ptr = c_null_ptr         ! null unless COM measurements
    integer     :: n_atoms = 0                   ! atoms in the selected space
    ! measurement definitions
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
    ! pre-allocated result arrays, (n_measure, n_frame) each
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
    type(s_error)       :: err
    type(s_trj_source)  :: source
    type(s_result_sink) :: sinks(TRJ_NSTREAM)
    type(s_trj_measure) :: measure
  end type t_trj_ctx

  private :: t_trj_ctx
  private :: set_com_inputs, measure_from_c
  private :: bind_result_rows, copy_int_array, release

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_c
  !> @brief        Trajectory analysis on an in-memory trajectory
  !! @authors      Claude Code
  !! @param[in]    s_trajes_c    : trajectories C structure
  !! @param[in]    ana_period    : analysis period
  !! @param[in]    dist_list_ptr : distance atom pairs (2, n_dist)
  !! @param[in]    angl_list_ptr : angle atom triplets (3, n_angl)
  !! @param[in]    tors_list_ptr : torsion atom quadruplets (4, n_tors)
  !! @param[in]    dist_ptr etc. : pre-allocated result arrays (n, nframe)
  !! @param[out]   nstru_out     : number of frames analyzed
  !! @param[out]   status, msg   : error status and message
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
    type(s_trajectories_c), intent(in), target :: s_trajes_c
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
    type(t_trj_ctx), target :: c

    c%trajes_ptr    = c_loc(s_trajes_c)
    c%n_atoms       = s_trajes_c%natom
    c%ana_period    = ana_period
    c%dist_list_ptr = dist_list_ptr
    c%n_dist        = n_dist
    c%angl_list_ptr = angl_list_ptr
    c%n_angl        = n_angl
    c%tors_list_ptr = tors_list_ptr
    c%n_tors        = n_tors
    c%dist_ptr      = dist_ptr
    c%dist_size     = dist_size
    c%angl_ptr      = angl_ptr
    c%angl_size     = angl_size
    c%tors_ptr      = tors_ptr
    c%tors_size     = tors_size

    call run_guarded(trj_analysis_body, c, c%err, status, msg, msglen)

    nstru_out = c%nstru_out
    call release(c)

  end subroutine trj_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_lazy_c
  !> @brief        Trajectory analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !! @param[in]    dcd_filename         : DCD file path (C string)
  !! @param[in]    filename_len         : length of filename
  !! @param[in]    trj_type             : trajectory type (1=COOR, 2=COOR+BOX)
  !! @param[in]    dcd_natom_expected   : physical atom count in the DCD file
  !! @param[in]    source_selection_ptr : selected atom indices (1-indexed DCD)
  !! @param[in]    n_source_selection   : number of selected atoms
  !! @param[in]    dist_list_ptr etc.   : measurement lists (selected space)
  !! @param[in]    n_atoms              : number of selected (logical) atoms
  !! @param[in]    ana_period           : analysis period
  !! @param[in]    n_frame              : number of result frames (columns)
  !! @param[in]    dist_ptr etc.        : pre-allocated result arrays
  !! @param[out]   nstru_out            : number of frames analyzed
  !! @param[out]   dcd_nframe_out       : total frames in the DCD
  !! @param[out]   dcd_natom_out        : atoms per frame in the DCD
  !! @param[out]   status, msg          : error status and message
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
    type(t_trj_ctx), target :: c

    call c_filename_to_fortran(dcd_filename, filename_len, c%filename_f)
    c%trj_type             = trj_type
    c%dcd_natom_expected   = dcd_natom_expected
    c%source_selection_ptr = source_selection_ptr
    c%n_source_selection   = n_source_selection
    c%n_atoms              = n_atoms
    c%ana_period           = ana_period
    c%n_frame              = n_frame
    c%dist_list_ptr        = dist_list_ptr
    c%n_dist               = n_dist
    c%angl_list_ptr        = angl_list_ptr
    c%n_angl               = n_angl
    c%tors_list_ptr        = tors_list_ptr
    c%n_tors               = n_tors
    c%dist_ptr             = dist_ptr
    c%dist_size            = dist_size
    c%angl_ptr             = angl_ptr
    c%angl_size            = angl_size
    c%tors_ptr             = tors_ptr
    c%tors_size            = tors_size

    call run_guarded(trj_analysis_body, c, c%err, status, msg, msglen)

    nstru_out      = c%nstru_out
    dcd_nframe_out = c%dcd_nframe_out
    dcd_natom_out  = c%dcd_natom_out
    call release(c)

  end subroutine trj_analysis_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_com_c
  !> @brief        Trajectory analysis with COM measurements (in-memory)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr, n_atoms  : masses of the selected atoms
  !! @param[in]    s_trajes_c         : trajectories C structure
  !! @param[in]    ana_period         : analysis period
  !! @param[in]    *_list_ptr         : atom-based measurement lists
  !! @param[in]    c*_atoms_ptr       : flat atom tables of the COM groups
  !! @param[in]    c*_offsets_ptr     : 0-based group offsets (n_groups + 1)
  !! @param[in]    c*_pairs/triplets/quads_ptr : 0-based group ids per measure
  !! @param[in]    *_ptr, *_size      : pre-allocated result arrays
  !! @param[out]   nstru_out          : number of frames analyzed
  !! @param[out]   status, msg        : error status and message
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

    ! Arguments
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in), target :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors
    type(c_ptr), value :: cdis_atoms_ptr
    integer(c_int), value :: n_cdis_atoms
    type(c_ptr), value :: cdis_offsets_ptr
    integer(c_int), value :: n_cdis_offsets
    type(c_ptr), value :: cdis_pairs_ptr
    integer(c_int), value :: n_cdis
    type(c_ptr), value :: cang_atoms_ptr
    integer(c_int), value :: n_cang_atoms
    type(c_ptr), value :: cang_offsets_ptr
    integer(c_int), value :: n_cang_offsets
    type(c_ptr), value :: cang_triplets_ptr
    integer(c_int), value :: n_cang
    type(c_ptr), value :: ctor_atoms_ptr
    integer(c_int), value :: n_ctor_atoms
    type(c_ptr), value :: ctor_offsets_ptr
    integer(c_int), value :: n_ctor_offsets
    type(c_ptr), value :: ctor_quads_ptr
    integer(c_int), value :: n_ctor
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
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_trj_ctx), target :: c

    c%trajes_ptr = c_loc(s_trajes_c)
    c%mass_ptr   = mass_ptr
    c%n_atoms    = n_atoms
    c%ana_period = ana_period
    call set_com_inputs(c, &
         dist_list_ptr, n_dist, angl_list_ptr, n_angl, tors_list_ptr, n_tors, &
         cdis_atoms_ptr, n_cdis_atoms, cdis_offsets_ptr, n_cdis_offsets, &
         cdis_pairs_ptr, n_cdis, &
         cang_atoms_ptr, n_cang_atoms, cang_offsets_ptr, n_cang_offsets, &
         cang_triplets_ptr, n_cang, &
         ctor_atoms_ptr, n_ctor_atoms, ctor_offsets_ptr, n_ctor_offsets, &
         ctor_quads_ptr, n_ctor, &
         dist_ptr, dist_size, angl_ptr, angl_size, tors_ptr, tors_size, &
         cdis_result_ptr, cdis_size, cang_result_ptr, cang_size, &
         ctor_result_ptr, ctor_size)

    call run_guarded(trj_analysis_body, c, c%err, status, msg, msglen)

    nstru_out = c%nstru_out
    call release(c)

  end subroutine trj_analysis_com_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_com_lazy_c
  !> @brief        Trajectory analysis with COM measurements, lazy DCD loading
  !! @authors      Claude Code
  !! @note         Arguments as in trj_analysis_lazy_c and trj_analysis_com_c
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

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    integer(c_int), value :: dcd_natom_expected
    type(c_ptr), value :: source_selection_ptr
    integer(c_int), value :: n_source_selection
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    integer(c_int), value :: n_frame
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors
    type(c_ptr), value :: cdis_atoms_ptr
    integer(c_int), value :: n_cdis_atoms
    type(c_ptr), value :: cdis_offsets_ptr
    integer(c_int), value :: n_cdis_offsets
    type(c_ptr), value :: cdis_pairs_ptr
    integer(c_int), value :: n_cdis
    type(c_ptr), value :: cang_atoms_ptr
    integer(c_int), value :: n_cang_atoms
    type(c_ptr), value :: cang_offsets_ptr
    integer(c_int), value :: n_cang_offsets
    type(c_ptr), value :: cang_triplets_ptr
    integer(c_int), value :: n_cang
    type(c_ptr), value :: ctor_atoms_ptr
    integer(c_int), value :: n_ctor_atoms
    type(c_ptr), value :: ctor_offsets_ptr
    integer(c_int), value :: n_ctor_offsets
    type(c_ptr), value :: ctor_quads_ptr
    integer(c_int), value :: n_ctor
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
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_trj_ctx), target :: c

    call c_filename_to_fortran(dcd_filename, filename_len, c%filename_f)
    c%trj_type             = trj_type
    c%dcd_natom_expected   = dcd_natom_expected
    c%source_selection_ptr = source_selection_ptr
    c%n_source_selection   = n_source_selection
    c%mass_ptr             = mass_ptr
    c%n_atoms              = n_atoms
    c%ana_period           = ana_period
    c%n_frame              = n_frame
    call set_com_inputs(c, &
         dist_list_ptr, n_dist, angl_list_ptr, n_angl, tors_list_ptr, n_tors, &
         cdis_atoms_ptr, n_cdis_atoms, cdis_offsets_ptr, n_cdis_offsets, &
         cdis_pairs_ptr, n_cdis, &
         cang_atoms_ptr, n_cang_atoms, cang_offsets_ptr, n_cang_offsets, &
         cang_triplets_ptr, n_cang, &
         ctor_atoms_ptr, n_ctor_atoms, ctor_offsets_ptr, n_ctor_offsets, &
         ctor_quads_ptr, n_ctor, &
         dist_ptr, dist_size, angl_ptr, angl_size, tors_ptr, tors_size, &
         cdis_result_ptr, cdis_size, cang_result_ptr, cang_size, &
         ctor_result_ptr, ctor_size)

    call run_guarded(trj_analysis_body, c, c%err, status, msg, msglen)

    nstru_out      = c%nstru_out
    dcd_nframe_out = c%dcd_nframe_out
    dcd_natom_out  = c%dcd_natom_out
    call release(c)

  end subroutine trj_analysis_com_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_com_inputs
  !> @brief        Copy the measurement and result arguments of the COM entry
  !!               points into the context
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_com_inputs(c, &
         dist_list_ptr, n_dist, angl_list_ptr, n_angl, tors_list_ptr, n_tors, &
         cdis_atoms_ptr, n_cdis_atoms, cdis_offsets_ptr, n_cdis_offsets, &
         cdis_pairs_ptr, n_cdis, &
         cang_atoms_ptr, n_cang_atoms, cang_offsets_ptr, n_cang_offsets, &
         cang_triplets_ptr, n_cang, &
         ctor_atoms_ptr, n_ctor_atoms, ctor_offsets_ptr, n_ctor_offsets, &
         ctor_quads_ptr, n_ctor, &
         dist_ptr, dist_size, angl_ptr, angl_size, tors_ptr, tors_size, &
         cdis_result_ptr, cdis_size, cang_result_ptr, cang_size, &
         ctor_result_ptr, ctor_size)
    implicit none
    type(t_trj_ctx), intent(inout) :: c
    type(c_ptr), value :: dist_list_ptr, angl_list_ptr, tors_list_ptr
    type(c_ptr), value :: cdis_atoms_ptr, cdis_offsets_ptr, cdis_pairs_ptr
    type(c_ptr), value :: cang_atoms_ptr, cang_offsets_ptr, cang_triplets_ptr
    type(c_ptr), value :: ctor_atoms_ptr, ctor_offsets_ptr, ctor_quads_ptr
    type(c_ptr), value :: dist_ptr, angl_ptr, tors_ptr
    type(c_ptr), value :: cdis_result_ptr, cang_result_ptr, ctor_result_ptr
    integer(c_int), value :: n_dist, n_angl, n_tors
    integer(c_int), value :: n_cdis_atoms, n_cdis_offsets, n_cdis
    integer(c_int), value :: n_cang_atoms, n_cang_offsets, n_cang
    integer(c_int), value :: n_ctor_atoms, n_ctor_offsets, n_ctor
    integer(c_int), value :: dist_size, angl_size, tors_size
    integer(c_int), value :: cdis_size, cang_size, ctor_size

    c%dist_list_ptr     = dist_list_ptr
    c%n_dist            = n_dist
    c%angl_list_ptr     = angl_list_ptr
    c%n_angl            = n_angl
    c%tors_list_ptr     = tors_list_ptr
    c%n_tors            = n_tors
    c%cdis_atoms_ptr    = cdis_atoms_ptr
    c%n_cdis_atoms      = n_cdis_atoms
    c%cdis_offsets_ptr  = cdis_offsets_ptr
    c%n_cdis_offsets    = n_cdis_offsets
    c%cdis_pairs_ptr    = cdis_pairs_ptr
    c%n_cdis            = n_cdis
    c%cang_atoms_ptr    = cang_atoms_ptr
    c%n_cang_atoms      = n_cang_atoms
    c%cang_offsets_ptr  = cang_offsets_ptr
    c%n_cang_offsets    = n_cang_offsets
    c%cang_triplets_ptr = cang_triplets_ptr
    c%n_cang            = n_cang
    c%ctor_atoms_ptr    = ctor_atoms_ptr
    c%n_ctor_atoms      = n_ctor_atoms
    c%ctor_offsets_ptr  = ctor_offsets_ptr
    c%n_ctor_offsets    = n_ctor_offsets
    c%ctor_quads_ptr    = ctor_quads_ptr
    c%n_ctor            = n_ctor
    c%dist_ptr          = dist_ptr
    c%dist_size         = dist_size
    c%angl_ptr          = angl_ptr
    c%angl_size         = angl_size
    c%tors_ptr          = tors_ptr
    c%tors_size         = tors_size
    c%cdis_result_ptr   = cdis_result_ptr
    c%cdis_size         = cdis_size
    c%cang_result_ptr   = cang_result_ptr
    c%cang_size         = cang_size
    c%ctor_result_ptr   = ctor_result_ptr
    c%ctor_size         = ctor_size

  end subroutine set_com_inputs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_body
  !> @brief        Guarded body shared by the four entry points
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_trj_ctx
  !
  !  Reports failures only through c%err; the wrapper converts them to the C
  !  status/message pair and releases the resources recorded in the context.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_body(ctx) bind(C, name="trj_analysis_body")
    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_trj_ctx), pointer :: c
    type(s_trajectories_c), pointer :: trajes
    integer, pointer :: source_selection_f(:)
    real(wp), pointer :: mass_f(:)
    real(wp), target  :: no_mass(0)
    integer :: init_status, nstru, n_frame

    call c_f_pointer(ctx, c)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Measurement definitions (validated against the selected atom space)
    call measure_from_c(c)
    if (error_has(c%err)) return

    ! Masses are only needed for COM measurements
    if (c_associated(c%mass_ptr)) then
      call c_f_pointer(c%mass_ptr, mass_f, [c%n_atoms])
    else
      mass_f => no_mass
    end if

    ! Trajectory source
    if (c_associated(c%trajes_ptr)) then

      write(MsgOut,'(A)') '[STEP1] Trajectory Analysis'
      write(MsgOut,'(A)') ' '

      call c_f_pointer(c%trajes_ptr, trajes)
      call init_source_memory(c%source, trajes%coords, trajes%pbc_boxes, &
                              trajes%natom, trajes%nframe, c%ana_period)
      n_frame = trajes%nframe / c%ana_period

    else

      if (c%n_frame <= 0) then
        call error_set(c%err, ERROR_INVALID_PARAM, &
                       "trj_analysis: n_frame must be positive")
        return
      end if
      if (c%n_source_selection /= c%n_atoms .or. &
          .not. c_associated(c%source_selection_ptr)) then
        call error_set(c%err, ERROR_INVALID_PARAM, &
                       "trj_analysis: invalid source selection")
        return
      end if
      call c_f_pointer(c%source_selection_ptr, source_selection_f, &
                       [c%n_source_selection])
      if (any(source_selection_f < 1) .or. &
          any(source_selection_f > c%dcd_natom_expected)) then
        call error_set(c%err, ERROR_INVALID_PARAM, &
                       "trj_analysis: source selection out of range")
        return
      end if

      write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for Trj Analysis'
      write(MsgOut,'(A)') ' '

      call init_source_lazy_dcd(c%source, trim(c%filename_f), c%trj_type, &
                                c%ana_period, source_selection_f, &
                                c%n_source_selection, init_status)
      if (init_status /= 0) then
        call error_set(c%err, init_status, &
                       "trj_analysis: unable to initialize DCD source")
        return
      end if

      ! Return DCD info (also on the atom-count error below)
      c%dcd_nframe_out = c%source%dcd_nframe
      c%dcd_natom_out  = c%source%dcd_natom
      if (c%source%dcd_natom /= c%dcd_natom_expected) then
        call error_set(c%err, ERROR_ATOM_COUNT, &
                       "trj_analysis: atom count mismatch")
        return
      end if
      n_frame = c%n_frame

    end if

    ! Result arrays -> sinks (streams without a buffer are left inactive)
    call bind_result_rows(c%sinks(TRJ_DIS),  c%dist_ptr, c%dist_size, &
                          c%n_dist, n_frame)
    call bind_result_rows(c%sinks(TRJ_ANG),  c%angl_ptr, c%angl_size, &
                          c%n_angl, n_frame)
    call bind_result_rows(c%sinks(TRJ_TOR),  c%tors_ptr, c%tors_size, &
                          c%n_tors, n_frame)
    call bind_result_rows(c%sinks(TRJ_CDIS), c%cdis_result_ptr, c%cdis_size, &
                          c%n_cdis, n_frame)
    call bind_result_rows(c%sinks(TRJ_CANG), c%cang_result_ptr, c%cang_size, &
                          c%n_cang, n_frame)
    call bind_result_rows(c%sinks(TRJ_CTOR), c%ctor_result_ptr, c%ctor_size, &
                          c%n_ctor, n_frame)

    write(MsgOut,'(A)') '[STEP2] Trajectory Analysis'
    write(MsgOut,'(A)') ' '

    call analyze_trj_unified(c%source, mass_f, c%measure, c%sinks, nstru)

    c%nstru_out = nstru

  end subroutine trj_analysis_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    measure_from_c
  !> @brief        Build s_trj_measure from the NumPy buffers of the context
  !! @authors      Claude Code
  !
  !  Python passes 1-based atom indices (selected atom space), 0-based group
  !  offsets and 0-based group ids; the latter two become 1-based here. All
  !  indices are validated so that a bad request raises instead of crashing.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine measure_from_c(c)
    implicit none
    type(t_trj_ctx), intent(inout) :: c

    integer, pointer :: list_f(:,:)

    ! atom-based measurements
    if (c%n_dist > 0 .and. c_associated(c%dist_list_ptr)) then
      call c_f_pointer(c%dist_list_ptr, list_f, [2, c%n_dist])
      c%measure%dist_list = list_f
      allocate(c%measure%dist_num(c%n_dist), c%measure%dist_weight(c%n_dist, 1))
      c%measure%dist_num    = 2
      c%measure%dist_weight = 1.0_wp
    end if
    if (c%n_angl > 0 .and. c_associated(c%angl_list_ptr)) then
      call c_f_pointer(c%angl_list_ptr, list_f, [3, c%n_angl])
      c%measure%angl_list = list_f
    end if
    if (c%n_tors > 0 .and. c_associated(c%tors_list_ptr)) then
      call c_f_pointer(c%tors_list_ptr, list_f, [4, c%n_tors])
      c%measure%tors_list = list_f
    end if

    ! COM measurements
    if (c%n_cdis > 0 .and. c%n_cdis_atoms > 0) then
      call copy_int_array(c%cdis_atoms_ptr,   c%n_cdis_atoms,   c%measure%cdis_atoms,   0)
      call copy_int_array(c%cdis_offsets_ptr, c%n_cdis_offsets, c%measure%cdis_offsets, 1)
      call copy_int_array(c%cdis_pairs_ptr,   2 * c%n_cdis,     c%measure%cdis_pairs,   1)
    end if
    if (c%n_cang > 0 .and. c%n_cang_atoms > 0) then
      call copy_int_array(c%cang_atoms_ptr,    c%n_cang_atoms,   c%measure%cang_atoms,    0)
      call copy_int_array(c%cang_offsets_ptr,  c%n_cang_offsets, c%measure%cang_offsets,  1)
      call copy_int_array(c%cang_triplets_ptr, 3 * c%n_cang,     c%measure%cang_triplets, 1)
    end if
    if (c%n_ctor > 0 .and. c%n_ctor_atoms > 0) then
      call copy_int_array(c%ctor_atoms_ptr,   c%n_ctor_atoms,   c%measure%ctor_atoms,   0)
      call copy_int_array(c%ctor_offsets_ptr, c%n_ctor_offsets, c%measure%ctor_offsets, 1)
      call copy_int_array(c%ctor_quads_ptr,   4 * c%n_ctor,     c%measure%ctor_quads,   1)
    end if

    ! Every atom index must stay within the selected atom range
    if (out_of_range(c%measure%dist_list, c%n_atoms) .or. &
        out_of_range(c%measure%angl_list, c%n_atoms) .or. &
        out_of_range(c%measure%tors_list, c%n_atoms)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis: measurement index out of range")
      return
    end if
    if (out_of_range1(c%measure%cdis_atoms, c%n_atoms) .or. &
        out_of_range1(c%measure%cang_atoms, c%n_atoms) .or. &
        out_of_range1(c%measure%ctor_atoms, c%n_atoms)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis: COM atom index out of range")
      return
    end if
    if (bad_groups(c%measure%cdis_pairs,    c%measure%cdis_offsets, c%n_cdis_atoms) .or. &
        bad_groups(c%measure%cang_triplets, c%measure%cang_offsets, c%n_cang_atoms) .or. &
        bad_groups(c%measure%ctor_quads,    c%measure%ctor_offsets, c%n_ctor_atoms)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "trj_analysis: COM group definition out of range")
      return
    end if

  contains

    logical function out_of_range(list, n)
      integer, allocatable, intent(in) :: list(:,:)
      integer, intent(in) :: n
      out_of_range = .false.
      if (allocated(list)) out_of_range = any(list < 1) .or. any(list > n)
    end function out_of_range

    logical function out_of_range1(list, n)
      integer, allocatable, intent(in) :: list(:)
      integer, intent(in) :: n
      out_of_range1 = .false.
      if (allocated(list)) out_of_range1 = any(list < 1) .or. any(list > n)
    end function out_of_range1

    logical function bad_groups(ids, offsets, natom)
      integer, allocatable, intent(in) :: ids(:), offsets(:)
      integer, intent(in) :: natom
      integer :: ngroup, i
      bad_groups = .false.
      if (.not. allocated(ids)) return
      if (.not. allocated(offsets)) then
        bad_groups = .true.
        return
      end if
      ngroup = size(offsets) - 1
      if (ngroup < 1) then
        bad_groups = .true.
        return
      end if
      bad_groups = any(ids < 1) .or. any(ids > ngroup) .or. &
                   any(offsets < 1) .or. any(offsets > natom + 1)
      do i = 1, ngroup
        if (offsets(i+1) < offsets(i)) bad_groups = .true.
      end do
    end function bad_groups

  end subroutine measure_from_c

  !> Copy an integer NumPy buffer into an allocatable array, adding `shift`
  !! (1 to turn 0-based ids/offsets into Fortran indices).
  subroutine copy_int_array(ptr, n, dest, shift)
    implicit none
    type(c_ptr), intent(in) :: ptr
    integer, intent(in) :: n, shift
    integer, allocatable, intent(out) :: dest(:)
    integer, pointer :: src(:)

    if (n <= 0 .or. .not. c_associated(ptr)) then
      allocate(dest(0))
      return
    end if
    call c_f_pointer(ptr, src, [n])
    allocate(dest(n))
    dest = src + shift
  end subroutine copy_int_array

  !> Point a row sink at a pre-allocated (n, n_frame) NumPy result array.
  subroutine bind_result_rows(sink, ptr, size_total, n, n_frame)
    implicit none
    type(s_result_sink), intent(inout) :: sink
    type(c_ptr), intent(in) :: ptr
    integer, intent(in) :: size_total, n, n_frame
    real(wp), pointer :: rows(:,:)

    if (n <= 0 .or. size_total <= 0 .or. .not. c_associated(ptr)) return
    call c_f_pointer(ptr, rows, [n, n_frame])
    call init_sink_array_rows(sink, rows, n, n_frame)
  end subroutine bind_result_rows

  !> Release what the body acquired (also after a longjmp).
  subroutine release(c)
    implicit none
    type(t_trj_ctx), intent(inout) :: c
    integer :: k

    do k = 1, TRJ_NSTREAM
      call finalize_sink(c%sinks(k))
    end do
    call finalize_source(c%source)
  end subroutine release

end module trj_c_mod
