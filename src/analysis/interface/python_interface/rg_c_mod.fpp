!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  rg_main
!! @brief   RG analysis
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use rg_analyze_mod           ! Use unified analysis from CLI module
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

  public :: rg_analysis_c
  public :: rg_analysis_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Type          t_rg_lazy_ctx
  !> @brief        State shared between rg_analysis_lazy_c and its guarded
  !!               body rg_analysis_lazy_body (see run_guarded in error_mod)
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  type :: t_rg_lazy_ctx
    ! inputs (copied from the bind(C) arguments)
    character(MaxFilename) :: filename_f = ''
    integer     :: trj_type = 0
    integer     :: dcd_natom_expected = 0
    type(c_ptr) :: source_selection_ptr = c_null_ptr
    integer     :: n_source_selection = 0
    type(c_ptr) :: mass_ptr = c_null_ptr
    integer     :: n_atoms = 0
    integer     :: ana_period = 1
    type(c_ptr) :: analysis_idx = c_null_ptr
    integer     :: n_analysis = 0
    integer     :: mass_weighted = 0
    type(c_ptr) :: result_ptr = c_null_ptr
    integer     :: result_size = 0
    ! outputs (copied back to the bind(C) arguments by the wrapper)
    integer     :: nstru_out = 0
    integer     :: dcd_nframe_out = 0
    integer     :: dcd_natom_out = 0
    ! error state and resources released by the wrapper
    type(s_error)        :: err
    type(s_trj_source)   :: source
    type(s_result_sink)  :: sink
    integer, allocatable :: idx_copy(:)
  end type t_rg_lazy_ctx

  private :: t_rg_lazy_ctx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rg_analysis_c
  !> @brief        RG analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr        : pointer to mass array (from Python NumPy)
  !! @param[in]    n_atoms         : number of atoms
  !! @param[in]    s_trajes_c      : trajectories C structure
  !! @param[in]    ana_period      : analysis period
  !! @param[in]    analysis_idx    : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis      : number of analysis atoms
  !! @param[in]    mass_weighted   : use mass weighting (0 or 1)
  !! @param[in]    result_ptr      : pointer to pre-allocated result array
  !! @param[in]    result_size     : size of result array
  !! @param[out]   nstru_out       : actual number of structures analyzed
  !! @param[out]   status          : error status
  !! @param[out]   msg             : error message
  !! @param[in]    msglen          : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rg_analysis_c(mass_ptr, n_atoms, s_trajes_c, ana_period, &
                           analysis_idx, n_analysis, mass_weighted, &
                           result_ptr, result_size, nstru_out, &
                           status, msg, msglen) &
        bind(C, name="rg_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: analysis_idx
    integer(c_int), value :: n_analysis
    integer(c_int), value :: mass_weighted
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    type(s_trj_source) :: source
    type(s_result_sink) :: sink
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: idx_f(:)
    integer, allocatable :: idx_copy(:)
    logical :: use_mass
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(analysis_idx, idx_f, [n_analysis])
    allocate(idx_copy(n_analysis))
    idx_copy(:) = idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize source (memory mode) and sink (array mode)
    call init_source_memory(source, s_trajes_c%coords, s_trajes_c%pbc_boxes, &
                            s_trajes_c%natom, s_trajes_c%nframe, ana_period)
    call init_sink_array(sink, result_f, result_size)

    ! Run unified RG analysis
    write(MsgOut,'(A)') '[STEP1] RG Analysis (unified)'
    write(MsgOut,'(A)') ' '

    call analyze_rg_unified(source, sink, mass_f, n_atoms, &
                            idx_copy, n_analysis, use_mass, nstru)

    nstru_out = nstru

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)
    deallocate(idx_copy)

  end subroutine rg_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rg_analysis_lazy_c
  !> @brief        RG analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !! @param[in]    dcd_filename    : DCD file path (C string)
  !! @param[in]    filename_len    : length of filename
  !! @param[in]    trj_type        : trajectory type (1=COOR, 2=COOR+BOX)
  !! @param[in]    mass_ptr        : pointer to mass array
  !! @param[in]    n_atoms         : number of atoms
  !! @param[in]    ana_period      : analysis period
  !! @param[in]    analysis_idx    : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis      : number of analysis atoms
  !! @param[in]    mass_weighted   : use mass weighting (0 or 1)
  !! @param[in]    result_ptr      : pointer to pre-allocated result array
  !! @param[in]    result_size     : size of result array
  !! @param[out]   nstru_out       : actual number of structures analyzed
  !! @param[out]   dcd_nframe_out  : total frames in DCD
  !! @param[out]   dcd_natom_out   : atoms per frame in DCD
  !! @param[out]   status          : error status
  !! @param[out]   msg             : error message
  !! @param[in]    msglen          : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rg_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                dcd_natom_expected, source_selection_ptr, &
                                n_source_selection, &
                                mass_ptr, n_atoms, ana_period, &
                                analysis_idx, n_analysis, mass_weighted, &
                                result_ptr, result_size, &
                                nstru_out, dcd_nframe_out, dcd_natom_out, &
                                status, msg, msglen) &
        bind(C, name="rg_analysis_lazy_c")
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
    type(c_ptr), value :: analysis_idx
    integer(c_int), value :: n_analysis
    integer(c_int), value :: mass_weighted
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_rg_lazy_ctx), target :: c

    ! Hand every input to the guarded body through the context
    call c_filename_to_fortran(dcd_filename, filename_len, c%filename_f)
    c%trj_type             = trj_type
    c%dcd_natom_expected   = dcd_natom_expected
    c%source_selection_ptr = source_selection_ptr
    c%n_source_selection   = n_source_selection
    c%mass_ptr             = mass_ptr
    c%n_atoms              = n_atoms
    c%ana_period           = ana_period
    c%analysis_idx         = analysis_idx
    c%n_analysis           = n_analysis
    c%mass_weighted        = mass_weighted
    c%result_ptr           = result_ptr
    c%result_size          = result_size

    ! Run the body under the library-mode error guard: reading the DCD may
    ! call error_msg -> exit(1) in CLI mode, which would kill the host Python
    ! process (see run_guarded in error_mod).
    call run_guarded(rg_analysis_lazy_body, c, c%err, status, msg, msglen)

    nstru_out      = c%nstru_out
    dcd_nframe_out = c%dcd_nframe_out
    dcd_natom_out  = c%dcd_natom_out

    ! Release what the body acquired; this also runs after a longjmp. The
    ! allocatable components of c are freed automatically on return.
    call finalize_sink(c%sink)
    call finalize_source(c%source)

  end subroutine rg_analysis_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rg_analysis_lazy_body
  !> @brief        Guarded body of rg_analysis_lazy_c
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_rg_lazy_ctx
  !
  !  Reports failures only through c%err; the wrapper converts them to the C
  !  status/message pair and releases the resources recorded in the context.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rg_analysis_lazy_body(ctx) bind(C, name="rg_analysis_lazy_body")
    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_rg_lazy_ctx), pointer :: c
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: idx_f(:)
    integer, pointer :: source_selection_f(:)
    logical :: use_mass
    integer :: nstru, init_status

    call c_f_pointer(ctx, c)

    ! Validate inputs
    if (c%n_analysis <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: n_analysis must be positive")
      return
    end if

    if (c%n_source_selection /= c%n_atoms .or. &
        .not. c_associated(c%source_selection_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: invalid source selection")
      return
    end if

    if (.not. c_associated(c%mass_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: mass_ptr is null")
      return
    end if

    if (.not. c_associated(c%result_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: result_ptr is null")
      return
    end if

    if (c%result_size <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: result_size must be positive")
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(c%mass_ptr, mass_f, [c%n_atoms])
    call C_F_POINTER(c%result_ptr, result_f, [c%result_size])
    call C_F_POINTER(c%source_selection_ptr, source_selection_f, &
                     [c%n_source_selection])
    if (any(source_selection_f < 1) .or. &
        any(source_selection_f > c%dcd_natom_expected)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: source selection out of range")
      return
    end if

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(c%analysis_idx, idx_f, [c%n_analysis])
    allocate(c%idx_copy(c%n_analysis))
    c%idx_copy(:) = idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (c%mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for RG'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(c%source, trim(c%filename_f), c%trj_type, &
                              c%ana_period, source_selection_f, &
                              c%n_source_selection, init_status)
    if (init_status /= 0) then
      call error_set(c%err, init_status, &
                     "rg_analysis_lazy_c: unable to initialize DCD source")
      return
    end if

    ! Return DCD info (also on the atom-count error below)
    c%dcd_nframe_out = c%source%dcd_nframe
    c%dcd_natom_out  = c%source%dcd_natom

    ! Check atom count
    if (c%source%dcd_natom /= c%dcd_natom_expected) then
      call error_set(c%err, ERROR_ATOM_COUNT, &
                     "rg_analysis_lazy_c: atom count mismatch")
      return
    end if

    ! Initialize sink (array mode)
    call init_sink_array(c%sink, result_f, c%result_size)

    ! Run unified RG analysis (lazy loading via source abstraction)
    write(MsgOut,'(A)') '[STEP2] RG Analysis (lazy loading, unified)'
    write(MsgOut,'(A)') ' '

    call analyze_rg_unified(c%source, c%sink, mass_f, c%n_atoms, &
                            c%idx_copy, c%n_analysis, use_mass, nstru)

    c%nstru_out = nstru

  end subroutine rg_analysis_lazy_body

end module rg_c_mod
