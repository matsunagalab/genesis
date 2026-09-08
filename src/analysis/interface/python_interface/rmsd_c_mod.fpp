!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  ra_main
!! @brief   RMSD analysis
!! @authors Takaharu Mori (TM), Yuji Sugita (YS), Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rmsd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use ra_analyze_mod           ! Use unified analysis from CLI module
  use trj_source_mod
  use result_sink_mod

  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use molecules_str_mod
  use string_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: rmsd_analysis_c
  public :: rmsd_analysis_fitting_c
  public :: rmsd_analysis_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Type          t_rmsd_lazy_ctx
  !> @brief        State shared between rmsd_analysis_lazy_c and its guarded
  !!               body rmsd_analysis_lazy_body
  !! @authors      Claude Code
  !
  !  The body runs under the library-mode error guard (fi_error_guard_run_ctx
  !  in error_mod) and therefore has to be a module procedure: taking
  !  c_funloc() of an internal procedure makes gfortran emit a trampoline on
  !  the stack, which requires an executable stack and is rejected by dlopen()
  !  on glibc >= 2.41. Everything the body needs is passed through this
  !  context instead of host association, and every resource the body may
  !  acquire lives here so that the wrapper can release it after the guard
  !  returns, including after a longjmp.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  type :: t_rmsd_lazy_ctx
    ! inputs (copied from the bind(C) arguments)
    character(MaxFilename) :: filename_f = ''
    integer     :: trj_type = 0
    integer     :: dcd_natom_expected = 0
    type(c_ptr) :: source_selection_ptr = c_null_ptr
    integer     :: n_source_selection = 0
    type(c_ptr) :: mass_ptr = c_null_ptr
    type(c_ptr) :: ref_coord_ptr = c_null_ptr
    integer     :: n_atoms = 0
    integer     :: ana_period = 1
    type(c_ptr) :: fitting_idx_ptr = c_null_ptr
    integer     :: n_fitting = 0
    type(c_ptr) :: analysis_idx_ptr = c_null_ptr
    integer     :: n_analysis = 0
    integer     :: fitting_method = 0
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
    integer, allocatable :: fitting_idx_copy(:)
    integer, allocatable :: analysis_idx_copy(:)
  end type t_rmsd_lazy_ctx

  private :: t_rmsd_lazy_ctx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_c
  !> @brief        RMSD analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr       : pointer to mass array (from Python NumPy)
  !! @param[in]    ref_coord_ptr  : pointer to reference coordinates (3, n_atoms)
  !! @param[in]    n_atoms        : number of atoms
  !! @param[in]    s_trajes_c     : trajectories C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    analysis_idx   : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    mass_weighted  : use mass weighting (0 or 1)
  !! @param[in]    result_ptr     : pointer to pre-allocated result array
  !! @param[in]    result_size    : size of pre-allocated result array
  !! @param[out]   nstru_out      : number of frames analyzed
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_c(mass_ptr, ref_coord_ptr, n_atoms, &
                             s_trajes_c, ana_period, &
                             analysis_idx, n_analysis, mass_weighted, &
                             result_ptr, result_size, nstru_out, &
                             status, msg, msglen) &
        bind(C, name="rmsd_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    type(c_ptr), value :: ref_coord_ptr
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
    real(wp), pointer :: ref_coord_f(:,:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: idx_f(:)
    integer, allocatable :: analysis_idx_copy(:)
    integer, allocatable :: dummy_fitting_idx(:)
    logical :: use_mass
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(ref_coord_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: ref_coord_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(ref_coord_ptr, ref_coord_f, [3, n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(analysis_idx, idx_f, [n_analysis])
    allocate(analysis_idx_copy(n_analysis))
    analysis_idx_copy(:) = idx_f(:)

    ! Create dummy fitting indices (no fitting for this function)
    allocate(dummy_fitting_idx(1))
    dummy_fitting_idx(1) = 1

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

    ! Run unified RMSD analysis (no fitting: n_fitting=0, method=FittingMethodNO)
    write(MsgOut,'(A)') '[STEP1] RMSD Analysis (no fitting, unified)'
    write(MsgOut,'(A)') ' '

    call analyze_rmsd_unified(source, sink, ref_coord_f, mass_f, n_atoms, &
                              dummy_fitting_idx, 0, &
                              analysis_idx_copy, n_analysis, &
                              FittingMethodNO, use_mass, nstru)

    nstru_out = nstru

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)
    deallocate(analysis_idx_copy)
    deallocate(dummy_fitting_idx)

  end subroutine rmsd_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_fitting_c
  !> @brief        RMSD analysis with fitting (zerocopy, pre-allocated)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr       : pointer to mass array (from Python NumPy)
  !! @param[in]    ref_coord_ptr  : pointer to reference coordinates (3, n_atoms)
  !! @param[in]    n_atoms        : number of atoms
  !! @param[in]    s_trajes_c     : trajectories C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    fitting_idx_ptr: pointer to fitting atom indices (1-indexed)
  !! @param[in]    n_fitting      : number of fitting atoms
  !! @param[in]    analysis_idx_ptr: pointer to analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    fitting_method : fitting method (1-6)
  !! @param[in]    mass_weighted  : use mass weighting (0 or 1)
  !! @param[in]    result_ptr     : pointer to pre-allocated result array
  !! @param[in]    result_size    : size of pre-allocated result array
  !! @param[out]   nstru_out      : number of frames analyzed
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_fitting_c(mass_ptr, ref_coord_ptr, n_atoms, &
                                     s_trajes_c, ana_period, &
                                     fitting_idx_ptr, n_fitting, &
                                     analysis_idx_ptr, n_analysis, &
                                     fitting_method, mass_weighted, &
                                     result_ptr, result_size, nstru_out, &
                                     status, msg, msglen) &
        bind(C, name="rmsd_analysis_fitting_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    type(c_ptr), value :: ref_coord_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: fitting_idx_ptr
    integer(c_int), value :: n_fitting
    type(c_ptr), value :: analysis_idx_ptr
    integer(c_int), value :: n_analysis
    integer(c_int), value :: fitting_method
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
    real(wp), pointer :: ref_coord_f(:,:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: fitting_idx_f(:)
    integer, pointer :: analysis_idx_f(:)
    integer, allocatable :: fitting_idx_copy(:)
    integer, allocatable :: analysis_idx_copy(:)
    logical :: use_mass
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_fitting <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: n_fitting must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(ref_coord_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: ref_coord_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(fitting_idx_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: fitting_idx_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(analysis_idx_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: analysis_idx_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(ref_coord_ptr, ref_coord_f, [3, n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert indices from C pointer to Fortran array
    call C_F_POINTER(fitting_idx_ptr, fitting_idx_f, [n_fitting])
    call C_F_POINTER(analysis_idx_ptr, analysis_idx_f, [n_analysis])

    allocate(fitting_idx_copy(n_fitting))
    allocate(analysis_idx_copy(n_analysis))
    fitting_idx_copy(:) = fitting_idx_f(:)
    analysis_idx_copy(:) = analysis_idx_f(:)

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

    ! Run unified RMSD analysis with fitting
    write(MsgOut,'(A)') '[STEP1] RMSD Analysis (with fitting, unified)'
    write(MsgOut,'(A)') ' '

    call analyze_rmsd_unified(source, sink, ref_coord_f, mass_f, n_atoms, &
                              fitting_idx_copy, n_fitting, &
                              analysis_idx_copy, n_analysis, &
                              fitting_method, use_mass, nstru)

    nstru_out = nstru

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)
    deallocate(fitting_idx_copy)
    deallocate(analysis_idx_copy)

  end subroutine rmsd_analysis_fitting_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_lazy_c
  !> @brief        RMSD analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !
  !  Thin bind(C) entry point: packs the arguments into a t_rmsd_lazy_ctx,
  !  runs rmsd_analysis_lazy_body under the library-mode error guard, reports
  !  the outcome to C and releases whatever the body acquired.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                  dcd_natom_expected, source_selection_ptr, &
                                  n_source_selection, &
                                  mass_ptr, ref_coord_ptr, n_atoms, &
                                  ana_period, &
                                  fitting_idx_ptr, n_fitting, &
                                  analysis_idx_ptr, n_analysis, &
                                  fitting_method, mass_weighted, &
                                  result_ptr, result_size, &
                                  nstru_out, dcd_nframe_out, dcd_natom_out, &
                                  status, msg, msglen) &
        bind(C, name="rmsd_analysis_lazy_c")
    implicit none

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    integer(c_int), value :: dcd_natom_expected
    type(c_ptr), value :: source_selection_ptr
    integer(c_int), value :: n_source_selection
    type(c_ptr), value :: mass_ptr
    type(c_ptr), value :: ref_coord_ptr
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    type(c_ptr), value :: fitting_idx_ptr
    integer(c_int), value :: n_fitting
    type(c_ptr), value :: analysis_idx_ptr
    integer(c_int), value :: n_analysis
    integer(c_int), value :: fitting_method
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
    type(t_rmsd_lazy_ctx), target :: c

    ! Hand every input to the guarded body through the context
    call c_filename_to_fortran(dcd_filename, filename_len, c%filename_f)
    c%trj_type             = trj_type
    c%dcd_natom_expected   = dcd_natom_expected
    c%source_selection_ptr = source_selection_ptr
    c%n_source_selection   = n_source_selection
    c%mass_ptr             = mass_ptr
    c%ref_coord_ptr        = ref_coord_ptr
    c%n_atoms              = n_atoms
    c%ana_period           = ana_period
    c%fitting_idx_ptr      = fitting_idx_ptr
    c%n_fitting            = n_fitting
    c%analysis_idx_ptr     = analysis_idx_ptr
    c%n_analysis           = n_analysis
    c%fitting_method       = fitting_method
    c%mass_weighted        = mass_weighted
    c%result_ptr           = result_ptr
    c%result_size          = result_size

    ! Run the body under the library-mode error guard: reading the DCD may
    ! call error_msg -> exit(1) in CLI mode, which would kill the host Python
    ! process. The body is a module procedure, so the callback needs neither
    ! a trampoline nor an executable stack (see run_guarded in error_mod).
    call run_guarded(rmsd_analysis_lazy_body, c, c%err, status, msg, msglen)

    nstru_out      = c%nstru_out
    dcd_nframe_out = c%dcd_nframe_out
    dcd_natom_out  = c%dcd_natom_out

    ! Release what the body acquired; this also runs after a longjmp, when
    ! the body never reached its own end. The allocatable components of c are
    ! freed automatically on return.
    call finalize_sink(c%sink)
    call finalize_source(c%source)

  end subroutine rmsd_analysis_lazy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_lazy_body
  !> @brief        Guarded body of rmsd_analysis_lazy_c
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_rmsd_lazy_ctx
  !
  !  Reports failures only through c%err; the wrapper converts them to the C
  !  status/message pair and releases the resources recorded in the context.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_lazy_body(ctx) bind(C, name="rmsd_analysis_lazy_body")
    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_rmsd_lazy_ctx), pointer :: c
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: ref_coord_f(:,:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: fitting_idx_f(:)
    integer, pointer :: analysis_idx_f(:)
    integer, pointer :: source_selection_f(:)
    logical :: use_mass
    integer :: nstru, n_fitting_use, fitting_method_use, init_status

    call c_f_pointer(ctx, c)

    ! Validate inputs
    if (c%n_analysis <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: n_analysis must be positive")
      return
    end if

    if (c%n_source_selection /= c%n_atoms .or. &
        .not. c_associated(c%source_selection_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: invalid source selection")
      return
    end if

    if (.not. c_associated(c%mass_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: mass_ptr is null")
      return
    end if

    if (.not. c_associated(c%ref_coord_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: ref_coord_ptr is null")
      return
    end if

    if (.not. c_associated(c%result_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: result_ptr is null")
      return
    end if

    if (c%result_size <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: result_size must be positive")
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(c%mass_ptr, mass_f, [c%n_atoms])
    call C_F_POINTER(c%ref_coord_ptr, ref_coord_f, [3, c%n_atoms])
    call C_F_POINTER(c%result_ptr, result_f, [c%result_size])
    call C_F_POINTER(c%source_selection_ptr, source_selection_f, &
                     [c%n_source_selection])
    if (any(source_selection_f < 1) .or. &
        any(source_selection_f > c%dcd_natom_expected)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: source selection out of range")
      return
    end if

    ! Convert analysis indices
    call C_F_POINTER(c%analysis_idx_ptr, analysis_idx_f, [c%n_analysis])
    allocate(c%analysis_idx_copy(c%n_analysis))
    c%analysis_idx_copy(:) = analysis_idx_f(:)

    ! Check for fitting and prepare fitting indices
    if (c%n_fitting > 0 .and. c%fitting_method > 0 .and. &
        c_associated(c%fitting_idx_ptr)) then
      call C_F_POINTER(c%fitting_idx_ptr, fitting_idx_f, [c%n_fitting])
      allocate(c%fitting_idx_copy(c%n_fitting))
      c%fitting_idx_copy(:) = fitting_idx_f(:)
      n_fitting_use = c%n_fitting
      fitting_method_use = c%fitting_method
    else
      ! No fitting - create dummy array
      allocate(c%fitting_idx_copy(1))
      c%fitting_idx_copy(1) = 1
      n_fitting_use = 0
      fitting_method_use = FittingMethodNO
    end if

    ! Convert mass_weighted to logical
    use_mass = (c%mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(c%source, trim(c%filename_f), c%trj_type, &
                              c%ana_period, source_selection_f, &
                              c%n_source_selection, init_status)
    if (init_status /= 0) then
      call error_set(c%err, init_status, &
                     "rmsd_analysis_lazy_c: unable to initialize DCD source")
      return
    end if

    ! Return DCD info (also on the atom-count error below)
    c%dcd_nframe_out = c%source%dcd_nframe
    c%dcd_natom_out  = c%source%dcd_natom

    ! Check atom count
    if (c%source%dcd_natom /= c%dcd_natom_expected) then
      call error_set(c%err, ERROR_ATOM_COUNT, &
                     "rmsd_analysis_lazy_c: atom count mismatch")
      return
    end if

    ! Initialize sink (array mode)
    call init_sink_array(c%sink, result_f, c%result_size)

    ! Run unified RMSD analysis (lazy loading via source abstraction)
    write(MsgOut,'(A)') '[STEP2] RMSD Analysis (lazy loading, unified)'
    write(MsgOut,'(A)') ' '

    call analyze_rmsd_unified(c%source, c%sink, ref_coord_f, mass_f, &
                              c%n_atoms, &
                              c%fitting_idx_copy, n_fitting_use, &
                              c%analysis_idx_copy, c%n_analysis, &
                              fitting_method_use, use_mass, nstru)

    c%nstru_out = nstru

  end subroutine rmsd_analysis_lazy_body

end module rmsd_c_mod
