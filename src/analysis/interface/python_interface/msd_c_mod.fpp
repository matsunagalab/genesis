!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   msd_c_mod
!> @brief   Python (bind(C)) entry point of msd_analysis
!! @authors Donatas Surblys (DS), Claude Code
!
!  The science lives in ma_analyze_mod (shared with the CLI program): the
!  entry point parses the control text and runs analyze_msd_unified on an
!  in-memory trajectory.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module msd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use ma_analyze_mod
  use ma_control_mod
  use ma_option_str_mod
  use trj_source_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: ma_analysis_c

  !> State shared between ma_analysis_c and its guarded body.
  type :: t_msd_ctx
    ! inputs
    type(c_ptr) :: molecule_ptr = c_null_ptr     ! caller's s_molecule_c
    type(c_ptr) :: trajes_ptr = c_null_ptr       ! caller's s_trajectories_c
    integer     :: ana_period = 1
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! output
    real(wp), allocatable :: msd(:,:)            ! (n_analysis_sets, delta)
    ! error state and resources released by the wrapper
    type(s_error)       :: err
    type(s_molecule)    :: f_molecule
    type(s_option)      :: option
    type(s_trj_source)  :: source
  end type t_msd_ctx

  private :: t_msd_ctx
  private :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    ma_analysis_c
  !> @brief        Mean square displacement of an in-memory trajectory
  !! @authors      Claude Code
  !! @param[in]    molecule          : molecule C structure
  !! @param[in]    s_trajes_c        : trajectories C structure
  !! @param[in]    ana_period        : analysis period
  !! @param[in]    ctrl_text         : control text ([SELECTION], [OPTION], ...)
  !! @param[out]   result_msd        : (num_delta, num_analysis_mols) C-order
  !!                                   array; the caller frees it with
  !!                                   deallocate_double2
  !! @param[out]   num_analysis_mols : number of analysis sets
  !! @param[out]   num_delta         : number of time lags
  !! @param[out]   status, msg       : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine ma_analysis_c(molecule, s_trajes_c, ana_period, &
                           ctrl_text, ctrl_len, &
                           result_msd, num_analysis_mols, num_delta, &
                           status, msg, msglen) &
        bind(C, name="ma_analysis_c")
    implicit none

    ! Arguments
    type(s_molecule_c), intent(in), target :: molecule
    type(s_trajectories_c), intent(in), target :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: result_msd
    integer, intent(out) :: num_analysis_mols
    integer, intent(out) :: num_delta
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    ! Local variables
    type(t_msd_ctx), target :: c
    real(wp), pointer :: msd_ptr(:,:)     ! handed to Python (deallocate_double2)

    result_msd        = c_null_ptr
    num_analysis_mols = 0
    num_delta         = 0

    c%molecule_ptr  = c_loc(molecule)
    c%trajes_ptr    = c_loc(s_trajes_c)
    c%ana_period    = ana_period
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    call run_guarded(ma_analysis_body, c, c%err, status, msg, msglen)

    if (.not. error_has(c%err) .and. allocated(c%msd)) then
      allocate(msd_ptr(size(c%msd, 1), size(c%msd, 2)))
      msd_ptr = c%msd
      num_analysis_mols = size(c%msd, 1)
      num_delta         = size(c%msd, 2)
      result_msd = c_loc(msd_ptr)
    end if

    call finalize_source(c%source)
    call dealloc_option(c%option)
    call dealloc_molecules_all(c%f_molecule)

  end subroutine ma_analysis_c

  !> Guarded body of ma_analysis_c (see run_guarded in error_mod).
  subroutine ma_analysis_body(ctx) bind(C, name="ma_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_msd_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule
    type(s_trajectories_c), pointer :: trajes
    character(kind=c_char), pointer :: ctrl_text(:)
    type(s_ctrl_data) :: ctrl_data
    type(s_output)    :: output

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%molecule_ptr, molecule)
    call c_f_pointer(c%trajes_ptr, trajes)
    call c_f_pointer(c%ctrl_text_ptr, ctrl_text, [c%ctrl_len])

    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    call c2f_s_molecule(molecule, c%f_molecule)

    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '
    call control_from_string(ctrl_text, c%ctrl_len, ctrl_data)

    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '
    call setup(c%f_molecule, ctrl_data, output, c%option)

    write(MsgOut,'(A)') '[STEP3] Analysis of trajectory files'
    write(MsgOut,'(A)') ' '
    call init_source_memory(c%source, trajes%coords, trajes%pbc_boxes, &
                            trajes%natom, trajes%nframe, c%ana_period)
    call analyze_msd_unified(c%f_molecule, c%source, c%option, c%msd)

  end subroutine ma_analysis_body

  !> Output, selection and option setup from the control data.
  subroutine setup(molecule, ctrl_data, output, option)
    use ma_option_mod
    use output_mod
    use select_mod
    use select_molecules_mod
    implicit none
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    type(s_selmols), dimension(:), allocatable      :: selmols
    type(s_one_molecule), dimension(:), allocatable :: allmols

    call setup_output(ctrl_data%out_info, output)
    call setup_selection(ctrl_data%sel_info, molecule)
    call setup_molselect(ctrl_data%molsel_info, ctrl_data%sel_info, &
      molecule, selmols, allmols=allmols)
    call setup_option(ctrl_data%opt_info, selmols, allmols, molecule, option)

    return

  end subroutine setup

end module msd_c_mod
