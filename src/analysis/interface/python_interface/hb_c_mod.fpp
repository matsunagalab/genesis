!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hb_c_mod
!> @brief   Python (bind(C)) entry point of hb_analysis
!! @authors Daisuke Matsuoka (DM), Claude Code
!
!  The science lives in hb_analyze_mod (shared with the CLI program): the
!  entry point parses the control text, runs analyze_hb_unified on an
!  in-memory trajectory and returns the result lines as one text block.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hb_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use hb_analyze_mod
  use hb_control_mod
  use hb_option_str_mod
  use trj_source_mod
  use result_sink_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: hb_analysis_c

  !> State shared between hb_analysis_c and its guarded body.
  type :: t_hb_ctx
    ! inputs
    type(c_ptr) :: molecule_ptr = c_null_ptr     ! caller's s_molecule_c
    type(c_ptr) :: trajes_ptr = c_null_ptr       ! caller's s_trajectories_c
    integer     :: ana_period = 1
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! error state and resources released by the wrapper
    type(s_error)       :: err
    type(s_molecule)    :: f_molecule
    type(s_option)      :: option
    type(s_trj_source)  :: source
    type(s_result_sink) :: out_sink        ! text sink: the result lines
    type(s_result_sink) :: hb_list_sink    ! inactive (no hb_listfile)
  end type t_hb_ctx

  private :: t_hb_ctx
  private :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hb_analysis_c
  !> @brief        Hydrogen bond analysis of an in-memory trajectory
  !! @authors      Claude Code
  !! @param[in]    molecule     : molecule C structure
  !! @param[in]    s_trajes_c   : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    ctrl_text    : control text ([SELECTION], [OPTION], ...)
  !! @param[out]   out_text_ptr : result lines (newline separated C string;
  !!                              the caller frees it with deallocate_c_string)
  !! @param[out]   status, msg  : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hb_analysis_c(molecule, s_trajes_c, ana_period, &
                           ctrl_text, ctrl_len, &
                           out_text_ptr, status, msg, msglen) &
        bind(C, name="hb_analysis_c")
    use conv_f_c_util
    implicit none

    ! Arguments
    type(s_molecule_c), intent(in), target :: molecule
    type(s_trajectories_c), intent(in), target :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: out_text_ptr
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    ! Local variables
    type(t_hb_ctx), target :: c
    character(kind=c_char), pointer :: out_text_c(:)

    out_text_ptr = c_null_ptr

    c%molecule_ptr  = c_loc(molecule)
    c%trajes_ptr    = c_loc(s_trajes_c)
    c%ana_period    = ana_period
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    call run_guarded(hb_analysis_body, c, c%err, status, msg, msglen)

    if (.not. error_has(c%err)) then
      if (c%out_sink%text_len > 0) then
        call f2c_string(c%out_sink%text(1:c%out_sink%text_len), out_text_c)
      else
        call f2c_string(' ', out_text_c)
      end if
      out_text_ptr = c_loc(out_text_c(1))
    end if

    call finalize_sink(c%out_sink)
    call finalize_source(c%source)
    call dealloc_option(c%option)
    call dealloc_molecules_all(c%f_molecule)

  end subroutine hb_analysis_c

  !> Guarded body of hb_analysis_c (see run_guarded in error_mod).
  subroutine hb_analysis_body(ctx) bind(C, name="hb_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_hb_ctx), pointer :: c
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

    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '
    call init_source_memory(c%source, trajes%coords, trajes%pbc_boxes, &
                            trajes%natom, trajes%nframe, c%ana_period)
    call init_sink_text(c%out_sink)
    call analyze_hb_unified(c%f_molecule, c%source, c%option, &
                            c%out_sink, c%hb_list_sink)

  end subroutine hb_analysis_body

  !> Selection, option and output setup from the control data.
  subroutine setup(molecule, ctrl_data, output, option)
    use hb_option_mod
    use output_mod
    use select_mod
    implicit none
    type(s_molecule),        intent(inout) :: molecule
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    call setup_selection(ctrl_data%sel_info, molecule)
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup

end module hb_c_mod
