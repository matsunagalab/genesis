!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   avecrd_c_mod
!> @brief   Python (bind(C)) entry point of avecrd_analysis
!! @authors Takaharu Mori (TM), Claude Code
!
!  The science lives in aa_analyze_mod (shared with the CLI program): the
!  entry point parses the control text, runs analyze_avecrd_unified on an
!  in-memory trajectory and returns the averaged structure as a PDB string.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module avecrd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use aa_analyze_mod
  use aa_control_mod
  use aa_option_str_mod
  use trj_source_mod
  use result_sink_mod
  use internal_file_type_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: aa_analysis_c

  !> State shared between aa_analysis_c and its guarded body.
  type :: t_avecrd_ctx
    ! inputs
    type(c_ptr) :: molecule_ptr = c_null_ptr     ! caller's s_molecule_c
    type(c_ptr) :: trajes_ptr = c_null_ptr       ! caller's s_trajectories_c
    integer     :: ana_period = 1
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! output
    character(len=:), allocatable :: out_pdb
    ! error state and resources released by the wrapper
    type(s_error)       :: err
    type(s_molecule)    :: f_molecule
    type(s_trj_source)  :: source
  end type t_avecrd_ctx

  private :: t_avecrd_ctx
  private :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    aa_analysis_c
  !> @brief        Average structure of an in-memory trajectory
  !! @authors      Claude Code
  !! @param[in]    molecule        : molecule C structure
  !! @param[in]    s_trajes_c      : trajectories C structure
  !! @param[in]    ana_period      : analysis period
  !! @param[in]    ctrl_text       : control text ([SELECTION], [FITTING], ...)
  !! @param[out]   out_pdb_ave_ptr : averaged structure as a PDB string
  !!                                 (C-allocated, freed by the Python side)
  !! @param[out]   status, msg     : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine aa_analysis_c(molecule, s_trajes_c, ana_period, &
                           ctrl_text, ctrl_len, &
                           out_pdb_ave_ptr, status, msg, msglen) &
        bind(C, name="aa_analysis_c")
    use conv_f_c_util
    implicit none

    ! Arguments
    type(s_molecule_c), intent(in), target :: molecule
    type(s_trajectories_c), intent(in), target :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: out_pdb_ave_ptr
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    ! Local variables
    type(t_avecrd_ctx), target :: c
    character(kind=c_char), pointer :: out_pdb_ave_c(:)

    out_pdb_ave_ptr = c_null_ptr
    c%molecule_ptr  = c_loc(molecule)
    c%trajes_ptr    = c_loc(s_trajes_c)
    c%ana_period    = ana_period
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    call run_guarded(aa_analysis_body, c, c%err, status, msg, msglen)

    call finalize_source(c%source)
    call dealloc_molecules_all(c%f_molecule)
    if (error_has(c%err)) return

    if (allocated(c%out_pdb)) then
      call f2c_string(c%out_pdb, out_pdb_ave_c)
      out_pdb_ave_ptr = c_loc(out_pdb_ave_c(1))
    end if

  end subroutine aa_analysis_c

  !> Guarded body of aa_analysis_c (see run_guarded in error_mod).
  subroutine aa_analysis_body(ctx) bind(C, name="aa_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_avecrd_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule
    type(s_trajectories_c), pointer :: trajes
    character(kind=c_char), pointer :: ctrl_text(:)
    type(s_ctrl_data)   :: ctrl_data
    type(s_fitting)     :: fitting
    type(s_output)      :: output
    type(s_option)      :: option
    type(s_result_sink) :: rms_sink       ! inactive: no RMSD file in Python
    type(s_pdb)         :: pdb_out

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
    call setup(ctrl_data, c%f_molecule, fitting, option, output)

    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '
    call init_source_memory(c%source, trajes%coords, trajes%pbc_boxes, &
                            trajes%natom, trajes%nframe, c%ana_period)
    call analyze_avecrd_unified(c%f_molecule, c%source, fitting, option, rms_sink)

    ! averaged structure as a PDB string
    call export_molecules(c%f_molecule, option%analysis_atom, pdb_out)
    call write_pdb_to_string(c%out_pdb, pdb_out, c%err)
    call dealloc_pdb_all(pdb_out)

  end subroutine aa_analysis_body

  !> Selection, fitting, option and output setup from the control data.
  subroutine setup(ctrl_data, molecule, fitting, option, output)
    use aa_option_mod
    use fitting_mod
    use output_mod
    use select_mod
    implicit none
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    call setup_selection(ctrl_data%sel_info, molecule)
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup

end module avecrd_c_mod
