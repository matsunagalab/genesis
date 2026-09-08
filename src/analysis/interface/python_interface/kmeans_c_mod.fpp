!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   kmeans_c_mod
!> @brief   Python (bind(C)) entry point of kmeans_clustering
!! @authors Takaharu Mori (TM), Claude Code
!
!  The science lives in kc_analyze_mod (shared with the CLI program): the
!  entry point parses the control text, runs analyze_kmeans_unified on an
!  in-memory trajectory and returns the cluster index of every structure and
!  the cluster centers as one PDB string (one MODEL per cluster).
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module kmeans_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use kc_analyze_mod
  use kc_control_mod
  use kc_option_str_mod
  use trj_source_mod
  use internal_file_type_mod
  use fitting_str_mod
  use trajectory_str_mod
  use input_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_trj_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: kc_analysis_c

  !> State shared between kc_analysis_c and its guarded body.
  type :: t_kmeans_ctx
    ! inputs
    type(c_ptr) :: molecule_ptr = c_null_ptr     ! caller's s_molecule_c
    type(c_ptr) :: trajes_ptr = c_null_ptr       ! caller's s_trajectories_c
    integer     :: ana_period = 1
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! outputs
    character(len=:), allocatable :: out_pdb
    integer,          allocatable :: cluster_index(:)
    ! error state and resources released by the wrapper
    type(s_error)       :: err
    type(s_molecule)    :: f_molecule
    type(s_option)      :: option
    type(s_trj_source)  :: source
  end type t_kmeans_ctx

  private :: t_kmeans_ctx
  private :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    kc_analysis_c
  !> @brief        k-means clustering of an in-memory trajectory
  !! @authors      Claude Code
  !! @param[in]    molecule           : molecule C structure
  !! @param[in]    s_trajes_c         : trajectories C structure
  !! @param[in]    ana_period         : analysis period
  !! @param[in]    ctrl_text          : control text
  !! @param[out]   out_pdb_ptr        : cluster centers as a PDB string
  !!                                    (C-allocated, freed by Python)
  !! @param[out]   cluster_index      : cluster of every structure
  !!                                    (C-allocated int array, freed by Python)
  !! @param[out]   size_cluster_index : number of structures
  !! @param[out]   status, msg        : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine kc_analysis_c(molecule, s_trajes_c, ana_period, &
          ctrl_text, ctrl_len, &
          out_pdb_ptr, cluster_index, size_cluster_index, status, msg, msglen)  &
        bind(C, name="kc_analysis_c")
    use conv_f_c_util
    implicit none

    ! Arguments
    type(s_molecule_c), intent(in), target :: molecule
    type(s_trajectories_c), intent(in), target :: s_trajes_c
    integer(c_int), intent(in) :: ana_period
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: out_pdb_ptr
    type(c_ptr), intent(out) :: cluster_index
    integer(c_int), intent(out) :: size_cluster_index
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    ! Local variables
    type(t_kmeans_ctx), target :: c
    character(kind=c_char), pointer :: out_pdb_c(:)

    out_pdb_ptr        = c_null_ptr
    cluster_index      = c_null_ptr
    size_cluster_index = 0

    c%molecule_ptr  = c_loc(molecule)
    c%trajes_ptr    = c_loc(s_trajes_c)
    c%ana_period    = ana_period
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    call run_guarded(kc_analysis_body, c, c%err, status, msg, msglen)

    if (.not. error_has(c%err)) then
      if (allocated(c%out_pdb)) then
        call f2c_string(c%out_pdb, out_pdb_c)
        out_pdb_ptr = c_loc(out_pdb_c(1))
      end if
      if (allocated(c%cluster_index)) then
        cluster_index = f2c_int_array(c%cluster_index)
        size_cluster_index = size(c%cluster_index)
      end if
    end if

    call finalize_source(c%source)
    call dealloc_option(c%option)
    call dealloc_molecules_all(c%f_molecule)

  end subroutine kc_analysis_c

  !> Guarded body of kc_analysis_c (see run_guarded in error_mod).
  subroutine kc_analysis_body(ctx) bind(C, name="kc_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_kmeans_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule
    type(s_trajectories_c), pointer :: trajes
    character(kind=c_char), pointer :: ctrl_text(:)
    type(s_ctrl_data)             :: ctrl_data
    type(s_input)                 :: input
    type(s_fitting)               :: fitting
    type(s_output)                :: output
    integer,          allocatable :: center_index(:)
    type(s_pdb),      allocatable :: center_pdb(:)
    type(s_trj_file), allocatable :: trj_out(:)    ! never allocated: no files
    integer :: iclst

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
    call setup(ctrl_data, c%f_molecule, input, fitting, c%option, output)

    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '
    call init_source_memory(c%source, trajes%coords, trajes%pbc_boxes, &
                            trajes%natom, trajes%nframe, c%ana_period)
    call analyze_kmeans_unified(c%f_molecule, input, c%source, fitting, &
                                c%option, c%cluster_index, center_index,   &
                                .true., center_pdb, trj_out)

    ! cluster centers as one PDB string, in cluster order
    write(MsgOut,'(A)') 'Analyze> output PDB of the cluster centers'
    write(MsgOut,'(A)') ''
    allocate(character(len=1024) :: c%out_pdb)
    c%out_pdb = ' '
    do iclst = 1, size(center_pdb)
      write(MsgOut,'(A,I10,A,I10)') '   cluster = ', iclst, &
                                    '  center structure = ', center_index(iclst)
      call append_pdb_to_string(c%out_pdb, center_pdb(iclst), c%err)
      if (error_has(c%err)) return
      call dealloc_pdb_all(center_pdb(iclst))
    end do

  end subroutine kc_analysis_body

  !> Input, selection, fitting, option and output setup from the control data.
  subroutine setup(ctrl_data, molecule, input, fitting, option, output)
    use kc_option_mod
    use fitting_mod
    use input_mod
    use output_mod
    use select_mod
    implicit none
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_input),           intent(inout) :: input
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    call setup_input(ctrl_data%inp_info, input)
    call setup_selection(ctrl_data%sel_info, molecule)
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup

end module kmeans_c_mod
