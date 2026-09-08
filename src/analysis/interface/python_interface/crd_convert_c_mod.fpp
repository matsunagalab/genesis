!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  cc_main
!! @brief   convert MD trajectory format
!! @authors Norio Takase (NT), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module crd_convert_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use crd_convert_impl_mod
  use conv_f_c_util

  use cc_control_mod
  use cc_option_str_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  implicit none

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Contexts shared between the bind(C) entry points below and their guarded
  !  bodies (see run_guarded in error_mod). Reference arguments of the entry
  !  points are stored as C pointers and re-materialised in the body.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  type :: t_crd_convert_ctx
    ! inputs
    type(c_ptr) :: molecule_ptr = c_null_ptr    ! caller's s_molecule_c
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! outputs (copied back to the bind(C) arguments by the wrapper)
    type(c_ptr)    :: s_trajes_c_array = c_null_ptr
    integer(c_int) :: num_trajs = 0
    type(c_ptr)    :: selected_atom_indices = c_null_ptr
    integer(c_int) :: num_selected_atoms = 0
    ! error state and resources released by the wrapper
    type(s_error)    :: err
    type(s_molecule) :: f_molecule
  end type t_crd_convert_ctx

  type :: t_crd_convert_info_ctx
    ! inputs
    type(c_ptr)    :: molecule_ptr = c_null_ptr
    type(c_ptr)    :: trj_filenames_ptr = c_null_ptr
    integer(c_int) :: n_trj_files = 0
    integer(c_int) :: filename_len = 0
    integer(c_int) :: trj_format = 0
    integer(c_int) :: trj_type = 0
    ! outputs (frame_counts is handed to Python on success)
    integer(c_int)          :: n_trajs = 0
    integer(c_int), pointer :: frame_counts(:) => null()
    ! error state and resources released by the wrapper
    type(s_error)    :: err
    type(s_molecule) :: f_molecule
  end type t_crd_convert_info_ctx

  type :: t_crd_convert_zerocopy_ctx
    ! inputs
    type(c_ptr)    :: molecule_ptr = c_null_ptr
    type(c_ptr)    :: trj_filenames_ptr = c_null_ptr
    integer(c_int) :: n_trj_files = 0
    integer(c_int) :: filename_len = 0
    integer(c_int) :: trj_format = 0
    integer(c_int) :: trj_type = 0
    type(c_ptr)    :: selected_indices = c_null_ptr
    integer(c_int) :: n_selected = 0
    integer(c_int) :: fitting_method = 0
    type(c_ptr)    :: fitting_indices = c_null_ptr
    integer(c_int) :: n_fitting = 0
    integer(c_int) :: mass_weighted = 0
    integer(c_int) :: do_centering = 0
    type(c_ptr)    :: centering_indices = c_null_ptr
    integer(c_int) :: n_centering = 0
    type(c_ptr)    :: center_coord = c_null_ptr
    integer(c_int) :: pbcc_mode = 0
    integer(c_int) :: ana_period = 1
    type(c_ptr)    :: frame_counts = c_null_ptr
    type(c_ptr)    :: coords_ptrs = c_null_ptr
    type(c_ptr)    :: pbc_box_ptrs = c_null_ptr
    ! error state and resources released by the wrapper
    type(s_error)    :: err
    type(s_molecule) :: f_molecule
  end type t_crd_convert_zerocopy_ctx

  private :: t_crd_convert_ctx, t_crd_convert_info_ctx, t_crd_convert_zerocopy_ctx
  private :: crd_convert_body, crd_convert_info_body, crd_convert_zerocopy_body

contains
  subroutine crd_convert_c( &
          molecule, ctrl_text, ctrl_len, s_trajes_c_array, num_trajs, &
          selected_atom_indices, num_selected_atoms, status, msg, msglen) &
          bind(C, name="crd_convert_c")
    implicit none
    type(s_molecule_c), intent(inout), target :: molecule
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: s_trajes_c_array
    integer(c_int), intent(out) :: num_trajs
    type(c_ptr), intent(out) :: selected_atom_indices
    integer(c_int), intent(out) :: num_selected_atoms
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    type(t_crd_convert_ctx), target :: c

    c%molecule_ptr  = c_loc(molecule)
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    ! Run under the library-mode error guard: a missing trajectory file calls
    ! error_msg -> exit(1) in CLI mode (see run_guarded in error_mod).
    call run_guarded(crd_convert_body, c, c%err, status, msg, msglen)

    ! The outputs keep their safe defaults (null / 0) when the body aborted,
    ! so the C side never reads uninitialised values.
    s_trajes_c_array      = c%s_trajes_c_array
    num_trajs             = c%num_trajs
    selected_atom_indices = c%selected_atom_indices
    num_selected_atoms    = c%num_selected_atoms

    call dealloc_molecules_all(c%f_molecule)
  end subroutine crd_convert_c

  !> Guarded body of crd_convert_c (see run_guarded in error_mod).
  subroutine crd_convert_body(ctx) bind(C, name="crd_convert_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_crd_convert_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule
    character(kind=c_char), pointer :: ctrl_text(:)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%molecule_ptr, molecule)
    call c_f_pointer(c%ctrl_text_ptr, ctrl_text, [c%ctrl_len])

    call c2f_s_molecule(molecule, c%f_molecule)
    call crd_convert_main(c%f_molecule, ctrl_text, c%ctrl_len, &
                          c%s_trajes_c_array, c%num_trajs, &
                          c%selected_atom_indices, c%num_selected_atoms, &
                          c%err)
  end subroutine crd_convert_body

  subroutine crd_convert_main(molecule, ctrl_text, ctrl_len, s_trajes_c_array, num_trajs, &
                              selected_atom_indices, num_selected_atoms, err)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer, intent(in) :: ctrl_len
    type(c_ptr), intent(out) :: s_trajes_c_array
    integer(c_int), intent(out) :: num_trajs
    type(c_ptr), intent(out) :: selected_atom_indices
    integer(c_int), intent(out) :: num_selected_atoms
    type(s_error),                   intent(inout) :: err
    type(s_ctrl_data)      :: ctrl_data
    type(s_trj_list)       :: trj_list
    type(s_trajectory)     :: trajectory
    type(s_fitting)        :: fitting
    type(s_option)         :: option
    type(s_output)         :: output


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.


    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Convert'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(ctrl_data,  &
               molecule,   &
               trj_list,   &
               trajectory, &
               fitting,    &
               option,     &
               output)


    ! [Step3] Convert trajectory files
    !
    write(MsgOut,'(A)') '[STEP3] Convert trajectory files'
    write(MsgOut,'(A)') ' '

    call convert(molecule,   &
                 trj_list,   &
                 trajectory, &
                 fitting,    &
                 option,     &
                 output,     &
                 s_trajes_c_array, &
                 num_trajs,  &
                 err)
    if (error_has(err)) return

    ! Extract selected atom indices from option%trjout_atom
    call extract_selected_atom_indices(option%trjout_atom, selected_atom_indices, num_selected_atoms)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_option(option)
    call dealloc_fitting(fitting)
    call dealloc_trajectory(trajectory)
    call dealloc_trj_list(trj_list)
  end subroutine crd_convert_main

  subroutine setup(ctrl_data,  &
                   molecule,   &
                   trj_list,   &
                   trajectory, &
                   fitting,    &
                   option,     &
                   output)
    use cc_control_mod
    use cc_option_mod
    use cc_option_str_mod
    use fitting_mod
    use fitting_str_mod
    use input_mod
    use output_mod
    use output_str_mod
    use trajectory_mod
    use trajectory_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_grocrd_mod
    use fileio_grotop_mod
    use fileio_ambcrd_mod
    use fileio_prmtop_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    implicit none

    ! formal arguments
    type(s_ctrl_data),  intent(in)    :: ctrl_data
    type(s_molecule),   intent(inout) :: molecule
    type(s_trj_list),   intent(inout) :: trj_list
    type(s_trajectory), intent(inout) :: trajectory
    type(s_fitting),    intent(inout) :: fitting
    type(s_option),     intent(inout) :: option
    type(s_output),     intent(inout) :: output

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: ref, ref_out
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd

    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)


    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, &
                          molecule, trj_list, trajectory)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup fitting
    !
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    ! export reference molecules
    !
    if (output%pdbfile /= '') then

      call export_molecules(molecule, option%trjout_atom, ref_out)
      call output_pdb(output%pdbfile, ref_out)
      call dealloc_pdb_all(ref_out)

    end if

    return

  end subroutine setup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    extract_selected_atom_indices
  !> @brief        extract selected atom indices from s_selatoms to C array
  !! @authors      Generated
  !! @param[in]    selatoms              : selected atoms structure
  !! @param[out]   selected_atom_indices : C pointer to integer array
  !! @param[out]   num_selected_atoms    : number of selected atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine extract_selected_atom_indices(selatoms, selected_atom_indices, num_selected_atoms)
    use, intrinsic :: iso_c_binding
    use select_atoms_str_mod
    use conv_f_c_util
    use messages_mod
    implicit none

    ! formal arguments
    type(s_selatoms), intent(in) :: selatoms
    type(c_ptr), intent(out) :: selected_atom_indices
    integer(c_int), intent(out) :: num_selected_atoms

    ! local variables
    integer :: i, nsel
    integer(c_int), pointer :: c_array(:)

    nsel = size(selatoms%idx)
    num_selected_atoms = nsel

    if (nsel > 0) then
      ! Allocate C array using conv_f_c_util
      selected_atom_indices = allocate_c_int_array(int(nsel, c_int))
      call c_f_pointer(selected_atom_indices, c_array, [nsel])
      do i = 1, nsel
        c_array(i) = int(selatoms%idx(i), c_int)
      end do
    else
      selected_atom_indices = c_null_ptr
    end if

  end subroutine extract_selected_atom_indices

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    crd_convert_info_c
  !> @brief        C interface to get trajectory info (frame counts)
  !! @param[in]    molecule_c       : molecule structure (C)
  !! @param[in]    trj_filenames    : packed trajectory filenames
  !! @param[in]    n_trj_files      : number of trajectory files
  !! @param[in]    filename_len     : max length per filename
  !! @param[in]    trj_format       : trajectory format
  !! @param[in]    trj_type         : trajectory type
  !! @param[out]   frame_counts_ptr : pointer to frame counts array
  !! @param[out]   n_trajs          : number of trajectories
  !! @param[out]   status           : error status
  !! @param[out]   msg              : error message
  !! @param[in]    msglen           : max message length
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine crd_convert_info_c( &
          molecule_c, &
          trj_filenames, n_trj_files, filename_len, &
          trj_format, trj_type, &
          frame_counts_ptr, n_trajs, &
          status, msg, msglen) &
          bind(C, name="crd_convert_info_c")
    implicit none

    type(s_molecule_c), intent(in), target :: molecule_c
    character(kind=c_char), intent(in), target :: trj_filenames(*)
    integer(c_int), value :: n_trj_files
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_format
    integer(c_int), value :: trj_type
    type(c_ptr), intent(out) :: frame_counts_ptr
    integer(c_int), intent(out) :: n_trajs
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    type(t_crd_convert_info_ctx), target :: c

    c%molecule_ptr      = c_loc(molecule_c)
    c%trj_filenames_ptr = c_loc(trj_filenames)
    c%n_trj_files       = n_trj_files
    c%filename_len      = filename_len
    c%trj_format        = trj_format
    c%trj_type          = trj_type

    ! Guard the file-reading section: a missing/unreadable trajectory calls
    ! error_msg -> exit(1) in CLI mode, which would kill the Python process.
    call run_guarded(crd_convert_info_body, c, c%err, status, msg, msglen)

    n_trajs = c%n_trajs
    if (error_has(c%err)) then
      if (associated(c%frame_counts)) deallocate(c%frame_counts)
      frame_counts_ptr = c_null_ptr
    else
      frame_counts_ptr = c_loc(c%frame_counts(1))
    end if

    call dealloc_molecules_all(c%f_molecule)
  end subroutine crd_convert_info_c

  !> Guarded body of crd_convert_info_c (see run_guarded in error_mod).
  subroutine crd_convert_info_body(ctx) bind(C, name="crd_convert_info_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_crd_convert_info_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule_c
    character(kind=c_char), pointer :: trj_filenames(:)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%molecule_ptr, molecule_c)
    call c_f_pointer(c%trj_filenames_ptr, trj_filenames, &
                     [c%n_trj_files * c%filename_len])

    call c2f_s_molecule(molecule_c, c%f_molecule)

    ! Allocate frame counts array (ownership passes to Python on success)
    allocate(c%frame_counts(c%n_trj_files))

    ! Get trajectory info
    call get_info(c%f_molecule, trj_filenames, c%n_trj_files, c%filename_len, &
                  c%trj_format, c%trj_type, c%frame_counts, c%n_trajs, c%err)
  end subroutine crd_convert_info_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    crd_convert_zerocopy_c
  !> @brief        C interface for zerocopy trajectory conversion
  !! @param[in]    molecule_c         : molecule structure (C)
  !! @param[in]    trj_filenames      : packed trajectory filenames
  !! @param[in]    n_trj_files        : number of trajectory files
  !! @param[in]    filename_len       : max length per filename
  !! @param[in]    trj_format         : trajectory format
  !! @param[in]    trj_type           : trajectory type
  !! @param[in]    selected_indices   : selected atom indices (1-indexed)
  !! @param[in]    n_selected         : number of selected atoms
  !! @param[in]    fitting_method     : fitting method
  !! @param[in]    fitting_indices    : fitting atom indices
  !! @param[in]    n_fitting          : number of fitting atoms
  !! @param[in]    mass_weighted      : use mass weighting
  !! @param[in]    do_centering       : enable centering
  !! @param[in]    centering_indices  : centering atom indices
  !! @param[in]    n_centering        : number of centering atoms
  !! @param[in]    center_coord       : target center coordinates
  !! @param[in]    pbcc_mode          : PBC correction mode
  !! @param[in]    ana_period         : analysis period
  !! @param[in]    frame_counts       : frame counts per trajectory
  !! @param[in]    coords_ptrs        : pre-allocated coords arrays
  !! @param[in]    pbc_box_ptrs       : pre-allocated pbc_box arrays
  !! @param[out]   status             : error status
  !! @param[out]   msg                : error message
  !! @param[in]    msglen             : max message length
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine crd_convert_zerocopy_c( &
          molecule_c, &
          trj_filenames, n_trj_files, filename_len, &
          trj_format, trj_type, &
          selected_indices, n_selected, &
          fitting_method, fitting_indices, n_fitting, mass_weighted, &
          do_centering, centering_indices, n_centering, center_coord, &
          pbcc_mode, ana_period, &
          frame_counts, &
          coords_ptrs, pbc_box_ptrs, &
          status, msg, msglen) &
          bind(C, name="crd_convert_zerocopy_c")
    implicit none

    type(s_molecule_c), intent(in), target :: molecule_c
    character(kind=c_char), intent(in), target :: trj_filenames(*)
    integer(c_int), value :: n_trj_files
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_format
    integer(c_int), value :: trj_type
    type(c_ptr), value :: selected_indices
    integer(c_int), value :: n_selected
    integer(c_int), value :: fitting_method
    type(c_ptr), value :: fitting_indices
    integer(c_int), value :: n_fitting
    integer(c_int), value :: mass_weighted
    integer(c_int), value :: do_centering
    type(c_ptr), value :: centering_indices
    integer(c_int), value :: n_centering
    type(c_ptr), value :: center_coord
    integer(c_int), value :: pbcc_mode
    integer(c_int), value :: ana_period
    type(c_ptr), value :: frame_counts
    type(c_ptr), value :: coords_ptrs
    type(c_ptr), value :: pbc_box_ptrs
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    type(t_crd_convert_zerocopy_ctx), target :: c

    c%molecule_ptr      = c_loc(molecule_c)
    c%trj_filenames_ptr = c_loc(trj_filenames)
    c%n_trj_files       = n_trj_files
    c%filename_len      = filename_len
    c%trj_format        = trj_format
    c%trj_type          = trj_type
    c%selected_indices  = selected_indices
    c%n_selected        = n_selected
    c%fitting_method    = fitting_method
    c%fitting_indices   = fitting_indices
    c%n_fitting         = n_fitting
    c%mass_weighted     = mass_weighted
    c%do_centering      = do_centering
    c%centering_indices = centering_indices
    c%n_centering       = n_centering
    c%center_coord      = center_coord
    c%pbcc_mode         = pbcc_mode
    c%ana_period        = ana_period
    c%frame_counts      = frame_counts
    c%coords_ptrs       = coords_ptrs
    c%pbc_box_ptrs      = pbc_box_ptrs

    ! Guard the conversion: a missing/unreadable trajectory calls error_msg ->
    ! exit(1) in CLI mode, which would kill the Python process.
    call run_guarded(crd_convert_zerocopy_body, c, c%err, status, msg, msglen)

    call dealloc_molecules_all(c%f_molecule)
  end subroutine crd_convert_zerocopy_c

  !> Guarded body of crd_convert_zerocopy_c (see run_guarded in error_mod).
  subroutine crd_convert_zerocopy_body(ctx) &
        bind(C, name="crd_convert_zerocopy_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_crd_convert_zerocopy_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule_c
    character(kind=c_char), pointer :: trj_filenames(:)
    integer(c_int), pointer :: sel_idx_f(:), fit_idx_f(:), cen_idx_f(:)
    integer(c_int), pointer :: frame_counts_f(:)
    real(c_double), pointer :: center_coord_f(:)
    type(c_ptr), pointer :: coords_ptrs_f(:), pbc_box_ptrs_f(:)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%molecule_ptr, molecule_c)
    call c_f_pointer(c%trj_filenames_ptr, trj_filenames, &
                     [c%n_trj_files * c%filename_len])

    call c2f_s_molecule(molecule_c, c%f_molecule)

    ! Get Fortran pointers to C arrays
    call c_f_pointer(c%selected_indices, sel_idx_f, [c%n_selected])
    call c_f_pointer(c%frame_counts, frame_counts_f, [c%n_trj_files])
    call c_f_pointer(c%coords_ptrs, coords_ptrs_f, [c%n_trj_files])
    call c_f_pointer(c%pbc_box_ptrs, pbc_box_ptrs_f, [c%n_trj_files])
    call c_f_pointer(c%center_coord, center_coord_f, [3])

    if (c%n_fitting > 0) then
      call c_f_pointer(c%fitting_indices, fit_idx_f, [c%n_fitting])
    else
      nullify(fit_idx_f)
    end if

    if (c%n_centering > 0) then
      call c_f_pointer(c%centering_indices, cen_idx_f, [c%n_centering])
    else
      nullify(cen_idx_f)
    end if

    ! Call the implementation
    call convert_zerocopy(c%f_molecule, &
                          trj_filenames, c%n_trj_files, c%filename_len, &
                          c%trj_format, c%trj_type, &
                          sel_idx_f, c%n_selected, &
                          c%fitting_method, fit_idx_f, c%n_fitting, &
                          c%mass_weighted, &
                          c%do_centering, cen_idx_f, c%n_centering, &
                          center_coord_f, &
                          c%pbcc_mode, c%ana_period, &
                          frame_counts_f, &
                          coords_ptrs_f, pbc_box_ptrs_f, &
                          c%err)
  end subroutine crd_convert_zerocopy_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_frame_counts_c
  !> @brief        Deallocate frame counts array
  !! @param[in]    frame_counts_ptr : pointer to frame counts array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_frame_counts_c(frame_counts_ptr) &
          bind(C, name="deallocate_frame_counts_c")
    implicit none

    type(c_ptr), value :: frame_counts_ptr
    integer(c_int), pointer :: frame_counts(:)

    if (c_associated(frame_counts_ptr)) then
      call c_f_pointer(frame_counts_ptr, frame_counts, [1])
      deallocate(frame_counts)
    end if

  end subroutine deallocate_frame_counts_c

end module crd_convert_c_mod
