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
  use cc_convert_mod
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

  private :: t_crd_convert_info_ctx, t_crd_convert_zerocopy_ctx

contains



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
    call count_trj_frames(trj_filenames, c%n_trj_files, c%filename_len, &
                          c%trj_format, c%trj_type, c%frame_counts)
    c%n_trajs = c%n_trj_files
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
    call convert_to_arrays(c%f_molecule, &
                           trj_filenames, c%n_trj_files, c%filename_len, &
                           c%trj_format, c%trj_type, &
                           sel_idx_f, c%n_selected, &
                           c%fitting_method, fit_idx_f, c%n_fitting, &
                           c%mass_weighted, &
                           c%do_centering, cen_idx_f, c%n_centering, &
                           center_coord_f, &
                           c%pbcc_mode, c%ana_period, &
                           frame_counts_f, &
                           coords_ptrs_f, pbc_box_ptrs_f)
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
