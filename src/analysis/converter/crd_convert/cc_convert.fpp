!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cc_convert_mod
!> @brief   convert trajectory files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module cc_convert_mod

  use, intrinsic :: iso_c_binding
  use cc_option_mod
  use cc_option_str_mod
  use pbc_correct_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_mod
  use molecules_str_mod
  use constants_mod
  use measure_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
  use string_mod
  use fileio_pdb_mod
 
  implicit none
  private

  ! subroutines
  public  :: convert
  public  :: count_trj_frames
  public  :: convert_to_arrays
  private :: convert_frame
  private :: centering
  private :: c2f_string_at_offset
  private :: output_split_trjpdb
  private :: get_filename

  real(wp), allocatable, save   :: tmp_coord(:,:)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert trajectory files
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine convert(molecule,   &
                    trj_list,   &
                    trajectory, &
                    fitting,    &
                    option,     &
                    output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_trj_file)         :: trj_in, trj_out
    integer                  :: nstru, irun, itrj
    integer                  :: rms_out, trr_out
    integer,     allocatable :: cen_idx(:)


    ! check-only
    if (option%check_only) &
      return

    ! open output files
    if (output%trjfile /= '' .and. .not. option%split_trjpdb) &
      call open_trj (trj_out,              &
                     output%trjfile,       &
                     option%trjout_format, &
                     option%trjout_type,   &
                     IOFileOutputNew)

    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)

    if (output%trrfile /= '') &
      call open_file(trr_out, output%trrfile, IOFileOutputNew)
    ! centering atoms (empty when centering is off)
    if (option%centering .and. allocated(option%centering_atom%idx)) then
      cen_idx = option%centering_atom%idx
    else
      allocate(cen_idx(0))
    end if

    nstru = 0
    do irun = 1, size(trj_list%md_steps)

      call open_trj(trj_in,                   &
                    trj_list%filenames(irun), &
                    trj_list%trj_format,      &
                    trj_list%trj_type,        &
                    IOFileInput)

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        call read_trj(trj_in, trajectory)

        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru          

          ! selection
          !
          call reselect_atom(molecule, &
                             option%trjout_atom_exp, &
                             trajectory%coord, &
                             option%trjout_atom, &
                             option%trjout_atom_trj)

          ! centering
          !
          call convert_frame(molecule, trajectory, fitting, &
                             option%centering .and. size(cen_idx) > 0, cen_idx, &
                             option%center_coord, option%pbcc_mode)

          ! write data
          !
          if (output%rmsfile /= '') &
            call out_rmsd (rms_out, nstru, fitting)

          if (output%trrfile /= '') &
            call out_trrot(trr_out, nstru, fitting)

          if (output%trjfile /= '') then
            if (option%split_trjpdb) then
              call output_split_trjpdb(nstru, molecule, trajectory, option, output)
            else
              call write_trj(trj_out, trajectory, option%trjout_atom_trj,molecule)
            end if
          end if

        end if

      end do

      call close_trj(trj_in)

    end do

    if (output%trrfile /= '') call close_file(trr_out)
    if (output%rmsfile /= '') call close_file(rms_out)
    if (output%trjfile /= '' .and. .not. option%split_trjpdb) &
                              call close_trj (trj_out)

    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    centering
  !> @brief        move the COM of the fitting target to the origin
  !! @authors      DM
  !! @param[in]    molecule   : molecule information
  !! @param[inout] coord      : atom coordinates
  !! @param[in]    fitting    : fitting information
  !! @param[in]    option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine centering(molecule, coord, option)

  ! formal arguments
  type(s_molecule), intent(in)    :: molecule
  real(wp),         intent(inout) :: coord(:,:)
  type(s_option),   intent(in)    :: option

  ! local variables
  integer  :: iatm, natm
  real(wp) :: com(3)

  if (.not. option%centering) return

  natm = size(molecule%atom_no)
  com  = compute_com(coord, molecule%mass, option%centering_atom%idx)

  do iatm = 1, natm
    coord(:,iatm) = coord(:,iatm) - com(:) + option%center_coord(:)
  end do

  return

  end subroutine centering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_split_trjpdb
  !> @brief        output split PDB files as trajectory
  !! @authors      TM
  !! @param[inout] nstru      : structure index
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_split_trjpdb(nstru, molecule, trajectory, option, output)

    ! formal arguments
    integer,                 intent(inout) :: nstru
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: tmp_pdb


    ! allocate tmp_coord to save molecule%atom_coord
    !
    if(.not. allocated(tmp_coord)) then
      allocate(tmp_coord(3,molecule%num_atoms))
    end if

    tmp_coord(:,:) = molecule%atom_coord(:,:)
    molecule%atom_coord(:,:) = trajectory%coord(:,:)

    call export_molecules(molecule, option%trjout_atom, tmp_pdb)

    if (option%trjout_type == TrjTypeCoorBox) then
      tmp_pdb%cryst_rec = .true.
      tmp_pdb%pbc_box(1:3,1:3) = trajectory%pbc_box(1:3,1:3)
    else
      tmp_pdb%cryst_rec = .false.
    end if

    tmp_pdb%model_rec = .false.

    call output_pdb(get_filename(output%trjfile,nstru), tmp_pdb)

    molecule%atom_coord(:,:) = tmp_coord(:,:)

    return

  end subroutine output_split_trjpdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_filename
  !> @brief        insert snapshot index into {} in the filename
  !! @authors      TM
  !! @param[in]    filename      : filename
  !! @param[in]    no            : index
  !! @note         this subroutine was originally made by NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_filename(filename, no)

    ! return
    character(Maxfilename)   :: get_filename

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br
    character(100)           :: fid


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Filename> {} is not found in the output trjfile name')

    write(fid,'(i0)') no
    get_filename = filename(:bl-1) // trim(fid) // filename(br+1:)

    return

  end function get_filename

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert_frame
  !> @brief        centering, PBC correction and fitting of one frame
  !!               (shared by the CLI and the Python interface)
  !! @authors      NT, Claude Code
  !! @param[inout] molecule      : molecule information
  !! @param[inout] trajectory    : trajectory information (modified in place)
  !! @param[inout] fitting       : fitting information (NO: no fitting)
  !! @param[in]    do_centering  : move the COM of centering_idx to center_coord
  !! @param[in]    centering_idx : atoms whose COM is centered
  !! @param[in]    center_coord  : target of the centering
  !! @param[in]    pbcc_mode     : PBC correction mode
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert_frame(molecule, trajectory, fitting, do_centering, &
                           centering_idx, center_coord, pbcc_mode)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    logical,                 intent(in)    :: do_centering
    integer,                 intent(in)    :: centering_idx(:)
    real(wp),                intent(in)    :: center_coord(3)
    integer,                 intent(in)    :: pbcc_mode

    ! local variables
    real(wp)                 :: com(3)
    integer                  :: iatm


    if (do_centering) then
      com = compute_com(trajectory%coord, molecule%mass, centering_idx)
      do iatm = 1, size(trajectory%coord, 2)
        trajectory%coord(:,iatm) = trajectory%coord(:,iatm) - com(:) + center_coord(:)
      end do
    end if

    call run_pbc_correct(pbcc_mode, molecule, trajectory)

    if (fitting%mass_weight) then
      call run_fitting(fitting, &
                       molecule%atom_coord, &
                       trajectory%coord, &
                       trajectory%coord, &
                       molecule%mass)
    else
      call run_fitting(fitting, &
                       molecule%atom_coord, &
                       trajectory%coord, &
                       trajectory%coord)
    end if

    return

  end subroutine convert_frame

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_trj_frames
  !> @brief        number of frames of every trajectory file
  !! @authors      Claude Code
  !! @param[in]    trj_filenames : file names, filename_len chars each
  !! @param[in]    n_trj_files   : number of files
  !! @param[in]    filename_len  : length of every name slot
  !! @param[in]    trj_format    : trajectory format
  !! @param[in]    trj_type      : trajectory type
  !! @param[out]   frame_counts  : frames per file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_trj_frames(trj_filenames, n_trj_files, filename_len, &
                              trj_format, trj_type, frame_counts)

    ! formal arguments
    character(kind=c_char),  intent(in)    :: trj_filenames(*)
    integer,                 intent(in)    :: n_trj_files
    integer,                 intent(in)    :: filename_len
    integer,                 intent(in)    :: trj_format
    integer,                 intent(in)    :: trj_type
    integer,                 intent(out)   :: frame_counts(n_trj_files)

    ! local variables
    integer                  :: i, offset
    character(MaxFilename)   :: filename_f


    do i = 1, n_trj_files
      offset = (i - 1) * filename_len
      call c2f_string_at_offset(trj_filenames, offset, filename_len, filename_f)
      frame_counts(i) = get_num_steps_trj(trim(filename_f), trj_format, trj_type)
      if (frame_counts(i) <= 0) &
        call error_msg('Count_Trj_Frames> Cannot read frame count from: ' &
                       // trim(filename_f), code=202)
    end do

    return

  end subroutine count_trj_frames

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert_to_arrays
  !> @brief        read trajectory files, process every frame with
  !!               convert_frame and store the selected atoms into
  !!               pre-allocated per-file arrays (Python zero-copy path)
  !! @authors      Claude Code
  !! @param[inout] molecule         : molecule information
  !! @param[in]    trj_filenames    : file names, filename_len chars each
  !! @param[in]    n_trj_files      : number of files
  !! @param[in]    filename_len     : length of every name slot
  !! @param[in]    trj_format       : trajectory format
  !! @param[in]    trj_type         : trajectory type
  !! @param[in]    selected_indices : atoms copied to the output (1-based)
  !! @param[in]    fitting_method   : fitting method (FittingMethodNO: none)
  !! @param[in]    fitting_indices  : fitting atoms (n_fitting may be 0)
  !! @param[in]    mass_weighted    : mass-weighted fitting (0/1)
  !! @param[in]    do_centering     : centering (0/1) of centering_indices
  !! @param[in]    center_coord     : target of the centering
  !! @param[in]    pbcc_mode        : PBC correction mode
  !! @param[in]    ana_period       : analysis period (frames p, 2p, ...)
  !! @param[in]    frame_counts     : frames of every file
  !! @param[in]    coords_ptrs      : per file, (3, n_selected, nframe) buffer
  !! @param[in]    pbc_box_ptrs     : per file, (3, 3, nframe) buffer
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert_to_arrays(molecule, &
                               trj_filenames, n_trj_files, filename_len, &
                               trj_format, trj_type, &
                               selected_indices, n_selected, &
                               fitting_method, fitting_indices, n_fitting, &
                               mass_weighted, &
                               do_centering, centering_indices, n_centering, &
                               center_coord, &
                               pbcc_mode, ana_period, &
                               frame_counts, &
                               coords_ptrs, pbc_box_ptrs)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    character(kind=c_char),  intent(in)    :: trj_filenames(*)
    integer,                 intent(in)    :: n_trj_files
    integer,                 intent(in)    :: filename_len
    integer,                 intent(in)    :: trj_format
    integer,                 intent(in)    :: trj_type
    integer,                 intent(in)    :: selected_indices(*)
    integer,                 intent(in)    :: n_selected
    integer,                 intent(in)    :: fitting_method
    integer,                 intent(in)    :: fitting_indices(*)
    integer,                 intent(in)    :: n_fitting
    integer,                 intent(in)    :: mass_weighted
    integer,                 intent(in)    :: do_centering
    integer,                 intent(in)    :: centering_indices(*)
    integer,                 intent(in)    :: n_centering
    real(wp),                intent(in)    :: center_coord(3)
    integer,                 intent(in)    :: pbcc_mode
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: frame_counts(n_trj_files)
    type(c_ptr),             intent(in)    :: coords_ptrs(n_trj_files)
    type(c_ptr),             intent(in)    :: pbc_box_ptrs(n_trj_files)

    ! local variables
    type(s_trj_file)         :: trj_in
    type(s_trajectory)       :: trajectory
    type(s_fitting)          :: fitting
    integer                  :: i, j, k, irun, itrj, frame_idx, offset
    logical                  :: centering
    character(MaxFilename)   :: filename_f
    real(wp),        pointer :: coords_f(:,:,:)   ! (3, n_selected, nframes)
    real(wp),        pointer :: pbc_box_f(:,:,:)  ! (3, 3, nframes)
    integer,     allocatable :: sel_idx(:), cen_idx(:)


    allocate(sel_idx(n_selected))
    do i = 1, n_selected
      sel_idx(i) = selected_indices(i)
    end do

    centering = (do_centering /= 0 .and. n_centering > 0)
    allocate(cen_idx(max(n_centering, 0)))
    do i = 1, n_centering
      cen_idx(i) = centering_indices(i)
    end do

    call alloc_trajectory(trajectory, molecule%num_atoms)

    fitting%fitting_method = FittingMethodNO
    if (fitting_method /= FittingMethodNO .and. n_fitting > 0) then
      call alloc_selatoms(fitting%fitting_atom, n_fitting)
      fitting%fitting_method = fitting_method
      fitting%fitting_atom%idx(1:n_fitting) = fitting_indices(1:n_fitting)
      fitting%mass_weight = (mass_weighted /= 0)
    end if

    do irun = 1, n_trj_files

      offset = (irun - 1) * filename_len
      call c2f_string_at_offset(trj_filenames, offset, filename_len, filename_f)
      call open_trj(trj_in, trim(filename_f), trj_format, trj_type, IOFileInput)

      call c_f_pointer(coords_ptrs(irun), coords_f, [3, n_selected, frame_counts(irun)])
      call c_f_pointer(pbc_box_ptrs(irun), pbc_box_f, [3, 3, frame_counts(irun)])

      frame_idx = 0
      do itrj = 1, frame_counts(irun)

        call read_trj(trj_in, trajectory)
        if (mod(itrj, ana_period) /= 0) cycle
        frame_idx = frame_idx + 1

        call convert_frame(molecule, trajectory, fitting, centering, &
                           cen_idx, center_coord, pbcc_mode)

        do j = 1, n_selected
          k = sel_idx(j)
          coords_f(:, j, frame_idx) = trajectory%coord(:, k)
        end do
        pbc_box_f(:, :, frame_idx) = trajectory%pbc_box(:, :)

      end do

      call close_trj(trj_in)

    end do

    if (fitting%fitting_method /= FittingMethodNO) &
      call dealloc_fitting(fitting)
    call dealloc_trajectory(trajectory)
    deallocate(sel_idx, cen_idx)

    return

  end subroutine convert_to_arrays

  !> Extract the i-th fixed-width C string of a packed name array.
  subroutine c2f_string_at_offset(c_str_array, offset, max_len, f_str)

    ! formal arguments
    character(kind=c_char),  intent(in)    :: c_str_array(*)
    integer,                 intent(in)    :: offset, max_len
    character(len=*),        intent(out)   :: f_str

    ! local variables
    integer                  :: i, pos

    f_str = ''
    do i = 1, min(max_len, len(f_str))
      pos = offset + i
      if (c_str_array(pos) == c_null_char) exit
      f_str(i:i) = c_str_array(pos)
    end do

    return

  end subroutine c2f_string_at_offset

end module cc_convert_mod
