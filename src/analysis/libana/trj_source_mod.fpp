!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   trj_source_mod
!> @brief   Abstract trajectory source for unified analysis
!! @authors Claude Code
!
!  (c) Copyright 2024 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module trj_source_mod

  use, intrinsic :: iso_c_binding
  use trajectory_str_mod
  use fileio_trj_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! Source type enumeration
  integer, parameter, public :: TRJ_SOURCE_FILE     = 1
  integer, parameter, public :: TRJ_SOURCE_MEMORY   = 2
  integer, parameter, public :: TRJ_SOURCE_LAZY_DCD = 3

  ! Abstract trajectory source type
  type, public :: s_trj_source
    integer :: source_type = TRJ_SOURCE_FILE

    ! File-based source (CLI mode)
    type(s_trj_file)          :: trj_file
    type(s_trj_list), pointer :: trj_list => null()
    integer                   :: current_file = 0
    integer                   :: current_step = 0
    integer                   :: file_natom = 0

    ! Memory-based source (Python mode)
    type(c_ptr)     :: coords_ptr = c_null_ptr      ! double[nframe][natom][3]
    type(c_ptr)     :: pbc_boxes_ptr = c_null_ptr   ! double[nframe][3][3]
    integer         :: mem_nframe = 0
    integer         :: mem_natom = 0
    integer         :: mem_current = 0

    ! Lazy DCD source (random access)
    integer                :: dcd_unit = -1
    integer                :: dcd_natom = 0
    integer                :: dcd_nframe = 0
    integer(8)             :: dcd_header_size = 0
    integer(8)             :: dcd_frame_size = 0
    logical                :: dcd_has_box = .false.
    logical                :: dcd_open = .false.
    integer                :: lazy_current = 0
    integer                :: dcd_selected_natom = 0
    integer, allocatable   :: dcd_selection(:)
    real(4), allocatable   :: dcd_x(:), dcd_y(:), dcd_z(:)

    ! Common
    integer :: ana_period = 1
    integer :: total_frames = 0           ! Total frames to process
    integer :: analyzed_count = 0         ! Frames actually analyzed
  end type s_trj_source

  ! Public subroutines
  public :: init_source_file
  public :: init_source_memory
  public :: init_source_lazy_dcd
  public :: get_next_frame
  public :: get_frame_by_index
  public :: get_total_frames
  public :: get_analyzed_frames
  public :: has_more_frames
  public :: reset_source
  public :: finalize_source
  public :: c_filename_to_fortran

contains

  subroutine c_filename_to_fortran(c_filename, filename_len, filename)

    character(kind=c_char), intent(in)  :: c_filename(*)
    integer,                intent(in)  :: filename_len
    character(*),           intent(out) :: filename
    integer :: i

    filename = ''
    do i = 1, min(filename_len, len(filename))
      if (c_filename(i) == c_null_char) exit
      filename(i:i) = c_filename(i)
    end do

  end subroutine c_filename_to_fortran

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_source_file
  !> @brief        Initialize file-based trajectory source (CLI mode)
  !! @authors      Claude Code
  !! @param[inout] source   : trajectory source
  !! @param[in]    trj_list : trajectory file list
  !! @param[in]    n_atoms  : number of atoms (for trajectory allocation)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_source_file(source, trj_list, n_atoms)

    ! formal arguments
    type(s_trj_source),       intent(inout) :: source
    type(s_trj_list), target, intent(in)    :: trj_list
    integer,                  intent(in)    :: n_atoms

    ! local variables
    integer :: ifile, total

    source%source_type = TRJ_SOURCE_FILE
    source%trj_list => trj_list
    source%current_file = 0
    source%current_step = 0
    source%file_natom = n_atoms
    source%analyzed_count = 0

    ! Calculate total frames to analyze
    total = 0
    do ifile = 1, size(trj_list%md_steps)
      total = total + (trj_list%md_steps(ifile) / trj_list%ana_periods(ifile))
    end do
    source%total_frames = total

    return

  end subroutine init_source_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_source_memory
  !> @brief        Initialize memory-based trajectory source (Python mode)
  !! @authors      Claude Code
  !! @param[inout] source     : trajectory source
  !! @param[in]    coords     : coordinates pointer (3, natom, nframe)
  !! @param[in]    pbc_boxes  : PBC boxes pointer (3, 3, nframe)
  !! @param[in]    natom      : number of atoms
  !! @param[in]    nframe     : number of frames
  !! @param[in]    ana_period : analysis period
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_source_memory(source, coords, pbc_boxes, natom, nframe, ana_period)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    type(c_ptr),        intent(in)    :: coords
    type(c_ptr),        intent(in)    :: pbc_boxes
    integer,            intent(in)    :: natom
    integer,            intent(in)    :: nframe
    integer,            intent(in)    :: ana_period

    source%source_type = TRJ_SOURCE_MEMORY
    source%coords_ptr = coords
    source%pbc_boxes_ptr = pbc_boxes
    source%mem_natom = natom
    source%mem_nframe = nframe
    source%mem_current = 0
    source%ana_period = ana_period
    source%total_frames = nframe / ana_period
    source%analyzed_count = 0

    return

  end subroutine init_source_memory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_source_lazy_dcd
  !> @brief        Initialize lazy DCD trajectory source (random access)
  !! @authors      Claude Code
  !! @param[inout] source     : trajectory source
  !! @param[in]    filename   : DCD file name
  !! @param[in]    trj_type   : trajectory type (COOR or COOR+BOX)
  !! @param[in]    ana_period : analysis period
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_source_lazy_dcd(source, filename, trj_type, ana_period, &
                                  selection_idx, n_selected, status)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    character(*),       intent(in)    :: filename
    integer,            intent(in)    :: trj_type
    integer,            intent(in)    :: ana_period
    integer,            intent(in)    :: selection_idx(:)
    integer,            intent(in)    :: n_selected
    integer,            intent(out)   :: status

    ! local variables
    integer :: unit_no, iostat
    integer :: header_int(20), ntitle
    integer(4) :: rec_size
    character(80) :: title
    character(4) :: signature
    integer :: i
    integer(8) :: file_size, payload_size
    integer(8) :: coor_frame_size, box_frame_size
    integer :: nframe_coor, nframe_box
    logical :: valid_coor, valid_box

    source%source_type = TRJ_SOURCE_LAZY_DCD
    status = 0
    if (ana_period <= 0 .or. n_selected <= 0 .or. &
        size(selection_idx) < n_selected) then
      status = 301
      return
    end if
    source%ana_period = ana_period
    source%lazy_current = 0
    source%analyzed_count = 0
    source%dcd_has_box = (trj_type == TrjTypeCoorBox)

    ! Open DCD file with stream access for fseek
    unit_no = get_unit_no()
    open(unit_no, file=trim(filename), status='old', &
         form='unformatted', access='stream', convert='little_endian', &
         iostat=iostat)
    if (iostat /= 0) then
      call free_unit_no(unit_no)
      status = 201
      return
    end if
    source%dcd_unit = unit_no
    source%dcd_open = .true.

    ! Read first record marker and reopen with the correct endian conversion.
    read(unit_no, iostat=iostat) rec_size
    if (iostat /= 0) then
      call close_file(unit_no)
      source%dcd_open = .false.
      status = 202
      return
    end if
    if (rec_size /= 84) then
      close(unit_no)
      open(unit_no, file=trim(filename), status='old', &
           form='unformatted', access='stream', convert='big_endian', &
           iostat=iostat)
      if (iostat /= 0) then
        call free_unit_no(unit_no)
        source%dcd_open = .false.
        status = 202
        return
      end if
      read(unit_no, iostat=iostat) rec_size
      if (iostat /= 0 .or. rec_size /= 84) then
        call close_file(unit_no)
        source%dcd_open = .false.
        status = 202
        return
      end if
    end if

    ! Read signature 'CORD' (4 bytes)
    read(unit_no, iostat=iostat) signature
    if (iostat /= 0 .or. &
        (signature /= 'CORD' .and. signature /= 'VELD')) then
      call close_file(unit_no)
      source%dcd_open = .false.
      status = 202
      return
    end if

    ! Read header integers (20 ints = 80 bytes)
    read(unit_no, iostat=iostat) header_int(1:20)
    if (iostat /= 0) goto 900
    read(unit_no, iostat=iostat) rec_size  ! end marker
    if (iostat /= 0 .or. rec_size /= 84) goto 900

    source%dcd_nframe = header_int(1)

    ! Read title section
    read(unit_no, iostat=iostat) rec_size  ! title block start marker
    if (iostat /= 0) goto 900
    read(unit_no, iostat=iostat) ntitle    ! number of titles
    if (iostat /= 0 .or. ntitle < 0) goto 900
    do i = 1, ntitle
      read(unit_no, iostat=iostat) title
      if (iostat /= 0) goto 900
    end do
    read(unit_no, iostat=iostat) rec_size  ! title block end marker
    if (iostat /= 0 .or. rec_size /= 4 + 80*ntitle) goto 900

    ! Read atom count
    read(unit_no, iostat=iostat) rec_size
    if (iostat /= 0 .or. rec_size /= 4) goto 900
    read(unit_no, iostat=iostat) source%dcd_natom
    if (iostat /= 0 .or. source%dcd_natom <= 0) goto 900
    read(unit_no, iostat=iostat) rec_size
    if (iostat /= 0 .or. rec_size /= 4) goto 900

    ! Calculate header size (bytes)
    ! First block: 4 + 4 + 80 + 4 = 92 (marker + signature + 20 ints + marker)
    ! Title block: 4 + 4 + 80*ntitle + 4 = 12 + 80*ntitle
    ! Atom block: 4 + 4 + 4 = 12
    source%dcd_header_size = 92 + 12 + 80*ntitle + 12

    ! Calculate frame size (bytes)
    ! Each coordinate block: 4 + 4*natom + 4 = 8 + 4*natom
    ! Three coordinate blocks (x, y, z): 3 * (8 + 4*natom)
    coor_frame_size = int(3 * (8 + 4*source%dcd_natom), 8)
    box_frame_size = coor_frame_size + 56_8

    inquire(unit_no, size=file_size, iostat=iostat)
    if (iostat /= 0 .or. file_size < source%dcd_header_size) goto 900
    payload_size = file_size - source%dcd_header_size
    valid_coor = (mod(payload_size, coor_frame_size) == 0)
    valid_box = (mod(payload_size, box_frame_size) == 0)
    nframe_coor = int(payload_size / coor_frame_size)
    nframe_box = int(payload_size / box_frame_size)

    if (valid_coor .and. valid_box .and. header_int(1) > 0) then
      if (nframe_box == header_int(1)) then
        valid_coor = .false.
      else if (nframe_coor == header_int(1)) then
        valid_box = .false.
      end if
    end if
    if (.not. valid_coor .and. .not. valid_box) goto 900
    if (valid_coor .neqv. valid_box) then
      source%dcd_has_box = valid_box
    end if
    if (source%dcd_has_box) then
      if (.not. valid_box) goto 900
      source%dcd_frame_size = box_frame_size
      source%dcd_nframe = nframe_box
    else
      if (.not. valid_coor) goto 900
      source%dcd_frame_size = coor_frame_size
      source%dcd_nframe = nframe_coor
    end if
    if (header_int(1) > 0 .and. source%dcd_nframe /= header_int(1)) goto 900

    source%dcd_selected_natom = n_selected
    allocate(source%dcd_selection(n_selected))
    source%dcd_selection(:) = selection_idx(1:n_selected)
    allocate(source%dcd_x(source%dcd_natom), source%dcd_y(source%dcd_natom), &
             source%dcd_z(source%dcd_natom))

    source%total_frames = source%dcd_nframe / ana_period

    return

900 continue
    call finalize_source(source)
    status = 202
    return

  end subroutine init_source_lazy_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_next_frame
  !> @brief        Get next frame from trajectory source
  !! @authors      Claude Code
  !! @param[inout] source     : trajectory source
  !! @param[inout] trajectory : output trajectory
  !! @param[out]   status     : 0 = success, 1 = no more frames
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_next_frame(source, trajectory, status)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    type(s_trajectory), intent(inout) :: trajectory
    integer,            intent(out)   :: status

    status = 0

    select case(source%source_type)

    case(TRJ_SOURCE_FILE)
      call get_next_frame_file(source, trajectory, status)

    case(TRJ_SOURCE_MEMORY)
      call get_next_frame_memory(source, trajectory, status)

    case(TRJ_SOURCE_LAZY_DCD)
      call get_next_frame_lazy_dcd(source, trajectory, status)

    case default
      status = 1

    end select

    if (status == 0) source%analyzed_count = source%analyzed_count + 1

    return

  end subroutine get_next_frame

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_next_frame_file
  !> @brief        Get next frame from file-based source
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_next_frame_file(source, trajectory, status)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    type(s_trajectory), intent(inout) :: trajectory
    integer,            intent(out)   :: status

    ! local variables
    integer :: num_files, md_steps, ana_period

    status = 1
    num_files = size(source%trj_list%md_steps)

    ! Ensure trajectory is allocated
    if (.not. allocated(trajectory%coord)) then
      allocate(trajectory%coord(3, source%file_natom))
    else if (size(trajectory%coord, 2) /= source%file_natom) then
      deallocate(trajectory%coord)
      allocate(trajectory%coord(3, source%file_natom))
    end if

    ! Loop through files and steps
    do while (.true.)

      ! Check if we need to open next file
      if (source%current_file == 0 .or. &
          source%current_step >= source%trj_list%md_steps(source%current_file)) then

        ! Close previous file if open
        if (source%current_file > 0) then
          call close_trj(source%trj_file)
        end if

        ! Move to next file
        source%current_file = source%current_file + 1
        source%current_step = 0

        ! Check if done
        if (source%current_file > num_files) then
          status = 1
          return
        end if

        ! Open next file
        call open_trj(source%trj_file, &
                      source%trj_list%filenames(source%current_file), &
                      source%trj_list%trj_format, &
                      source%trj_list%trj_type, &
                      IOFileInput)
      end if

      ! Read next step
      source%current_step = source%current_step + 1
      call read_trj(source%trj_file, trajectory)

      ! Check if this step should be analyzed
      ana_period = source%trj_list%ana_periods(source%current_file)
      if (mod(source%current_step, ana_period) == 0) then
        status = 0
        return
      end if

    end do

    return

  end subroutine get_next_frame_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_next_frame_memory
  !> @brief        Get next frame from memory-based source
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_next_frame_memory(source, trajectory, status)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    type(s_trajectory), intent(inout) :: trajectory
    integer,            intent(out)   :: status

    ! local variables
    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)

    status = 1

    ! Find next frame to analyze
    do while (source%mem_current < source%mem_nframe)
      source%mem_current = source%mem_current + 1

      if (mod(source%mem_current, source%ana_period) == 0) then
        ! Allocate trajectory if needed
        if (.not. allocated(trajectory%coord)) then
          allocate(trajectory%coord(3, source%mem_natom))
        else if (size(trajectory%coord, 2) /= source%mem_natom) then
          deallocate(trajectory%coord)
          allocate(trajectory%coord(3, source%mem_natom))
        end if

        ! Get data via C_F_POINTER
        call C_F_POINTER(source%coords_ptr, coords, &
                         [3, source%mem_natom, source%mem_nframe])
        call C_F_POINTER(source%pbc_boxes_ptr, boxes, &
                         [3, 3, source%mem_nframe])

        trajectory%coord(:,:) = coords(:,:,source%mem_current)
        trajectory%pbc_box(:,:) = boxes(:,:,source%mem_current)

        status = 0
        return
      end if
    end do

    return

  end subroutine get_next_frame_memory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_next_frame_lazy_dcd
  !> @brief        Get next frame from lazy DCD source (sequential)
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_next_frame_lazy_dcd(source, trajectory, status)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    type(s_trajectory), intent(inout) :: trajectory
    integer,            intent(out)   :: status

    status = 1

    ! Find next frame to analyze
    do while (source%lazy_current < source%dcd_nframe)
      source%lazy_current = source%lazy_current + 1

      if (mod(source%lazy_current, source%ana_period) == 0) then
        call get_frame_by_index(source, source%lazy_current, trajectory, status)
        return
      end if
    end do

    return

  end subroutine get_next_frame_lazy_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_frame_by_index
  !> @brief        Get specific frame by index (1-indexed, random access)
  !! @authors      Claude Code
  !! @param[inout] source     : trajectory source
  !! @param[in]    frame_idx  : frame index (1-indexed)
  !! @param[inout] trajectory : output trajectory
  !! @param[out]   status     : 0 = success, 1 = error
  !! @note         Only supported for MEMORY and LAZY_DCD sources
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_frame_by_index(source, frame_idx, trajectory, status)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source
    integer,            intent(in)    :: frame_idx
    type(s_trajectory), intent(inout) :: trajectory
    integer,            intent(out)   :: status

    ! local variables
    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)
    integer(8) :: byte_offset
    real(8)    :: box_data(6), pbcbox(6)
    integer(4) :: rec_size
    integer    :: i, natom, selected_natom, iostat

    status = 1

    select case(source%source_type)

    case(TRJ_SOURCE_MEMORY)
      if (frame_idx < 1 .or. frame_idx > source%mem_nframe) return

      ! Allocate trajectory if needed
      if (.not. allocated(trajectory%coord)) then
        allocate(trajectory%coord(3, source%mem_natom))
      else if (size(trajectory%coord, 2) /= source%mem_natom) then
        deallocate(trajectory%coord)
        allocate(trajectory%coord(3, source%mem_natom))
      end if

      call C_F_POINTER(source%coords_ptr, coords, &
                       [3, source%mem_natom, source%mem_nframe])
      call C_F_POINTER(source%pbc_boxes_ptr, boxes, &
                       [3, 3, source%mem_nframe])

      trajectory%coord(:,:) = coords(:,:,frame_idx)
      trajectory%pbc_box(:,:) = boxes(:,:,frame_idx)
      status = 0

    case(TRJ_SOURCE_LAZY_DCD)
      if (frame_idx < 1 .or. frame_idx > source%dcd_nframe) return
      if (.not. source%dcd_open) return

      natom = source%dcd_natom
      selected_natom = source%dcd_selected_natom

      ! Allocate trajectory if needed
      if (.not. allocated(trajectory%coord)) then
        allocate(trajectory%coord(3, selected_natom))
      else if (size(trajectory%coord, 2) /= selected_natom) then
        deallocate(trajectory%coord)
        allocate(trajectory%coord(3, selected_natom))
      end if

      ! Calculate byte offset for frame (1-indexed position for stream access)
      byte_offset = source%dcd_header_size + &
                    int(frame_idx - 1, 8) * source%dcd_frame_size + 1

      ! Read box if present
      if (source%dcd_has_box) then
        read(source%dcd_unit, pos=byte_offset, iostat=iostat) rec_size
        if (iostat /= 0 .or. rec_size /= 48) goto 910
        read(source%dcd_unit, iostat=iostat) box_data(1:6)
        if (iostat /= 0) goto 910
        read(source%dcd_unit, iostat=iostat) rec_size
        if (iostat /= 0 .or. rec_size /= 48) goto 910
        call to_pbcbox(box_data, pbcbox)
        trajectory%pbc_box(:,:) = 0.0_wp
        trajectory%pbc_box(1,1) = pbcbox(1)
        trajectory%pbc_box(2,2) = pbcbox(2)
        trajectory%pbc_box(3,3) = pbcbox(3)
        read(source%dcd_unit, iostat=iostat) rec_size
        if (iostat /= 0) goto 910
      else
        read(source%dcd_unit, pos=byte_offset, iostat=iostat) rec_size
        if (iostat /= 0) goto 910
        trajectory%pbc_box(:,:) = 0.0_wp
      end if

      ! Read X coordinates
      if (rec_size /= 4*natom) goto 910
      read(source%dcd_unit, iostat=iostat) source%dcd_x(1:natom)
      if (iostat /= 0) goto 910
      read(source%dcd_unit, iostat=iostat) rec_size
      if (iostat /= 0 .or. rec_size /= 4*natom) goto 910

      ! Read Y coordinates
      read(source%dcd_unit, iostat=iostat) rec_size
      if (iostat /= 0 .or. rec_size /= 4*natom) goto 910
      read(source%dcd_unit, iostat=iostat) source%dcd_y(1:natom)
      if (iostat /= 0) goto 910
      read(source%dcd_unit, iostat=iostat) rec_size
      if (iostat /= 0 .or. rec_size /= 4*natom) goto 910

      ! Read Z coordinates
      read(source%dcd_unit, iostat=iostat) rec_size
      if (iostat /= 0 .or. rec_size /= 4*natom) goto 910
      read(source%dcd_unit, iostat=iostat) source%dcd_z(1:natom)
      if (iostat /= 0) goto 910
      read(source%dcd_unit, iostat=iostat) rec_size
      if (iostat /= 0 .or. rec_size /= 4*natom) goto 910

      ! Copy the selected atoms to the logical trajectory view.
      do i = 1, selected_natom
        trajectory%coord(1, i) = real( &
          source%dcd_x(source%dcd_selection(i)), wp)
        trajectory%coord(2, i) = real( &
          source%dcd_y(source%dcd_selection(i)), wp)
        trajectory%coord(3, i) = real( &
          source%dcd_z(source%dcd_selection(i)), wp)
      end do

      status = 0

    case default
      ! FILE source does not support random access
      status = 1

    end select

    return

910 continue
    call finalize_source(source)
    call error_msg('get_frame_by_index> Invalid or truncated DCD frame', code=202)
    status = 1
    return

  end subroutine get_frame_by_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_total_frames
  !> @brief        Get total number of frames to be analyzed
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_total_frames(source) result(total)

    ! formal arguments
    type(s_trj_source), intent(in) :: source
    integer :: total

    total = source%total_frames

    return

  end function get_total_frames

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_analyzed_frames
  !> @brief        Get number of frames analyzed so far
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_analyzed_frames(source) result(count)

    ! formal arguments
    type(s_trj_source), intent(in) :: source
    integer :: count

    count = source%analyzed_count

    return

  end function get_analyzed_frames

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      has_more_frames
  !> @brief        Check if more frames are available
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function has_more_frames(source) result(has_more)

    ! formal arguments
    type(s_trj_source), intent(in) :: source
    logical :: has_more

    select case(source%source_type)

    case(TRJ_SOURCE_FILE)
      has_more = source%current_file <= size(source%trj_list%md_steps)

    case(TRJ_SOURCE_MEMORY)
      has_more = source%mem_current < source%mem_nframe

    case(TRJ_SOURCE_LAZY_DCD)
      has_more = source%lazy_current < source%dcd_nframe

    case default
      has_more = .false.

    end select

    return

  end function has_more_frames

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reset_source
  !> @brief        Reset source to beginning
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reset_source(source)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source

    select case(source%source_type)

    case(TRJ_SOURCE_FILE)
      if (source%current_file > 0) then
        call close_trj(source%trj_file)
      end if
      source%current_file = 0
      source%current_step = 0

    case(TRJ_SOURCE_MEMORY)
      source%mem_current = 0

    case(TRJ_SOURCE_LAZY_DCD)
      source%lazy_current = 0

    end select

    source%analyzed_count = 0

    return

  end subroutine reset_source

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    finalize_source
  !> @brief        Finalize and cleanup trajectory source
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine finalize_source(source)

    ! formal arguments
    type(s_trj_source), intent(inout) :: source

    select case(source%source_type)

    case(TRJ_SOURCE_FILE)
      ! Nested so that a default-initialised source (trj_list => null()) can
      ! be finalized safely: Fortran does not guarantee short-circuit .and.
      if (source%current_file > 0 .and. associated(source%trj_list)) then
        if (source%current_file <= size(source%trj_list%md_steps)) then
          call close_trj(source%trj_file)
        end if
      end if
      source%current_file = 0
      source%trj_list => null()

    case(TRJ_SOURCE_MEMORY)
      source%coords_ptr = c_null_ptr
      source%pbc_boxes_ptr = c_null_ptr
      source%mem_nframe = 0
      source%mem_natom = 0

    case(TRJ_SOURCE_LAZY_DCD)
      if (source%dcd_open) then
        call close_file(source%dcd_unit)
        source%dcd_open = .false.
      end if
      source%dcd_unit = -1
      source%dcd_nframe = 0
      source%dcd_natom = 0
      source%dcd_selected_natom = 0
      if (allocated(source%dcd_selection)) deallocate(source%dcd_selection)
      if (allocated(source%dcd_x)) deallocate(source%dcd_x)
      if (allocated(source%dcd_y)) deallocate(source%dcd_y)
      if (allocated(source%dcd_z)) deallocate(source%dcd_z)

    end select

    source%analyzed_count = 0

    return

  end subroutine finalize_source

end module trj_source_mod
