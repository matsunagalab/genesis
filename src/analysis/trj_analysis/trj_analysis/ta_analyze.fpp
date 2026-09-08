!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT), Takaharu Mori (TM), Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_analyze_mod

  use ta_option_str_mod
  use trj_source_mod
  use result_sink_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private

  ! Measurement definitions shared by the CLI (built from s_option) and the
  ! Python interface (built from NumPy arrays). Atom indices are 1-based.
  !
  ! A distance is a weighted sum of one or more atom pairs:
  !   dist(i) = sum_j dist_weight(i,j) * |x(dist_list(2j-1,i)) - x(dist_list(2j,i))|
  !   for j = 1 .. dist_num(i)/2   (plain pairs: dist_num = 2, weight = 1)
  ! COM measurements use flat atom tables: the atoms of group g are
  ! atoms(offsets(g) : offsets(g+1)-1), and pairs/triplets/quads hold the
  ! 1-based group ids of each measurement, 2/3/4 consecutive entries each.
  type, public :: s_trj_measure
    integer,  allocatable :: dist_list(:,:)
    integer,  allocatable :: dist_num(:)
    real(wp), allocatable :: dist_weight(:,:)
    integer,  allocatable :: angl_list(:,:)      ! (3, n_angl)
    integer,  allocatable :: tors_list(:,:)      ! (4, n_tors)
    integer,  allocatable :: cdis_atoms(:)
    integer,  allocatable :: cdis_offsets(:)
    integer,  allocatable :: cdis_pairs(:)       ! (2 * n_cdis)
    integer,  allocatable :: cang_atoms(:)
    integer,  allocatable :: cang_offsets(:)
    integer,  allocatable :: cang_triplets(:)    ! (3 * n_cang)
    integer,  allocatable :: ctor_atoms(:)
    integer,  allocatable :: ctor_offsets(:)
    integer,  allocatable :: ctor_quads(:)       ! (4 * n_ctor)
  end type s_trj_measure

  ! Result streams of analyze_trj_unified (indices into its sinks array)
  integer, public, parameter :: TRJ_DIS     = 1
  integer, public, parameter :: TRJ_ANG     = 2
  integer, public, parameter :: TRJ_TOR     = 3
  integer, public, parameter :: TRJ_CDIS    = 4
  integer, public, parameter :: TRJ_CANG    = 5
  integer, public, parameter :: TRJ_CTOR    = 6
  integer, public, parameter :: TRJ_NSTREAM = 6

  ! subroutines
  public  :: analyze
  public  :: analyze_trj_unified
  public  :: measure_count
  private :: measure_from_option
  private :: measure_frame
  private :: print_output_info

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories (CLI driver)
  !! @authors      NT, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information (kept for the caller;
  !!                            frames are read through trj_source_mod)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule
    type(s_trj_list), target, intent(in)    :: trj_list
    type(s_output),           intent(in)    :: output
    type(s_option),           intent(inout) :: option
    type(s_trajectory),       intent(inout) :: trajectory

    ! local variables
    type(s_trj_source)       :: source
    type(s_result_sink)      :: sinks(TRJ_NSTREAM)
    type(s_trj_measure)      :: measure
    integer                  :: nstru, k


    if (option%check_only) &
      return

    call measure_from_option(option, measure)

    ! open output files
    !
    if (option%out_dis)  &
      call init_sink_file_rows(sinks(TRJ_DIS),  output%disfile,    &
                               measure_count(measure, TRJ_DIS))
    if (option%out_ang)  &
      call init_sink_file_rows(sinks(TRJ_ANG),  output%angfile,    &
                               measure_count(measure, TRJ_ANG))
    if (option%out_tor)  &
      call init_sink_file_rows(sinks(TRJ_TOR),  output%torfile,    &
                               measure_count(measure, TRJ_TOR))
    if (option%out_cdis) &
      call init_sink_file_rows(sinks(TRJ_CDIS), output%comdisfile, &
                               measure_count(measure, TRJ_CDIS))
    if (option%out_cang) &
      call init_sink_file_rows(sinks(TRJ_CANG), output%comangfile, &
                               measure_count(measure, TRJ_CANG))
    if (option%out_ctor) &
      call init_sink_file_rows(sinks(TRJ_CTOR), output%comtorfile, &
                               measure_count(measure, TRJ_CTOR))

    ! analysis loop
    !
    call init_source_file(source, trj_list, size(molecule%mass))
    call analyze_trj_unified(source, molecule%mass, measure, sinks, nstru)
    call finalize_source(source)

    ! close output files
    !
    do k = TRJ_NSTREAM, 1, -1
      call finalize_sink(sinks(k))
    end do

    ! Output summary
    !
    call print_output_info(output, option)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_trj_unified
  !> @brief        analysis loop shared by the CLI and the Python interface
  !! @authors      Claude Code
  !! @param[inout] source    : trajectory source (file, memory or lazy DCD)
  !! @param[in]    mass      : atomic masses (only used by COM measurements)
  !! @param[in]    measure   : measurement definitions
  !! @param[inout] sinks     : one result sink per stream (TRJ_DIS ... TRJ_CTOR);
  !!                           streams without measurements are skipped
  !! @param[out]   nstru_out : number of analyzed structures
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_trj_unified(source, mass, measure, sinks, nstru_out)

    ! formal arguments
    type(s_trj_source),  intent(inout) :: source
    real(wp),            intent(in)    :: mass(:)
    type(s_trj_measure), intent(in)    :: measure
    type(s_result_sink), intent(inout) :: sinks(TRJ_NSTREAM)
    integer,             intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory)       :: trajectory
    real(wp), allocatable    :: dis(:), ang(:), tor(:)
    real(wp), allocatable    :: cdis(:), cang(:), ctor(:)
    integer                  :: n(TRJ_NSTREAM)
    integer                  :: nstru, k, frame_status


    do k = 1, TRJ_NSTREAM
      n(k) = measure_count(measure, k)
    end do
    allocate(dis(n(TRJ_DIS)), ang(n(TRJ_ANG)), tor(n(TRJ_TOR)), &
             cdis(n(TRJ_CDIS)), cang(n(TRJ_CANG)), ctor(n(TRJ_CTOR)))

    nstru = 0

    do while (has_more_frames(source))

      call get_next_frame(source, trajectory, frame_status)
      if (frame_status /= 0) exit

      nstru = nstru + 1
      write(MsgOut,*) '      number of structures = ', nstru

      call measure_frame(trajectory%coord, mass, measure, &
                         dis, ang, tor, cdis, cang, ctor)

      if (n(TRJ_DIS)  > 0) call write_result_row(sinks(TRJ_DIS),  dis)
      if (n(TRJ_ANG)  > 0) call write_result_row(sinks(TRJ_ANG),  ang)
      if (n(TRJ_TOR)  > 0) call write_result_row(sinks(TRJ_TOR),  tor)
      if (n(TRJ_CDIS) > 0) call write_result_row(sinks(TRJ_CDIS), cdis)
      if (n(TRJ_CANG) > 0) call write_result_row(sinks(TRJ_CANG), cang)
      if (n(TRJ_CTOR) > 0) call write_result_row(sinks(TRJ_CTOR), ctor)

    end do

    nstru_out = nstru

    return

  end subroutine analyze_trj_unified

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      measure_count
  !> @brief        number of measurements of one stream
  !! @authors      Claude Code
  !! @param[in]    measure : measurement definitions
  !! @param[in]    stream  : TRJ_DIS ... TRJ_CTOR
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function measure_count(measure, stream) result(n)

    ! formal arguments
    type(s_trj_measure), intent(in) :: measure
    integer,             intent(in) :: stream

    ! return value
    integer :: n

    n = 0
    select case(stream)
    case(TRJ_DIS)
      if (allocated(measure%dist_num))      n = size(measure%dist_num)
    case(TRJ_ANG)
      if (allocated(measure%angl_list))     n = size(measure%angl_list, 2)
    case(TRJ_TOR)
      if (allocated(measure%tors_list))     n = size(measure%tors_list, 2)
    case(TRJ_CDIS)
      if (allocated(measure%cdis_pairs))    n = size(measure%cdis_pairs) / 2
    case(TRJ_CANG)
      if (allocated(measure%cang_triplets)) n = size(measure%cang_triplets) / 3
    case(TRJ_CTOR)
      if (allocated(measure%ctor_quads))    n = size(measure%ctor_quads) / 4
    end select

    return

  end function measure_count

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    measure_from_option
  !> @brief        build the measurement definitions from the CLI options
  !! @authors      Claude Code
  !! @param[in]    option  : option information
  !! @param[out]   measure : measurement definitions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine measure_from_option(option, measure)

    ! formal arguments
    type(s_option),      intent(in)  :: option
    type(s_trj_measure), intent(out) :: measure

    ! local variables
    integer,  allocatable    :: atoms(:), offsets(:)
    integer                  :: ngroup, natom, g, i, k


    if (option%out_dis) then
      measure%dist_list   = option%dist_list
      measure%dist_num    = option%dist_num
      measure%dist_weight = option%dist_weight
    end if

    if (option%out_ang) measure%angl_list = option%angl_list
    if (option%out_tor) measure%tors_list = option%tors_list

    ! COM measurements: flatten the selection groups into one atom table
    !
    if ((option%out_cdis .or. option%out_cang .or. option%out_ctor) .and. &
        allocated(option%selatoms)) then

      ngroup = size(option%selatoms)
      natom  = 0
      do g = 1, ngroup
        natom = natom + size(option%selatoms(g)%idx)
      end do

      allocate(atoms(natom), offsets(ngroup + 1))
      k = 0
      do g = 1, ngroup
        offsets(g) = k + 1
        do i = 1, size(option%selatoms(g)%idx)
          k = k + 1
          atoms(k) = option%selatoms(g)%idx(i)
        end do
      end do
      offsets(ngroup + 1) = k + 1

      if (option%out_cdis) then
        measure%cdis_atoms   = atoms
        measure%cdis_offsets = offsets
        measure%cdis_pairs   = reshape(option%cdist_group, &
                                       [size(option%cdist_group)])
      end if
      if (option%out_cang) then
        measure%cang_atoms    = atoms
        measure%cang_offsets  = offsets
        measure%cang_triplets = reshape(option%cangl_group, &
                                        [size(option%cangl_group)])
      end if
      if (option%out_ctor) then
        measure%ctor_atoms   = atoms
        measure%ctor_offsets = offsets
        measure%ctor_quads   = reshape(option%ctor_group, &
                                       [size(option%ctor_group)])
      end if

    end if

    return

  end subroutine measure_from_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    measure_frame
  !> @brief        compute every measurement for one frame
  !! @authors      NT, TM, SI, Claude Code
  !! @param[in]    coord   : coordinates (3, natom)
  !! @param[in]    mass    : atomic masses (COM measurements only)
  !! @param[in]    measure : measurement definitions
  !! @param[inout] dis, ang, tor, cdis, cang, ctor : results of each stream
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine measure_frame(coord, mass, measure, dis, ang, tor, cdis, cang, ctor)

    ! formal arguments
    real(wp),            intent(in)    :: coord(:,:)
    real(wp),            intent(in)    :: mass(:)
    type(s_trj_measure), intent(in)    :: measure
    real(wp),            intent(inout) :: dis(:), ang(:), tor(:)
    real(wp),            intent(inout) :: cdis(:), cang(:), ctor(:)

    ! local variables
    real(wp)                 :: c(3,4)
    integer                  :: i, j, g


    ! distances (weighted sums of atom pairs)
    do i = 1, size(dis)
      dis(i) = 0.0_wp
      do j = 1, measure%dist_num(i)/2
        dis(i) = dis(i) + measure%dist_weight(i,j) &
               * compute_dis(coord(:,measure%dist_list(2*j-1,i)), &
                             coord(:,measure%dist_list(2*j,  i)))
      end do
    end do

    ! angles
    do i = 1, size(ang)
      ang(i) = compute_ang(coord(:,measure%angl_list(1,i)), &
                           coord(:,measure%angl_list(2,i)), &
                           coord(:,measure%angl_list(3,i)))
    end do

    ! torsions
    do i = 1, size(tor)
      tor(i) = compute_dih(coord(:,measure%tors_list(1,i)), &
                           coord(:,measure%tors_list(2,i)), &
                           coord(:,measure%tors_list(3,i)), &
                           coord(:,measure%tors_list(4,i)))
    end do

    ! COM distances
    do i = 1, size(cdis)
      do j = 1, 2
        g = measure%cdis_pairs(2*(i-1) + j)
        c(:,j) = compute_com(coord, mass, measure%cdis_atoms( &
                   measure%cdis_offsets(g):measure%cdis_offsets(g+1)-1))
      end do
      cdis(i) = compute_dis(c(:,1), c(:,2))
    end do

    ! COM angles
    do i = 1, size(cang)
      do j = 1, 3
        g = measure%cang_triplets(3*(i-1) + j)
        c(:,j) = compute_com(coord, mass, measure%cang_atoms( &
                   measure%cang_offsets(g):measure%cang_offsets(g+1)-1))
      end do
      cang(i) = compute_ang(c(:,1), c(:,2), c(:,3))
    end do

    ! COM torsions
    do i = 1, size(ctor)
      do j = 1, 4
        g = measure%ctor_quads(4*(i-1) + j)
        c(:,j) = compute_com(coord, mass, measure%ctor_atoms( &
                   measure%ctor_offsets(g):measure%ctor_offsets(g+1)-1))
      end do
      ctor(i) = compute_dih(c(:,1), c(:,2), c(:,3), c(:,4))
    end do

    return

  end subroutine measure_frame

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_output_info
  !> @brief        print detailed output information
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_output_info(output, option)

    ! formal arguments
    type(s_output),          intent(in) :: output
    type(s_option),          intent(in) :: option

    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''

    if (option%out_dis)  then
      write(MsgOut,'(A)') '  [disfile] ' // trim(output%disfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Distance (angstrom)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Distance for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_ang)  then
      write(MsgOut,'(A)') '  [angfile] ' // trim(output%angfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Angle (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_tor)  then
      write(MsgOut,'(A)') '  [torfile] ' // trim(output%torfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Torsion angle (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Torsion angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_cdis)  then
      write(MsgOut,'(A)') '  [comdisfile] ' // trim(output%comdisfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Distance between the centers of mass (angstrom)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Distance for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_cang)  then
      write(MsgOut,'(A)') '  [comangfile] ' // trim(output%comangfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Angle between the centers of mass (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_ctor)  then
      write(MsgOut,'(A)') '  [comtorfile] ' // trim(output%comtorfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Torsion angle between the centers of mass (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Torsion angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine print_output_info

end module ta_analyze_mod
