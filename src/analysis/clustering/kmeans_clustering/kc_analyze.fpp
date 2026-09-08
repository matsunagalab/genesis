!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   kc_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module kc_analyze_mod

  use kc_option_str_mod
  use fitting_mod
  use trj_source_mod
  use fileio_trj_mod
  use fitting_str_mod
  use trajectory_str_mod
  use input_str_mod
  use output_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use string_mod
  use random_mod
  use select_atoms_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_kmeans_unified
  private :: assign_mass
  private :: get_replicate_name1

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM
  !! @param[inout] molecule   : molecule information
  !! @param[in]    input      : input information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, input, trj_list, trajectory, fitting, option, output)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_input),            intent(in)    :: input
    type(s_trj_list), target, intent(inout) :: trj_list
    type(s_trajectory),       intent(inout) :: trajectory
    type(s_fitting),          intent(inout) :: fitting
    type(s_option),           intent(inout) :: option
    type(s_output),           intent(inout) :: output

    ! local variables
    type(s_trj_source)            :: source
    integer,          allocatable :: cluster_index(:)
    integer,          allocatable :: center_index(:)
    type(s_pdb),      allocatable :: center_pdb(:)
    type(s_trj_file), allocatable :: trj_out(:)
    integer                       :: nclst, iclst, istru, idx_out


    if (option%check_only) &
      return

    nclst = option%num_clusters

    ! trajectory files of the clusters are written during the last pass
    !
    if (output%trjfile /= '') then
      allocate(trj_out(nclst))
      do iclst = 1, nclst
        call open_trj(trj_out(iclst),                             &
                      get_replicate_name1(output%trjfile, iclst), &
                      option%trjout_format, &
                      option%trjout_type,   &
                      IOFileOutputNew)
      end do
    end if

    ! clustering (shared with the Python interface)
    !
    call init_source_file(source, trj_list, molecule%num_atoms)
    call analyze_kmeans_unified(molecule, input, source, fitting, option, &
                                cluster_index, center_index,             &
                                output%pdbfile /= '', center_pdb, trj_out)
    call finalize_source(source)

    ! cluster index of every structure
    !
    if (output%indexfile /= '') then
      call open_file(idx_out, output%indexfile, IOFileOutputNew)
      do istru = 1, size(cluster_index)
        write(idx_out,'(I10,1X,I10)') istru, cluster_index(istru)
      end do
      call close_file(idx_out)
    end if

    ! PDB files of the cluster centers
    !
    if (output%pdbfile /= '') then
      write(MsgOut,'(A)') 'Analyze> output PDB files of the cluster centers'
      write(MsgOut,'(A)') ''
      do iclst = 1, nclst
        write(MsgOut,'(A,I10,A)') '   structure = ', center_index(iclst), '  >  ' // &
                                   trim(get_replicate_name1(output%pdbfile,iclst))
        call output_pdb(get_replicate_name1(output%pdbfile,iclst), center_pdb(iclst))
        call dealloc_pdb_all(center_pdb(iclst))
      end do
    end if

    write(MsgOut,'(A)') ''

    if (output%trjfile /= '') then
      do iclst = 1, nclst
        call close_trj (trj_out(iclst))
      end do
      deallocate(trj_out)
    end if

    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [indexfile] ' // trim(output%indexfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Index of the cluster to which the structure belongs'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_kmeans_unified
  !> @brief        k-means clustering (shared by the CLI and the Python
  !!               interface)
  !! @authors      TM, Claude Code
  !! @param[inout] molecule        : molecule information (coordinates are
  !!                                 overwritten during the last pass)
  !! @param[in]    input           : input information (initial index file)
  !! @param[inout] source          : trajectory source, read several times
  !! @param[inout] fitting         : fitting information
  !! @param[inout] option          : option information
  !! @param[out]   cluster_index   : cluster of every analyzed structure
  !! @param[out]   center_index    : structure number of every cluster center
  !! @param[in]    want_center_pdb : export the cluster centers into center_pdb
  !! @param[out]   center_pdb      : (num_clusters) PDB data of the centers
  !!                                 (allocated when want_center_pdb)
  !! @param[inout] trj_out         : if allocated, every structure is written,
  !!                                 fitted to its center, to trj_out(cluster)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_kmeans_unified(molecule, input, source, fitting, option, &
                                    cluster_index, center_index,             &
                                    want_center_pdb, center_pdb, trj_out)

    ! formal arguments
    type(s_molecule),              intent(inout) :: molecule
    type(s_input),                 intent(in)    :: input
    type(s_trj_source),            intent(inout) :: source
    type(s_fitting),               intent(inout) :: fitting
    type(s_option),                intent(inout) :: option
    integer,          allocatable, intent(out)   :: cluster_index(:)
    integer,          allocatable, intent(out)   :: center_index(:)
    logical,                       intent(in)    :: want_center_pdb
    type(s_pdb),      allocatable, intent(out)   :: center_pdb(:)
    type(s_trj_file), allocatable, intent(inout) :: trj_out(:)

    ! local variables
    type(s_trajectory)       :: trajectory
    real(wp)                 :: rmsd, min_rmsd, convergency, min_convergency
    real(wp)                 :: diff_coord(3)
    integer                  :: i, j, k, idx, nclst, iclst, iseed
    integer                  :: idx_in
    integer                  :: natom, nstru
    integer                  :: iatom, iiter, istru, frame_status
    integer                  :: alloc_stat
    logical                  :: converged
    character(MaxLine)       :: linein

    real(wp),         allocatable :: av0_coord_tmp(:,:), av0_coord_tmp2(:,:)
    real(wp),         allocatable :: ave_coord_tmp(:,:)
    real(wp),         allocatable :: av0_coord(:,:,:)
    real(wp),         allocatable :: ave_coord(:,:,:)
    real(wp),         allocatable :: cnt_coord(:,:,:)
    real(wp),         allocatable :: trj_coord(:,:)
    real(wp),         allocatable :: sqrt_mass(:)
    real(wp),         allocatable :: min_rmsd_clst(:)
    real(wp),         allocatable :: sum_rmsd(:)
    integer,          allocatable :: diff1(:), diff2(:)
    integer,          allocatable :: cluster_index_old(:)
    integer,          allocatable :: ndata(:)
    logical,          allocatable :: init_cluster(:)


    natom   = molecule%num_atoms
    nclst   = option%num_clusters
    iseed   = option%iseed

    allocate(sqrt_mass(natom),          &
             av0_coord_tmp(3,natom),    &
             ave_coord_tmp(3,natom),    &
             av0_coord_tmp2(3,natom),   &
             av0_coord (nclst,3,natom), &
             ave_coord (nclst,3,natom), &
             cnt_coord (nclst,3,natom), &
             trj_coord(3,natom),        &
             init_cluster(nclst),       &
             sum_rmsd(nclst),           &
             center_index(nclst),       &
             min_rmsd_clst(nclst),      &
             diff1(nclst),              &
             diff2(nclst),              &
             ndata(nclst), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! check mass
    !
    if (fitting%mass_weight) then
      call assign_mass(molecule)
    end if

    if (fitting%mass_weight) then
      do iatom = 1, natom
        sqrt_mass(iatom) = sqrt(molecule%mass(iatom))
      end do
    else
      do iatom = 1, natom
        sqrt_mass(iatom) = 1.0_wp
      end do
    end if

    ! number of structures to be analyzed
    !
    nstru = get_total_frames(source)
    allocate(cluster_index(nstru),cluster_index_old(nstru))
    cluster_index(1:nstru) = 0
    center_index(1:nclst) = 0

    ! initial cluster index
    !
    if (input%indexfile /= '') then
      call open_file(idx_in, input%indexfile, IOFileInput)
      do while (.true.)
        read (idx_in,'(A)',end=10) linein
        read (linein,'(I10)') istru
        if (istru <= nstru) then
          read (linein,'(I10,1x,I10)') istru, cluster_index(istru)
        end if
      end do
10    continue
      call close_file(idx_in)

      do i = 1, nstru
        if (cluster_index(i) <= 0 .or. cluster_index(i) > nclst) then
          iclst = int(random_get_legacy(iseed)*nclst) + 1
          if (iclst <= 0    ) iclst = 1
          if (iclst >  nclst) iclst = nclst
          cluster_index(i) = iclst
        end if
      end do
    else
      do i = 1, nstru
        iclst = int(random_get_legacy(iseed)*nclst) + 1
        if (iclst <= 0    ) iclst = 1
        if (iclst >  nclst) iclst = nclst
        cluster_index(i) = iclst
      end do
    end if

    init_cluster(1:nclst) = .false.

    write(MsgOut,'(A)') 'Analyze> initial cluster index'
    do i = 1, nstru
      write(MsgOut,'(5x,I10,1x,I10)') i, cluster_index(i)
    end do
    write(MsgOut,'(A)') ' '

    ! k-means iterations
    !
    do iiter = 1, option%max_iteration

      write(MsgOut,'(A,i10)') 'Analyze> k-means iteration = ', iiter

      cluster_index_old(1:nstru) = cluster_index(1:nstru)

      ! average structure of every cluster
      !
      do k = 1, option%num_iterations

        do j = 1, nclst
          do iatom = 1, natom
            ave_coord(j,1:3,iatom) = 0.0_wp
          end do
        end do

        istru = 0
        ndata(1:nclst) = 0
        call reset_source(source)
        do while (has_more_frames(source))
          call get_next_frame(source, trajectory, frame_status)
          if (frame_status /= 0) exit
          istru = istru + 1
          iclst = cluster_index(istru)
          if (.not. init_cluster(iclst)) then
            do iatom = 1, natom
              av0_coord(iclst,1:3,iatom) = trajectory%coord(1:3,iatom)
              cnt_coord(iclst,1:3,iatom) = trajectory%coord(1:3,iatom)
            end do
            center_index(iclst)  = istru
            init_cluster(iclst) = .true.
          end if
          do iatom = 1, natom
            av0_coord_tmp(1:3,iatom) = av0_coord(iclst,1:3,iatom)
          end do
          call run_fitting(fitting, av0_coord_tmp, trajectory%coord, trajectory%coord)
          do iatom = 1, natom
            av0_coord_tmp2(1:3,iatom) = av0_coord(iclst,1:3,iatom) *sqrt_mass(iatom)
            trj_coord(1:3,iatom)      = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
          end do
          call run_fitting(fitting, av0_coord_tmp2, trj_coord, trj_coord)
          ndata(iclst) = ndata(iclst) + 1
          do iatom = 1, natom
            ave_coord(iclst,1:3,iatom) = ave_coord(iclst,1:3,iatom) + trj_coord(1:3,iatom)
          end do
        end do

        do iclst = 1, nclst
          do iatom = 1, natom
            ave_coord(iclst,1:3,iatom) = ave_coord(iclst,1:3,iatom) / real(ndata(iclst), wp)
          end do
          do iatom = 1, natom
            av0_coord_tmp2(1:3,iatom) = av0_coord(iclst,1:3,iatom)*sqrt_mass(iatom)
            ave_coord_tmp (1:3,iatom) = ave_coord(iclst,1:3,iatom)
          end do
          call run_fitting(fitting, av0_coord_tmp2, ave_coord_tmp, ave_coord_tmp)
          do iatom = 1, natom
            av0_coord(iclst,1:3,iatom) = ave_coord_tmp(1:3,iatom)
            av0_coord(iclst,1:3,iatom) = av0_coord(iclst,1:3,iatom)/sqrt_mass(iatom)
          end do
        end do

      end do

      ! assign every structure to the closest cluster
      !
      istru                  = 0
      sum_rmsd     (1:nclst) = 0.0_wp
      ndata        (1:nclst) = 0
      min_rmsd_clst(1:nclst) = 999999999.9_wp

      call reset_source(source)
      do while (has_more_frames(source))
        call get_next_frame(source, trajectory, frame_status)
        if (frame_status /= 0) exit
        istru = istru + 1
        min_rmsd = 999999999.99_wp
        do iclst = 1, nclst
          do iatom = 1, natom
            av0_coord_tmp(1:3,iatom) = av0_coord(iclst,1:3,iatom)
          end do
          call run_fitting(fitting, av0_coord_tmp, trajectory%coord, trajectory%coord)
          do iatom = 1, natom
            av0_coord_tmp2(1:3,iatom) = av0_coord(iclst,1:3,iatom) *sqrt_mass(iatom)
            trj_coord(1:3,iatom)      = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
          end do
          call run_fitting(fitting, av0_coord_tmp2, trj_coord, trj_coord)
          rmsd = 0.0_wp
          do iatom = 1, size(option%analysis_atom%idx)
            idx = option%analysis_atom%idx(iatom)
            diff_coord(1:3) = av0_coord_tmp2(1:3, idx) - trj_coord(1:3, idx)
            rmsd = rmsd + dot_product(diff_coord, diff_coord)
          end do
          rmsd = sqrt(rmsd / real(size(option%analysis_atom%idx), wp))
          if (rmsd <= min_rmsd) then
            min_rmsd = rmsd
            cluster_index(istru) = iclst
          end if
        end do
        if (min_rmsd <= min_rmsd_clst(cluster_index(istru))) then
          min_rmsd_clst(cluster_index(istru)) = min_rmsd
          center_index(cluster_index(istru))  = istru
          iclst = cluster_index(istru)
          do iatom = 1, natom
            cnt_coord(iclst,1:3,iatom) = trajectory%coord(1:3,iatom)
          end do
        end if
        sum_rmsd(cluster_index(istru)) = sum_rmsd(cluster_index(istru)) + min_rmsd
        ndata   (cluster_index(istru)) = ndata   (cluster_index(istru)) + 1
      end do

      write(MsgOut,'(A)') '   Cluster index    Cluster center   # of structures    Cluster radius'
      do iclst = 1, nclst
        if (ndata(iclst) /= 0) then
          write(MsgOut,'(I16,I18,I18,F18.5)') iclst, center_index(iclst), ndata(iclst), &
                                              sum_rmsd(iclst) / real(ndata(iclst), wp)
        else
          write(MsgOut,'(I16,I18,I18,F18.5)') iclst, 0, ndata(iclst), 0.0_wp
        end if
      end do
      write(MsgOut,'(A)') ''

      do iclst = 1, nclst
        if (ndata(iclst) == 0) then
          write(MsgOut,'(A)') ' Warning: Empty cluster was generated'
          write(MsgOut,'(A)') ' Closest structure to the old center was selected as a new center'
          write(MsgOut,'(A)') ''
          do iatom = 1, natom
            av0_coord(iclst,1:3,iatom) = cnt_coord(iclst,1:3,iatom)
          end do
          cluster_index(center_index(iclst)) = iclst
        end if
      end do

      ! convergence check
      !
      if (iiter == 1) then
        diff1(1:nclst) = 0
        do istru = 1, nstru
          if (cluster_index_old(istru) /= cluster_index(istru)) then
            diff1(cluster_index(istru)) = diff1(cluster_index(istru)) + 1
          end if
        end do
        converged = .false.
      else if (iiter == 2) then
        diff2(1:nclst) = 0
        do istru = 1, nstru
          if (cluster_index_old(istru) /= cluster_index(istru)) then
            diff2(cluster_index(istru)) = diff2(cluster_index(istru)) + 1
          end if
        end do
        converged = .false.
      else
        diff1(1:nclst) = diff2(1:nclst)
        diff2(1:nclst) = 0
        do istru = 1, nstru
          if (cluster_index_old(istru) /= cluster_index(istru)) then
            diff2(cluster_index(istru)) = diff2(cluster_index(istru)) + 1
          end if
        end do
        min_convergency = 999999999
        do iclst = 1, nclst
          convergency = 100.0_wp - 100.0_wp * abs(diff2(iclst) - diff1(iclst))/ndata(iclst)
          if (convergency < min_convergency) then
            min_convergency = convergency
          end if
        end do
        converged = .true.
        if (min_convergency < option%stop_threshold) then
          converged = .false.
        end if
        write(MsgOut,'(A,F10.5,A)') '   Convergency = ', min_convergency,' %'
        write(MsgOut,'(A)') ''
      end if

      if (converged) exit

    end do

    ! last pass: cluster centers and (optionally) fitted trajectories
    !
    if (want_center_pdb .or. allocated(trj_out)) then

      if (want_center_pdb) allocate(center_pdb(nclst))

      istru = 0
      call reset_source(source)
      do while (has_more_frames(source))
        call get_next_frame(source, trajectory, frame_status)
        if (frame_status /= 0) exit
        istru = istru + 1
        do iclst = 1, nclst
          if (want_center_pdb) then
            if (istru == center_index(iclst)) then
              molecule%atom_coord(1:3,1:natom) = cnt_coord(iclst,1:3,1:natom)
              call run_fitting(fitting,             &
                               molecule%atom_coord, &
                               trajectory%coord,    &
                               trajectory%coord)
              molecule%atom_coord(1:3,1:natom) = trajectory%coord(1:3,1:natom)
              call export_molecules(molecule, option%trjout_atom, center_pdb(iclst))
            end if
          end if
          if (allocated(trj_out)) then
            if (iclst == cluster_index(istru)) then
              molecule%atom_coord(1:3,1:natom) = cnt_coord(iclst,1:3,1:natom)
              call run_fitting(fitting,             &
                               molecule%atom_coord, &
                               trajectory%coord,    &
                               trajectory%coord)
              call write_trj(trj_out(iclst), trajectory, option%trjout_atom, molecule)
            end if
          end if
        end do
      end do

    end if

    deallocate(sqrt_mass, av0_coord_tmp, ave_coord_tmp, av0_coord_tmp2, &
               av0_coord, ave_coord, cnt_coord, trj_coord,              &
               init_cluster, sum_rmsd, min_rmsd_clst,                   &
               diff1, diff2, ndata, cluster_index_old)

    return

  end subroutine analyze_kmeans_unified

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_mass
  !> @brief        assign mass
  !! @authors      NT, TM
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_mass(molecule)

    ! parameters
    real(wp),                parameter     :: MassH   =  1.008000_wp
    real(wp),                parameter     :: MassC   = 12.011000_wp
    real(wp),                parameter     :: MassN   = 14.007000_wp
    real(wp),                parameter     :: MassO   = 15.999000_wp
    real(wp),                parameter     :: MassS   = 32.060000_wp
    real(wp),                parameter     :: MassP   = 30.974000_wp
    real(wp),                parameter     :: MassMG  = 24.305000_wp
    real(wp),                parameter     :: MassZN  = 65.370000_wp

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i


    write(MsgOut,'(a)'),      'WARNING: atom mass is not assigned.'
    write(MsgOut,'(a)'),      '   uses default mass.'
    write(MsgOut,'(a,f9.6)')  '      1) H   : ', MassH
    write(MsgOut,'(a,f9.6)')  '      2) C   : ', MassC
    write(MsgOut,'(a,f9.6)')  '      3) N   : ', MassN
    write(MsgOut,'(a,f9.6)')  '      4) O   : ', MassO
    write(MsgOut,'(a,f9.6)')  '      5) S   : ', MassS
    write(MsgOut,'(a,f9.6)')  '      6) P   : ', MassP
    write(MsgOut,'(a,f9.6)')  '      7) MG  : ', MassMG
    write(MsgOut,'(a,f9.6)')  '      8) ZN  : ', MassZN


    do i = 1, molecule%num_atoms

      if (molecule%atom_name(i)(1:1) == 'H') then
        molecule%mass(i) = MassH
      else if (molecule%atom_name(i)(1:1) == 'C') then
        molecule%mass(i) = MassC
      else if (molecule%atom_name(i)(1:1) == 'N') then
        molecule%mass(i) = MassN
      else if (molecule%atom_name(i)(1:1) == 'O') then
        molecule%mass(i) = MassO
      else if (molecule%atom_name(i)(1:1) == 'S') then
        molecule%mass(i) = MassS
      else if (molecule%atom_name(i)(1:1) == 'P') then
        molecule%mass(i) = MassP
      else if (molecule%atom_name(i)(1:2) == 'MG') then
        molecule%mass(i) = MassMG
      else if (molecule%atom_name(i)(1:2) == 'ZN') then
        molecule%mass(i) = MassZN
      else
        write(MsgOut,'(a,a)') 'Assign_Mass> Unknown atom :', &
             molecule%atom_name(i)
      end if

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine assign_mass

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
         filename(:bl-1),no,filename(br+1:)

    return

  end function get_replicate_name1

end module kc_analyze_mod
