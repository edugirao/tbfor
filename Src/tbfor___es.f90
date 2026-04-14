      SUBROUTINE tbfor(local_comm_ent,inputfile_ent)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This file and the containing routines are part of the TBfor 
      ! (Tight-Binding: Fortran Operational Resource) Project for the 
      ! calculation of electronic structure of nanoscaled systems using
      ! a tight-binding +Hubbard approach. This program has been 
      ! written by Eduardo Costa Girao (Universidade Federal do Piaui) and
      ! Vincent Meunier (Rensselaer Polytechnic Institute) in Nov/21/2012
      ! and its use is limited to the authorization of the authors.
      ! Contact: edu@ufpi.edu.br ; meuniv@rpi.edu
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE tbfor_var
      USE tbfor_read
      USE tbfor_init
      USE tbfor_kvec      
      USE tbfor_rcov
      USE tbfor_task
      USE tbfor_self
      USE tbfor_post      
      IMPLICIT none
      INTEGER:: local_comm_ent,ierr
      CHARACTER*14 inputfile_ent
!      INCLUDE 'mpif.h'
      inputfile=inputfile_ent
      local_comm=local_comm_ent
      CALL print_version
      bcalc=.false.
      onlyham=.false.
      onlyjoin=.false. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading the input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL reading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Performing pre-calculation tasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL initials
      CALL k_space      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recovering previous unfinished run (if it is the case)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL recover
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Self-consistency cycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL self_consistency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop Along the path on the reciprocal space to calculate the bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Post-processing tasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL post_processing

      END SUBROUTINE tbfor