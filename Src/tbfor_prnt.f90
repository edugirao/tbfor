      MODULE tbfor_prnt
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
      IMPLICIT none
      CHARACTER*66:: exclamation_bar      
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE print1
        USE tbfor_var
        IMPLICIT none
        
        IF(my_rank.eq.0) WRITE(*,*) 'Self-consistency reached!!!'
        IF(my_rank.eq.0) WRITE(*,'(A)') exclamation_bar

        END SUBROUTINE print1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE print2
        USE tbfor_var
        IMPLICIT none

        IF(my_rank.eq.0) WRITE(*,'(A)') exclamation_bar        
        
        END SUBROUTINE print2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE print3(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER::i
        
        WRITE(*,*) 'Iteration step ',i
        WRITE(*,*) 'System= ',TRIM(label)
              
        END SUBROUTINE print3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE print4
        USE tbfor_var
        IMPLICIT none

        WRITE(*,*) 'Ef=',ef_scc

        END SUBROUTINE print4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE print5(xx)
        USE tbfor_var
        IMPLICIT none
        REAL(KIND=8):: xx

        WRITE(*,*)
        WRITE(*,*) 'N_el= ',xx

        END SUBROUTINE print5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE print6
        USE tbfor_var
        IMPLICIT none

        WRITE(*,*)
        WRITE(*,*) 'Maximal difference=',dn_max 
        WRITE(*,'(A)') exclamation_bar
        
        END SUBROUTINE print6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SUBROUTINE print
!        USE tbfor_var
!        IMPLICIT none
!
!        
!        END SUBROUTINE print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_prnt
