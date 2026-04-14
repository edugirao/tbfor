      MODULE tbfor_allo
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
      IMPLICIT none
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation1(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> ev_tmp,ec_tmp'
          STOP
        END IF      

        END SUBROUTINE allocation1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation2(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> S_up2,S_dw2'
          STOP
        END IF           

        END SUBROUTINE allocation2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation3(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> n_up_new,n_dw_new'
          STOP
        END IF

        END SUBROUTINE allocation3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation4(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> pn_up_new,pn_dw_new'
          STOP
        END IF        

        END SUBROUTINE allocation4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation5(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> Rpul,pop_up,pop_dw'
          STOP
        END IF          

        END SUBROUTINE allocation5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation6(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> n_up_k,n_dw_k'
          STOP
        END IF           

        END SUBROUTINE allocation6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation7(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> pn_up_k,pn_dw_k'
          STOP
        END IF  

        END SUBROUTINE allocation7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation8(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> ec,ev'
          STOP
        END IF         

        END SUBROUTINE allocation8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation9(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> n_pop_tmp_up,n_pop_tmp_dw'
          STOP
        END IF    

        END SUBROUTINE allocation9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation10(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> pn_pop_tmp_up,pn_pop_tmp_dw'
          STOP
        END IF   

        END SUBROUTINE allocation10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation11(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> Rmat,alpha,ipiv'
          STOP
        END IF   

        END SUBROUTINE allocation11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation12(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> dif'
          STOP
        END IF

        END SUBROUTINE allocation12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation13(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> E_tmp_up,E_tmp_dw'
          STOP
        END IF

        END SUBROUTINE allocation13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation14(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to allocate -> E_tmp,spin,pos'
          STOP
        END IF

        END SUBROUTINE allocation14
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation15(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation16(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation16
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation17(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation17
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation18(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation18
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation19(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation19
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation20(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation20
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation21(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation21
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation22(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation23(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation23
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation24(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation24
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation25(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation25
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation26(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation26
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation27(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation27
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation28(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation28
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation29(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation29
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE allocation30(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE allocation30
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation1(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0) PRINT*,'Failed to deallocate -> ec,ec_tmp,ev,ev_tmp'
        IF(i.ne.0) STOP        

        END SUBROUTINE deallocation1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation2(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> S_up2,S_dw2'
          STOP
        END IF

        END SUBROUTINE deallocation2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation3(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> n_up_new,n_dw_new'
          STOP
        END IF

        END SUBROUTINE deallocation3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation4(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> n_up_k,n_dw_k'
          STOP
        END IF

        END SUBROUTINE deallocation4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation5(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> Rpul,pop_up,pop_dw'
          STOP
        END IF

        END SUBROUTINE deallocation5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation6(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp,spin,pos'
          STOP
        END IF

        END SUBROUTINE deallocation6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation7(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> H_up,H_dw'
          STOP
        END IF

        END SUBROUTINE deallocation7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation8(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> S_up,S_dw'
          STOP
        END IF

        END SUBROUTINE deallocation8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation9(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> Rmat,alpha,ipiv'
          STOP
        END IF   

        END SUBROUTINE deallocation9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation10(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> n_pop_tmp_up,n_pop_tmp_dw'
          STOP
        END IF   

        END SUBROUTINE deallocation10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation11(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> pn_pop_tmp_up,pn_pop_tmp_dw'
          STOP
        END IF   

        END SUBROUTINE deallocation11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation12(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        IF(i.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp_up,E_tmp_dw'
          STOP
        END IF        

        END SUBROUTINE deallocation12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation13(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation14(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation14
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation15(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation16(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation16
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation17(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation17
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation18(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation18
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation19(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation19
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation20(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation20
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation21(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation21
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation22(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation23(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation23
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation24(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation24
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation25(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation25
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation26(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation26
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation27(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation27
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation28(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation28
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation29(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation29
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE deallocation30(i)
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i

        

        END SUBROUTINE deallocation30
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_allo
