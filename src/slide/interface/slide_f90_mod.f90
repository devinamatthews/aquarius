! Copyright (c) 2013, Devin Matthews
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following
! conditions are met:
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!      * Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
! OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
! SUCH DAMAGE.

module slide

    implicit none

    integer, parameter :: SLIDE_OP_E = 0
    integer, parameter :: SLIDE_OP_C2Z = 1
    integer, parameter :: SLIDE_OP_C2Y = 2
    integer, parameter :: SLIDE_OP_C2X = 3
    integer, parameter :: SLIDE_OP_I = 4
    integer, parameter :: SLIDE_OP_SXY = 5
    integer, parameter :: SLIDE_OP_SXZ = 6
    integer, parameter :: SLIDE_OP_SYZ = 7

    integer, parameter :: SLIDE_GROUP_C1 = 0
    integer, parameter :: SLIDE_GROUP_C2 = 1
    integer, parameter :: SLIDE_GROUP_CS = 2
    integer, parameter :: SLIDE_GROUP_CI = 3
    integer, parameter :: SLIDE_GROUP_C2V = 4
    integer, parameter :: SLIDE_GROUP_D2 = 5
    integer, parameter :: SLIDE_GROUP_C2H = 6
    integer, parameter :: SLIDE_GROUP_D2H = 7

    integer, parameter :: SLIDE_IRREP_TOT_SYM = 0

    integer, parameter :: SLIDE_IRREP_C1_A = 0

    integer, parameter :: SLIDE_IRREP_C2_A = 0
    integer, parameter :: SLIDE_IRREP_C2_B = 1

    integer, parameter :: SLIDE_IRREP_CS_AP = 0
    integer, parameter :: SLIDE_IRREP_CS_APP = 1

    integer, parameter :: SLIDE_IRREP_CI_AG = 0
    integer, parameter :: SLIDE_IRREP_CI_AU = 1

    integer, parameter :: SLIDE_IRREP_C2V_A1 = 0
    integer, parameter :: SLIDE_IRREP_C2V_A2 = 1
    integer, parameter :: SLIDE_IRREP_C2V_B1 = 2
    integer, parameter :: SLIDE_IRREP_C2V_B2 = 3

    integer, parameter :: SLIDE_IRREP_D2_A = 0
    integer, parameter :: SLIDE_IRREP_D2_B1 = 1
    integer, parameter :: SLIDE_IRREP_D2_B2 = 2
    integer, parameter :: SLIDE_IRREP_D2_B3 = 3

    integer, parameter :: SLIDE_IRREP_C2H_AG = 0
    integer, parameter :: SLIDE_IRREP_C2H_AU = 1
    integer, parameter :: SLIDE_IRREP_C2H_BG = 2
    integer, parameter :: SLIDE_IRREP_C2H_BU = 3

    integer, parameter :: SLIDE_IRREP_D2H_AG = 0
    integer, parameter :: SLIDE_IRREP_D2H_B1G = 1
    integer, parameter :: SLIDE_IRREP_D2H_B2G = 2
    integer, parameter :: SLIDE_IRREP_D2H_B3G = 3
    integer, parameter :: SLIDE_IRREP_D2H_AU = 4
    integer, parameter :: SLIDE_IRREP_D2H_B1U = 5
    integer, parameter :: SLIDE_IRREP_D2H_B2U = 6
    integer, parameter :: SLIDE_IRREP_D2H_B3U = 7

    integer, parameter :: SLIDE_ORDER_ISCF = 0
    integer, parameter :: SLIDE_ORDER_ISFC = 1
    integer, parameter :: SLIDE_ORDER_SICF = 2
    integer, parameter :: SLIDE_ORDER_SIFC = 3
    integer, parameter :: SLIDE_ORDER_SCIF = 4
    integer, parameter :: SLIDE_ORDER_SCFI = 5
    integer, parameter :: SLIDE_ORDER_SFIC = 6
    integer, parameter :: SLIDE_ORDER_SFCI = 7

    integer, parameter :: SLIDE_SHELL_ORDER_FCD = 0
    integer, parameter :: SLIDE_SHELL_ORDER_FDC = 1
    integer, parameter :: SLIDE_SHELL_ORDER_CFD = 2
    integer, parameter :: SLIDE_SHELL_ORDER_CDF = 3
    integer, parameter :: SLIDE_SHELL_ORDER_DFC = 4
    integer, parameter :: SLIDE_SHELL_ORDER_DCF = 5

    integer, parameter :: SLIDE_MAX_L = 8

    interface

        subroutine SLIDE_element_symbol(element, symbol)

            integer*8 :: element
            character*2 :: symbol

        end subroutine SLIDE_element_symbol

        integer function SLIDE_element_nucleon(element)

            integer*8 :: element

        end function SLIDE_element_nucleon

        integer function SLIDE_element_spin(element)

            integer*8 :: element

        end function SLIDE_element_spin

        double precision function SLIDE_element_mass(element)

            integer*8 :: element

        end function SLIDE_element_mass

        double precision function SLIDE_element_charge(element)

            integer*8 :: element

        end function SLIDE_element_charge

        integer*8 function SLIDE_copy_element(element, charge)

            integer*8 :: element
            double precision :: charge

        end function SLIDE_copy_element

        subroutine SLIDE_free_element(element)

            integer*8 :: element

        end subroutine SLIDE_free_element

        integer*8 function SLIDE_get_element(symbol)

            character*2 :: symbol

        end function SLIDE_get_element

        integer function SLIDE_center_degeneracy(center)

            integer*8 :: center

        end function SLIDE_center_degeneracy

        integer function SLIDE_center_element(center)

            integer*8 :: center

        end function SLIDE_center_element

        integer function SLIDE_center_stabilizer(center)

            integer*8 :: center

        end function SLIDE_center_stabilizer

        subroutine SLIDE_center(center, degen, pos)

            integer :: degen
            integer*8 :: center
            double precision, dimension(3) :: pos

        end subroutine SLIDE_center

        subroutine SLIDE_center_after_op(center, op, pos)

            integer :: op
            integer*8 :: center
            double precision, dimension(3) :: pos

        end subroutine SLIDE_center_after_op

        integer*8 function SLIDE_new_center(pos, element)

            double precision, dimension(3) :: pos
            integer*8 :: element

        end function SLIDE_new_center

        subroutine SLIDE_free_center(center)

            integer*8 :: center

        end subroutine SLIDE_free_center

        integer function SLIDE_shell_L(shell)

            integer*8 :: shell

        end function SLIDE_shell_L

        integer function SLIDE_shell_nprim(shell)

            integer*8 :: shell

        end function SLIDE_shell_nprim

        integer function SLIDE_shell_ncontr(shell)

            integer*8 :: shell

        end function SLIDE_shell_ncontr

        integer function SLIDE_shell_nfunc(shell)

            integer*8 :: shell

        end function SLIDE_shell_nfunc

        subroutine SLIDE_shell_nfunc_irrep(shell, nfunc_irrep)

            integer*8 :: shell
            integer, dimension(8) :: nfunc_irrep

        end subroutine SLIDE_shell_nfunc_irrep

        integer function SLIDE_shell_irrep(shell, func, degen)

            integer*8 :: shell
            integer :: func, degen

        end function SLIDE_shell_irrep

        subroutine SLIDE_shell_irreps(shell, irreps)

            integer*8 :: shell
            integer, dimension(8,*) :: irreps

        end subroutine SLIDE_shell_irreps

        subroutine SLIDE_shell_idx(shell, idx)

            integer*8 :: shell
            integer, dimension(8) :: idx
            logical :: set

        end subroutine SLIDE_shell_idx

        subroutine SLIDE_shell_set_idx(shell, idx)

            integer*8 :: shell
            integer, dimension(8) :: idx
            logical :: set

        end subroutine SLIDE_shell_set_idx

        logical function SLIDE_shell_spherical(shell)

            integer*8 :: shell

        end function SLIDE_shell_spherical

        logical function SLIDE_shell_contaminants(shell)

            integer*8 :: shell

        end function SLIDE_shell_contaminants

        integer*8 function SLIDE_shell_pos(shell)

            integer*8 :: shell

        end function SLIDE_shell_pos

        integer*8 function SLIDE_new_shell(pos, L, nprim, ncontr, spherical, contaminants, &
                                           exponents, coefficients, idx)

            integer*8 :: pos
            integer :: L, nprim, ncontr
            logical :: spherical, contaminants
            double precision, dimension(nprim) :: exponents
            double precision, dimension(nprim,ncontr) :: coefficients
            integer, dimension(8) :: idx

        end function SLIDE_new_shell

        subroutine SLIDE_shell_ao_to_so(shell, ao_ordering, aoso, ld)

            integer*8 :: shell
            integer :: ao_ordering, ld
            double precision, dimension(ld,*) :: aoso

        end subroutine SLIDE_shell_ao_to_so

        integer*8 function SLIDE_copy_shell(shell, spherical, contaminants)

            integer*8 :: shell
            logical :: spherical, contaminants

        end function SLIDE_copy_shell

        subroutine SLIDE_free_shell(shell)

            integer*8 :: shell

        end subroutine SLIDE_free_shell

        integer*8 function SLIDE_context_num_integrals(context)

            integer*8 :: context

        end function SLIDE_context_num_integrals

        integer*8 function SLIDE_context_num_processed(context)

            integer*8 :: context

        end function SLIDE_context_num_processed

        subroutine SLIDE_context_set_num_processed(context, num_processed)

            integer*8 :: context, num_processed

        end subroutine SLIDE_context_set_num_processed

        integer*8 function SLIDE_new_context()
        end function SLIDE_new_context

        subroutine SLIDE_free_context(context)

            integer*8 :: context

        end subroutine SLIDE_free_context

        integer function SLIDE_set_cartesian_ordering(L, ordering)

            integer :: L, ordering(*)

        end function SLIDE_set_cartesian_ordering

        subroutine SLIDE_get_cartesian_ordering(L, ordering)

            integer :: L, ordering(*)

        end subroutine SLIDE_get_cartesian_ordering

        integer function SLIDE_set_spherical_ordering(L, ordering)

            integer :: L, ordering(*)

        end function SLIDE_set_spherical_ordering

        subroutine SLIDE_get_spherical_ordering(L, ordering)

            integer :: L, ordering(*)

        end subroutine SLIDE_get_spherical_ordering

        integer function SLIDE_set_group(group)

            integer :: group

        end function SLIDE_set_group

        integer function SLIDE_get_group()
        end function SLIDE_get_group

        subroutine SLIDE_set_ordering(ordering)

            integer :: ordering

        end subroutine SLIDE_set_ordering

        logical function SLIDE_get_ordering()
        end function SLIDE_get_ordering

        subroutine SLIDE_set_geometry_tolerance(geometry_tolerance)

            double precision :: geometry_tolerance

        end subroutine SLIDE_set_geometry_tolerance

        double precision function SLIDE_get_geometry_tolerance()
        end function SLIDE_get_geometry_tolerance

        integer function SLIDE_init()
        end function SLIDE_init

        integer function SLIDE_finish()
        end function SLIDE_finish

        integer*8 function SLIDE_process_2e_ints(context, nprocess, integrals, indices, cutoff)

            integer*8 :: context
            integer :: nprocess
            double precision :: cutoff
            double precision, dimension(nprocess) :: integrals
            integer*8, dimension(nprocess) :: indices

        end function SLIDE_process_2e_ints

        integer*8 function SLIDE_process_1e_ints(context, nprocess, integrals, indices, cutoff)

            integer*8 :: context
            integer :: nprocess
            double precision :: cutoff
            double precision, dimension(nprocess) :: integrals
            integer*4, dimension(nprocess) :: indices

        end function SLIDE_process_1e_ints

        integer function SLIDE_calc_eri(context, alpha, beta, a, b, c, d)

            double precision :: alpha, beta
            integer*8 :: context, a, b, c, d

        end function SLIDE_calc_eri

        integer function SLIDE_calc_ovi(context, alpha, beta, a, b)

            double precision :: alpha, beta
            integer*8 :: context, a, b

        end function SLIDE_calc_ovi

        integer function SLIDE_calc_kei(context, alpha, beta, a, b)

            double precision :: alpha, beta
            integer*8 :: context, a, b

        end function SLIDE_calc_kei

        integer function SLIDE_calc_nai(context, alpha, beta, a, b, centers, ncenters)

            double precision :: alpha, beta
            integer*8 :: context, a, b
            integer :: ncenters
            integer*8, dimension(ncenters) :: centers

        end function SLIDE_calc_nai

        integer function SLIDE_irrep(idx, ijkl)

            integer*8 :: idx
            integer :: ijkl

        end function SLIDE_irrep

        integer function SLIDE_orbital(idx, ijkl)

            integer*8 :: idx
            integer :: ijkl

        end function SLIDE_orbital

        integer function SLIDE_index_iscf(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_iscf

        integer function SLIDE_index_isfc(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_isfc

        integer function SLIDE_index_sicf(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_sicf

        integer function SLIDE_index_sifc(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_sifc

        integer function SLIDE_index_scif(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_scif

        integer function SLIDE_index_sfic(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_sfic

        integer function SLIDE_index_scfi(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_scfi

        integer function SLIDE_index_sfci(shell, func, contr, degen, irrep)

            integer :: func, contr, degen, irrep
            integer*8 :: shell

        end function SLIDE_index_sfci

        integer function SLIDE_func_cart(x, y, z)

            integer :: x, y, z

        end function SLIDE_func_cart

        integer function SLIDE_func_spher(n, l, m)

            integer :: n, l, m

        end function SLIDE_func_spher

        subroutine SLIDE_group_label(group, label)

            integer :: group
            character*3 :: label

        end subroutine SLIDE_group_label

        subroutine SLIDE_L_label(L, label)

            integer :: L
            character*1 :: label

        end subroutine SLIDE_L_label

        subroutine SLIDE_op_label(op, label)

            integer :: op
            character*3 :: label

        end subroutine SLIDE_op_label

        subroutine SLIDE_irrep_label(irrep, label)

            integer :: irrep
            character*3 :: label

        end subroutine SLIDE_irrep_label

        subroutine SLIDE_func_label(L, spherical, func, label)

            logical :: spherical
            integer :: L, func
            character*4 :: label

        end subroutine SLIDE_func_label

        subroutine SLIDE_shell_func_label(shell, func, label)

            integer*8 :: shell
            integer :: func
            character*4 :: label

        end subroutine SLIDE_shell_func_label

        subroutine SLIDE_stabilizer(center, stab)

            integer*8 :: center
            integer :: stab(8)

        end subroutine SLIDE_stabilizer

        integer function SLIDE_group_order()
        end function SLIDE_group_order

        subroutine SLIDE_shell_func_centers(shell, func, irrep, proj)

            integer*8 :: shell
            integer :: func, irrep
            integer, dimension(8) :: proj

        end subroutine SLIDE_shell_func_centers

        double precision function SLIDE_nuclear_repulsion(centers, ncenters)

            integer :: ncenters
            integer*8 :: centers(ncenters)

        end function SLIDE_nuclear_repulsion

    end interface

end module slide
