!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
! #VC# VERSION: 1.0  (23 January 2018)
!
!##############################################################################
module globals

    use toolbox

    implicit none

    ! number of years the household lives
    integer, parameter :: JJ = 12
!~     integer, parameter :: JJ = 80
    
    ! number of years the household retires
    integer, parameter :: JR = 10
!~     integer, parameter :: JR = 40

    ! number of productivity type
    integer, parameter :: NP = 2 ! ip = 1,2

    ! number of technological shock process values
    integer, parameter :: NS = 5 ! is = 1,...,5

    ! number of points on the risky asset grid
    integer, parameter :: NK = 10
!~     integer, parameter :: NK = 100
    
    ! number of points on the risk-free asset grid
    integer, parameter :: NB = 20
!~     integer, parameter :: NB = 100

    ! household preference parameters 
    ! sigma and eta parameter value from Heer & Maussner 2012. Neeed gamma's value
    real*8, parameter :: gamma = 0.5 ! entropic risk aversion in Value function
    real*8, parameter :: beta  = 0.998
    real*8, parameter :: sigma = 0.5 !2d0
    real*8, parameter :: eta =  0.6 ! 7
    real*8, parameter :: phiu = 0.26 !for phi in utility. Value from HM(2012)

    ! household risk process
    real*8, parameter :: sigma_theta = 0.23d0
    
    ! transaction cost for changes in risk asset holdings
    real*8, parameter :: varphi1 = 0.06
    real*8, parameter :: varphi2 = 1.6 ! as in Cwik (2015), but too low
    
    ! production parameters
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 1d0-(1d0-0.0823d0)
    real*8, parameter :: Omega = 1.60d0 ! Need shock process
    
    ! aggregate risk process [just use parameter from Fehr's HH risk process. Need to adjust later]
    real*8, parameter :: rho         = 0.98d0
    real*8, parameter :: sigma_eps   = 0.05d0

    ! size of the risky asset grid
    real*8, parameter :: k_l    = 0.0d0
    real*8, parameter :: k_u    = 40d0
    real*8, parameter :: k_grow = 0.05d0
    
    ! size of the risk-free asset grid
    real*8, parameter :: b_l    = 0d0
    real*8, parameter :: b_u    = 40d0
    real*8, parameter :: b_grow = 0.05d0

    ! demographic parameters
    real*8, parameter :: n_p   = 0.01d0

    ! simulation parameters
    real*8, parameter :: damp    = 0.30d0
    real*8, parameter :: sig     = 1d-4
    integer, parameter :: itermax = 50

    ! counter variables
    integer :: iter

    ! macroeconomic variables
    real*8 :: rb, rbn, rk, rkn,  w, wn, p
    real*8 :: KK, AA, BB, LL, HH
    real*8 :: YY, CC, II, GG, INC

    ! government variables
    real*8 :: tauc, tauw, taurb, taurk, taup, kappa
    real*8 :: gy, by, pen(JJ), PP, taxrev(4)
    integer :: tax
    logical :: ageing_on

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), l_coh(JJ), k_coh(JJ), b_coh(JJ), v_coh(JJ)

    ! the shock process
    real*8 :: dist_theta(NP), theta(NP)
    real*8 :: pi(NS, NS), EtaShock(NS)
!~     real*8 :: eta(NS)
    integer :: is_initial = 3


    ! demographic and other model parameters
    real*8 :: m(JJ)
    real*8 :: eff(JJ)

    ! individual variables
    real*8 :: k(0:NK), kplus(JJ, 0:NK, 0:NB, NP, NS)
    real*8 :: b(0:NB), bplus(JJ, 0:NK, 0:NB, NP, NS)
    real*8 :: c(JJ, 0:NK, 0:NB, NP, NS), l(JJ, 0:NK, 0:NB, NP, NS)
    real*8 :: phi(JJ, 0:NK, 0:NB, NP, NS), V(JJ, 0:NK, 0:NB, NP, NS) = 0d0

    ! numerical variables
    real*8 :: RHS(JJ, 0:NK, 0:NB, NP, NS), EV(JJ, 0:NK, 0:NB, NP, NS)
    real*8 :: RHS_b(JJ, 0:NK, 0:NB, NP, NS), RHS_k(JJ, 0:NK, 0:NB, NP, NS)
    real*8 :: RHSN_b(JJ, 0:NK, 0:NB, NP, NS), RHSN_k(JJ, 0:NK, 0:NB, NP, NS)
    real*8 :: RHSD_b(JJ, 0:NK, 0:NB, NP, NS), RHSD_k(JJ, 0:NK, 0:NB, NP, NS)
    integer :: i, ij_com, ik_com, ib_com, ip_com, is_com, it_com
    real*8 :: k_in_next, b_in_next, k_temp, b_temp
    real*8 :: k_in, b_in
    real*8 :: cons_com, lab_com, DIFF, INC_init

contains

    ! calculates implicit function for labour
    function implicitl(l_in)
    
        implicit none
        
        real*8, intent(in) :: l_in
        real*8 :: lab_com, rk
        real*8 :: implicitl, wage
        
        ! calculate labour
        lab_com = l_in

        ! calculate the wage rate
        wage = wn*eff(ij_com)*theta(ip_com)
        
        ! calculate implicit value of labour
        implicitl = ((1d0/phiu)*(1d0-l_in)**eta*(1d0-tauw)*wage)**(1/sigma) &
        + bplus(ij_com, ik_com, ib_com, ip_com, is_com)+kplus(ij_com, ik_com, ib_com, ip_com, is_com) &
        - (1d0+rb)*b(ib_com) - ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) &
        -(1-tauw)*wage*l_in-varphi1*(abs(kplus(ij_com, ik_com, ib_com, ip_com, is_com)-k(ik_com)))**varphi2
        
    end function
    
!~     ! calculates implicit function for labour
!~     function implicitlb(l_in)
    
!~         implicit none
        
!~         real*8, intent(in) :: l_in
!~         real*8 :: lab_com, rk, bplus
!~         real*8 :: implicitlb, wage
        
!~         ! calculate labour
!~         lab_com = l_in
!~         bplus = b_in
 
!~         ! calculate the wage rate
!~         wage = wn*eff(ij_com)*theta(ip_com)
        
!~         ! calculate implicit value of labour
!~         implicitlb = ((1d0/phiu)*(1d0-l_in)**eta*(1d0-tauw)*wage)**(1/sigma) &
!~         + bplus + kplus(ij_com, ik_com, ib_com, ip_com, is_com) &
!~         - (1d0+rb)*b(ib_com) - ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) &
!~         -(1-tauw)*wage*l_in-varphi1*(abs(kplus(ij_com, ik_com, ib_com, ip_com, is_com)-k(ik_com)))**varphi2
        
!~     end function
    
!~     ! calculates implicit function for labour
!~     function implicitlk(l_in)
    
!~         implicit none
        
!~         real*8, intent(in) :: l_in
!~         real*8 :: lab_com, rk, kplus
!~         real*8 :: implicitlk, wage
        
!~         ! calculate labour
!~         lab_com = l_in
!~         kplus = k_in
 
!~         ! calculate the wage rate
!~         wage = wn*eff(ij_com)*theta(ip_com)
        
!~         ! calculate implicit value of labour
!~         implicitlk = ((1d0/phiu)*(1d0-l_in)**eta*(1d0-tauw)*wage)**(1/sigma) &
!~         + bplus(ij_com, ik_com, ib_com, ip_com, is_com) + kplus &
!~         - (1d0+rb)*b(ib_com) - ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) &
!~         -(1-tauw)*wage*l_in-varphi1*(abs(kplus-k(ik_com)))**varphi2
        
!~     end function

    ! calculated the first and second FOCs vary in k_plus (eq 20 & 21)
    function FOCK(k_in)
    
        implicit none
        real*8, intent(in) :: k_in
        real*8 :: fock, kplus, l_in, varphi, tomorrow_k, wage
        real*8 :: available
        integer :: ikl, ikr
        logical :: check

        ! calculate tomorrows assets
        kplus = k_in

        ! calculate the wage rate
        wage = wn*eff(ij_com)*theta(ip_com)

        ! calculate available resources
        available = ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) + (1d0+rbn)*b(ib_com) + pen(ij_com)

        ! determine labour
        if (ij_com < JR) then
            ! solve lab_com from the implicit FOC using root finding
            call fzero(l_in, implicitl, check)
            lab_com = l_in
        else
            lab_com = 0d0
        endif

        if (lab_com < 0) then
            lab_com = 0d0
        end if
        
        ! calculate consumption
        cons_com = max((available + (1-tauw)*wage*lab_com - kplus &
         - bplus(ij_com, ik_com, ib_com, ip_com, is_com)- varphi1*(abs(kplus-k(ik_com)))**varphi2), 1d-10)

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(kplus, k_l, k_u, k_grow, NK, ikl, ikr, varphi)

        tomorrow_k = varphi*RHS_k(ij_com+1, ikl, ib_com, ip_com, is_com) + &
                  (1d0-varphi)*RHS_k(ij_com+1, ikr, ib_com,  ip_com, is_com) 

        ! calculate the first order condition for consumption
        fock = abs(cons_com**(-sigma) - tomorrow_k)

!~         ! print for debug ***********************************************
!~         if (ij_com < JR) then
!~             print*,'ij=', ij_com, 'ik=', ik_com,'ib=', ib_com,'ip=', ip_com, &
!~             'is=', is_com, 'wage=', wage, 'kplus=', kplus, 'lab_com=', lab_com, 'cons_com=', &
!~              cons_com, 'tomorrow_k=', tomorrow_k, 'fock=', fock, 'bplus =', bplus(ij_com, ik_com, ib_com, ip_com, is_com)
!~         end if
        
      
    end function
    
    ! calculated the first and second FOCs vary in bplus (eq 20 & 21)
    function FOCB(b_in)
    
        implicit none
        real*8, intent(in) :: b_in
        real*8 :: focb, l_in, bplus, varphi, tomorrow_b, wage
        real*8 :: available
        integer :: ibl, ibr
        logical :: check

        ! calculate tomorrows assets
        bplus = b_in

        wage = wn*eff(ij_com)*theta(ip_com)
        
        ! calculate available resources
        available = ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) + (1d0+rbn)*b(ib_com) + pen(ij_com)

        ! determine labour
        if (ij_com < JR) then
            ! solve lab_com from the implicit FOC using root finding
            call fzero(l_in, implicitl, check)
            lab_com = l_in
        else
            lab_com = 0d0
        endif

        if (lab_com < 0) then
            lab_com = 0d0
        end if

        ! calculate consumption
        cons_com = max((available + (1-tauw)*wage*lab_com - kplus(ij_com, ik_com, ib_com, ip_com, is_com)  & 
        - bplus - varphi1*(abs(kplus(ij_com, ik_com, ib_com, ip_com, is_com)-k(ik_com)))**varphi2) ,1d-10)

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(bplus, b_l, b_u, b_grow, NB, ibl, ibr, varphi)

        tomorrow_b = varphi*RHS_b(ij_com+1, ik_com, ibl, ip_com, is_com) + &
                  (1d0-varphi)*RHS_b(ij_com+1, ik_com, ibr, ip_com, is_com)


        ! calculate the first order condition for consumption
        focb = abs(cons_com**(-sigma) - tomorrow_b)
        
        ! print for debug ***********************************************
!~         if (ij_com < JR) then
!~             print*,'ij=', ij_com, 'ik=', ik_com,'ib=', ib_com,'ip=', ip_com, &
!~             'is=', is_com, 'wage=', wage, 'bplus=', bplus, 'lab_com=', lab_com, 'cons_com=', &
!~              cons_com, 'tomorrow_b=', tomorrow_b, 'focb=', focb, 'kplus =', kplus(ij_com, ik_com, ib_com, ip_com, is_com)
!~         end if

    end function
    

    ! calculates the value function
    function valuefunc(kplus, bplus, cons, lab, ij, ip, is)

        implicit none
        integer, intent(in) :: ij, ip, is
        real*8, intent(in) :: kplus, bplus, cons, lab
        real*8 :: valuefunc, varphik, varphib, c_help, l_help, egam
        integer :: ikl, ikr, ibl, ibr

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)
        l_help = min(max(lab, 0d0),1d0-1d-10)

        ! get tomorrows utility using bilinear interbolation
        call linint_Grow(kplus, k_l, k_u, k_grow, NK, ikl, ikr, varphik)
        call linint_Grow(bplus, b_l, b_u, b_grow, NK, ibl, ibr, varphib)

        ! calculate tomorrow's part of the value function 
        valuefunc = 1d0 
        if(ij < JJ)then
            valuefunc = max(varphik*varphib*EV(ij+1, ikl, ibl, ip, is) &
                       + (1d0-varphik)*varphib*EV(ij+1, ikr, ibl, ip, is) &
                       + varphik*(1d0-varphib)*EV(ij+1, ikl, ibr, ip, is) &
                       + (1d0-varphik)*(1d0-varphib)*EV(ij+1, ikr, ibr, ip, is) &
                       , 1d-10)
        endif

        ! add todays part and discount
        valuefunc = (c_help**(1d0-sigma))/(1d0 - sigma) & 
        + phiu*((1d0-l_help)**(1d0-eta))/(1d0 - eta) & 
        - (beta/gamma)*dlog(valuefunc) 

    end function

end module
