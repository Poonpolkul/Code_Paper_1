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
    
    ! number of years the household retires
    integer, parameter :: JR = 10

    ! number of white noise (zeta) shocks
    integer, parameter :: NW = 1 !7

    ! number of rate of return (vtheta) shocks 
    ![in our paper = aggregate productivity shock]
    integer, parameter :: NR = 7
    
    ! number of eta shocks
    integer, parameter :: NE = 1 !5

    ! number of points on the risky asset grid
    integer, parameter :: NK = 10
    
    ! number of points on the risk-free asset grid
    integer, parameter :: NB = 20

    ! household preference parameters 
    ! sigma and eta parameter value from Heer & Maussner 2012. Neeed gamma's value
    real*8, parameter :: gamma = 0.1 !0.5d0 ! entropic risk aversion in Value function
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta  = 0.98d0**5

    ! household risk process
    real*8, parameter :: sigma_zeta   = 5d0*0.0738d0
    real*8, parameter :: sigma_eps    = sqrt(20d0)*0.007d0 ! from HM
    real*8, parameter :: sigma_vtheta = 0.007d0*sqrt(20d0) !from HM !sqrt(5d0)*(0.157d0**2d0) !0.00016d0!HM!#####################
    real*8, parameter :: rho         = 0.95d0 !from HM #####################
        
    ! size of the risky asset grid
    real*8, parameter :: k_l    = 0.0d0
    real*8, parameter :: k_u    = 40d0
    real*8, parameter :: k_grow = 0.05d0
    
    ! size of the risk-free asset grid
    real*8, parameter :: b_l    = -40d0 !0d0
    real*8, parameter :: b_u    = 40d0
    real*8, parameter :: b_grow = 0.05d0

    ! initial risk premium
    real*8, parameter :: mu_r = 0.2d0 !discrete ##########

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW)
    real*8 :: pi_eta(NE, NE), eta(NE)
    real*8 :: pi_TProd(NR, NR), TProd(NR)
    integer :: iq_initial = 1, iv_initial = 1
!~     integer :: iq_initial = (NE+1)/2, iv_initial = (NR+1)/2

    ! production parameters
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 1d0-(1d0-0.0823d0)**5d0 !delta = 1d0-(1d0-0.0823d0)
    real*8, parameter :: TProd_bar = 1.60d0 ! Need shock process
    
    ! demographic parameters
    real*8, parameter :: n_p   = 0.01d0

    ! simulation parameters
    real*8, parameter :: damp    = 0.30d0
    real*8, parameter :: sig     = 1d-4
    integer, parameter :: itermax = 50

    ! counter variables
    integer :: iter

    ! macroeconomic variables
    real*8 :: rb, rk(NR)
    real*8 :: KK, BB, LL
    real*8 :: YY, CC, II, INC
    
    ! wages, transfer payments (old-age) and survival probabilities
    real*8 :: w(NR), wn(NR), eff(JJ), psi(JJ+1)

    ! demographic and other model parameters
    real*8 :: m(JJ)
    
    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! government variables
    real*8 :: tauw
    real*8 :: pen(JJ, NR), taxrev
    real*8 :: total_pen, total_INC

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), l_coh(JJ), v_coh(JJ)
    real*8 :: k_coh(JJ), b_coh(JJ)

    ! individual variables
    real*8 :: k(0:NK), kplus(JJ, 0:NK, 0:NB, NE, NW, NR), kplus_opt
    real*8 :: b(0:NB), bplus(JJ, 0:NK, 0:NB, NE, NW, NR), bplus_opt
    real*8 :: c(JJ, 0:NK, 0:NB, NE, NW, NR)
    real*8 :: phi(JJ, 0:NK, 0:NB, NE, NW, NR), V(JJ, 0:NK, 0:NB, NE, NW, NR) = 0d0
    real*8 :: phij_eta(JJ, NR), phij_Tprod(JJ, NR)
    real*8 :: a_bor(JJ)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NK, 0:NB, NE, NW, NR), EV(JJ, 0:NK, 0:NB, NE, NW, NR)
    real*8 :: RHS_b(JJ, 0:NK, 0:NB, NE, NW, NR), RHS_k(JJ, 0:NK, 0:NB, NE, NW, NR)
    integer :: i, ij_com, ik_com, ib_com, iq_com, ig_com, iv_com
    real*8 :: k_in, b_in
    real*8 :: cons_com, DIFF, INC_init

contains

    
    function FOC(a_in)
        
        implicit none
        real*8, intent(in) :: a_in(:)
        real*8 :: foc(size(a_in,1))
        real*8 :: focb, fock, l_in, varphik, varphib
        real*8 :: tomorrow_k, tomorrow_b, wage
        real*8 :: available, Root, Earnings, R_port
        integer :: ikl, ikr, ibl, ibr
        logical :: check
        
        kplus_opt = a_in(1)
        bplus_opt = a_in(2)
        
        !initialize root value
        root = 1d0 !?????
        
        ! calculate effective labour earnings
        Earnings  = wn(iv_com)*eff(ij_com)*exp(eta(iq_com) + zeta(ig_com)) &
                    + pen(ij_com, iv_com)

        ! calculate return asset return
        R_port = (1d0+rk(iv_com)-delta)*k(ik_com) + (1d0+rb)*b(ib_com)

        ! calculate available resources
        available = earnings + R_port

        ! calculate consumption
        cons_com = max((available - kplus_opt - bplus_opt) ,1d-10)

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(kplus_opt, k_l, k_u, k_grow, NK, ikl, ikr, varphik)
        call linint_Grow(bplus_opt, b_l, b_u, b_grow, NB, ibl, ibr, varphib)
        
        tomorrow_k = varphik*RHS_k(ij_com+1, ikl, ib_com, iq_com, ig_com, iv_com) + &
        (1d0-varphik)*RHS_k(ij_com+1, ikr, ib_com,  iq_com, ig_com, iv_com) 
        
        tomorrow_b = varphib*RHS_b(ij_com+1, ik_com, ibl, iq_com, ig_com, iv_com) + &
        (1d0-varphib)*RHS_b(ij_com+1, ik_com, ibr, iq_com, ig_com, iv_com)

        ! calculate the first order conditions for consumption
        foc(1) = margu(cons_com)**(-gamma) - tomorrow_k
        foc(2) = margu(cons_com)**(-gamma) - tomorrow_b

    end function

    ! calculates marginal utility of consumption
    function margu(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu

        margu = cons**(-1/gamma)

    end function
    
    ! calculates the value function
    function valuefunc(kplus, bplus, cons, ij, iq, ig, iv)

        implicit none
        integer, intent(in) :: ij, iq, ig, iv
        real*8, intent(in) :: kplus, bplus, cons
        real*8 :: valuefunc, varphik, varphib, c_help
        integer :: ikl, ikr, ibl, ibr

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! get tomorrows utility using bilinear interbolation
        call linint_Grow(kplus, k_l, k_u, k_grow, NK, ikl, ikr, varphik)
        call linint_Grow(bplus, b_l, b_u, b_grow, NB, ibl, ibr, varphib)

        ! calculate tomorrow's part of the value function 
        valuefunc = 1d0 
        if(ij < JJ)then
            valuefunc = max(varphik*varphib*EV(ij+1, ikl, ibl, iq, ig, iv) &
                       + (1d0-varphik)*varphib*EV(ij+1, ikr, ibl, iq, ig, iv) &
                       + varphik*(1d0-varphib)*EV(ij+1, ikl, ibr, iq, ig, iv) &
                       + (1d0-varphik)*(1d0-varphib)*EV(ij+1, ikr, ibr, iq, ig, iv) &
                       , 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = c_help**egam/egam + beta*psi(ij+1)*valuefunc

    end function

end module
