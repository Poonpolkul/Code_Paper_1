!##############################################################################
! MODULE globals
!
! Paper 1 Option 7
! Build on the code of Hans Fehr and Fabian Kindermann
!
!##############################################################################
module globals

    use toolbox

    implicit none

    ! number of years the household lives
    integer, parameter :: JJ = 45
    
    ! number of years the household retires
    integer, parameter :: JR = 80

    ! number of white noise (zeta) shocks
    integer, parameter :: NW = 7

!~     ! number of transitory (epsilon) shocks
!~     integer, parameter :: NS = 7

    ! number of rate of return (vtheta) shocks 
    ![in our paper = aggregate productivity shock]
    integer, parameter :: NR = 7
    
!~     ! number of eps-vtheta shocks
!~     integer, parameter :: NSR = NS*NR
    
    ! number of eta shocks
    integer, parameter :: NE = 5

    ! number of points on the asset grid
    integer, parameter :: NA = 200

    ! number of points on the risky share grid
    integer, parameter :: NO = 200

!~     ! number of points on the cash-on-hand-grid
!~     integer, parameter :: NX = 200

    ! household preference parameters 
    real*8, parameter :: gamma = 0.10d0 
    real*8, parameter :: beta  = 0.96d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    
    ! household risk process
    real*8, parameter :: sigma_zeta   = 0.0738d0
    real*8, parameter :: sigma_eta    = 0.05d0
    real*8, parameter :: sigma_vtheta = 0.157d0**2d0
    real*8, parameter :: rho         = 0.98d0
    
    ! size of the asset grid
    real*8, parameter :: a_l    = -500.0d0
    real*8, parameter :: a_u    = 500d0
    real*8, parameter :: a_grow = 0.04d0

    ! size of risky share grid
    real*8, parameter :: omega_l    = 0.0d0
    real*8, parameter :: omega_u    = 1.0d0
    real*8, parameter :: omega_grow = 0.04d0

!~     ! lower and upper bound for the eta-grid
!~     real*8 :: eta_l(JJ), eta_u(JJ)

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW)
    real*8 :: dist_vtheta(NSR), vtheta(NSR)
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 3
    
!~     real*8 :: dist_epsvtheta(NSR), eps(NSR), vtheta(NSR)

!~     ! should cohort averages and variance be calculated analytically
!~     logical, parameter :: analytical = .true.
!~     !logical, parameter :: analytical = .false.     ! for quintiles

    ! production parameters
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 1d0-(1d0-0.0823d0)
    real*8, parameter :: Omega_bar = 1.60d0
    real*8 :: Omega(0:NR)
    
    ! demographic parameters 
    real*8, parameter :: n_p   = 0.01d0 

    ! simulation parameters
    real*8, parameter :: damp    = 0.30d0
    real*8, parameter :: sig     = 1d-4
    integer, parameter :: itermax = 50

    ! counter variables
    integer :: iter

    ! macroeconomic variables
    real*8 :: rb, rk(0:NR)
    real*8 :: AA, KK, BB, LL, HH
    real*8 :: YY, CC, II, INC
    
    ! wages, transfer payments (old-age) and survival probabilities
    real*8 :: w(0:NR), wn(0:NR), eff(JJ), pen(JJ, 0:NR), psi(JJ+1)

    ! demographic and other model parameters
    real*8 :: m(JJ)
    
    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! government variables
    real*8 :: tauw
    real*8 :: pen(JJ, 0:NR), taxrev

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), a_coh(JJ), omega_coh(JJ), l_coh(JJ), v_coh(JJ)
    real*8 :: k_coh(JJ), b_coh(JJ)
    real*8 :: cv_coh(JJ), yv_coh(JJ), av_coh(JJ), ov_coh(JJ), lv_coh(JJ), vv_coh(JJ)
    real*8 :: kv_coh(JJ), bv_coh(JJ)

    ! individual variables
    real*8 :: a(0:NA), omega(0:NO)
    real*8 :: a_bor(JJ)
    real*8 :: omega_plus(JJ, 0:NA, 0:NS), Q(JJ, 0:NA, 0:NS)
    real*8 :: c(JJ, 0:NA, 0:NO, 0:NS, 0:NW, 0:NSR), a_plus(JJ, 0:NA, 0:NO, 0:NS, 0:NW, 0:NSR)
    real*8 :: V(JJ, 0:NA, 0:NO, 0:NS, 0:NW, 0:NSR) = 0d0
    real*8 :: phi(JJ, 0:NA, 0:NO, 0:NS, 0:NW, 0:NSR)
    real*8 :: phi_a(0:NA), phi_aoe(0:NA, 0:NO, 0:NS)
    
    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, 0:NS)
    integer :: i, ij_com, ia_com, io_com, ig_com, iv_com, iq_com
    real*8 :: cons_com, lab_com, DIFF, INC_init ! INC_init??????

contains

    ! the first order condition regarding consumption
    function foc_cons(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_cons, a_plus, varphi, tomorrow, earnings
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus  = x_in

        ! calculate effective labour earnings
        earnings  = (1-tau)*w(iv)*eff(ij_com+1)*exp(eta(iq_com) + zeta(ig_com)) &
                    + pen(ij_com, iv_com)
                
        ! calculate return on next-period portfolio, R
        R_port = 1d0 + rb + omega(io_com)*(rk(iv_com) - rb)

        ! calculate current consumption 
        cons_com = R_port*a(ia_com) + earnings - a_plus
        
        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(a_plus, a_l, a_u, a_grow, &
                    NA, ial, iar, varphi)
        
        tomorrow = varphi*RHS(ij_com, ial, iq_com) + (1d0-varphi)*RHS(ij_com, iar, iq_com)
        
        ! calculate first order condition for consumption
        foc_cons = cons_com - tomorrow
    end function


!##############################################################################

    ! the first order condition regarding portfolio choice
    function foc_port(p)

        implicit none
        real*8, intent(in) :: p
        real*8 :: foc_port, omega_p, c_p, varphi, dist
        integer :: ioml, iomr, ig, iv, iq

        ! store portfolio share
        omega_p  = p

        foc_port = 0d0
           
        ! find the portfolio FOC
        do ig = 1, NW
            do iv = 1, NR
                do iq = 1, NE
                    ! derive interpolation weights
                    call linint_Grow(omega_p, omega_l, omega_u, omega_grow, &
                    NO, ioml, iomr, varphi)

                    ! get distributional weight
                    dist = dist_zeta(ig)*dist_vtheta(iv)*pi(iq_com, iq)

                    ! calculate consumption and FOC
                    c_p = varphi      *c(ij_com+1, ia_com, ioml, ig, iv, iq) + &
                          (1d0-varphi)*c(ij_com+1, ia_com, iomr, ig, iv, iq)
                    c_p = max(c_p, 1d-10)
                    foc_port = foc_port + dist*(rk(iv)-rb))*a(ia_com)*margu(c_p)
                enddo
            enddo
        enddo

    end function

!##############################################################################

    ! calculates marginal utility of consumption
    function margu(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu

        margu = cons**(-1/gamma)

    end function

!##############################################################################

    ! calculates the value function
    function valuefunc(a_plus, cons, ij, iq)

        implicit none
        integer, intent(in) :: ij
        real*8, intent(in) :: a_plus, cons
        real*8 :: valuefunc, varphi, c_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij < JJ)then
            valuefunc = max(varphi*Q(ij, ial, iq) + (1d0-varphi)*Q(ij, iar, iq), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = c_help**egam/egam + beta*psi(ij+1)*valuefunc

    end function

end module
