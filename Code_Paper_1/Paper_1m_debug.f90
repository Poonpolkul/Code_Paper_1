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
    integer, public, parameter :: JJ = 12
!~     integer, parameter :: JJ = 80
    
    ! number of years the household retires
    integer, public, parameter :: JR = 10
!~     integer, parameter :: JR = 40

    ! number of productivity type
    integer, public, parameter :: NP = 2 ! ip = 1,2

    ! number of technological shock process values
    integer, public, parameter :: NS = 5 ! is = 1,...,5

    ! number of points on the risky asset grid
    integer, public, parameter :: NK = 10
!~     integer, parameter :: NK = 100
    
    ! number of points on the risk-free asset grid
    integer, public, parameter :: NB = 20
!~     integer, parameter :: NB = 100

    ! household preference parameters 
    ! sigma and eta parameter value from Heer & Maussner 2012. Neeed gamma's value
    real*8, public, parameter :: gamma = 0.5 ! entropic risk aversion in Value function
    real*8, public, parameter :: beta  = 0.998
    real*8, public, parameter :: sigma = 0.5 !2d0
    real*8, public, parameter :: eta =  0.6 ! 7
    real*8, public, parameter :: phiu = 0.26 !for phi in utility. Value from HM(2012)

    ! household risk process
    real*8, public, parameter :: sigma_theta = 0.23d0
    
    ! transaction cost for changes in risk asset holdings
    real*8, public, parameter :: varphi1 = 0.06
    real*8, public, parameter :: varphi2 = 1.6 ! as in Cwik (2015), but too low
    
    ! production parameters
    real*8, public, parameter :: alpha = 0.36d0
    real*8, public, parameter :: delta = 1d0-(1d0-0.0823d0) !delta = 1d0-(1d0-0.0823d0)
    real*8, public, parameter :: Omega = 1.60d0 ! Need shock process
    
    ! aggregate risk process [just use parameter from Fehr's HH risk process. Need to adjust later]
    real*8, public, parameter :: rho         = 0.98d0
    real*8, public, parameter :: sigma_eps   = 0.05d0

    ! size of the risky asset grid
    real*8, public, parameter :: k_l    = 0.0d0
    real*8, public, parameter :: k_u    = 40d0
    real*8, public, parameter :: k_grow = 0.05d0
    
    ! size of the risk-free asset grid
    real*8, public, parameter :: b_l    = 0d0
    real*8, public, parameter :: b_u    = 40d0
    real*8, public, parameter :: b_grow = 0.05d0

    ! demographic parameters
    real*8, public, parameter :: n_p   = 0.01d0

    ! simulation parameters
    real*8, public, parameter :: damp    = 0.30d0
    real*8, public, parameter :: sig     = 1d-4
    integer, public, parameter :: itermax = 50

    ! counter variables
    integer :: iter

    ! macroeconomic variables
    real*8, public :: rb, rbn, rk, rkn,  w, wn, p
    real*8, public :: KK, AA, BB, LL, HH
    real*8, public :: YY, CC, II, GG, INC

    ! government variables
    real*8, public :: tauc, tauw, taurb, taurk, taup, kappa
    real*8, public :: gy, by, pen(JJ), PP, taxrev(4)
    integer, public :: tax
    logical, public :: ageing_on

    ! cohort aggregate variables
    real*8, public :: c_coh(JJ), y_coh(JJ), l_coh(JJ), k_coh(JJ), b_coh(JJ), v_coh(JJ)

    ! the shock process
    real*8, public :: dist_theta(NP), theta(NP)
    real*8, public :: pi(NS, NS), EtaShock(NS)
!~     real*8 :: eta(NS)
    integer, public :: is_initial = 3


    ! demographic and other model parameters
    real*8, public :: m(JJ)
    real*8, public :: eff(JJ)

    ! individual variables
    real*8, public :: k(0:NK), kplus(JJ, 0:NK, 0:NB, NP, NS), kplus_opt
    real*8, public :: b(0:NB), bplus(JJ, 0:NK, 0:NB, NP, NS), bplus_opt
    real*8, public :: c(JJ, 0:NK, 0:NB, NP, NS), l(JJ, 0:NK, 0:NB, NP, NS)
    real*8, public :: phi(JJ, 0:NK, 0:NB, NP, NS), V(JJ, 0:NK, 0:NB, NP, NS) = 0d0

    ! numerical variables
    real*8, public :: RHS(JJ, 0:NK, 0:NB, NP, NS), EV(JJ, 0:NK, 0:NB, NP, NS)
    real*8, public :: RHS_b(JJ, 0:NK, 0:NB, NP, NS), RHS_k(JJ, 0:NK, 0:NB, NP, NS)
    real*8, public :: RHSN_b(JJ, 0:NK, 0:NB, NP, NS), RHSN_k(JJ, 0:NK, 0:NB, NP, NS)
    real*8, public :: RHSD_b(JJ, 0:NK, 0:NB, NP, NS), RHSD_k(JJ, 0:NK, 0:NB, NP, NS)
    integer, public :: i, ij_com, ik_com, ib_com, ip_com, is_com, it_com
    real*8, public :: k_in_next, b_in_next, k_temp, b_temp
    real*8, public :: k_in, b_in
    real*8, public :: cons_com, lab_com, DIFF, INC_init
    
    ! make functions/subroutines public for use in 'Paper_1.f90'
    public :: valuefunc, foc
!~     public :: broydn
    
!##############################################################################
!##############################################################################
! Declaration of Variables
!##############################################################################
!##############################################################################

    ! declare everything as private by default
    private

    ! Level of tolerance for all routines
    real*8,  private  :: tbox_gftol_root = 1d-8

    ! Maximum number of iterations for broydn
    integer, private  :: itermax_root = 200


!##############################################################################
!##############################################################################
! Define public access points
!##############################################################################
!##############################################################################

    ! matrix operations
    public :: lu_solve, lu_invert, lu_dec
    public :: cholesky

    ! rootfinding
    public :: settol_root
    public :: setiter_root
    public :: fzero

    ! optimization
    public :: settol_min
    public :: setiter_min
    public :: fminsearch

    ! linear programming with constraints
    public :: solve_lin

    ! integration methods
    public :: legendre

    ! discretization of normal distributions
    public :: normal_discrete
    public :: log_normal_discrete

    ! discretization of AR(1) process
    public :: discretize_AR, discretize_log_AR

    ! probabilities and distributions
    public :: uniformPDF, uniformCDF, uniformCDF_Inv
    public :: normalPDF, normalCDF, normalCDF_Inv
    public :: log_normalPDF, log_normalCDF, log_normalCDF_Inv
    public :: GammaPDF, GammaCDF
    public :: betaPDF, betaCDF
    public :: bernoulliPDF, bernoulliCDF
    public :: binomialPDF, binomialCDF
    public :: binomial_coefficient

    ! simulation
    public :: simulate_uniform
    public :: simulate_normal
    public :: simulate_log_normal
    public :: simulate_Gamma
    public :: simulate_beta
    public :: simulate_bernoulli
    public :: simulate_binomial
    public :: simulate_AR

    ! discretization of continuous intervals
    public :: grid_Cons_Equi, grid_Val_Equi, grid_Inv_Equi
    public :: grid_Cons_Cheb, grid_Val_Cheb, grid_Inv_Cheb
    public :: grid_Cons_Grow, grid_Val_Grow, grid_Inv_Grow

    ! interpolation
    public :: poly_interpol
    public :: linint_Equi, linint_Cheb, linint_Grow, linint_Gen
    public :: spline_interp, spline_eval, spline

    ! plotting
    public :: plot
    public :: plot3d
    public :: plot_hist
    public :: execplot

    ! sorting
    public :: sort

    ! the clock
    public :: tic, toc


!~ !##############################################################################
!~ ! INTERFACE assert_eq
!~ !
!~ ! Interface for equality assertions by assert_eqx functions.
!~ !##############################################################################
    interface assert_eq

        module procedure assert_eq2, assert_eq3, assert_eq4, assert_eq5, &
            assert_eqn

    end interface


contains

    Subroutine bisection(f,x1,x2,eps,Root,check)

        implicit none
        real*8 :: f, x1, x2, eps, Root
        real*8 :: a, b, c
        integer :: i
        logical :: check
        integer, parameter:: iter=200
        
        !* check the bisection condition
        if(f(x1)>0.0 .AND. f(x2)>0.0) then
          Root = 1
          check = .false.
          return
        else if (f(x1)<0.0 .AND. f(x2)<0.0) then
          Root = 0
          check = .false.
          return
        end if
        
        !* initialize calculations
        a=x1
        b=x2

        !* Iterative refining the solution 
        do i=1,iter
!~         print *, i
          c=(b+a)/2.0
          if(f(a)*f(c).le.0.0) then
              b = c
            else
              a=c
          end if
        ! condition(s) to stop iterations)
          if(abs(b-a)<= eps) exit  
        end do
        Root=(b+a)/2.0

    end subroutine bisection
    
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
        + bplus_opt+kplus_opt &
        - (1d0+rb)*b(ib_com) - ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) &
        -(1-tauw)*wage*l_in-varphi1*(abs(kplus_opt-k(ik_com)))**varphi2
        
    end function


    function FOC(a_in)
        
        implicit none
        real*8, intent(in) :: a_in(:)
        real*8 :: foc(size(a_in,1))
        real*8 :: focb, fock, l_in, varphik, varphib
        real*8 :: tomorrow_k, tomorrow_b, wage
        real*8 :: x1, x2, available, eps, Root
        integer :: ikl, ikr, ibl, ibr
        logical :: check
        
        kplus_opt = a_in(1)
        bplus_opt = a_in(2)
        
        !initialize root value
        root = 1d0
        
        ! calculate the wage rate
        wage = wn*eff(ij_com)*theta(ip_com)

        ! calculate available resources
        available = ((1+(YY-wn*LL)/KK)-delta)*k(ik_com) + (1d0+rbn)*b(ib_com) + pen(ij_com)

        ! determine labour
        x1=0d0
        x2=1d0
        eps=1.0e-6
        if (ij_com < JR) then
            call bisection(implicitl, x1, x2, eps, Root, check)
            lab_com = Root
!~             print*, kplus_opt, bplus_opt, root
        else
            lab_com = 0d0
        endif
        
        ! calculate consumption
        cons_com = max((available + (1-tauw)*wage*lab_com - kplus_opt  & 
        - bplus_opt - varphi1*(abs(kplus_opt-k(ik_com)))**varphi2) ,1d-10)

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(kplus_opt, k_l, k_u, k_grow, NK, ikl, ikr, varphik)
        call linint_Grow(bplus_opt, b_l, b_u, b_grow, NB, ibl, ibr, varphib)
        
        tomorrow_k = varphik*RHS_k(ij_com+1, ikl, ib_com, ip_com, is_com) + &
        (1d0-varphik)*RHS_k(ij_com+1, ikr, ib_com,  ip_com, is_com) 
        
        tomorrow_b = varphib*RHS_b(ij_com+1, ik_com, ibl, ip_com, is_com) + &
        (1d0-varphib)*RHS_b(ij_com+1, ik_com, ibr, ip_com, is_com)

        ! calculate the first order conditions for consumption
        foc(1) = abs(cons_com**(-sigma) - tomorrow_k)
        foc(2) = abs(cons_com**(-sigma) - tomorrow_b)

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

    !##############################################################################
    ! SUBROUTINE broydn
    !
    ! Find root of multidimensional function of.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press. 
    !##############################################################################
    subroutine broydn(x, funcv, check_return)
    
        implicit none
     
     
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! initial guess and root of the function
        real*8, intent(inout) :: x(:)
     
        ! check is true if broydn converged to local minimum or can make no
        !     further progress
        logical, intent(out), optional :: check_return
     
     
        !##### OTHER VARIABLES ####################################################
     
        real*8, parameter :: eps = epsilon(x)
        real*8 :: tolf, tolmin
        real*8, parameter :: tolx = eps
        real*8, parameter :: stpmx = 100d0
        integer :: i, j, k, its, n
        real*8 :: f, fold, stpmax
        real*8, dimension(size(x)) :: c, d, fvcold, g, p, s, t, w, xold
        real*8, dimension(size(x),size(x)) :: qt, r
        logical :: restrt, sing, check
        real*8 :: fvec(size(x,1))
     
     
        !##### INTERFACES #########################################################
     
        ! interface for the function
        interface
            function funcv(p)
                implicit none
                real*8, intent(in) :: p(:)
                real*8 :: funcv(size(p, 1))
            end function funcv
        end interface
     
     
        !##### ROUTINE CODE #######################################################
     
        ! set tolerance levels
        tolf = tbox_gftol_root
        tolmin = tbox_gftol_root
        if(present(check_return))check_return = .false.
     
        ! get size of x
        n = size(x)   
     
        ! calculate function euklidean norm at starting point
        f = fmin(x, funcv)
     
        ! check if root has been found
        if (maxval(abs(fvec(:))) < 0.01d0*tolf) then
            if(present(check_return))check_return = .false.
            return
        endif
     
        stpmax = stpmx*max(sqrt(dot_product(x(:),x(:))),dble(n))
        restrt = .true.
     
        ! iterate broydn steps
        do its=1,itermax_root
     
            ! If restart then calculate jacobian of function
            if (restrt) then
     
                ! calculate jacobian of func at x
                call fdjac(x, fvec, r, funcv)
     
                ! make q-r-decomposition of jacobian
                call qrdcmp(r, c, d, sing)
     
                ! throw error if jacobian is singular
                if(sing)then
                    call warning('fzero', 'singular jacobian')
                    if(present(check_return))check_return = .true.
                    return
                endif
     
                ! create unity matrix
                qt(:,:) = 0d0
                    do j = 1, n
                            qt(j, j) = 1d0
                    enddo
     
                ! for Q^T explicitly
                do k = 1, n-1
                    if (abs(c(k)) >= 1d-100) then
                        qt(k:n,:) = qt(k:n, :)-outerprod(r(k:n, k), &
                            matmul(r(k:n, k), qt(k:n, :)))/c(k)
                    endif
                enddo
                where(lower_triangle(n,n))r(:, :) = 0d0
     
                ! puts diagonal elements of R matrix to r
                do j = 1, n
                    r(j, j) = d(j)
                enddo
     
            ! else do Broydn update step
            else
     
                ! set up s as delta x
                s(:) = x(:)-xold(:)
     
                ! t = R*delta x
                do i = 1, n
                    t(i) = dot_product(r(i,i:n), s(i:n))
                enddo
     
                ! w = delta f - B*s = delta f - R*s*Q^T
                w(:) = fvec(:)-fvcold(:)-matmul(t(:), qt(:,:))
     
                ! if w entries are small enough, set them to zero
                where(abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) w(:) = 0d0
     
                ! update for non-noisy components of w
                if(any(abs(w(:)) >= 1d-100))then
     
                    ! update t and s
                    t(:) = matmul(qt(:,:),w(:))
                    s(:)=s(:)/dot_product(s,s)
     
                    ! update R and Q^T
                    call qrupdt(r,qt,t,s)
     
                    ! get diagonal of matrix r
                    do j = 1, size(r,1)
                        d(j) = r(j,j)
                    enddo
     
                    ! if any diagonal value of r is 0, then jacobian is singular
                    if(any(abs(d(:)) <= 1d-100))then
                        call warning('fzero', 'singular jacobian')
                        if(present(check_return))check_return = .true.
                        return  
                    endif
                endif
            endif
     
            ! perform the newton step by inverting jacobian
            p(:) = -matmul(qt(:,:), fvec(:))
            do i = 1, n
                g(i) = -dot_product(r(1:i,i), p(1:i))
            enddo
     
            ! store old x, function value and function norm
            xold(:) = x(:)
            fvcold(:) = fvec(:)
            fold = f
     
            ! solve linear equation with upper triangular matrix r
            call rsolv(r, d, p)
     
            ! searches along the new gradient direction for new x and f
            call lnsrch(xold, fold, g, p, x, f, stpmax, check, funcv)
     
            ! check whether root was found
            if(maxval(abs(fvec(:))) < tolf)then
                if(present(check_return))check_return = .false.
                return
            endif
     
            ! if check is true
            if(check)then
     
                ! check if improvement can be made, if not, return
                if(restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
                        1d0)/max(f, 0.5d0*n)) < tolmin)then
                    if(present(check_return))check_return = check
                    return
                endif
     
                ! else calculate new jacobian
                restrt=.true.
     
            ! if check is false
            else
     
                ! do broydn step
                restrt=.false.
     
                ! check for convergence
                if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
                    1.0d0)) < tolx)then
                    if(present(check_return))check_return = check
                    return
                endif
            endif
        enddo
     
        ! throw warning if broydn didn't converge
        if(present(check_return))check_return = .true.
     
     
    !##### SUBROUTINES AND FUNCTIONS ##########################################
     
    contains
     
     
        !##########################################################################
        ! FUNCTION fdjac
        !
        ! Calculates finite difference jacobian.
        !##########################################################################
        subroutine fdjac(x, fvec, df, funcv)
     
            implicit none
     
     
            !##### INPUT/OUTPUT VARIABLES #########################################
     
            ! value where to calculate finite difference jacobian
            real*8, intent(inout) :: x(:)
     
            ! function value at x
            real*8, intent(in) :: fvec(:)
     
            ! resulting finite difference jacobian
            real*8, intent(out) :: df(:, :)
     
     
            !##### OTHER VARIABLES ################################################
     
            real*8, parameter :: eps = 1.0e-6
            integer :: j, n
            real*8, dimension(size(x)) :: xsav, xph, h
     
     
            !##### INTERFACES #####################################################
     
            ! interface for the function
            interface
                function funcv(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: funcv(size(p, 1))
                end function funcv
            end interface
     
     
            !##### ROUTINE CODE ###################################################
     
            ! check equality of sizes
            n = assert_eq(size(x), size(fvec), size(df,1), size(df,2), 'fdjac')
     
            ! store old x
            xsav = x
     
            ! calculate difference
            h = eps*abs(xsav)
            where(abs(h) <= 1d-100)h = EPS
     
            ! calculate x + h
            xph = xsav + h
            h = xph - xsav
     
            ! itertate over dimensions and calculate difference
            do j = 1, n
                x(j) = xph(j)
                df(:,j) = (funcv(x)-fvec(:))/h(j)
                x(j) = xsav(j)
            enddo
     
        end subroutine fdjac
     
     
        !##########################################################################
        ! FUNCTION lnsrch
        !
        ! Finds point along a line, given function value and gradient, where
        !     function has decreased sufficiently (for one dimensional function).
        !##########################################################################
        subroutine lnsrch(xold, fold, g, p, x, f, stpmax, check, funcv)
     
            implicit none
     
     
            !##### INPUT/OUTPUT VARIABLES #########################################
     
            ! point where to start line search
            real*8, intent(in) :: xold(:)
     
            ! the old function value
            real*8, intent(in) :: fold
     
            ! gradient at this point
            real*8, intent(in) :: g(:)
     
            ! a line search direction
            real*8, intent(inout) :: p(:)
     
            ! new value along the search line
            real*8, intent(out) :: x(:)
     
            ! function value at new x
            real*8, intent(out) :: f
     
            ! maximum size of steps such that lnsrch does not search un undefined
            !     areas
            real*8, intent(in) :: stpmax
     
            ! is true if x is too close at xold
            logical, intent(out) :: check
     
     
            !##### OTHER VARIABLES ################################################
     
            real*8, parameter :: alf = 1.0e-4
            real*8, parameter :: tolx = epsilon(x)
            integer :: ndum
            real*8 :: a, alam, alam2=0d0, alamin, b, disc, f2=0d0, pabs, rhs1, rhs2, &
                slope, tmplam
     
     
            !##### INTERFACES #####################################################
     
            ! interface for the function
            interface
                function funcv(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: funcv(size(p, 1))
                end function funcv
            end interface
     
     
            !##### ROUTINE CODE ###################################################
     
            ! assert sizes or arrays
            ndum = assert_eq(size(g), size(p), size(x), size(xold), 'lnsrch')
            ndum = ndum
     
            ! set check's default value
            check=.false.
     
            ! calculate norm of p
            pabs = sqrt(dot_product(p, p))
     
            ! restrict p to maximum stepsize
            if(pabs > stpmax)p(:) = p(:)*stpmax/pabs
     
            ! calculate slope
            slope = dot_product(g, p)
     
            ! throw error if you would go uphill
            if(slope >= 0d0)then
                call warning('lnsrch', 'roundoff problem, I cannot go uphill')
                return
            endif
     
            ! calculate newton stepsize
            alamin = tolx/maxval(abs(p(:))/max(abs(xold(:)),1d0))
            alam = 1d0
     
            ! start iteration
            do
                ! calculate calculate new x
                x(:) = xold(:)+alam*p(:)
     
                ! calculate new function value at x
                f = fmin(x, funcv)
     
                ! if new x is not away enough return with check=true
                if(alam < alamin)then
                    x(:) = xold(:)
                    check = .true.
                    return
     
                ! if optimal value found return with false
                elseif(f <= fold+alf*alam*slope)then
                    return
     
                ! else do backtracking
                else
                    if(abs(alam -1d0) <= 1d-100)then
                        tmplam = -slope/(2d0*(f-fold-slope))
                    else
                        rhs1 = f-fold-alam*slope
                        rhs2 = f2-fold-alam2*slope
                        a = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                        b = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
                        if(abs(a) <= 1d-100)then
                            tmplam = -slope/(2d0*b)
                        else
                            disc = b*b-3d0*a*slope
                            if(disc < 0d0)then
                                tmplam = 0.5d0*alam
                            elseif(b <= 0d0)then
                                tmplam = (-b+sqrt(disc))/(3d0*a)
                            else
                                tmplam = -slope/(b+sqrt(disc))
                            endif
                        endif
                        if(tmplam > 0.5d0*alam)tmplam = 0.5d0*alam
                    endif
                endif
                alam2 = alam
                f2 = f
                alam = max(tmplam,0.1d0*alam)
            enddo
     
        end subroutine lnsrch
     
     
        !##########################################################################
        ! FUNCTION fmin
        !
        ! Calculates vector norm of multidimensional function.
        !##########################################################################
        function fmin(x, funcv)
     
            implicit none
     
     
            !##### INPUT/OUTPUT VARIABLES #########################################
     
            ! value where to evaluate function
            real*8, intent(in) :: x(:)
     
            ! euklidean square norm of function at x
            real*8 :: fmin
     
     
            !##### INTERFACES #####################################################
     
            ! interface for the function
            interface
                function funcv(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: funcv(size(p, 1))
                end function funcv
            end interface
     
            ! calculate function value
            fvec = funcv(x)
     
            ! calculate squared norm
            fmin = 0.5d0*dot_product(fvec, fvec)
     
        end function fmin
     
     
        !##########################################################################
        ! SUBROUTINE qrdcmp
        !
        ! Calculates QR decomposition of a matrix.
        !##########################################################################
        subroutine qrdcmp(a, c, d, sing)
     
            implicit none
     
            real*8, intent(inout) :: a(:, :)
            real*8, intent(out) :: c(:), d(:)
            logical, intent(out) :: sing
            integer :: k, n
            real*8 :: scale, sigma
     
            n = assert_eq(size(a,1), size(a,2), size(c), size(d), 'qrdcmp')
            sing = .false.
            do k = 1, n-1
                scale = maxval(abs(a(k:n, k)))
                if(abs(scale) <= 1d-100)then
                        sing = .true.
                        c(k) = 0d0
                        d(k) = 0d0
                else
                        a(k:n, k) = a(k:n, k)/scale
                        sigma = sign(sqrt(dot_product(a(k:n, k),a(k:n, k))),a(k, k))
                        a(k,k) = a(k, k)+sigma
                        c(k) = sigma*a(k, k)
                        d(k) = -scale*sigma
                        a(k:n, k+1:n) = a(k:n, k+1:n)-outerprod(a(k:n, k),&
                                matmul(a(k:n, k),a(k:n, k+1:n)))/c(k)
                endif
            enddo
            d(n) = a(n, n)
            if (abs(d(n)) <= 1d-100) sing = .true.
     
        end subroutine qrdcmp
     
     
        !##########################################################################
        ! SUBROUTINE qrupdt
        !
        ! Updates qr-matrices.
        !##########################################################################
        subroutine qrupdt(r,qt,u,v)
     
            implicit none
     
            real*8, intent(inout) :: r(:, :), qt(:, :)
            real*8, intent(inout) :: u(:)
            real*8, intent(in) :: v(:)
            integer :: i, k, n
     
            n = assert_eq((/ size(r,1), size(r,2), size(qt,1), size(qt,2), &
                size(u), size(v)/), 'qrupdt')
            k = n+1-ifirstloc(abs(u(n:1:-1)) >= 1d-100)
            if(k < 1)k=1
            do i = k-1, 1, -1
                call rotate(r,qt,i,u(i),-u(i+1))
                u(i) = pythag(u(i),u(i+1))
            enddo
            r(1,:) = r(1,:)+u(1)*v
            do i = 1,k-1
                call rotate(r,qt,i,r(i,i),-r(i+1,i))
            enddo
        end subroutine qrupdt
     
        !##########################################################################
        ! SUBROUTINE rsolv
        !
        ! Solves upper diagonal system.
        !##########################################################################
        subroutine rsolv(a, d, b)
     
            implicit none
     
            real*8, intent(in) :: a(:, :), d(:)
            real*8, intent(inout) :: b(:)
            integer :: i, n
     
            n = assert_eq(size(a,1), size(a,2), size(b), size(d), 'rsolv')
            b(n) = b(n)/d(n)
            do i = n-1, 1, -1
                    b(i) =( b(i)-dot_product(a(i, i+1:n),b(i+1:n)))/d(i)
            enddo
     
        end subroutine rsolv
     
     
        subroutine rotate(r, qt, i, a, b)
     
            implicit none
     
            real*8, intent(inout) :: r(:, :), qt(:, :)
            integer, intent(in) :: i
            real*8, intent(in) :: a, b
            integer :: n
            real*8 :: c, fact, s, temp(size(r,1))
     
            n = assert_eq(size(r,1), size(r,2), size(qt,1), size(qt,2), 'rotate')
            if(abs(a) <= 1d-100)then
                c = 0d0
                s = sign(1d0, b)
            elseif(abs(a) > abs(b))then
                fact = b/a
                c = sign(1d0/sqrt(1d0+fact**2), a)
                s = fact*c
            else
                fact = a/b
                s = sign(1d0/sqrt(1d0+fact**2), b)
                c=fact*s
            endif
            temp(i:n) = r(i, i:n)
            r(i, i:n) = c*temp(i:n)-s*r(i+1, i:n)
            r(i+1, i:n) = s*temp(i:n)+c*r(i+1, i:n)
            temp = qt(i, :)
            qt(i, :) = c*temp-s*qt(i+1, :)
            qt(i+1, :) = s*temp+c*qt(i+1, :)
     
        end subroutine rotate
     
     
        function pythag(a, b)
     
            implicit none
     
            real*8, intent(in) :: a, b
            real*8 :: pythag
            real*8 :: absa, absb
     
            absa = abs(a)
            absb = abs(b)
            if(absa > absb)then
                pythag = absa*sqrt(1d0+(absb/absa)**2)
            else
                if(abs(absb) <= 1d-100)then
                    pythag = 0d0
                else
                    pythag = absb*sqrt(1d0+(absa/absb)**2)
                endif
            endif
     
        end function pythag
     
     
        function ifirstloc(mask)
     
            logical, intent(in) :: mask(:)
            integer :: ifirstloc, loca(1)
     
            loca = maxloc(merge(1, 0, mask))
            ifirstloc = loca(1)
            if(.not. mask(ifirstloc))ifirstloc = size(mask)+1
     
        end function ifirstloc
     
     
        function lower_triangle(j, k, extra)
     
            integer, intent(in) :: j, k
            integer, intent(in), optional :: extra
            logical :: lower_triangle(j, k)
            integer :: n
     
            n = 0
            if(present(extra))n = extra
     
            lower_triangle = (outerdiff(arth_i(1, 1, j), arth_i(1, 1, k)) > -n)
     
        end function lower_triangle
     
     
        function outerdiff(a, b)
     
            integer, intent(in) :: a(:), b(:)
            integer :: outerdiff(size(a, 1),size(b, 1))
     
            outerdiff = spread(a, dim=2, ncopies=size(b, 1)) - &
                spread(b, dim=1, ncopies=size(a, 1))
     
        end function outerdiff
     
     
        function outerprod(a, b)
     
            real*8, intent(in) :: a(:), b(:)
            real*8 :: outerprod(size(a, 1),size(b, 1))
     
            outerprod = spread(a, dim=2, ncopies=size(b, 1)) * &
                spread(b, dim=1, ncopies=size(a, 1))
     
        end function outerprod
     
     
        function arth_i(first, increment, n)
     
            integer, intent(in) :: first, increment, n
            integer, parameter :: npar_arth = 16
            integer, parameter :: npar2_arth = 8
            integer :: arth_i(n)
            integer :: k, k2, temp
     
            if(n > 0)arth_i(1) = first
            if(n <= npar_arth) then
                do k = 2, n
                    arth_i(k) = arth_i(k-1) + increment
                enddo
            else
                do k = 2, npar2_arth
                    arth_i(k) = arth_i(k-1) + increment
                enddo
                temp = increment*npar2_arth
                k = npar2_arth
                do
                    if(k >= n)exit
                    k2 = k+k
                    arth_i(k+1:min(k2,n)) = temp+arth_i(1:min(k,n-k))
                    temp = temp + temp
                    k = k2
                enddo
            endif
        end function arth_i
    
    end subroutine broydn

!############################################################################## 
!##############################################################################
! MODULE assertions
!##############################################################################
!##############################################################################
 
 
    !##############################################################################
    ! FUNCTION assert_eq2
    !
    ! Checks equality for two integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq2(n1, n2, string)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! integers to compare for equality
        integer, intent(in) :: n1, n2
        
        ! routine from which error should be thrown
        character(len=*), intent(in) :: string
        
        ! return value
        integer :: assert_eq2
        
        
        !##### ROUTINE CODE #######################################################
        
        ! if equality, set return value to n1
        if (n1 == n2)then
            assert_eq2 = n1
        
        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq2')
        end if
    
    end function assert_eq2
    
    
    !##############################################################################
    ! FUNCTION assert_eq3
    !
    ! Checks equality for three integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq3(n1, n2, n3, string)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! integers to compare for equality
        integer, intent(in) :: n1, n2, n3
        
        ! routine from which error should be thrown
        character(len=*), intent(in) :: string
        
        ! return value
        integer :: assert_eq3
        
        
        !##### ROUTINE CODE #######################################################
        
        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3)then
            assert_eq3 = n1
        
        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq3')
        end if
    
    end function assert_eq3
    
    
    !##############################################################################
    ! FUNCTION assert_eq4
    !
    ! Checks equality for four integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq4(n1, n2, n3, n4, string)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! integers to compare for equality
        integer, intent(in) :: n1, n2, n3, n4
        
        ! routine from which error should be thrown
        character(len=*), intent(in) :: string
        
        ! return value
        integer :: assert_eq4
        
        
        !##### ROUTINE CODE #######################################################
        
        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4)then
            assert_eq4 = n1
        
        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq4')
        end if
    
    end function assert_eq4
    
    
    !##############################################################################
    ! FUNCTION assert_eq5
    !
    ! Checks equality for five integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq5(n1, n2, n3, n4, n5, string)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! integers to compare for equality
        integer, intent(in) :: n1, n2, n3, n4, n5
        
        ! routine from which error should be thrown
        character(len=*), intent(in) :: string
        
        ! return value
        integer :: assert_eq5
        
        
        !##### ROUTINE CODE #######################################################
        
        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4 .and. n4 == n5)then
            assert_eq5 = n1
        
        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq5')
        end if
    
    end function assert_eq5

    
    
    !##############################################################################
    ! FUNCTION assert_eqn
    !
    ! Checks equality for n integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eqn(nn, string)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! integers to compare for equality
        integer, intent(in) :: nn(:)
        
        ! routine from which error should be thrown
        character(len=*), intent(in) :: string
        
        ! return value
        integer :: assert_eqn
        
        
        !##### ROUTINE CODE #######################################################
        
        ! if equality, set return value to n1
        if (all(nn(2:) == nn(1)))then
            assert_eqn = nn(1)
        
        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eqn')
        end if
    
    end function assert_eqn
 




end module
