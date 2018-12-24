!##############################################################################
! PROGRAM PortfolioChoice
!
! ## Portfolio choice in the life cycle model (with different aggregations)
! with engodenous rates of returns on bond and capital
!
! This code is adjusted from Hans Fehr and Fabian Kindermann
!
! #VC# VERSION: 7.2  (20 Dec 2019)
!
!##############################################################################
include "Paper_1m_V7.2.f90"

program PortfolioChoice

    use globals

    implicit none

    ! initialize remaining variables
    call initialize()

    ! start the clock
    call tic()

    ! iterate until value function converges
    do iter = 1, itermax

        iterb = 0
        do while (abs(BB)>= sig .or. iterb == 0)
            iterb = iterb + 1

Print*, ' ################### iteration (outer, inner) =', iter, iterb, '############################'

            ! derive prices
            call prices()

            ! solve the household problem
            call solve_household()

            ! calculate the distribution of households over state space
            call get_distribution()

            ! aggregate individual decisions
            call aggregation()

            ! update bond market return
            call bond_return()

            ! determine the government parameters
            call government()

            if(abs(DIFF/YY)*100d0 < sig .and. abs(BB) < sig)then
                call toc
                call output()
                return
            endif
        enddo
    enddo
    
    ! stop the clock
    call toc()
    call output()
        
    write(*,'(a/)')'ATTENTION: NO CONVERGENCE !!!'

    ! close files
    close(21)

contains


    ! initializes all remaining variables
    subroutine initialize

        implicit none
        integer*8 :: ij

        ! set up population structure
        do ij = 1, JJ
            m(ij) = (1d0+n_p)**(1d0-ij)  
        enddo

        ! set survival probabilities
        psi = (/1.00000d0, 0.97455d0, 0.94909d0, 0.92364d0, 0.89818d0, &
                0.87273d0, 0.84727d0, 0.82182d0, 0.79636d0, 0.77091d0, &
                0.74545d0, 0.72000d0, 0d0/)

        ! initialize age earnings process
        eff(1) = 1.0000d0
        eff(2) = 1.3527d0
        eff(3) = 1.6952d0
        eff(4) = 1.8279d0
        eff(5) = 1.9606d0
        eff(6) = 1.9692d0
        eff(7) = 1.9692d0
        eff(8) = 1.9392d0
        eff(9) = 1.9007d0
        eff(JR:JJ) = 0d0

        ! discretize zeta shocks
        call normal_discrete(zeta, dist_zeta, 0d0, sigma_zeta)
        zeta = exp(zeta)

        ! discretize eps shocks        
        call normal_discrete(eps, dist_eps, 0d0, sigma_eps)
        eps = exp(eps)

        ! discretize TProd shocks        
        call normal_discrete(TProd, dist_TProd, TProd_bar, sigma_vtheta)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)


        ! initialize tax and transfers
        tauw  = 0.0d0

        ! initial guesses for macro variables
        KK = 1d0
        BB = 0.5d0
        LL =  1d0
        YY = 1d0
        II = (n_p+delta)*KK

        ! open files
        open(21, file='output.out')



    end subroutine

!##############################################################################

    ! subroutine for calculating prices
    subroutine prices()

        implicit none
        integer :: ir
        
        ! calculate wage rate
        do ir = 1, NR
            w(ir) = TProd(ir)*(1d0-alpha)*(KK/LL)**alpha
print*, 'w', w(ir)
        enddo

        ! calculate return on risky asset
        do ir = 1, NR
            rk(ir) = (TProd(ir)*KK**alpha*LL**(1d0-alpha)-w(ir)*LL)/KK
print*, 'rk', rk(ir)
        enddo

        if (iter==1 .and. iterb == 1) then
            rb = 0
            do ir = 1, NR
                rb = rb + (rk(ir)-mu_r)/NR

print*, 'rb', rb
            enddo
        endif

print*, 'rb_end', rb       
        ! calculate after-tax wage rate
        do ir = 1, NR
            wn(ir) = w(ir)*(1d0-tauw)
!~ print*, 'wn(iv)', wn(ir)
        enddo
        
        ! old-age transfers
        pen(:,:) = 0d0
        do ir = 1, NR
            pen(JR:JJ, ir) = max(kappa*w(ir)*eff(JR-1), 1d-10)
!~ print*, 'pension iv', pen(JJ, ir)
        enddo
        
        ! endogenous lower and upper bound of cash-on-hand grid
        X_l = min((1d0-tauw)*w(1)*minval(eff(1:JR-1))*minval(eps(:))*zeta(1), pen(JR, 1))
        X_u = (1d0 + rk(NR))*a_u + &
                    w(NR)*maxval(eff(1:JR-1))*maxval(eps(:))*zeta(NW)
        call grid_Cons_Grow(X, X_l, X_u, X_grow)
!~ Print*, 'X', X

    end subroutine

!##############################################################################

    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ix, ia

        ! get decisions in last period of life
        omega_plus(JJ, :) = 0d0
        do ix = 0, NX
            a_plus(JJ, ix) = 0d0
            c(JJ, ix) = X(ix)
            V(JJ, ix) = valuefunc(0d0, c(JJ, ix), JJ)
        enddo

        do ij = JJ-1, 1, -1

            ! determine optimal portfolio choice for all others
            do ia = 1, NA
                call solve_portfolio(ij, ia)
            enddo

            ! set omega for zero savings consistent with next gridpoint
            omega_plus(ij, 0) = omega_plus(ij, 1)

            ! interpolate individual RHS and value function
            call interpolate(ij)

            ! determine consumption-savings solution
            do ix = 0, NX
                call solve_consumption(ij, ix)
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

        enddo

    end subroutine


    ! solve the household's portfolio decision
    subroutine solve_portfolio(ij, ia)

        implicit none
        integer, intent(in) :: ij, ia
        real*8 :: x_in, port0, port1, tolerance
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia

        ! check for corner solutions
        port0 = foc_port(0d0)
        port1 = foc_port(1d0)

!##############################################################################
        ! use intermediate value theorem
        if(port0*port1 > 0d0 .and. abs(port0) <= abs(port1))then
                omega_plus(ij, ia) = 0d0
            return
        else
!##############################################################################

            ! get order of magnitude of foc
            tolerance = 1d-5*abs(port0-port1)

            tolerance = min(tolerance, 1d-8)

            call settol_root(tolerance)

            ! get best guess for the root of foc_port
            x_in = -port0/(port1-port0)
            check = .false.
            ! solve the household problem using rootfinding
            call fzero(x_in, foc_port, check)

            ! write screen output in case of a problem
            if(check)write(*,'(a, 2i4)')'ERROR IN ROOTFINDING PORT : ', ij, ia

            omega_plus(ij, ia) = x_in

            ! reset tolerance level to original value
            call settol_root(1d-8)
        endif

    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption(ij, ix)

        implicit none
        integer, intent(in) :: ij, ix
        real*8 :: x_in
        logical :: check

        ! determine decision for zero cash-on-hand
        if(X(ix) < 1d-10)then
            a_plus(ij, ix) = 0d0
            c(ij, ix) = 0d0
            V(ij, ix) = valuefunc(0d0, 0d0, ij)
            return
        endif

        ! set up communication variables
        ij_com = ij
        ix_com = ix

        ! get best initial guess from future period
        x_in = a_plus(ij+1, ix)
        check = .false.

        ! solve the household problem using rootfinding
        call fzero(x_in, foc_cons, check)

        ! write screen output in case of a problem
        if(check)write(*,'(a, 2i4)')'ERROR IN ROOTFINDING CONS : ', ij, ix

        ! check for borrowing constraint
        if(x_in < 0d0)then
            x_in = 0d0
            cons_com = X(ix)
        endif

        ! copy decisions
        a_plus(ij, ix) = x_in
        c(ij, ix) = cons_com
        V(ij, ix) = valuefunc(x_in, cons_com, ij)

    end subroutine


    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, iw, is, ir !, isr
        real*8 :: X_p, c_p, varphi, dist, EV, R_port
        integer :: ixl, ixr

        RHS(ij, :) = 0d0
        Q(ij, :) = 0d0

        do ia = 0, NA

            ! case agent is retired tomorrow
            if(ij >= JR-1)then

                do ir = 1, NR

                    ! get return on the portfolio
                    R_port = 1d0 + rb + omega_plus(ij, ia)*(rk(ir) - rb)
        

                    ! get tomorrow's cash-on-hand (epsilon^+ = 0)
                    X_p = R_port*a(ia) + pen(ij+1, ir)

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_TProd(ir)

                    ! get future consumption value
                    c_p = max(varphi*c(ij+1, ixl) + (1d0-varphi)*c(ij+1, ixr), 1d-10)

                    ! get tomorrow's value function
                    EV = varphi      *(egam*V(ij+1, ixl))**(1d0/egam) + &
                         (1d0-varphi)*(egam*V(ij+1, ixr))**(1d0/egam)

                    ! get RHS of foc and Q
                    RHS(ij, ia) = RHS(ij, ia) + dist*R_port*margu(c_p)
                    Q(ij, ia)   = Q(ij, ia) + dist*EV**egam/egam

                enddo

            ! agent is working
            else
                do iw = 1, NW
                    do ir = 1, NR
                        do is = 1, NS

                            ! get return on the portfolio

                            R_port = 1d0 + rb + omega_plus(ij, ia)*(rk(ir) - rb)
                            
                            ! get tomorrow's cash on hand
                            X_p = R_port*a(ia)/eps(is) + wn(ir)*eff(ij+1)*zeta(iw)

                            ! derive interpolation weights
                            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                            ! get distributional weight
                            dist = dist_zeta(iw)*dist_eps(is)*dist_TProd(ir)

                            ! get future consumption value
                            c_p = max(varphi*c(ij+1, ixl) + (1d0-varphi)*c(ij+1, ixr), 1d-10)

                            ! get tomorrow's value function
                            EV = varphi      *(egam*V(ij+1, ixl))**(1d0/egam) + &
                                 (1d0-varphi)*(egam*V(ij+1, ixr))**(1d0/egam)

                            ! get RHS of foc and Q
                            RHS(ij, ia) = RHS(ij, ia) + dist*R_port*margu(eps(is)*c_p)

                            Q(ij, ia)   = Q(ij, ia) + dist*(eps(is)*EV)**egam/egam
                        enddo
                    enddo
                enddo
            endif

            RHS(ij, ia) = (beta*psi(ij+1)*RHS(ij, ia))**(-gamma)

            Q(ij, ia)   = (egam*Q(ij, ia))**(1d0/egam)

        enddo

    end subroutine


    ! determines the invariant distribution over state space
    subroutine get_distribution()

        implicit none
        integer :: ij

        ! set distributions to zero
        phi_X(:, :) = 0d0
        phi_a(:, :) = 0d0

        do ij = 1, JJ

            ! get distribution on cash-on-hand grid
            call get_distribution_X(ij)

            ! get distribution on asset grid
            call get_distribution_a(ij)
        enddo

    end subroutine


    ! to calculate distribution on cash-on-hand grid
    subroutine get_distribution_X(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, iw, is, ir, ixl, ixr !isr
        real*8 :: varphi, X_p, R_port, dist

        if(ij == 1)then

            ! get initial distribution at age 1 of cash-on-hand
            do iw = 1, NW
                do ir = 1, NR

                    ! get initial cash-on-hand
                    X_p = wn(ir)*eff(1)*zeta(iw)

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_zeta(iw)*dist_TProd(ir)

                    ! initialize the distribution
                    phi_X(1, ixl) = phi_X(1, ixl) + dist*varphi
                    phi_X(1, ixr) = phi_X(1, ixr) + dist*(1d0-varphi)
                enddo
            enddo

        elseif(ij <= JR-1)then

            ! iterate over yesterdays asset distribution
            do ia = 0, NA

                ! iterate over current shocks
                do iw = 1, NW
                    do ir = 1, NR
                        do is = 1, NS

                            ! get today's cash-on-hand

                            R_port = 1d0 + rb + omega_plus(ij-1, ia)*(rk(ir) - rb)
                            
                            X_p = R_port*a(ia)/eps(is) + wn(ir)*eff(ij)*zeta(iw)

                            ! derive interpolation weights
                            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                            ! get distributional weight
                            dist = dist_zeta(iw)*dist_eps(is)*dist_TProd(ir)

                            phi_X(ij, ixl) = phi_X(ij, ixl) + dist*varphi*phi_a(ij-1, ia)
                            phi_X(ij, ixr) = phi_X(ij, ixr) + dist*(1d0-varphi)*phi_a(ij-1, ia)
                        enddo
                    enddo
                enddo
            enddo

        else

            ! iterate over yesterdays asset distribution
            do ia = 0, NA

                ! iterate over current shocks
                do ir = 1, NR

                    ! get today's cash-on-hand

                    R_port = 1d0 + rb + omega_plus(ij-1, ia)*(rk(ir) - rb)
                    
                    X_p = R_port*a(ia) + pen(ij, ir)

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_TProd(ir)

                    phi_X(ij, ixl) = phi_X(ij, ixl) + dist*varphi*phi_a(ij-1, ia)
                    phi_X(ij, ixr) = phi_X(ij, ixr) + dist*(1d0-varphi)*phi_a(ij-1, ia)
                enddo
            enddo
        endif

    end subroutine


    ! to calculate the end of period asset distribution
    subroutine get_distribution_a(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ix, ial, iar
        real*8 :: varphi

        ! iterate over todays cash on hand
        do ix = 0, NX

            ! interpolate asset decision
            call linint_Grow(a_plus(ij, ix), a_l, a_u, a_grow, NA, ial, iar, varphi)

            ! restrict values to grid just in case
            ial = min(ial, NA)
            iar = min(iar, NA)
            varphi = min(varphi, 1d0)

            ! get end of period asset distribution
            phi_a(ij, ial) = phi_a(ij, ial) + varphi*phi_X(ij, ix)
            phi_a(ij, iar) = phi_a(ij, iar) + (1d0-varphi)*phi_X(ij, ix)
        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ix, ia, ie, iw, ir
        real*8 :: sigma_eta(JJ), mu_exp(JJ), sigma_exp(JJ)
        real*8 :: var_c(JJ), var_a(JJ), var_y(JJ), var_o(JJ)
        real*8 :: workpop, LL_old, KK_old, BB_old

        ! generate eta distribution if not analytical calculation
        if(.not. analytical)call generate_eta()

        ! Store current aggregate quantities
        LL_old = LL 
        KK_old = KK

        ! calculate cohort averages
        c_coh(:) = 0d0
        y_coh(:) = 0d0
        a_coh(:) = 0d0
        o_coh(:) = 0d0
        l_coh(:) = 0d0
        k_coh(:) = 0d0
        b_coh(:) = 0d0

        ! analytical approach or not
        if(analytical)then

            do ij = 1, JJ

                do iw = 1, NW
                    do ir = 1, NR
                        y_coh(ij) = y_coh(ij) + (w(ir)*eff(ij)*zeta(iw)&
                                    + pen(ij, ir))*dist_zeta(iw)*dist_TProd(ir)
                    enddo
                enddo

                do iw = 1, NW
                    l_coh(ij) = l_coh(ij) + eff(ij)*zeta(iw)*dist_zeta(iw)
                enddo

                do ix = 0, NX
                    c_coh(ij) = c_coh(ij) + c(ij, ix)*phi_X(ij, ix)
                enddo

                if(ij > 1)then
                    do ia = 0, NA
                        a_coh(ij) = a_coh(ij) + a(ia)*phi_a(ij-1, ia)
                        o_coh(ij) = o_coh(ij) + omega_plus(ij-1, ia)*phi_a(ij-1, ia)
                    enddo
                endif
                
                if(ij > 1)then
                    k_coh(ij) = o_coh(ij)*a_coh(ij)
                    b_coh(ij) = (1-o_coh(ij))*a_coh(ij)
print*, 'ij:', ij, 'o:', o_coh(ij),'a:', a_coh(ij),'k:', k_coh(ij),'b:', b_coh(ij),'c:', c_coh(ij)
                endif

            enddo

        else

            do ij = 1, JJ
                do ie = 0, NE

                    do iw = 1, NW
                        do ir = 1, NR
                            y_coh(ij) = y_coh(ij) + (w(ir)*eff(ij)*zeta(iw) + pen(ij, ir))*eta(ij, ie) &
                                             *dist_zeta(iw)*phi_e(ij, ie)
                        enddo
                    enddo

                    do ix = 0, NX
                        c_coh(ij) = c_coh(ij) + c(ij, ix)*eta(ij, ie)*phi_X(ij, ix)*phi_e(ij, ie)
                    enddo

                    if(ij > 1)then
                        do ia = 0, NA
                            a_coh(ij) = a_coh(ij) + a(ia)*eta(ij, ie)*phi_a(ij-1, ia)*phi_e(ij, ie)
                            o_coh(ij) = o_coh(ij) + omega_plus(ij-1, ia)*phi_a(ij-1, ia)*phi_e(ij, ie)
                        enddo
                    endif
                enddo
            enddo

        endif

        ! calculate aggregate quantities
        BB = 0d0
        do ij = 1, JJ
            BB = BB + b_coh(ij)*m(ij)
        enddo
        
        if (abs(BB) < sig) then
            CC = 0d0
            LL = 0d0
            AA = 0d0
            KK = 0d0
            BB = 0d0
            workpop = 0d0
            do ij = 1, JJ
                CC = CC + c_coh(ij)*m(ij)
                LL = LL + l_coh(ij)*m(ij)
                AA = AA + a_coh(ij)*m(ij)
                KK = KK + k_coh(ij)*m(ij)
                BB = BB + b_coh(ij)*m(ij)
                if(ij < JR) workpop = workpop + m(ij)
            enddo

            ! damping and other quantities [damping acording to Gauss-Seidel procedure]
            KK = damp*(KK) + (1d0-damp)*KK_old 
!~             LL = damp*LL + (1d0-damp)*LL_old
            II = (n_p)*KK !(n_p+delta)*KK
            YY = TProd_bar * KK ** alpha * LL ** (1d0-alpha)

            ! get difference on goods market
            DIFF = YY-CC-II
        endif
        
print*, 'KK:', KK, 'BB', BB, 'LL:', LL, 'YY:', YY, 'DIFF:', DIFF

        ! calculate variances
        var_c = 0d0
        var_y = 0d0
        var_a = 0d0
        var_o = 0d0

        ! analytical approach or not
        if(analytical)then

            do ij = 1, JJ

                do iw = 1, NW
                    do ir = 1, NR
                        var_y(ij) = var_y(ij) + ((w(ir)*eff(ij)*zeta(iw) + pen(ij, ir)))**2*dist_zeta(iw)
                    enddo
                enddo

                do ix = 0, NX
                    var_c(ij) = var_c(ij) + c(ij, ix)**2*phi_X(ij, ix)
                enddo

                do ia = 0, NA
                    if(ij > 1)then
                        var_a(ij) = var_a(ij) + a(ia)**2*phi_a(ij-1, ia)
                        var_o(ij) = var_o(ij) + omega_plus(ij-1, ia)**2*phi_a(ij-1, ia)
                    endif
                enddo
            enddo

        else

            do ij = 1, JJ
                do ie = 0, NE

                    do iw = 1, NW
                        do ir = 1, NR
                            var_y(ij) = var_y(ij) + ((w(ir)*eff(ij)*zeta(iw) + pen(ij, ir))*eta(ij, ie))**2 &
                                             *dist_zeta(iw)*phi_e(ij, ie)
                        enddo
                    enddo

                    do ix = 0, NX
                        var_c(ij) = var_c(ij) + (c(ij, ix)*eta(ij, ie))**2*phi_X(ij, ix)*phi_e(ij, ie)
                    enddo

                    if(ij > 1)then
                        do ia = 0, NA
                            var_a(ij) = var_a(ij) + (a(ia)*eta(ij, ie))**2*phi_a(ij-1, ia)*phi_e(ij, ie)
                            var_o(ij) = var_o(ij) + omega_plus(ij-1, ia)**2*phi_a(ij-1, ia)*phi_e(ij, ie)
                        enddo
                    endif
                enddo
            enddo
        endif

        var_c = var_c - c_coh**2
        var_y = var_y - y_coh**2
        var_a = var_a - a_coh**2
        var_o = var_o - o_coh**2

        ! in case of analytical, add extra eta effect
        if(analytical)then

            ! get age dependent variance of eta
            sigma_eta = sigma_eps*(/(dble(min(ij, JR-1)-1), ij=1,JJ)/)

            ! calculate age specific expectations and variance of exp(eta)
            mu_exp = exp(0.5d0*sigma_eta)
            sigma_exp = exp(sigma_eta)*(exp(sigma_eta)-1d0)

            ! add variance effects
            var_y = mu_exp**2*var_y + sigma_exp*y_coh**2 + sigma_exp*var_y
            var_c = mu_exp**2*var_c + sigma_exp*c_coh**2 + sigma_exp*var_c
            var_a = mu_exp**2*var_a + sigma_exp*a_coh**2 + sigma_exp*var_a

            ! add level effect to averages
            y_coh = mu_exp*y_coh
            c_coh = mu_exp*c_coh
            a_coh = mu_exp*a_coh
        endif

        ! calculate coefficients of variation
        cv_y = sqrt(max(var_y, 0d0))/max(y_coh, 1d-10)
        cv_c = sqrt(max(var_c, 0d0))/max(c_coh, 1d-10)
        cv_a = sqrt(max(var_a, 0d0))/max(a_coh, 1d-10)
        cv_o = sqrt(max(var_o, 0d0))/max(o_coh, 1d-10)

        if(.not. analytical)call calculate_quantiles

    end subroutine


    ! generates the eta distribution
    subroutine generate_eta()

        implicit none
        integer :: ij, iel, ier, ie, is, ir!, isr
        real*8 :: varphi, eta_temp

        ! set bounds and grid for working ages
        eta_l(1)  = 0d0
        eta_u(1)  = 0d0
        eta(1, :) = 0d0

        do ij = 2, JR-1
            eta_l(ij)  = (ij-1)*(minval(log(eps)))
            eta_u(ij)  = (ij-1)*(maxval(log(eps)))
            call grid_Cons_Equi(eta(ij, :), eta_l(ij), eta_u(ij))
        enddo

        ! no innovations throughout retirement
        do ij = JR, JJ
            eta_l(ij)  = eta_l(JR-1)
            eta_u(ij)  = eta_u(JR-1)
            eta(ij, :) = eta(JR-1, :)
        enddo

        phi_e = 0d0

        ! initial distribution at age 1
        phi_e(1, :) = 1d0/dble(NE+1)

        ! iterate over working different years
        do ij = 2, JR-1

            ! last period's etas
            do ie = 0, NE

                ! new innovations
                do is = 1, NR

                    ! distribute on new grid
                    eta_temp = eta(ij-1, ie) + log(eps(is))
                    call linint_Equi(eta_temp, eta_l(ij), eta_u(ij), NE, iel, ier, varphi)
                    phi_e(ij, iel) = phi_e(ij, iel) + dist_eps(is)*varphi*phi_e(ij-1, ie)
                    phi_e(ij, ier) = phi_e(ij, ier) + dist_eps(is)*(1d0-varphi)*phi_e(ij-1, ie)
                enddo
            enddo
        enddo

        ! no innovations throughout retirement
        do ij = JR, JJ
            phi_e(ij, :) = phi_e(JR-1, :)
        enddo

        ! take exponentials
        eta = exp(eta)

    end subroutine


    ! subroutine to calculate age specific quantiles of the distribution
    subroutine calculate_quantiles()

        use toolbox

        implicit none
        integer :: ij, ie, ia, ic, it, NC
        real*8 :: a_sort((NA+1)*(NE+1)), a_dist((NA+1)*(NE+1)), a_cdist((NA+1)*(NE+1))
        integer :: iorder((NA+1)*(NE+1))
        real*8 :: thresholds(5), quantiles(size(thresholds, 1), JJ), slope, ages(JJ)

        ! define quantile thresholds
        thresholds = (/0.05d0, 0.25d0, 0.50d0, 0.75d0, 0.95d0/)
        quantiles = 0d0

        ! iterate over ages
        do ij = 2, JJ

            a_sort = 0d0
            a_dist = 0d0

            ! copy savings into one-dimensional array
            ic = 1
            do ie = 0, NE
                do ia = 0, NA
                    if(phi_a(ij-1, ia)*phi_e(ij, ie) > 1d-12)then
                        a_sort(ic) = a(ia)*eta(ij, ie)
                        a_dist(ic) = phi_a(ij-1, ia)*phi_e(ij, ie)
                        ic = ic + 1
                    endif
                enddo
            enddo
            NC = ic -1

            ! sort array and distribution
            call sort(a_sort(1:NC), iorder(1:NC))

            ! calculate cumulative distribution (attention ordering)
            a_cdist(1) = a_dist(iorder(1))
            do ic = 2, NC
                a_cdist(ic) = a_cdist(ic-1) + a_dist(iorder(ic))
            enddo

            ! get quantiles
            do it = 1, size(thresholds, 1)
                if(thresholds(it) <= a_cdist(1))then
                    quantiles(it, ij) = a_sort(1)
                else
                    do ic = 2, NC
                        if(thresholds(it) < a_cdist(ic))then
                            slope = (a_sort(ic)-a_sort(ic-1))/(a_cdist(ic)-a_cdist(ic-1))
                            quantiles(it, ij) = a_sort(ic-1) + slope*(thresholds(it)-a_cdist(ic-1))
                            exit
                        elseif(ic == NC)then
                            quantiles(it, ij) = a_sort(NC)
                        endif
                    enddo
                endif
            enddo

        enddo

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        call plot(ages, quantiles(1, :), legend='5%')
        call plot(ages, quantiles(2, :), legend='25%')
        call plot(ages, quantiles(3, :), legend='50%')
        call plot(ages, quantiles(4, :), legend='75%')
        call plot(ages, quantiles(5, :), legend='95%')

        call execplot(title='Quantiles of Asset Distributions', &
            xlabel='Age j', ylabel='Assets')

    end subroutine

    ! subroutine for updating bond return
    subroutine bond_return()
        implicit none
        real*8 :: rb_temp

print*, 'rb1', rb
        if (iterb == 1 .and. BB > 0d0 .and. abs(BB)>= sig) then
            rb_a = rb
            rb_b = rb-0.5
print*, 'rb2', rb
!~             rb_b = rb*(1d0-Damp_rb)
print*, 'rb_b', rb_b

            BB_a = BB
            rb = rb_b
        elseif (iterb == 1 .and. BB < 0d0 .and. abs(BB)>= sig) then
            rb_a = rb
            rb_b = rb+0.5
!~             rb_b = rb*(1d0+Damp_rb)
            BB_a = BB
            rb = rb_b
        elseif (iterb == 2) then

            BB_b = BB
            if(BB_a*BB_b >= 0d0)then
                stop 'Error: There is no root in [rb_old, rb_new]'
            endif
            rb_c = (rb_a+rb_b)/2d0
            rb = rb_c
        else
            BB_c = BB
            rb_c = rb
            if (BB_a*BB_c < 0d0) then
                rb_b = rb_c
                BB_b = BB_c
                rb = (rb_a+rb_b)/2d0
            else
                rb_a = rb_c
                BB_a = BB_c
                rb = (rb_a+rb_b)/2d0
            endif
        endif
   
print*, 'rb_a:', rb_a, 'rb_b:', rb_b, 'rb_c:', rb_c, 'rb', rb
print*, 'BB_a:', BB_a, 'BB_b:', BB_b, 'BB_c:', BB_c, 'BB', BB

    end subroutine

    ! subroutine for calculating government parameters
    subroutine government()

        implicit none
        integer :: ij, ir, iq, ig, iw
        
        ! calculate pension
        pen(:,:) = 0d0
        do ir = 1, NR
            pen(JR:JJ, ir) = kappa*w(ir)*eff(JR-1)
        enddo

        ! get total pension spending
        do ij = 1, JJ
            do ir = 1, NR
                total_pen = total_pen + m(ij)*dist_Tprod(ir)*pen(ij, ir)
            enddo
        enddo

        ! calculate total working income
        do ij = 1, JJ
            do iw = 1, NW
                do ir = 1, NR
                    total_INC = total_INC + (w(ir)*eff(ij)*zeta(iw))&
                                *m(ij)*dist_zeta(iw)*dist_Tprod(ir)
                enddo
            enddo
        enddo
        

        ! calculate budget-balancing income tax rate
        tauw = total_pen/total_INC
print*, 'tauw', tauw
    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ixmax(JJ), iamax(JJ), ages(JJ)

        ! check for the maximium grid points used
        call check_grid_X(ixmax)
        call check_grid_a(iamax)

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)*5

        write(21, '(a,a)')' IJ      CONS    INCOME    ASSETS', &
            '     OMEGA     Capital    Bond    CV(C)     CV(Y)     CV(A)     CV(O)     IXMAX     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,10f10.3,2i10)')ages(ij), c_coh(ij), y_coh(ij), a_coh(ij), o_coh(ij), k_coh(ij), &
                    b_coh(ij), cv_c(ij), cv_y(ij), cv_a(ij), cv_o(ij), ixmax(ij), iamax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption (Mean)')
        call plot(dble(ages), y_coh, legend='Labor Income (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income', ylim=(/0d0, 4d0/))

        call plot(dble(ages(2:JJ)), o_coh(2:JJ))
        call execplot(xlabel='Age j', ylabel='Portfolio Share', ylim=(/0d0, 1d0/))

        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

    end subroutine


    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                if(phi_a(ij, ia) > 1d-8)iamax(ij) = ia
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_X(ixmax)

        implicit none
        integer :: ixmax(JJ), ij, ix

        ixmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                if(phi_X(ij, ix) > 1d-8)ixmax(ij) = ix
            enddo
        enddo

    end subroutine

end program
