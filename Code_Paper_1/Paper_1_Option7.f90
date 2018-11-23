!##############################################################################
! PROGRAM SOLG_LR
!
! Paper 1 Option 7
! Build on the code of Hans Fehr and Fabian Kindermann
!
!##############################################################################
include "Paper_1m_Option7.f90"

program SOLG_LR

    use globals

    implicit none

    ! initialize variables
    call initialize()

    open(21, file='output.out')

    ! calculate initial equilibrium
    call get_SteadyState()

    close(21)

contains

!##############################################################################

    ! computes the initial steady state of the economy
    subroutine get_SteadyState()

        implicit none
        integer :: iter

        ! start timer
        call tic()

        ! iterate until value function converges
        do iter = 1, itermax

            ! derive prices
            call prices()

            ! solve the household problem
            call solve_household()

            ! calculate the distribution of households over state space
            call get_distribution()

            ! aggregate individual decisions over cohorts
            call aggregation()

            ! update bond market return
            call bond_return()

            ! determine the government parameters
            call government()

            write(*,'(i4,6f8.2,f12.5)')iter, (/5d0*KK, CC, II/)/YY*100d0, &
                                       rb, rk, w, DIFF/YY*100d0
            if(abs(DIFF/YY)*100d0 < sig)then
                call toc
                call output()
                return
            endif
        enddo

        call toc
        call output()

        write(*,'(a/)')'ATTENTION: NO CONVERGENCE !!!'

    end subroutine

!##############################################################################

    ! initializes the remaining model parameters and variables
    subroutine initialize

        implicit none
        integer :: ij, ip, is, ib, ik, iv

        write(*,'(/a/)')'INITIAL EQUILIBRIUM'
        write(*,'(a)')'ITER     K/Y     C/Y     I/Y       rb       rk       w        DIFF'

        ! set up population structure
        do ij = 1, JJ
            m(ij) = (1d0+n_p)**(1d0-ij)  
        enddo

        ! set survival probabilities
        psi = (/1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, &
                0.99906d0, 0.99908d0, 0.99906d0, 0.99907d0, 0.99901d0, &
                0.99899d0, 0.99896d0, 0.99893d0, 0.99890d0, 0.99887d0, &
                0.99886d0, 0.99878d0, 0.99871d0, 0.99862d0, 0.99853d0, &
                0.99841d0, 0.99835d0, 0.99819d0, 0.99801d0, 0.99785d0, &
                0.99757d0, 0.99735d0, 0.99701d0, 0.99676d0, 0.99650d0, &
                0.99614d0, 0.99581d0, 0.99555d0, 0.99503d0, 0.99471d0, &
                0.99435d0, 0.99393d0, 0.99343d0, 0.99294d0, 0.99237d0, &
                0.99190d0, 0.99137d0, 0.99085d0, 0.99000d0, 0.98871d0, &
                0.98871d0, 0.98721d0, 0.98612d0, 0.98462d0, 0.98376d0, &
                0.98226d0, 0.98062d0, 0.97908d0, 0.97682d0, 0.97514d0, &
                0.97250d0, 0.96925d0, 0.96710d0, 0.96330d0, 0.95965d0, &
                0.95619d0, 0.95115d0, 0.94677d0, 0.93987d0, 0.93445d0, &
                0.92717d0, 0.91872d0, 0.91006d0, 0.90036d0, 0.88744d0, &
                0.87539d0, 0.85936d0, 0.84996d0, 0.82889d0, 0.81469d0, &
                0.79705d0, 0.78081d0, 0.76174d0, 0.74195d0, 0.72155d0, &
                0.00000d0/)

        ! initialize age earnings process
        eff(1:JR-1) = &
            (/1.0000d0, 1.0719d0, 1.1438d0, 1.2158d0, 1.2842d0, 1.3527d0, &
              1.4212d0, 1.4897d0, 1.5582d0, 1.6267d0, 1.6952d0, 1.7217d0, &
              1.7438d0, 1.7748d0, 1.8014d0, 1.8279d0, 1.8545d0, 1.8810d0, &
              1.9075d0, 1.9341d0, 1.9606d0, 1.9623d0, 1.9640d0, 1.9658d0, &
              1.9675d0, 1.9692d0, 1.9709d0, 1.9726d0, 1.9743d0, 1.9760d0, &
              1.9777d0, 1.9700d0, 1.9623d0, 1.9546d0, 1.9469d0, 1.9392d0, &
              1.9315d0, 1.9238d0, 1.9161d0, 1.9084d0, 1.9007d0, 1.8354d0, &
              1.7701d0, 1.7048d0/)
        eff(JR:JJ) = 0d0

        ! discretize zeta shocks
        call normal_discrete(zeta, dist_zeta, 0d0, sigma_zeta)
        zeta = exp(zeta)

        ! calculate the shock process for aggregate productivity (TProd)
        call discretize_AR(rho, TProd_bar, sigma_vtheta, TProd, pi_TProd)
        TProd = exp(TProd) ! Do we need exp()????

        ! calculate the shock process for eta
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi_eta)
        eta = exp(eta)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)
        call grid_Cons_Grow(omega, omega_l, omega_u, omega_grow)
        
        ! initialize tax and transfers
        tauw  = 0.0d0

        ! initial guesses for macro variables
        KK = 1d0
        BB = 0d0
        LL = 1d0
        YY = 1d0
        II = (n_p+delta)*KK

        ! open files
        open(21, file='output.out')

    end subroutine

!##############################################################################

    ! subroutine for calculating prices
    subroutine prices()

        implicit none
        integer :: iv, ij, ijj
        real*8 :: abor_temp
        
        do iv = 1, NR
            rk(iv) = TProd(iv)*alpha*(KK/LL)**(alpha-1d0)-delta 
        enddo
        
        if (iter==1) then
            rb = 0
            do iv = 1, NR
                rb = rb + (rk(iv)-delta)/NR
            enddo
        endif
        
        ! calculate wage rate
        do iv = 1, NR
            w(iv) = TProd(iv)*(1d0-alpha)*(KK/LL)**alpha
        enddo
        
        ! calculate after-tax wage rate
        do iv = 1, NR
            wn(iv) = w(iv)*(1d0-tauw)
        enddo
        
        ! old-age transfers
        pen(:,:) = 0d0
        do iv = 1, NR
            pen(JR:JJ, iv) = kappa*w(iv)*eff(JR-1)
        enddo
        
        ! endogenous lower bound for asset
        do ij = 1, JJ
            ! calculate natural borrowing limit
            abor_temp = 0d0
            do ijj = JJ, ij+1, -1
                abor_temp = abor_temp/(1d0+rb) + wn(1)*eff(ijj) &
                            *exp(eta(1) + zeta(1))+ pen(ijj, 1)
            enddo
                                
            ! set maximum of natural and exogenous borrowing constraint
            a_bor(ij) = -abor_temp/(1d0+rb)
        enddo

    end subroutine

!##############################################################################

    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ia, io, iq, ig, iv
        real*8 :: k_in, b_in, wage, available
        real*8 :: a_in(2)
        logical :: check
        
        ! get decision in the last period of life
        do ia = 1, NA
            do io = 1, NO
                do iq = 1, 1
                    do ig = 1, 1
                        do iv = 1, NR
                            a_plus(JJ, ia, io, iq, ig, iv) = 0d0
                            omega_plus(JJ, ia, iq, iv) = 0d0
                            c(JJ, ia, io, iq, ig, iv) = (1d0 + rb + omega(io)*(rk(iv)-rb))*a(ia) &
                            + pen(JJ, iv)
                            V(JJ, ia, io, iq, ig, iv) = valuefunc(0d0, c(JJ, ia, io, iq, ig, iv), JJ, iq, iv)
                        enddo
                    enddo                
                enddo
            enddo
        enddo
       
        ! get decision for other cohorts
        do ij = JJ-1, 1, -1

            ! determine optimal portfolio choice for all others
            do ia = 1, NA
                do iq = 1, NE
                    do iv = 1, NR
                        ! assign omega = 0 for asset <= 0 (borrow from bond)
                        if (a(ia) <= 0) then
                            omega_plus(ij, ia, iq, iv) = 0.0d0
                        ! solve for 0 < omega <= 1
                        else
                            call solve_portfolio(ij, ia, iq, iv)
                        endif
                    enddo
                enddo
            enddo

            ! interpolate individual RHS and value function
            do iq = 1, NE
                do iv = 1, NR
                    call interpolate(ij, iq, iv)
                enddo
            enddo


            ! determine consumption-savings solution
            do ia = 0, NA
                do io = 0, NO
                    do iq = 0, NE
                        do ig = 0, NW
                            do iv = 0, NR
                                call solve_consumption(ij, ia, io, iq, ig, iv)
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine

!##############################################################################

    ! solve the household's portfolio decision
    subroutine solve_portfolio(ij, ia, iq, iv)

        implicit none
        integer, intent(in) :: ij, ia, iq, iv
        real*8 :: x_in, port0, port1, tolerance
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia
        iq_com = iq
        iv_com = iv

        ! check for corner solutions
        port0 = foc_port(0d0)
        port1 = foc_port(1d0)

        ! use intermediate value theorem
        if(port0*port1 > 0d0)then
            if(abs(port0) > abs(port1))then
                omega_plus(ij, ia, iq, iv) = 1d0
            else
                omega_plus(ij, ia, iq, iv) = 0d0
            endif
            return
        else

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
            if(check)write(*,'(a, 2i4)')'ERROR IN ROOTFINDING PORT : ', ij, ia, iq

            omega_plus(ij, ia, iq, iv) = x_in

            ! reset tolerance level to original value
            call settol_root(1d-8)
        endif

    end subroutine

!##############################################################################

    ! solve the household's consumption-savings decision
    subroutine solve_consumption(ij, ia, io, iq, ig, iv)

        implicit none
        integer, intent(in) :: ij, ia, io, iq, ig, iv
        real*8 :: x_in
        logical :: check

!~         ! determine decision for zero cash-on-hand
!~         if(X(ix) < 1d-10)then
!~             a_plus(ij, ix) = 0d0
!~             c(ij, ix) = 0d0
!~             V(ij, ix) = valuefunc(0d0, 0d0, ij)
!~             return
!~         endif

        ! set up communication variables
        ij_com = ij
        ia_com = ia
        io_com = io
        iq_com = iq
        ig_com = ig
        iv_com = iv

        ! get best initial guess from future period
        x_in = a_plus(ij+1, ia, io, iq, ig, iv)
        check = .false.

        ! solve the household problem using rootfinding
        call fzero(x_in, foc_cons, check)

        ! write screen output in case of a problem
        if(check)write(*,'(a, 2i4)')'ERROR IN ROOTFINDING CONS : ', ij, ia

        ! check for borrowing constraint
        if(x_in < a_bor(ij_com))then
            cons_com = cons_com + x_in - a_bor(ij_com)
            x_in = a_bor(ij_com)
        endif

        ! copy decisions
        a_plus(ij, ia, io, iq, ig, iv) = x_in
        c(ij, ia, io, iq, ig, iv) = cons_com
        V(ij, ia, io, iq, ig, iv) = valuefunc(x_in, cons_com, ij, iq, iv)

    end subroutine


!##############################################################################

    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij, iq, iv)

        implicit none
        integer, intent(in) :: ij, iq, iv
        integer :: ia, ig, iv_next, iq_next
        real*8 :: X_p, c_p, varphi, dist, EV, R_port
        integer :: ioml, iomr

        RHS(ij, :, :, :) = 0d0
        Q(ij, :, :, :) = 0d0
        iq_com = iq
        iv_com = iv

        do ia = 0, NA
            do ig = 1, NW
                do iv_next = 1, NR
                    do iq_next = 1, NE
                        ! get return on the portfolio
                        R_port = 1d0 + rb + omega_plus(ij, ia, iq_next, iv_next)*(rk(iv_next) - rb)
                        
                        ! derive interpolation weights
                        call linint_Grow(omega_plus(ij, ia, iq_next, iv_next), omega_l, &
                        omega_u, omega_grow, NO, ioml, iomr, varphi)
                        
                        ! calculate next-period consumption
                        c_p = varphi*c(ij_com+1, ia_com, ioml, iq_next, ig, iv_next) + &
                          (1d0-varphi)*c(ij_com+1, ia_com, iomr, iq_next, ig, iv_next)
                    
                        ! get distributional weight
                        dist = dist_zeta(ig)*pi_Tprod(iv_com, iv_next)*pi_eta(iq_com, iq_next)
                        
                        ! get RHS of foc and Q
                        RHS(ij, ia, iq_com, iv_com) = RHS(ij, ia, iq_com, iv_com) + &
                                            dist*R_port*margu(c_p)
                        Q(ij, ia, iq_com, iv_com)   = Q(ij, ia, iq_com, iv_com) + &
                                            dist*(varphi*V(ij+1, ia, ioml, iq_next, ig, iv_next)+&
                                            (1-varphi)*V(ij+1, ia, iomr, iq_next, ig, iv_next))**egam/egam
                        !?? why **egam/egam in the line above and (egam*Q(ij, ia))**(1d0/egam) below?
                    enddo
                enddo
            enddo
            
            RHS(ij, ia, iq_com, iv_com) = (beta*psi(ij+1)*RHS(ij, ia, iq_com, iv_com))**(-1/gamma)
            Q(ij, ia, iq_com, iv_com)   = (egam*Q(ij, ia, iq_com, iv_com))**(1d0/egam)

        enddo
    end subroutine                    

!##############################################################################

    ! determines the invariant distribution over state space
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, io, iq, ig, iv
        real*8 :: varphi

        ! set distributions to zero
        phi_ij(:, :, :, :, :, :) = 0d0
        phi_aplus(:, :) = 0d0
        phi_aoep(:, :, :, :, :) = 0d0
        phi_eta(:, :) = 0d0
        phi_Tprod(:, :) = 0d0
        

        ! get initial distribution in age 1
        do ig = 1, NW
            phi_ij(1, 0, 0, iq_initial, ig, iv_initial) = dist_zeta(ig)
        enddo

        do ij = 1, JJ-1

            ! get distribution of next period asset
            call get_distribution_a(ij)
            
            ! get joint-distribution of next-period asset, next-period omega, eta, and Omega
            call get_distribution_aoeO(ij)
            
            ! get distribution cohort ij+1 across all states z(t+1)
            call get_distribution_ij(ij)
        enddo

    end subroutine

!##############################################################################

    ! to calculate distribution of a_plus for cohort ij
    subroutine get_distribution_a(ij)
        implicit none
        integer :: ij, ia, io, iq, ig, iv
        integer :: ial, iar
        real*8 :: varphi_a
        
        ! derive interpolation weights for a+ and assign prob to discrete phi_aplus(ij, ia)
        do ia = 0, NA
            do io = 0, NO
                do iq = 0, NE
                    do ig = 0, NW
                        do iv = 0, NR
                            call linint_Grow(a_plus(ij, ia, io, iq, ig, iv), a_l, a_u, a_grow, &
                                 NA, ial, iar, varphi_a)
                                 
                            phi_aplus(ij, ial) = phi_aplus(ij, ial) + &
                                                 varphi_a*phi_ij(ij, ia, io, iq, ig, iv)
                                                 
                            phi_aplus(ij, iar) = phi_aplus(ij, iar) + &
                                                 (1-varphi_a)*phi_ij(ij, ia, io, iq, ig, iv)
                        enddo
                    enddo
                enddo
            enddo
        enddo
                        
    end subroutine

!##############################################################################

    ! to calculate the joint distribution of a_plus, eta, TProd of cohort ij
    subroutine get_distribution_aoeO(ij)
        implicit none
        integer :: ij, ia, io, iq, ig, iv
        integer ::ioml, iomr
        real*8 :: varphi_o
        
        ! derive distribution of discrete eta and TProd
        do ia = 0, NA
            do io = 0, NO
                do iq = 0, NE
                    do ig = 0, NW
                        do iv = 0, NR
                            phi_eta(ij, iq) = phi_eta(ij, iq) + phi_ij(ij, ia, io, iq, ig, iv)
                            phi_Tprod(ij, iv) = phi_Tprod(ij, iv) + phi_ij(ij, ia, io, iq, ig, iv)
                        enddo
                    enddo
                enddo
            enddo
        enddo    
     
        ! derive interpolation weights for omega_plus and joint-distribute
        do ia = 0, NA
            do iq = 0, NE
                do iv = 0, NR
                    call linint_Grow(omega_plus(ij, ia, iq, iv), omega_l, &
                        omega_u, omega_grow, NA, ioml, iomr, varphi_o)                
                                 
                    phi_aoep(ij, ia, ioml, iq, iv) = phi_aoep(ij, ia, ioml, iq, iv) + &
                                                     varphi_o*phi_eta(ij, iq)*phi_Tprod(ij, iv)&
                                                     *phi_ij(ij, ia, io, iq, ig, iv)
                                         
                    phi_aoep(ij, ia, iomr, iq, iv) = phi_aoep(ij, ia, ioml, iq, iv) + &
                                                     (1-varphi_o)*phi_eta(ij, iq)*phi_Tprod(ij, iv)&
                                                     *phi_ij(ij, ia, io, iq, ig, iv)
                enddo
            enddo
        enddo
        
    end subroutine

!##############################################################################
        
    ! to calculate the distribution of cohort ij+1 across all states z(t+1)
    subroutine get_distribution_ij(ij)
        implicit none
        integer :: ij, ia, io, iq, ig, iv
        
        !distribute from phi_aoep to next period cohort across all states
        do ia = 0, NA
            do io = 0, NO
                do iq = 0, NE
                    do iq_com = 0, NE
                        do ig = 0, NW
                            do iv = 0, NR
                                do iv_com = 0, NR
                                    phi_ij(ij+1, ia, io, iq, ig, iv) = phi_ij(ij+1, ia, io, iq, ig, iv) +&
                                        pi_eta(iq_com, iq)*dist_zeta(ig)*pi_TProd(iv_com, iv)*&
                                        phi_aoep(ij, ia, io, iq_com, iv_com)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        
    end subroutine

!##############################################################################

    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ia, io, iq, ig, iv
        real*8 :: workpop, LL_old, KK_old, BB_old

        ! Store current aggregate quantities
        LL_old = LL 
        KK_old = KK

        ! calculate cohort aggregates
        c_coh(:) = 0d0
        l_coh(:) = 0d0
        y_coh(:) = 0d0
        k_coh(:) = 0d0
        b_coh(:) = 0d0
        v_coh(:) = 0d0

        do ij = 1, JJ
            do ia = 0, NA
                do io = 0, NO
                    do iq = 0, NE
                        do ig = 0, NW
                            do iv = 0, NR
                                c_coh(ij) = c_coh(ij) + phi_ij(ij, ia, io, iq, ig, iv)*&
                                    c(ij, ia, io, iq, ig, iv)
                                l_coh(ij) = l_coh(ij) + phi_ij(ij, ia, io, iq, ig, iv)*&
                                    eff(ij)*exp(eta(iq) + zeta(ig))
                                a_coh(ij) = a_coh(ij) + phi_ij(ij, ia, io, iq, ig, iv)*&
                                    a(ia)
                                b_coh(ij) = b_coh(ij) + phi_ij(ij, ia, io, iq, ig, iv)*&
                                    a(ia)*(1-omega_plus(ij, ia, iq, iv))
                                k_coh(ij) = k_coh(ij) + phi_ij(ij, ia, io, iq, ig, iv)*&
                                    a(ia)*omega_plus(ij, ia, iq, iv)
                                o_coh(ij) = o_coh(ij) + phi_ij(ij, ia, io, iq, ig, iv)*&
                                    omega(io)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! calculate aggregate quantities
        CC = 0d0
        LL = 0d0
        AA = 0d0
        KK = 0d0
        BB = 0d0
        workpop = 0d0
        do ij = 1, JJ
            CC = CC + c_coh(ij)*m(ij)
            LL = LL + l_coh(ij)*m(ij)
            AA = AA + y_coh(ij)*m(ij)
            KK = KK + k_coh(ij)*m(ij)
            BB = BB + b_coh(ij)*m(ij)
            if(ij < JR) workpop = workpop + m(ij)
        enddo

        ! damping and other quantities [damping acording to Gauss-Seidel procedure]
        KK = damp*(KK) + (1d0-damp)*KK_old 
        LL = damp*LL + (1d0-damp)*LL_old
        II = (n_p+delta)*KK
        YY = TProd_bar * KK ** alpha * LL ** (1d0-alpha)

        ! get average income and average working hours
        do ij = 1, JJ
            do ia = 0, NA
                do io = 0, NO
                    do iq = 0, NE
                        do ig = 0, NW
                            do iv = 0, NR
                                INC = INC + w(iv)*phi_ij(ij, ia, io, iq, ig, iv)*&
                                    eff(ij)*exp(eta(iq) + zeta(ig))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        INC = INC/workpop ! average income
!~         HH  = HH/workpop 

        ! get difference on goods market
        DIFF = YY-CC-II
    
    end subroutine

!##############################################################################

    ! subroutine for updating bond return
    subroutine bond_return()
        implicit none
        
        if (BB > 0d0) then
            rb = rb*0.9d0
        elseif (BB < 0d0) then
            rb = rb*1.1d0
        endif
        
    end subroutine

!##############################################################################

    ! subroutine for calculating government parameters
    subroutine government()

        implicit none
        integer :: ij, iv, iq, ig
        
        ! calculate pension
        pen(:,:) = 0d0
        do iv = 1, NR
            pen(JR:JJ, iv) = kappa*w(iv)*eff(JR-1)
        enddo

        ! get total pension spending
        do ij = 1, JJ
            do iv = 1, NR
                total_pen = total_pen + m(ij)*pen(ij, iv)*phi_Tprod(ij, iv)
            enddo
        enddo
        
        ! calculate total working income
        do ij = 1, JJ
            do iq = 0, NE
                do ig = 0, NW
                    do iv = 0, NR
                        total_INC = total_INC + m(ij)*phi_Tprod(ij, iv)*phi_eta(ij, iq)&
                                    *dist_zeta(ig)*w(iv)*eff(ij)*exp(eta(iq) + zeta(ig))
                    enddo
                enddo
            enddo
        enddo

        ! calculate budget-balancing income tax rate
        tauw = total_pen/total_INC

    end subroutine

!##############################################################################

    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, iamax(JJ), ages(JJ)

        ! check for the maximium grid points used
        call check_grid_a(iamax)

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        write(21, '(a,a)')' IJ      CONS    INCOME    ASSETS', &
            '     OMEGA     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,4f10.3,i10)')ages(ij), c_coh(ij), y_coh(ij), a_coh(ij), o_coh(ij), &
                    iamax(ij)
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

!##############################################################################

    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                if(phi_aplus(ij, ia) > 1d-8)iamax(ij) = ia
            enddo
        enddo

    end subroutine

end program
