!##############################################################################
! PROGRAM SOLG_LR
!
! ## Long-run equilibria in the stochastic OLG model
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

include "Paper_1m_Option5.f90"

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

            write(*,'(i4, 6f8.2)')iter, (/5d0*KK, CC, II/)/YY*100d0, &
                                       rb, rk(4), w(4)!, DIFF/YY*100d0

!~             write(*,'(i4,6f8.2,f12.5)')iter, (/5d0*KK, CC, II/)/YY*100d0, &
!~                                        rb, rk, w, DIFF/YY*100d0

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
        integer :: ij, ib, ik, iq, ig, iv

        write(*,'(/a/)')'INITIAL EQUILIBRIUM'
        write(*,'(a)')'ITER     K/Y     C/Y     I/Y       rb       rk       w        DIFF'

        ! set up population structure
        do ij = 1, JJ
            m(ij) = (1d0+n_p)**(1d0-ij) 
        enddo

        ! initialize asset grid
        call grid_Cons_Grow(k, k_l, k_u, k_grow)
        call grid_Cons_Grow(b, b_l, b_u, b_grow)

        ! get initial guess for savings decision
        do ij = 1, JJ
            do ik = 0, NK
                do ib = 0, NB
                    do iq = 1,NE
                        do ig = 1,NW
                            do iv = 1, NR
                                kplus(ij, :, ib, iq, ig, iv) = max(k(:)/2d0, k(1)/2d0) 
                                !!why include k(1)/2d0? isnt it already in k(:)/2d0? 
                                bplus(ij, ik, :, iq, ig, iv) = max(b(:)/2d0, b(1)/2d0)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! set survival probabilities (interpolate linearly from 1.00 and 0.72)
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



        ! calculate the shock process
!~         call normal_discrete(zeta, dist_zeta, 0d0, sigma_zeta)
!~         zeta = exp(zeta)
        zeta(1) = 0d0
        dist_zeta(1) = 1d0
        
        ! calculate the shock process for aggregate productivity (TProd)
        call discretize_AR(rho, TProd_bar, sigma_vtheta, TProd, pi_TProd)
        TProd=exp(TProd)

!~         ! calculate the shock process for eta
!~         call discretize_AR(rho, 0d0, sigma_eps, eta, pi_eta)
!~         eta = exp(eta)
        eta(1) = 0d0
        pi_eta(1,1) = 1d0
        
        ! initialize tax and transfers
        tauw  = 0.0d0

        ! initial guesses for macro variables
        KK = 1d0
        BB = 0d0
        LL = 1d0
        YY = 1d0
        II = (n_p+delta)*KK



    end subroutine

!##############################################################################

    ! subroutine for calculating prices
    subroutine prices()

                implicit none
        integer :: iv, ij, ijj, ig, iq
        real*8 :: abor_temp

        ! calculate wage rate
        do iv = 1, NR
            w(iv) = TProd(iv)*(1d0-alpha)*(KK/LL)**alpha
        enddo

        ! calculate return on risky asset
        do iv = 1, NR
            rk(iv) = (TProd(iv)*KK**alpha*LL**(1d0-alpha)-w(iv)*LL)/KK
        enddo
       
        if (iter==1) then
            rb = 0
            do iv = 1, NR
                rb = rb + (rk(iv)-mu_r)/NR
            enddo
        endif
       
        ! calculate after-tax wage rate
        do iv = 1, NR
            wn(iv) = w(iv)*(1d0-tauw)
        enddo
        
        ! old-age transfers
        pen(:,:) = 0d0
        do iv = 1, NR
            pen(JR:JJ, iv) = max(kappa*w(iv)*eff(JR-1), 1d-10)
        enddo
        
        ! endogenous lower bound for asset borrowing
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
        integer :: ij, ik, ib, ig, ig_max, iq, iq_max, iv
        real*8 :: k_in, b_in, wage, available
        real*8 :: a_in(2)
        logical :: check
        
        ! get decision in the last period of life
        do ik = 0, NK
            do ib = 0, NB 
                do iv = 1, NR
                    kplus(JJ, ik, ib, :, :, :) = 0d0
                    bplus(JJ, ik, ib, :, :, :) = 0d0
                    c(JJ, ik, ib, :, :, :) = max((1d0+rk(iv))*k(ik) + (1d0+rb)*b(ib) + pen(JJ, iv), 1d-10)
                    V(JJ, ik, ib, :, :, :) = valuefunc(0d0, 0d0, c(JJ, ik, ib, 1, 1, iv), JJ, 1, 1, iv)
                enddo
            enddo
        enddo

        ! interpolate individual RHS
        call interpolate_b(JJ)
        call interpolate_k(JJ)

        ! get decision for other cohorts
        do ij = JJ-1, 1, -1

            ! check about how many is to iterate
            if(ij >= JR)then
                iq_max = 1
                ig_max = 1
            else
                iq_max = NE
                ig_max = NW
            endif

            do ik = 0, NK
                do ib = 0, NB

                    ! determine decision for zero assets at retirement without pension
                    if(ij >= JR .and. ik == 0 .and. ib == 0 .and. kappa <= 1d-10)then
                        kplus(ij, ik, ib, :, :, :) = 0d0
                        bplus(ij, ik, ib, :, :, :) = 0d0
                        c(ij, ik, ib, :, :, :) = 0d0
                        V(ij, ik, ib, :, :, :) = valuefunc(0d0, 0d0, 0d0, ij, 1, 1, iv)
                        cycle
                    endif

                    do iq = 1, iq_max
                        do ig = 1, ig_max
                            do iv = 1 , NR
                                ! get initial guess for the individual choices
                                k_in = kplus(ij+1, ik, ib, iq, ig, iv)
                                b_in = bplus(ij+1, ik, ib, iq, ig, iv)
                                
                                a_in(:) = (/k_in, b_in/)
                                
                                ! set up communication variables
                                ij_com = ij
                                ik_com = ik
                                ib_com = ib
                                iq_com = iq
                                ig_com = ig
                                iv_com = iv

                                call fzero(a_in, foc, check)
                                
                                ! write screen output in case of a problem
                                if(check)write(*,'(a, 6i4)')'ERROR IN ROOTFINDING : ', ij, ik, ib, iq, ig, iv

                                ! check for borrowing constraint
                                if(k_in < 0d0)then
                                    k_in = 0d0
                                endif

                                ! copy decisions
                                kplus(ij, ik, ib, iq, ig, iv) = k_in
                                bplus(ij, ik, ib, iq, ig, iv) = b_in
                                c(ij, ik, ib, iq, ig, iv) = cons_com
                                V(ij, ik, ib, iq, ig, iv) = valuefunc(k_in, b_in, cons_com, ij, iq, ig, iv)
                            enddo
                        
                            ! copy decision in retirement age
                            if(ij >= JR)then
                                do iv = 1, NR
                                    kplus(ij, ik, ib, :, :, iv) = kplus(ij, ik, ib, 1, 1, iv)
                                    bplus(ij, ik, ib, :, :, iv) = bplus(ij, ik, ib, 1, 1, iv)
                                    c(ij, ik, ib, :, :, iv) = c(ij, ik, ib, 1, 1, iv)
                                    V(ij, ik, ib, :, :, iv) = V(ij, ik, ib, 1, 1, iv)
                                enddo
                            endif
                        
                        enddo
                    enddo
                enddo
            enddo

            ! interpolate individual RHS
            call interpolate_b(ij)
            call interpolate_k(ij)
        enddo

    end subroutine

!##############################################################################

    ! for calculating the rhs of the first order condition for bond at age ij
    subroutine interpolate_b(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ik, ib, iq, iq_p, ig, iv, iv_p
        real*8 :: chelp

        do ik = 0, NK
            do ib = 0, NB
                do iq = 1, NE
                    do ig = 1, NW
                        do iv = 1, NR

                            ! calculate the RHS of the first order condition of b+
                            RHS_b(ij, ik, ib, iq, ig, iv) = 0d0
                            EV(ij, ik, ib, iq, ig, iv) = 0d0
                            do iq_p = 1, NE
                                do iv_p = 1, NR
                                    
                                    chelp = max(c(ij, ik, ib, iq_p, ig, iv_p),1d-10)
                                    RHS_b(ij, ik, ib, iq, ig, iv) = RHS_b(ij, ik, ib, iq, ig, iv) + &
                                        pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*margu(chelp)
                                    EV(ij, ik, ib, iq, ig, iv)  = EV(ij, ik, ib, iq, ig, iv) + &
                                        pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*V(ij, ik, ib, iq_p, ig, iv_p)
                                enddo
                            enddo
                            
                            RHS_b(ij, ik, ib, iq, ig, iv) = (1d0+rb)*beta*RHS_b(ij, ik, ib, iq, ig, iv)**(-gamma)
                            EV(ij, ik, ib, iq, ig, iv) = (egam*EV(ij, ik, ib, iq, ig, iv))**(1d0/egam)
                            ! we transform EV for interpolation accuracy (see Fehr's p414)
                        
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine

!##############################################################################

! for calculating the rhs of the first order condition for capital at age ij
    subroutine interpolate_k(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ik, ib, iq, iq_p, ig, iv, iv_p
        real*8 :: chelp

        do ik = 0, NK
            do ib = 0, NB
                do iq = 1, NE
                    do ig = 1, NW
                        do iv = 1, NR

                            ! calculate the RHS of the 2nd first order condition
                            ! RHSN = numerator term, RHSD = denominator term
                            RHS_k(ij, ik, ib, iq, ig, iv) = 0d0
                            do iq_p = 1, NE
                                do iv_p = 1, NR
                                
                                    chelp = max(c(ij, ik, ib, iq_p, ig, iv_p),1d-10)
                                    RHS_k(ij, ik, ib, iq, ig, iv) = RHS_k(ij, ik, ib, iq, ig, iv) + &
                                        pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*margu(chelp)*(1+rk(iv_p)-delta)
                                
                                enddo
                            enddo
                            
                            RHS_k(ij, ik, ib, iq, ig, iv) = beta*RHS_k(ij, ik, ib, iq, ig, iv)**(-gamma)
                        
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine

!##############################################################################

    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ik, ib, iq, iq_p, ig, iv, iv_p, ikl, ikr, ibl, ibr
        real*8 :: varphib, varphik

        ! set distribution to zero
        phi(:, :, :, :, :, :) = 0d0

        ! get initial distribution in age 1
        do ig = 1, NW
            phi(1, 0, 0, iq_initial, ig, iv_initial) = dist_zeta(ig)
        enddo
        
        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ik = 0, NK
                do ib = 0, NB
                    do iq = 1, NE
                        do ig = 1, NW
                            do iv = 1, NR
                        
                                ! interpolate yesterday's savings decision
                                call linint_Grow(kplus(ij, ik, ib, iq, ig, iv), k_l, k_u, k_grow, NK, ikl, ikr, varphik)
                                call linint_Grow(bplus(ij, ik, ib, iq, ig, iv), b_l, b_u, b_grow, NB, ibl, ibr, varphib)

                                ! restrict values to grid just in case
                                ikl = min(ikl, NK)
                                ikr = min(ikr, NK)
                                ibl = min(ibl, NB)
                                ibr = min(ibr, NB)
                                varphik = min(varphik, 1d0)
                                varphib = min(varphib, 1d0)

                                ! redistribute households
                                do iq_p = 1, NE
                                    do iv_p = 1, NR
                                        phi(ij, ikl, ibl, iq_p, ig, iv_p) &
                                            = phi(ij, ikl, ibl, iq_p, ig, iv_p) +&
                                            pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*varphik*varphib*&
                                            phi(ij-1, ik, ib, iq, ig, iv)
                                            
                                        phi(ij, ikl, ibr, iq_p, ig, iv_p) &
                                            = phi(ij, ikl, ibr, iq_p, ig, iv_p) +&
                                            pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*varphik*(1-varphib)*&
                                            phi(ij-1, ik, ib, iq, ig, iv)
                                            
                                        phi(ij, ikr, ibl, iq_p, ig, iv_p) & 
                                            = phi(ij, ikr, ibl, iq_p, ig, iv_p) +&
                                            pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*(1d0-varphik)*varphib*&
                                            phi(ij-1, ik, ib, iq, ig, iv)
                                            
                                        phi(ij, ikr, ibr, iq_p, ig, iv_p) &
                                            = phi(ij, ikr, ibr, iq_p, ig, iv_p) +&
                                            pi_eta(iq, iq_p)*pi_TProd(iv, iv_p)*(1-varphik)*(1-varphib)*&
                                            phi(ij-1, ik, ib, iq, ig, iv)
                                    enddo
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
        integer :: ij, ik, ib, iq, ig, iv
        real*8 :: workpop, LL_old, KK_old, BB_old

        ! Store current aggregate quantities
        LL_old = LL 
        KK_old = KK
        BB_old = BB

        ! calculate cohort aggregates
        c_coh(:) = 0d0
        l_coh(:) = 0d0
        y_coh(:) = 0d0
        k_coh(:) = 0d0
        b_coh(:) = 0d0
        v_coh(:) = 0d0

        do ij = 1, JJ
            do ik = 0, NK
                do ib = 0, NB
                    do iq = 1, NE
                        do ig = 1, NW
                            do iv = 1, NR
                            
                            c_coh(ij) = c_coh(ij) + c(ij, ik, ib, iq, ig, iv)*phi(ij, ik, ib, iq, ig, iv)
                            l_coh(ij) = l_coh(ij) + eff(ij)*exp(eta(iq) + zeta(ig))*phi(ij, ik, ib, iq, ig, iv)
                            k_coh(ij) = k_coh(ij) + k(ik)*phi(ij, ik, ib, iq, ig, iv)
                            b_coh(ij) = b_coh(ij) + b(ib)*phi(ij, ik, ib, iq, ig, iv)
                            v_coh(ij) = v_coh(ij) + V(ij, ik, ib, iq, ig, iv)*phi(ij, ik, ib, iq, ig, iv)
                        
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! calculate aggregate quantities
        CC = 0d0
        LL = 0d0
        KK = 0d0
        BB = 0d0
        workpop = 0d0
        do ij = 1, JJ
            CC = CC + c_coh(ij)*m(ij)
            LL = LL + l_coh(ij)*m(ij)
            KK = KK + k_coh(ij)*m(ij)
            BB = BB + b_coh(ij)*m(ij)
            if(ij < JR) workpop = workpop + m(ij)
        enddo

        ! damping and other quantities [damping acording to Gauss-Seidel procedure]
        KK = damp*(KK) + (1d0-damp)*KK_old !check this
        LL = damp*LL + (1d0-damp)*LL_old
        II = (n_p+delta)*KK
        YY = TProd_bar * KK ** alpha * LL ** (1d0-alpha)

        ! get average income and average working hours
        do ij = 1, JJ
            do ik = 0, NK
                do ib = 0, NB
                    do iq = 1, NE
                        do ig = 1, NW
                            do iv = 1, NR
                                INC = INC + w(iv)*phi(ij, ik, ib, iq, ig, iv)*&
                                    eff(ij)*exp(eta(iq) + zeta(ig))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo


        ! get average income and average working hours
        INC = INC/workpop ! average labour earning

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
        integer :: ik, ib, ij, iv, iq, ig
        
        ! calculate pension
        pen(:,:) = 0d0
        do iv = 1, NR
            pen(JR:JJ, iv) = kappa*w(iv)*eff(JR-1)
        enddo

        ! derive distribution of eta and TProd for age j
        do ij = 1, JJ
            do ik = 0, NK
                do ib = 0, NB
                    do iq = 1, NE
                        do ig = 1, NW
                            do iv = 1, NR
                                phij_eta(ij, iq) = phij_eta(ij, iq) + phi(ij, ik, ib, iq, ig, iv)
                                phij_Tprod(ij, iv) = phij_Tprod(ij, iv) + phi(ij, ik, ib, iq, ig, iv)
                            enddo
                        enddo
                    enddo
                enddo
            enddo  
        enddo

        ! get total pension spending
        do ij = 1, JJ
            do iv = 1, NR
                total_pen = total_pen + m(ij)*pen(ij, iv)*phij_Tprod(ij, iv)
            enddo
        enddo
        
        ! calculate total working income
        do ij = 1, JJ
            do iq = 1, NE
                do ig = 1, NW
                    do iv = 1, NR
                        total_INC = total_INC + m(ij)*phij_Tprod(ij, iv)*phij_eta(ij, iq)&
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

!~         ! check for the maximium grid points used
!~         call check_grid_a(iamax)

        ! set up age variable
        ages = 20 + (/(ij*5, ij=1,JJ)/) !need to adjust this

        write(21, '(a,a)')' IJ      CONS    labour    ASSETS', &
            '     BONDS' !add income y later
        do ij = 1, JJ
            write(21,'(i3,4f10.3,i10)')ages(ij), c_coh(ij), l_coh(ij), k_coh(ij), b_coh(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption (Mean)')
        call plot(dble(ages), l_coh, legend='Labor (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income', ylim=(/0d0, 4d0/))

        call plot(dble(ages), k_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

        call plot(dble(ages), b_coh)
        call execplot(xlabel='Age j', ylabel='Bonds')

    end subroutine

end program
