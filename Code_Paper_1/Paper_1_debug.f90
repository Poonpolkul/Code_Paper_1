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
include "Paper_1m_debug.f90"

program SOLG_LR

    use globals

    implicit none

    ! initialize variables
    call initialize()

    open(21, file='output_initial.out')

    ! calculate initial equilibrium
    ageing_on = .false.
    call get_SteadyState()

    close(21)

    open(21, file='output_final.out')

    ! set reform variables
    ageing_on = .true.

    ! set reform values (Table 11.3, different rows)
    tax = 1    ;    tauw = 0d0    ;    taurb = 0d0  ; taurk = 0d0 ! Simulation (1)
    !tax = 3    ;    taur = 0d0                     ! Simulation (2)
    !kappa = 0d0                                    ! Simulation (3)

    ! calculate final equilibrium
    call get_SteadyState()

    close(21)

contains


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


    ! initializes the remaining model parameters and variables
    subroutine initialize

        implicit none
        integer :: ij, ip, is, ib, ik

        write(*,'(/a/)')'INITIAL EQUILIBRIUM'
        write(*,'(a)')'ITER     K/Y     C/Y     I/Y       rb       rk       w        DIFF'

        ! set up population structure
        do ij = 1, JJ
            m(ij) = (1d0+n_p)**(1d0-ij) 
            ! m = relative population share to cohort age ij=1 (see P506) 
        enddo

        ! initialize asset grid
        call grid_Cons_Grow(k, k_l, k_u, k_grow)
        call grid_Cons_Grow(b, b_l, b_u, b_grow)

        ! get initial guess for savings decision
        do ij = 1, JJ
            do ip = 1, NP
                do is = 1, NS
                    do ib = 1,NB
                        do ik = 1,NK
                            kplus(ij, :, ib, ip, is) = max(k(:)/2d0, k(1)/2d0) !!why include k(1)/2d0? isnt it already in k(:)/2d0? 
                            bplus(ij, ik, :, ip, is) = max(b(:)/2d0, b(1)/2d0)
                        enddo
                    enddo
                enddo
            enddo
        enddo

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

        ! initialize fixed effect
        dist_theta = 1d0/dble(NP)
        theta(1)   = -sqrt(sigma_theta)
        theta(2)   = sqrt(sigma_theta)
        theta = exp(theta)


        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, EtaShock, pi)
        EtaShock = exp(EtaShock)

        ! tax and transfers
        tax   = 2
        tauc  = 0.075d0
        tauw  = 0.0d0
        taurk = 0.0d0
        taurb  = 0.0d0
        taup  = 0.1d0
        kappa = 0.5d0
        gy    = 0.19d0
!~         by    = 0.60d0/5d0

        ! initial guesses for macro variables
        KK = 1d0
        BB = 0d0 ! because of the bond market clearing condition ?
        LL = 1d0
        YY = 1d0
        II = (n_p+delta)*KK

        GG = gy*YY
        BB = by*YY

        pen = 0d0
        pen(JR:JJ) = kappa

    end subroutine


    ! subroutine for calculating prices
    subroutine prices()

        implicit none

        rk = Omega*alpha*(KK/LL)**(alpha-1d0)-delta !  ???
        rb = rk-delta 
        w = Omega*(1d0-alpha)*(KK/LL)**alpha
        rkn = rk*(1d0-taurk)
        rbn = rb*(1d0-taurb)
        wn = w*(1d0-tauw-taup)
        p = 1d0 + tauc

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ik, ib, ip, ip_max, is, is_max
        real*8 :: k_in, b_in, wage, available
        real*8 :: a_in(2)
        logical :: check
        
        ! get decision in the last period of life
        do ik = 0, NK
            do ib = 0, NB
                kplus(JJ, ik, ib, :, :) = 0d0
                bplus(JJ, ik, ib, :, :) = 0d0
                c(JJ, ik, ib, :, :) = ((1d0+rkn)*k(ik) + (1d0+rbn)*b(ib) + pen(JJ))/p
                l(JJ, ik, ib, :, :) = 0d0
                V(JJ, ik, ib, :, :) = valuefunc(0d0, 0d0, c(JJ, ik, ib, 1, 1), l(JJ, ik, ib, 1, 1), JJ, 1, 1)
            end do
        enddo

        ! interpolate individual RHS
        call interpolate_b(JJ)
        call interpolate_k(JJ)

        do ij = JJ-1, 1, -1

            ! check about how many is to iterate
            if(ij >= JR)then
                ip_max = 1
                is_max = 1
            else
                ip_max = NP
                is_max = NS
            endif

            do ik = 0, NK
                do ib = 0, NB

                    ! determine decision for zero assets at retirement without pension
                    if(ij >= JR .and. ik == 0 .and. ib == 0 .and. kappa <= 1d-10)then
                        kplus(ij, ik, ib, :, :) = 0d0
                        bplus(ij, ik, ib, :, :) = 0d0
                        c(ij, ik, ib, :, :) = 0d0
                        l(ij, ik, ib, :, :) = 0d0
                        V(ij, ik, ib, :, :) = valuefunc(0d0, 0d0, 0d0, 0d0, ij, 1, 1)
                        cycle
                    endif

                    do ip = 1, ip_max
                        do is = 1, is_max
                            ! get initial guess for the individual choices
                            k_in = kplus(ij+1, ik, ib, ip, is) ! initial guess using value from the next period's kplus in the same states
                            b_in = bplus(ij+1, ik, ib, ip, is)
                            
                            a_in(:) = (/k_in, b_in/)
                            
                            ! set up communication variables
                            ij_com = ij
                            ik_com = ik
                            ib_com = ib
                            ip_com = ip
                            is_com = is

                            call broydn(a_in, foc, check)
!~                             call fzero(b_in, focb, check)
!~                             call fzero(k_in, fock, check)
                            
                            ! write screen output in case of a problem
                            if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ik, ib, ip, is
   
  

                            ! check for borrowing constraint
                            if(k_in < 0d0)then
                                k_in = 0d0
!~                                 wage = wn*eff(ij)*theta(ip)*eta(is)
!~                                 available = (1d0+rn)*a(ia) + pen(ij)
!~                                 if(ij < JR)then
!~                                     lab_com = min( max(nu-(1d0-nu)*available/wage , 0d0) , 1d0-1d-10)
!~                                 else
!~                                     lab_com = 0d0
!~                                 endif
!~                                 cons_com = max( (available + wage*lab_com)/p , 1d-10)
                            endif

                            ! copy decisions
                            kplus(ij, ik, ib, ip, is) = k_in
                                                        
                            bplus(ij, ik, ib, ip, is) = b_in
                            c(ij, ik, ib, ip, is) = cons_com
                            l(ij, ik, ib, ip, is) = lab_com
                            V(ij, ik, ib, ip, is) = valuefunc(k_in, b_in, cons_com, lab_com, ij, ip, is)
                        enddo
                        ! copy decision in retirement age
                        if(ij >= JR)then
                            kplus(ij, ik, ib, :, :) = kplus(ij, ik, ib, 1, 1)
                            bplus(ij, ik, ib, :, :) = bplus(ij, ik, ib, 1, 1)
                            c(ij, ik, ib, :, :) = c(ij, ik, ib, 1, 1)
                            l(ij, ik, ib, :, :) = l(ij, ik, ib, 1, 1)
                            V(ij, ik, ib, :, :) = V(ij, ik, ib, 1, 1)
                        endif
                    enddo
                enddo
            enddo

            ! interpolate individual RHS
            call interpolate_b(ij)
            call interpolate_k(ij)
        enddo

    end subroutine


    ! for calculating the rhs of the first order condition for bond at age ij
    subroutine interpolate_b(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ik, ib, ip, is, is_p
        real*8 :: chelp, lhelp

        do ik = 0, NK
            do ib = 0, NB
                do ip = 1, NP
                    do is = 1, NS

                        ! calculate the RHS of the 1st first order condition and EV
                        ! RHSN = numerator term, RHSD = denominator term
                        RHSN_b(ij, ik, ib, ip, is) = 0d0
                        EV(ij, ik, ib, ip, is) = 0d0
                        do is_p = 1, NS
                            chelp = max(c(ij, ik, ib, ip, is_p),1d-10)
                            lhelp = max(l(ij, ik, ib, ip, is_p),1d-10)
                            RHSN_b(ij, ik, ib, ip, is) = RHSN_b(ij, ik, ib, ip, is) + &
                                pi(is, is_p)*exp(-gamma*V(ij, ik, ib, ip, is_p))*chelp**(-sigma)*(1+rb)
                            RHSD_b(ij, ik, ib, ip, is) = RHSD_b(ij, ik, ib, ip, is) + &
                                pi(is, is_p)*exp(-gamma*V(ij, ik, ib, ip, is_p))
                            EV(ij, ik, ib, ip, is)  = EV(ij, ik, ib, ip, is) + pi(is, is_p)*exp(-gamma*V(ij, ik, ib, ip, is_p))
                        enddo
                        
                        RHS_b(ij, ik, ib, ip, is) = beta*RHSN_b(ij, ik, ib, ip, is)*(1/RHSD_b(ij, ik, ib, ip, is))
!~                         EV(ij, ik, ib, ip, is) = EV(ij, ik, ib, ip, is))
                    enddo
                enddo
            enddo
        enddo

    end subroutine

! for calculating the rhs of the first order condition for capital at age ij
    subroutine interpolate_k(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ik, ib, ip, is, is_p
        real*8 :: chelp, lhelp

        do ik = 0, NK
            do ib = 0, NB
                do ip = 1, NP
                    do is = 1, NS

                        ! calculate the RHS of the 2nd first order condition
                        ! RHSN = numerator term, RHSD = denominator term
                        RHSN_k(ij, ik, ib, ip, is) = 0d0
!~                         EV(ij, ik, ib, ip, is) = 0d0
                        do is_p = 1, NS
                            chelp = max(c(ij, ik, ib, ip, is_p),1d-10)
                            lhelp = max(l(ij, ik, ib, ip, is_p),1d-10)
                            RHSN_k(ij, ik, ib, ip, is) = RHSN_k(ij, ik, ib, ip, is) + &
                                pi(is, is_p)*exp(-gamma*V(ij, ik, ib, ip, is_p))*chelp**(-sigma)*(1+(YY-wn*LL)/KK-delta)
                            RHSD_k(ij, ik, ib, ip, is) = RHSD_k(ij, ik, ib, ip, is) + &
                                pi(is, is_p)*exp(-gamma*V(ij, ik, ib, ip, is_p))
                        enddo
                        
                        RHS_k(ij, ik, ib, ip, is) = beta*RHSN_k(ij, ik, ib, ip, is)*(1/RHSD_k(ij, ik, ib, ip, is))
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ik, ib, ip, is, is_p, ikl, ikr, ibl, ibr
        real*8 :: varphib, varphik

        ! set distribution to zero
        phi(:, :, :, :, :) = 0d0

        ! get initial distribution in age 1
        do ip = 1, NP
            phi(1, 0, 0, ip, is_initial) = dist_theta(ip)
        enddo

print*, 'get distribution', 1d0-(1d0-0.0823d0)

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ik = 0, NK
                do ib = 0, NB
                    do ip = 1 , NP
                        do is = 1, NS
                        
                            ! interpolate yesterday's savings decision
                            call linint_Grow(kplus(ij, ik, ib, ip, is), k_l, k_u, k_grow, NK, ikl, ikr, varphik)
                            call linint_Grow(bplus(ij, ik, ib, ip, is), b_l, b_u, b_grow, NB, ibl, ibr, varphib)

                            ! restrict values to grid just in case
                            ikl = min(ikl, NK)
                            ikr = min(ikr, NK)
                            ibl = min(ibl, NB)
                            ibr = min(ibr, NB)
                            varphik = min(varphik, 1d0)
                            varphib = min(varphib, 1d0)

!~ print*, kplus(ij-1, ik, ib, ip, is), varphik, bplus(ij-1, ik, ib, ip, is), varphib

                            ! redistribute households
                            do is_p = 1, NS
                                phi(ij, ikl, ibl, ip, is_p) = phi(ij, ikl, ibl, ip, is_p) + &
                                                            pi(is, is_p)*varphik*varphib*phi(ij-1, ik, ib, ip, is)
                                phi(ij, ikl, ibr, ip, is_p) = phi(ij, ikl, ibl, ip, is_p) + &
                                                            pi(is, is_p)*varphik*(1d0-varphib)*phi(ij-1, ik, ib, ip, is)
                                phi(ij, ikr, ibl, ip, is_p) = phi(ij, ikr, ib, ip, is_p) + &
                                                            pi(is, is_p)*(1d0-varphik)*varphib*phi(ij-1, ik, ib, ip, is)
                                phi(ij, ikl, ibr, ip, is_p) = phi(ij, ikl, ibl, ip, is_p) + &
                                                            pi(is, is_p)*(1-varphik)*(1-varphib)*phi(ij-1, ik, ib, ip, is)
!~                                 print*, 'phi', phi(ij, ikl, ibl, ip, is_p), phi(ij, ikl, ibr, ip, is_p), &
!~                                 phi(ij, ikr, ibl, ip, is_p), phi(ij, ikl, ibr, ip, is_p)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
print*, 'End of get distribution'
    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ik, ib, ip, is
        real*8 :: workpop, LL_old, KK_old, BB_old

        LL_old = LL ! store current labour supply to LL_old
        KK_old = KK
        BB_old = BB

        ! calculate cohort aggregates
        c_coh(:) = 0d0
        l_coh(:) = 0d0
        y_coh(:) = 0d0
        k_coh(:) = 0d0
        b_coh(:) = 0d0
        v_coh(:) = 0d0

print*, 'Aggregation1'

        do ij = 1, JJ
            do ik = 0, NK
                do ib = 0, NB
                    do ip = 1, NP
                        do is = 1, NS
                            c_coh(ij) = c_coh(ij) + c(ij, ik, ib, ip, is)*phi(ij, ik, ib, ip, is)
                            l_coh(ij) = l_coh(ij) + l(ij, ik, ib, ip, is)*phi(ij, ik, ib, ip, is)
                            y_coh(ij) = y_coh(ij) + eff(ij)*theta(ip)*l(ij, ik, ib, ip, is)*phi(ij, ik, ib, ip, is)
                            k_coh(ij) = k_coh(ij) + k(ik)*phi(ij, ik, ib, ip, is)
!~                             print*, k_coh(ij)
                            b_coh(ij) = b_coh(ij) + b(ib)*phi(ij, ik, ib, ip, is)
                            v_coh(ij) = v_coh(ij) + V(ij, ik, ib, ip, is)*phi(ij, ik, ib, ip, is)
                        enddo
                    enddo
                enddo
            enddo
        enddo

print*, 'Aggregation2'

        ! calculate aggregate quantities
        CC = 0d0
        LL = 0d0
        HH = 0d0
        KK = 0d0
        BB = 0d0
        workpop = 0d0
        do ij = 1, JJ
            CC = CC + c_coh(ij)*m(ij)
            LL = LL + l_coh(ij)*m(ij)
            HH = HH + y_coh(ij)*m(ij)
            KK = KK + k_coh(ij)*m(ij)
            BB = BB + b_coh(ij)*m(ij)
            if(ij < JR) workpop = workpop + m(ij)
        enddo

print*, 'Aggregation3'

        ! damping and other quantities [damping acording to Gauss-Seidel procedure]
        KK = damp*(KK) + (1d0-damp)*KK_old !check this
print*, 'KK =' , KK       
        BB = damp*(BB) + (1d0-damp)*BB_old !check this
print*, 'BB =', BB        
        LL = damp*LL + (1d0-damp)*LL_old
print*, 'LL =', LL        
        II = (n_p+delta)*KK
print*, 'II =', II    
        YY = Omega * KK ** alpha * LL ** (1d0-alpha)
print*, 'YY =', YY 

        ! get average income and average working hours
        INC = w*LL/workpop ! average labour earning
        HH  = HH/workpop

        ! get difference on goods market
        DIFF = YY-CC-II-GG
print*, 'Aggregation5'
    end subroutine


    ! subroutine for calculating government parameters
    subroutine government()

        implicit none
        integer :: ij
        real*8 :: expend

        ! set government quantities and pension payments
!~         if(.not. reform_on)then
!~             GG = gy*YY
!~             BB = by*YY
!~         endif

        ! calculate government expenditure
!~         expend = GG + (1d0+r)*BB - (1d0+n_p)*BB

!~         ! get budget balancing tax rate
!~         if(tax == 1)then
!~             tauc = (expend - (tauw*w*LL + taur*r*AA))/CC
!~             p    = 1d0 + tauc
!~         elseif(tax == 2)then
!~             tauw = (expend - tauc*CC)/(w*LL + r*AA)
!~             taur = tauw
!~         elseif(tax == 3)then
!~             tauw = (expend - (tauc*CC + taur*r*AA))/(w*LL)
!~         else
!~             taur = (expend - (tauc*CC + tauw*w*LL))/(r*AA)
!~         endif

!~         taxrev(1) = tauc*CC
!~         taxrev(2) = tauw*w*LL
!~         taxrev(3) = taur*r*AA
!~         taxrev(4) = sum(taxrev(1:3))

        ! get budget balancing social security contribution
        pen(JR:JJ) = kappa*INC
        PP = 0d0
        do ij = JR, JJ
            PP = PP + pen(ij)*m(ij)
        enddo

        taup = PP/(w*LL)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ik, ib, ip, is, iamax(JJ)
        real*8 :: temp
        real*8 :: exp_c(JJ), exp_l(JJ), exp_y(JJ)
        real*8 :: var_c(JJ), var_l(JJ), var_y(JJ)
        real*8 :: mas_c(JJ), mas_l(JJ), mas_y(JJ)
        

        ! calculate cohort specific variances of logs
        exp_c = 0d0 ; var_c = 0d0 ; mas_c = 0d0
        exp_l = 0d0 ; var_l = 0d0 ; mas_l = 0d0
        exp_y = 0d0 ; var_y = 0d0 ; mas_y = 0d0
        do ij = 1, JJ
            do ik = 0, NK
                do ib = 0, NB
                    do ip = 1, NP
                        do is = 1, NS

                            ! consumption
                            if(c(ij, ik, ib, ip, is) > 0d0)then
                                temp = log(c(ij, ik, ib, ip, is))
                                exp_c(ij) = exp_c(ij) + temp*phi(ij, ik, ib, ip, is)
                                var_c(ij) = var_c(ij) + temp**2*phi(ij, ik, ib, ip, is)
                                mas_c(ij) = mas_c(ij) + phi(ij, ik, ib, ip, is)
                            endif

                            if(l(ij, ik, ib, ip, is) > 0.01d0)then

                                ! hours
                                temp = log(l(ij, ik, ib, ip, is))
                                exp_l(ij) = exp_l(ij) + temp*phi(ij, ik, ib, ip, is)
                                var_l(ij) = var_l(ij) + temp**2*phi(ij, ik, ib, ip, is)
                                mas_l(ij) = mas_l(ij) + phi(ij, ik, ib, ip, is)

                                ! earnings
                                temp = log(w*eff(ij)*theta(ip)*l(ij, ik, ib, ip, is))
                                exp_y(ij) = exp_y(ij) + temp*phi(ij, ik, ib, ip, is)
                                var_y(ij) = var_y(ij) + temp**2*phi(ij, ik, ib, ip, is)
                                mas_y(ij) = mas_y(ij) + phi(ij, ik, ib, ip, is)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        exp_c = exp_c/max(mas_c, 1d-4) ; var_c = var_c/max(mas_c, 1d-4)
        exp_l = exp_l/max(mas_l, 1d-4) ; var_l = var_l/max(mas_l, 1d-4)
        exp_y = exp_y/max(mas_y, 1d-4) ; var_y = var_y/max(mas_y, 1d-4)
        var_c = var_c - exp_c**2
        var_l = var_l - exp_l**2
        var_y = var_y - exp_y**2

        ! save initial equilibrium average income if no reform
        if(.not. ageing_on)INC_init = INC

        ! Output
        write(21,'(a/)')'STEADY STATE EQUILIBRIUM'
        write(21,'(a)')'CAPITAL        K       B       rk      rb  rk p.a.  rb p.a.'
        write(21,'(8x,5f8.2)')KK, BB, rk, rb , ((1d0+rk)**(1d0/5d0)-1d0)*100d0, ((1d0+rb)**(1d0/5d0)-1d0)*100d0 !check
        write(21,'(a,3f8.2/)')'(in %)  ',(/KK, BB/)/YY*500d0

        write(21,'(a)')'LABOR          L      HH     INC       w'
        write(21,'(8x,4f8.2/)')LL, HH*100d0, INC, w

        write(21,'(a)')'GOODS          Y       C       I       G    DIFF'
        write(21,'(8x,4f8.2,f8.3)')YY,CC,II,GG,diff
        write(21,'(a,4f8.2,f8.3/)')'(in %)  ',(/YY, CC, II, GG, diff/)/YY*100d0

!~         write(21,'(a)')'GOV         TAUC    TAUW    TAUR   TOTAL       G       B'
!~         write(21,'(8x,6f8.2)')taxrev(1:4),GG,BB
!~         write(21,'(a,6f8.2)')'(in %)  ',(/taxrev(1:4), GG, BB*5d0/)/YY*100d0
!~         write(21,'(a,3f8.2/)')'(rate)  ',(/tauc, tauw, taur/)*100d0

        write(21,'(a)')'PENS        TAUP     PEN      PP'
        write(21,'(8x,6f8.2)')taup*w*LL, pen(JR), PP
        write(21,'(a,3f8.2/)')'(in %)  ',(/taup, kappa, PP/YY/)*100d0

        ! check for the maximium grid point used
        call check_grid(iamax)

        write(21, '(a,a)')' IJ      CONS     LABOR  EARNINGS    INCOME    INCTAX      PENS    ASSETS', &
            '    BONDs    VAR(C)    VAR(L)    VAR(Y)     VALUE     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,11f10.3,i10)')ij, c_coh(ij)/INC_init, l_coh(ij), (/w*y_coh(ij), wn*y_coh(ij)+rkn*k_coh(ij), &
                    tauw*w*y_coh(ij)+taurk*rk*k_coh(ij), pen(ij)-taup*w*y_coh(ij), 5d0*k_coh(ij), 5d0*b_coh(ij)/)/INC_init, &
                    var_c(ij), var_l(ij), var_y(ij), v_coh(ij), iamax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid(iamax)

        implicit none
        integer :: iamax(JJ), ij, ik, ib, ip, is

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ik = 0, NK
                do ib = 0, NB
                    do ip = 1, NP
                        do is = 1, NS
                            if(phi(ij, ik, ib, ip, is) > 1d-8)iamax(ij) = ik
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine



end program
