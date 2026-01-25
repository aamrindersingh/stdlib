! Test least squares solver
module test_linalg_least_squares
    use testdrive, only: error_type, check, new_unittest, unittest_type
    use stdlib_linalg_constants
    use stdlib_linalg, only: lstsq, solve_lstsq, weighted_lstsq
    use stdlib_linalg_state, only: linalg_state_type

    implicit none (type,external)
    private
    
    public :: test_least_squares 

    contains

    !> Solve sample least squares problems
    subroutine test_least_squares(tests)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: tests(:)
        
        allocate(tests(0))
        
        call add_test(tests,new_unittest("issue_823",test_issue_823))

        call add_test(tests,new_unittest("least_squares_s",test_lstsq_one_s))
        call add_test(tests,new_unittest("least_squares_randm_s",test_lstsq_random_s))
        call add_test(tests,new_unittest("weighted_lstsq_s",test_weighted_lstsq_s))
        call add_test(tests,new_unittest("weighted_lstsq_effect_s",test_weighted_lstsq_effect_s))
        call add_test(tests,new_unittest("weighted_lstsq_negative_s",test_weighted_lstsq_negative_s))
        call add_test(tests,new_unittest("least_squares_d",test_lstsq_one_d))
        call add_test(tests,new_unittest("least_squares_randm_d",test_lstsq_random_d))
        call add_test(tests,new_unittest("weighted_lstsq_d",test_weighted_lstsq_d))
        call add_test(tests,new_unittest("weighted_lstsq_effect_d",test_weighted_lstsq_effect_d))
        call add_test(tests,new_unittest("weighted_lstsq_negative_d",test_weighted_lstsq_negative_d))

    end subroutine test_least_squares
    
    !> Simple polynomial fit
    subroutine test_lstsq_one_s(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp) :: rank

        !> Example scattered data
        real(sp), parameter :: x(*)  = real([1.0, 2.5, 3.5, 4.0, 5.0, 7.0, 8.5], sp)
        real(sp), parameter :: y(*)  = real([0.3, 1.1, 1.5, 2.0, 3.2, 6.6, 8.6], sp)
        real(sp), parameter :: ab(*) = real([0.20925829,  0.12013861], sp)

        real(sp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        call check(error, all(abs(p-ab)<1.0e-4_sp), 'data converged')
        if (allocated(error)) return
        
        call check(error, rank==2, 'matrix rank == 2')
        if (allocated(error)) return

    end subroutine test_lstsq_one_s
    
    !> Fit from random array
    subroutine test_lstsq_random_s(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp), parameter :: n = 12, m = 3
        real :: Arnd(n,m),xrnd(m)
        real(sp), allocatable :: x(:)
        real(sp) :: xsol(m),y(n),A(n,m)

        ! Random coefficient matrix and solution
        call random_number(Arnd)
        call random_number(xrnd)
        
        ! Compute rhs
        A    = real(Arnd,sp)
        xsol = real(xrnd,sp)
        y    = matmul(A,xsol)

        ! Find polynomial
        x = lstsq(A,y,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        ! Check size
        call check(error,size(x)==m)
        if (allocated(error)) return
        
        call check(error, all(abs(x-xsol)<1.0e-4_sp), 'data converged')
        if (allocated(error)) return
        
    end subroutine test_lstsq_random_s    
    
    !> Simple polynomial fit
    subroutine test_lstsq_one_d(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp) :: rank

        !> Example scattered data
        real(dp), parameter :: x(*)  = real([1.0, 2.5, 3.5, 4.0, 5.0, 7.0, 8.5], dp)
        real(dp), parameter :: y(*)  = real([0.3, 1.1, 1.5, 2.0, 3.2, 6.6, 8.6], dp)
        real(dp), parameter :: ab(*) = real([0.20925829,  0.12013861], dp)

        real(dp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        call check(error, all(abs(p-ab)<1.0e-4_dp), 'data converged')
        if (allocated(error)) return
        
        call check(error, rank==2, 'matrix rank == 2')
        if (allocated(error)) return

    end subroutine test_lstsq_one_d
    
    !> Fit from random array
    subroutine test_lstsq_random_d(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp), parameter :: n = 12, m = 3
        real :: Arnd(n,m),xrnd(m)
        real(dp), allocatable :: x(:)
        real(dp) :: xsol(m),y(n),A(n,m)

        ! Random coefficient matrix and solution
        call random_number(Arnd)
        call random_number(xrnd)
        
        ! Compute rhs
        A    = real(Arnd,dp)
        xsol = real(xrnd,dp)
        y    = matmul(A,xsol)

        ! Find polynomial
        x = lstsq(A,y,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        ! Check size
        call check(error,size(x)==m)
        if (allocated(error)) return
        
        call check(error, all(abs(x-xsol)<1.0e-4_dp), 'data converged')
        if (allocated(error)) return
        
    end subroutine test_lstsq_random_d    
    

    !-------------------------------------------------------------
    !-----     Weighted Least-Squares Tests                  -----
    !-------------------------------------------------------------

    !> Test basic weighted least-squares
    subroutine test_weighted_lstsq_s(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp), parameter :: m = 4, n = 2
        real(sp), parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n), b(m)
        real(sp) :: w(m)
        real(sp), allocatable :: x(:)

        ! Simple test case
        A(:,1) = 1.0_sp
        A(:,2) = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
        b = [2.0_sp, 4.0_sp, 5.0_sp, 4.0_sp]
        w = [1.0_sp, 1.0_sp, 1.0_sp, 1.0_sp]  ! Uniform weights = OLS

        x = weighted_lstsq(w, A, b, err=state)

        call check(error, state%ok(), 'weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error, size(x)==n, 'weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! Verify residual is small
        call check(error, norm2(matmul(A,x) - b) < 2.0_sp, 'weighted_lstsq: residual too large')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_s

    !> Test that non-uniform weights change the solution
    subroutine test_weighted_lstsq_effect_s(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp), parameter :: m = 4, n = 2
        real(sp), parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n), b(m)
        real(sp) :: w_uniform(m), w_nonuniform(m)
        real(sp), allocatable :: x_uniform(:), x_weighted(:)

        ! Setup problem
        A(:,1) = 1.0_sp
        A(:,2) = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
        b = [1.0_sp, 3.0_sp, 2.0_sp, 5.0_sp]
        
        w_uniform = 1.0_sp
        w_nonuniform = [10.0_sp, 1.0_sp, 1.0_sp, 10.0_sp]  ! Weight first and last more

        x_uniform = weighted_lstsq(w_uniform, A, b, err=state)
        call check(error, state%ok(), 'uniform weighted_lstsq failed')
        if (allocated(error)) return

        x_weighted = weighted_lstsq(w_nonuniform, A, b, err=state)
        call check(error, state%ok(), 'non-uniform weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error, any(abs(x_uniform - x_weighted) > tol), &
                   'weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_effect_s

    !> Test error on negative weights
    subroutine test_weighted_lstsq_negative_s(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        real(sp) :: A(3,2), b(3)
        real(sp) :: w(3)
        real(sp), allocatable :: x(:)

        A = 1.0_sp
        b = 1.0_sp
        w = [-1.0_sp, 1.0_sp, 1.0_sp]  ! Invalid: negative weight!

        x = weighted_lstsq(w, A, b, err=state)

        call check(error, state%error(), 'weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_negative_s

    !> Test basic weighted least-squares
    subroutine test_weighted_lstsq_d(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp), parameter :: m = 4, n = 2
        real(dp), parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n), b(m)
        real(dp) :: w(m)
        real(dp), allocatable :: x(:)

        ! Simple test case
        A(:,1) = 1.0_dp
        A(:,2) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        b = [2.0_dp, 4.0_dp, 5.0_dp, 4.0_dp]
        w = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]  ! Uniform weights = OLS

        x = weighted_lstsq(w, A, b, err=state)

        call check(error, state%ok(), 'weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error, size(x)==n, 'weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! Verify residual is small
        call check(error, norm2(matmul(A,x) - b) < 2.0_dp, 'weighted_lstsq: residual too large')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_d

    !> Test that non-uniform weights change the solution
    subroutine test_weighted_lstsq_effect_d(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        integer(ilp), parameter :: m = 4, n = 2
        real(dp), parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n), b(m)
        real(dp) :: w_uniform(m), w_nonuniform(m)
        real(dp), allocatable :: x_uniform(:), x_weighted(:)

        ! Setup problem
        A(:,1) = 1.0_dp
        A(:,2) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        b = [1.0_dp, 3.0_dp, 2.0_dp, 5.0_dp]
        
        w_uniform = 1.0_dp
        w_nonuniform = [10.0_dp, 1.0_dp, 1.0_dp, 10.0_dp]  ! Weight first and last more

        x_uniform = weighted_lstsq(w_uniform, A, b, err=state)
        call check(error, state%ok(), 'uniform weighted_lstsq failed')
        if (allocated(error)) return

        x_weighted = weighted_lstsq(w_nonuniform, A, b, err=state)
        call check(error, state%ok(), 'non-uniform weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error, any(abs(x_uniform - x_weighted) > tol), &
                   'weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_effect_d

    !> Test error on negative weights
    subroutine test_weighted_lstsq_negative_d(error)
        type(error_type), allocatable, intent(out) :: error

        type(linalg_state_type) :: state
        real(dp) :: A(3,2), b(3)
        real(dp) :: w(3)
        real(dp), allocatable :: x(:)

        A = 1.0_dp
        b = 1.0_dp
        w = [-1.0_dp, 1.0_dp, 1.0_dp]  ! Invalid: negative weight!

        x = weighted_lstsq(w, A, b, err=state)

        call check(error, state%error(), 'weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_negative_d

    
    ! Test issue #823
    subroutine test_issue_823(error)
        type(error_type), allocatable, intent(out) :: error
        
        ! Dimension of the problem.
        integer(ilp), parameter :: n = 42
        ! Data for the least-squares problem.
        complex(dp) :: A(n+1, n), b(n+1), x_true(n), x_lstsq(n)
        ! Internal variables.
        real(dp), allocatable :: tmp(:, :, :), tmp_vec(:, :)
        ! Error handler
        type(linalg_state_type) :: state

        ! Zero-out data.
        A = 0.0_dp
        b = 0.0_dp
        x_lstsq = 0.0_dp
        allocate(tmp(n+1, n, 2), tmp_vec(n, 2), source=0.0_dp)

        ! Generate a random complex least-squares problem of size (n+1, n).
        call random_number(tmp)
        call random_number(tmp_vec)
        
        A      = cmplx(tmp(:, :, 1), tmp(:, :, 2), kind=dp)
        x_true = cmplx(tmp_vec(:, 1), tmp_vec(:, 2), kind=dp)
        b      = matmul(A, x_true)

        ! Solve the lstsq problem.
        call solve_lstsq(A, b, x_lstsq, err=state)
          
        ! Check that no segfault occurred
        call check(error,state%ok(),'issue 823 returned '//state%print())
        if (allocated(error)) return

        ! Check that least squares are verified
        call check(error,all(abs(x_true-x_lstsq)<sqrt(epsilon(0.0_dp))),'issue 823 results')
        if (allocated(error)) return

    end subroutine test_issue_823

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type), allocatable, intent(inout) :: tests(:)    
        type(unittest_type), intent(in) :: new_test
        
        integer :: n
        type(unittest_type), allocatable :: new_tests(:)
        
        if (allocated(tests)) then 
            n = size(tests)
        else
            n = 0
        end if
        
        allocate(new_tests(n+1))
        if (n>0) new_tests(1:n) = tests(1:n)
                 new_tests(1+n) = new_test
        call move_alloc(from=new_tests,to=tests)        
        
    end subroutine add_test

end module test_linalg_least_squares

program test_lstsq
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_linalg_least_squares, only : test_least_squares
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_least_squares", test_least_squares) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program test_lstsq
