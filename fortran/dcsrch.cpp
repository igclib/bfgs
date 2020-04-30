#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

void
dcstep(...)
{
  throw std::runtime_error(
    "Missing function implementation: dcstep");
}

using fem::common;

void
dcsrch(
  double& stp,
  double const& f,
  double const& g,
  double const& ftol,
  double const& gtol,
  double const& xtol,
  str_ref task,
  double const& stpmin,
  double const& stpmax,
  arr_ref<int> isave,
  arr_ref<double> dsave)
{
  isave(dimension(2));
  dsave(dimension(13));
  const double zero = 0.0e0;
  bool brackt = fem::bool0;
  int stage = fem::int0;
  double finit = fem::double0;
  double ginit = fem::double0;
  double gtest = fem::double0;
  double width = fem::double0;
  const double p5 = 0.5e0;
  double width1 = fem::double0;
  double stx = fem::double0;
  double fx = fem::double0;
  double gx = fem::double0;
  double sty = fem::double0;
  double fy = fem::double0;
  double gy = fem::double0;
  double stmin = fem::double0;
  const double xtrapu = 4.0e0;
  double stmax = fem::double0;
  double ftest = fem::double0;
  double fm = fem::double0;
  double fxm = fem::double0;
  double fym = fem::double0;
  double gm = fem::double0;
  double gxm = fem::double0;
  double gym = fem::double0;
  const double p66 = 0.66e0;
  const double xtrapl = 1.1e0;
  //C     **********
  //C
  //C     Subroutine dcsrch
  //C
  //C     This subroutine finds a step that satisfies a sufficient
  //C     decrease condition and a curvature condition.
  //C
  //C     Each call of the subroutine updates an interval with
  //C     endpoints stx and sty. The interval is initially chosen
  //C     so that it contains a minimizer of the modified function
  //C
  //C           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
  //C
  //C     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
  //C     interval is chosen so that it contains a minimizer of f.
  //C
  //C     The algorithm is designed to find a step that satisfies
  //C     the sufficient decrease condition
  //C
  //C           f(stp) <= f(0) + ftol*stp*f'(0),
  //C
  //C     and the curvature condition
  //C
  //C           abs(f'(stp)) <= gtol*abs(f'(0)).
  //C
  //C     If ftol is less than gtol and if, for example, the function
  //C     is bounded below, then there is always a step which satisfies
  //C     both conditions.
  //C
  //C     If no step can be found that satisfies both conditions, then
  //C     the algorithm stops with a warning. In this case stp only
  //C     satisfies the sufficient decrease condition.
  //C
  //C     A typical invocation of dcsrch has the following outline:
  //C
  //C     Evaluate the function at stp = 0.0d0; store in f.
  //C     Evaluate the gradient at stp = 0.0d0; store in g.
  //C     Choose a starting step stp.
  //C
  //C     task = 'START'
  //C  10 continue
  //C        call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
  //C    +               isave,dsave)
  //C        if (task .eq. 'FG') then
  //C           Evaluate the function and the gradient at stp
  //C           go to 10
  //C           end if
  //C
  //C     NOTE: The user must not alter work arrays between calls.
  //C
  //C     The subroutine statement is
  //C
  //C       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
  //C                         task,isave,dsave)
  //C     where
  //C
  //C       stp is a double precision variable.
  //C         On entry stp is the current estimate of a satisfactory
  //C            step. On initial entry, a positive initial estimate
  //C            must be provided.
  //C         On exit stp is the current estimate of a satisfactory step
  //C            if task = 'FG'. If task = 'CONV' then stp satisfies
  //C            the sufficient decrease and curvature condition.
  //C
  //C       f is a double precision variable.
  //C         On initial entry f is the value of the function at 0.
  //C            On subsequent entries f is the value of the
  //C            function at stp.
  //C         On exit f is the value of the function at stp.
  //C
  //C       g is a double precision variable.
  //C         On initial entry g is the derivative of the function at 0.
  //C            On subsequent entries g is the derivative of the
  //C            function at stp.
  //C         On exit g is the derivative of the function at stp.
  //C
  //C       ftol is a double precision variable.
  //C         On entry ftol specifies a nonnegative tolerance for the
  //C            sufficient decrease condition.
  //C         On exit ftol is unchanged.
  //C
  //C       gtol is a double precision variable.
  //C         On entry gtol specifies a nonnegative tolerance for the
  //C            curvature condition.
  //C         On exit gtol is unchanged.
  //C
  //C       xtol is a double precision variable.
  //C         On entry xtol specifies a nonnegative relative tolerance
  //C            for an acceptable step. The subroutine exits with a
  //C            warning if the relative difference between sty and stx
  //C            is less than xtol.
  //C         On exit xtol is unchanged.
  //C
  //C       task is a character variable of length at least 60.
  //C         On initial entry task must be set to 'START'.
  //C         On exit task indicates the required action:
  //C
  //C            If task(1:2) = 'FG' then evaluate the function and
  //C            derivative at stp and call dcsrch again.
  //C
  //C            If task(1:4) = 'CONV' then the search is successful.
  //C
  //C            If task(1:4) = 'WARN' then the subroutine is not able
  //C            to satisfy the convergence conditions. The exit value of
  //C            stp contains the best point found during the search.
  //C
  //C            If task(1:5) = 'ERROR' then there is an error in the
  //C            input arguments.
  //C
  //C         On exit with convergence, a warning or an error, the
  //C            variable task contains additional information.
  //C
  //C       stpmin is a double precision variable.
  //C         On entry stpmin is a nonnegative lower bound for the step.
  //C         On exit stpmin is unchanged.
  //C
  //C       stpmax is a double precision variable.
  //C         On entry stpmax is a nonnegative upper bound for the step.
  //C         On exit stpmax is unchanged.
  //C
  //C       isave is an integer work array of dimension 2.
  //C
  //C       dsave is a double precision work array of dimension 13.
  //C
  //C     Subprograms called
  //C
  //C       MINPACK-2 ... dcstep
  //C
  //C     MINPACK-1 Project. June 1983.
  //C     Argonne National Laboratory.
  //C     Jorge J. More' and David J. Thuente.
  //C
  //C     MINPACK-2 Project. November 1993.
  //C     Argonne National Laboratory and University of Minnesota.
  //C     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
  //C
  //C     **********
  //C
  //C     Initialization block.
  //C
  if (task(1, 5) == "START") {
    //C
    //C        Check the input arguments for errors.
    //C
    if (stp < stpmin) {
      task = "ERROR: STP .LT. STPMIN";
    }
    if (stp > stpmax) {
      task = "ERROR: STP .GT. STPMAX";
    }
    if (g >= zero) {
      task = "ERROR: INITIAL G .GE. ZERO";
    }
    if (ftol < zero) {
      task = "ERROR: FTOL .LT. ZERO";
    }
    if (gtol < zero) {
      task = "ERROR: GTOL .LT. ZERO";
    }
    if (xtol < zero) {
      task = "ERROR: XTOL .LT. ZERO";
    }
    if (stpmin < zero) {
      task = "ERROR: STPMIN .LT. ZERO";
    }
    if (stpmax < stpmin) {
      task = "ERROR: STPMAX .LT. STPMIN";
    }
    //C
    //C        Exit if there are errors on input.
    //C
    if (task(1, 5) == "ERROR") {
      return;
    }
    //C
    //C        Initialize local variables.
    //C
    brackt = false;
    stage = 1;
    finit = f;
    ginit = g;
    gtest = ftol * ginit;
    width = stpmax - stpmin;
    width1 = width / p5;
    //C
    //C        The variables stx, fx, gx contain the values of the step,
    //C        function, and derivative at the best step.
    //C        The variables sty, fy, gy contain the value of the step,
    //C        function, and derivative at sty.
    //C        The variables stp, f, g contain the values of the step,
    //C        function, and derivative at stp.
    //C
    stx = zero;
    fx = finit;
    gx = ginit;
    sty = zero;
    fy = finit;
    gy = ginit;
    stmin = zero;
    stmax = stp + xtrapu * stp;
    task = "FG";
    //C
    goto statement_10;
    //C
  }
  else {
    //C
    //C        Restore local variables.
    //C
    if (isave(1) == 1) {
      brackt = true;
    }
    else {
      brackt = false;
    }
    stage = isave(2);
    ginit = dsave(1);
    gtest = dsave(2);
    gx = dsave(3);
    gy = dsave(4);
    finit = dsave(5);
    fx = dsave(6);
    fy = dsave(7);
    stx = dsave(8);
    sty = dsave(9);
    stmin = dsave(10);
    stmax = dsave(11);
    width = dsave(12);
    width1 = dsave(13);
    //C
  }
  //C
  //C     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
  //C     algorithm enters the second stage.
  //C
  ftest = finit + stp * gtest;
  if (stage == 1 && f <= ftest && g >= zero) {
    stage = 2;
  }
  //C
  //C     Test for warnings.
  //C
  if (brackt && (stp <= stmin || stp >= stmax)) {
    task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS";
  }
  if (brackt && stmax - stmin <= xtol * stmax) {
    task = "WARNING: XTOL TEST SATISFIED";
  }
  if (stp == stpmax && f <= ftest && g <= gtest) {
    task = "WARNING: STP = STPMAX";
  }
  if (stp == stpmin && (f > ftest || g >= gtest)) {
    task = "WARNING: STP = STPMIN";
  }
  //C
  //C     Test for convergence.
  //C
  if (f <= ftest && fem::abs(g) <= gtol * (-ginit)) {
    task = "CONVERGENCE";
  }
  //C
  //C     Test for termination.
  //C
  if (task(1, 4) == "WARN" || task(1, 4) == "CONV") {
    goto statement_10;
  }
  //C
  //C     A modified function is used to predict the step during the
  //C     first stage if a lower function value has been obtained but
  //C     the decrease is not sufficient.
  //C
  if (stage == 1 && f <= fx && f > ftest) {
    //C
    //C        Define the modified function and derivative values.
    //C
    fm = f - stp * gtest;
    fxm = fx - stx * gtest;
    fym = fy - sty * gtest;
    gm = g - gtest;
    gxm = gx - gtest;
    gym = gy - gtest;
    //C
    //C        Call dcstep to update stx, sty, and to compute the new step.
    //C
    dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
    //C
    //C        Reset the function and derivative values for f.
    //C
    fx = fxm + stx * gtest;
    fy = fym + sty * gtest;
    gx = gxm + gtest;
    gy = gym + gtest;
    //C
  }
  else {
    //C
    //C       Call dcstep to update stx, sty, and to compute the new step.
    //C
    dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
    //C
  }
  //C
  //C     Decide if a bisection step is needed.
  //C
  if (brackt) {
    if (fem::abs(sty - stx) >= p66 * width1) {
      stp = stx + p5 * (sty - stx);
    }
    width1 = width;
    width = fem::abs(sty - stx);
  }
  //C
  //C     Set the minimum and maximum steps allowed for stp.
  //C
  if (brackt) {
    stmin = fem::min(stx, sty);
    stmax = fem::max(stx, sty);
  }
  else {
    stmin = stp + xtrapl * (stp - stx);
    stmax = stp + xtrapu * (stp - stx);
  }
  //C
  //C     Force the step to be within the bounds stpmax and stpmin.
  //C
  stp = fem::max(stp, stpmin);
  stp = fem::min(stp, stpmax);
  //C
  //C     If further progress is not possible, let stp be the best
  //C     point obtained during the search.
  //C
  if (brackt && (stp <= stmin || stp >= stmax) || (brackt && stmax -
      stmin <= xtol * stmax)) {
    stp = stx;
  }
  //C
  //C     Obtain another function and derivative.
  //C
  task = "FG";
  //C
  statement_10:
  //C
  //C     Save local variables.
  //C
  if (brackt) {
    isave(1) = 1;
  }
  else {
    isave(1) = 0;
  }
  isave(2) = stage;
  dsave(1) = ginit;
  dsave(2) = gtest;
  dsave(3) = gx;
  dsave(4) = gy;
  dsave(5) = finit;
  dsave(6) = fx;
  dsave(7) = fy;
  dsave(8) = stx;
  dsave(9) = sty;
  dsave(10) = stmin;
  dsave(11) = stmax;
  dsave(12) = width;
  dsave(13) = width1;
  //C
}

} // namespace placeholder_please_replace
