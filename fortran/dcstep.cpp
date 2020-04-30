#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

using fem::common;

void
dcstep(
  double& stx,
  double& fx,
  double& dx,
  double& sty,
  double& fy,
  double& dy,
  double& stp,
  double const& fp,
  double const& dp,
  bool& brackt,
  double const& stpmin,
  double const& stpmax)
{
  //C     **********
  //C
  //C     Subroutine dcstep
  //C
  //C     This subroutine computes a safeguarded step for a search
  //C     procedure and updates an interval that contains a step that
  //C     satisfies a sufficient decrease and a curvature condition.
  //C
  //C     The parameter stx contains the step with the least function
  //C     value. If brackt is set to .true. then a minimizer has
  //C     been bracketed in an interval with endpoints stx and sty.
  //C     The parameter stp contains the current step.
  //C     The subroutine assumes that if brackt is set to .true. then
  //C
  //C           min(stx,sty) < stp < max(stx,sty),
  //C
  //C     and that the derivative at stx is negative in the direction
  //C     of the step.
  //C
  //C     The subroutine statement is
  //C
  //C       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
  //C                         stpmin,stpmax)
  //C
  //C     where
  //C
  //C       stx is a double precision variable.
  //C         On entry stx is the best step obtained so far and is an
  //C            endpoint of the interval that contains the minimizer.
  //C         On exit stx is the updated best step.
  //C
  //C       fx is a double precision variable.
  //C         On entry fx is the function at stx.
  //C         On exit fx is the function at stx.
  //C
  //C       dx is a double precision variable.
  //C         On entry dx is the derivative of the function at
  //C            stx. The derivative must be negative in the direction of
  //C            the step, that is, dx and stp - stx must have opposite
  //C            signs.
  //C         On exit dx is the derivative of the function at stx.
  //C
  //C       sty is a double precision variable.
  //C         On entry sty is the second endpoint of the interval that
  //C            contains the minimizer.
  //C         On exit sty is the updated endpoint of the interval that
  //C            contains the minimizer.
  //C
  //C       fy is a double precision variable.
  //C         On entry fy is the function at sty.
  //C         On exit fy is the function at sty.
  //C
  //C       dy is a double precision variable.
  //C         On entry dy is the derivative of the function at sty.
  //C         On exit dy is the derivative of the function at the exit sty.
  //C
  //C       stp is a double precision variable.
  //C         On entry stp is the current step. If brackt is set to .true.
  //C            then on input stp must be between stx and sty.
  //C         On exit stp is a new trial step.
  //C
  //C       fp is a double precision variable.
  //C         On entry fp is the function at stp
  //C         On exit fp is unchanged.
  //C
  //C       dp is a double precision variable.
  //C         On entry dp is the the derivative of the function at stp.
  //C         On exit dp is unchanged.
  //C
  //C       brackt is an logical variable.
  //C         On entry brackt specifies if a minimizer has been bracketed.
  //C            Initially brackt must be set to .false.
  //C         On exit brackt specifies if a minimizer has been bracketed.
  //C            When a minimizer is bracketed brackt is set to .true.
  //C
  //C       stpmin is a double precision variable.
  //C         On entry stpmin is a lower bound for the step.
  //C         On exit stpmin is unchanged.
  //C
  //C       stpmax is a double precision variable.
  //C         On entry stpmax is an upper bound for the step.
  //C         On exit stpmax is unchanged.
  //C
  //C     MINPACK-1 Project. June 1983
  //C     Argonne National Laboratory.
  //C     Jorge J. More' and David J. Thuente.
  //C
  //C     MINPACK-2 Project. November 1993.
  //C     Argonne National Laboratory and University of Minnesota.
  //C     Brett M. Averick and Jorge J. More'.
  //C
  //C     **********
  //C
  double sgnd = dp * (dx / fem::abs(dx));
  //C
  //C     First case: A higher function value. The minimum is bracketed.
  //C     If the cubic step is closer to stx than the quadratic step, the
  //C     cubic step is taken, otherwise the average of the cubic and
  //C     quadratic steps is taken.
  //C
  const double three = 3.0e0;
  double theta = fem::double0;
  double s = fem::double0;
  double gamma = fem::double0;
  double p = fem::double0;
  double q = fem::double0;
  double r = fem::double0;
  double stpc = fem::double0;
  const double two = 2.0e0;
  double stpq = fem::double0;
  double stpf = fem::double0;
  const double zero = 0.0e0;
  const double p66 = 0.66e0;
  if (fp > fx) {
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = fem::max(fem::abs(theta), fem::abs(dx), fem::abs(dp));
    gamma = s * fem::sqrt(fem::pow2((theta / s)) - (dx / s) * (dp / s));
    if (stp < stx) {
      gamma = -gamma;
    }
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p / q;
    stpc = stx + r * (stp - stx);
    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx);
    if (fem::abs(stpc - stx) < fem::abs(stpq - stx)) {
      stpf = stpc;
    }
    else {
      stpf = stpc + (stpq - stpc) / two;
    }
    brackt = true;
    //C
    //C     Second case: A lower function value and derivatives of opposite
    //C     sign. The minimum is bracketed. If the cubic step is farther from
    //C     stp than the secant step, the cubic step is taken, otherwise the
    //C     secant step is taken.
    //C
  }
  else if (sgnd < zero) {
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = fem::max(fem::abs(theta), fem::abs(dx), fem::abs(dp));
    gamma = s * fem::sqrt(fem::pow2((theta / s)) - (dx / s) * (dp / s));
    if (stp > stx) {
      gamma = -gamma;
    }
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p / q;
    stpc = stp + r * (stx - stp);
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (fem::abs(stpc - stp) > fem::abs(stpq - stp)) {
      stpf = stpc;
    }
    else {
      stpf = stpq;
    }
    brackt = true;
    //C
    //C     Third case: A lower function value, derivatives of the same sign,
    //C     and the magnitude of the derivative decreases.
    //C
  }
  else if (fem::abs(dp) < fem::abs(dx)) {
    //C
    //C        The cubic step is computed only if the cubic tends to infinity
    //C        in the direction of the step or if the minimum of the cubic
    //C        is beyond stp. Otherwise the cubic step is defined to be the
    //C        secant step.
    //C
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = fem::max(fem::abs(theta), fem::abs(dx), fem::abs(dp));
    //C
    //C        The case gamma = 0 only arises if the cubic does not tend
    //C        to infinity in the direction of the step.
    //C
    gamma = s * fem::sqrt(fem::max(zero, fem::pow2((theta / s)) - (
      dx / s) * (dp / s)));
    if (stp > stx) {
      gamma = -gamma;
    }
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p / q;
    if (r < zero && gamma != zero) {
      stpc = stp + r * (stx - stp);
    }
    else if (stp > stx) {
      stpc = stpmax;
    }
    else {
      stpc = stpmin;
    }
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    //C
    if (brackt) {
      //C
      //C           A minimizer has been bracketed. If the cubic step is
      //C           closer to stp than the secant step, the cubic step is
      //C           taken, otherwise the secant step is taken.
      //C
      if (fem::abs(stpc - stp) < fem::abs(stpq - stp)) {
        stpf = stpc;
      }
      else {
        stpf = stpq;
      }
      if (stp > stx) {
        stpf = fem::min(stp + p66 * (sty - stp), stpf);
      }
      else {
        stpf = fem::max(stp + p66 * (sty - stp), stpf);
      }
    }
    else {
      //C
      //C           A minimizer has not been bracketed. If the cubic step is
      //C           farther from stp than the secant step, the cubic step is
      //C           taken, otherwise the secant step is taken.
      //C
      if (fem::abs(stpc - stp) > fem::abs(stpq - stp)) {
        stpf = stpc;
      }
      else {
        stpf = stpq;
      }
      stpf = fem::min(stpmax, stpf);
      stpf = fem::max(stpmin, stpf);
    }
    //C
    //C     Fourth case: A lower function value, derivatives of the same sign,
    //C     and the magnitude of the derivative does not decrease. If the
    //C     minimum is not bracketed, the step is either stpmin or stpmax,
    //C     otherwise the cubic step is taken.
    //C
  }
  else {
    if (brackt) {
      theta = three * (fp - fy) / (sty - stp) + dy + dp;
      s = fem::max(fem::abs(theta), fem::abs(dy), fem::abs(dp));
      gamma = s * fem::sqrt(fem::pow2((theta / s)) - (dy / s) * (dp / s));
      if (stp > sty) {
        gamma = -gamma;
      }
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p / q;
      stpc = stp + r * (sty - stp);
      stpf = stpc;
    }
    else if (stp > stx) {
      stpf = stpmax;
    }
    else {
      stpf = stpmin;
    }
  }
  //C
  //C     Update the interval which contains a minimizer.
  //C
  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  }
  else {
    if (sgnd < zero) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }
  //C
  //C     Compute the new step.
  //C
  stp = stpf;
  //C
}

} // namespace placeholder_please_replace
