# Matter and local-gravity gates

## Linear matter response

The reconstructed action has `alpha_M=alpha_T=0` and `c_s^2=1`.  In the
subhorizon quasistatic limit its scalar response is fixed, rather than fitted:

\[
\mu_\infty(a)=1+\frac{\alpha_B^2}{2D c_s^2},\qquad
\eta_\infty(a)=1,\qquad \Sigma_\infty(a)=\mu_\infty(a).
\]

`matter_qs.py` integrates

\[
D_{+,NN}+\left(2+\frac{E_N}{E}\right)D_{+,N}
-\frac32\Omega_m(a)\mu_\infty(a)D_+=0
\]

from `a=10^-3` with the same matter-era initial condition as the constant-R
LCDM comparison.  This is an action-derived late-time prediction.  It is not
an angular CMB calculation and is not presented as one.

## Solar-system screening

At the present point, the cubic action coefficient is

\[
\Lambda_3^3=\frac{M_{\rm Pl}H_0^2}{|\hat C_0|}.
\]

The linear coupling is defined by the exact quasistatic normalization

\[
2\beta_{\rm eff}^2=\mu_\infty-1
=\frac{\alpha_{B0}^2}{2D_0}.
\]

The point-source decoupling-limit Vainshtein scale is

\[
r_V^3=4\pi\beta_{\rm eff}|\hat C_0|\frac{r_s}{(H_0/c)^2}.
\]

Inside this radius, the conservative cubic-Galileon envelope is

\[
|\gamma-1|\le4\beta_{\rm eff}^2\left(\frac{r}{r_V}\right)^{3/2}.
\]

`ppn_screening.py` evaluates this expression at a Cassini impact parameter
`1.6 R_sun` and compares it against `|gamma-1| < 2.3e-5`.

The result is an analytic screening consistency gate for the exact action;
an ephemeris re-fit is a different, data-intensive likelihood calculation and
is deliberately not mislabelled as completed here.
