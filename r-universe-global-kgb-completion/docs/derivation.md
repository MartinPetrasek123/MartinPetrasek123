# Global R-alpha KGB completion

## Defined covariant theory

In reduced Planck units, the action is

\[
S=\int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R
+A(\phi)X+B(\phi)X^2-V(\phi)-C(\phi)X\Box\phi\right]+S_m[g,\psi_m],
\qquad X=-\frac12(\nabla\phi)^2.
\]

Matter is minimally coupled and the Ricci coefficient is constant. Therefore

\[
Q_T=M_{\rm Pl}^2,\qquad c_T^2=1,\qquad
d_L^{\rm GW}/d_L^{\rm EM}=1.
\]

The field coordinate is a definition of the action functions, not a gauge
condition imposed off shell:

\[
\phi/M_{\rm Pl}=N=\ln a,\qquad a=e^{\phi/M_{\rm Pl}}.
\]

For every real `phi`, `scripts/ru_kgb.py` defines the unique expanding root
`E(phi)` and then the four functions `A(phi), B(phi), C(phi), V(phi)`.
The CSV is only an audit output; it is not the definition of the theory.

## Global finite-window background

For every observed epoch, `a <= 1`, this is exactly the fitted R-alpha law:

\[
E^2=\Omega_m a^{-3}+\Omega_r a^{-4}+r,
\qquad r=\Omega_{R0}E^{\alpha a}.
\]

The old formula has no prescribed infinite-future completion.  Define a smooth
function `a_eff(a)` that equals `a` for `a <= 1`, equals `a_sat=2` for
`a >= 2`, and joins the two intervals with a compact C-infinity step.  The
global equation is

\[
r=\Omega_{R0}E^{\alpha a_{\rm eff}(a)}.
\]

It leaves all past-light-cone likelihood predictions unchanged.  Because
`alpha a_sat < 2`, it tends to the unique de Sitter root

\[
E_\infty=\Omega_{R0}^{1/(2-\alpha a_{\rm sat})}.
\]

Writing `x=alpha a_eff`, differentiation gives

\[
E_N=\frac{r x_N\ln E-3\Omega_m a^{-3}-4\Omega_r a^{-4}}
{2E-rx/E},\qquad
r_N=r\left(x_N\ln E+x\frac{E_N}{E}\right).
\]

The denominator is checked numerically on the entire released domain.

## Stable scalar prescription

The KGB braiding is fixed without an additional fitted cosmological parameter:

\[
\alpha_B=-b_{\rm rad}+(1+b_{\rm rad})
\frac{x}{1+x}\Omega_R,
\qquad b_{\rm rad}=\alpha\Omega_{r0}.
\]

The small radiation floor is determined by already present R-alpha quantities.
Its negative sign is fixed by the scalar gradient condition in radiation
domination, where a positive constant braiding would give a negative gradient
numerator. It vanishes in the exact `alpha=0` LCDM limit, but keeps the scalar
kinetic coefficient nonzero instead of leaving a vanishing-Qs branch.

For luminal KGB coupled to conserved matter and radiation, the standard
Bellini--Sawicki coefficients, with `alpha_M=alpha_T=0`, are

\[
Q_s/M_{\rm Pl}^2=\frac{2D}{(2-\alpha_B)^2},
\]

\[
\mathcal N_s=(2-\alpha_B)
\left(-\frac{E_N}{E}+\frac{\alpha_B}{2}\right)
+\alpha_{B,N}-(3\Omega_m+4\Omega_r),
\qquad c_s^2=\frac{\mathcal N_s}{D}.
\]

The matter/radiation enthalpy is necessary: in the `alpha_B=0` canonical
scalar limit it returns `c_s^2=1`, rather than incorrectly counting fluid
enthalpy as scalar gradient energy.

The scalar sector is closed by choosing

\[
D=\frac{\mathcal N_s}{c_\star^2},\qquad
\alpha_K=D-\frac32\alpha_B^2,\qquad c_\star^2=1.
\]

Thus `N_s>0`, `D>0`, `Q_s>0`, and `c_s^2=1` are identities once the remaining
background gate `N_s>0` is checked.  `validate_ru_kgb.py` checks that gate
over `10^-8 <= a <= 10^3`.

## Reconstruction of action functions

Let `X=E^2/2` and `C=alpha_B/E^2`.  Along `phi/Mpl=N`, `C_phi=C_N`.
The target density and pressure are

\[
\rho_\phi=3r,\qquad
p_\phi=-(2EE_N+3E^2)-\Omega_r a^{-4}.
\]

Define

\[
R_\rho=\rho_\phi-6E^2XC+2C_\phi X^2,
\qquad
R_p=p_\phi+2C_\phi X^2+2CXEE_N.
\]

The exact action functions are

\[
B=\frac{\alpha_K E^2-(R_\rho+R_p)+8C_\phi X^2-12E^2XC}{8X^2},
\]

\[
A=\frac{R_\rho+R_p-4BX^2}{2X},\qquad
V=\frac{R_\rho-R_p}{2}-BX^2.
\]

Substitution gives the density and pressure of the stated covariant action
exactly.  Since `phi_dot=H Mpl` never vanishes on the branch, the verified
continuity equation is equivalent to the homogeneous scalar equation of
motion by the Bianchi identity.

## Gates that this package closes, and what a likelihood is

The package closes the mathematical background, global-future, action-density,
scalar no-ghost, scalar no-gradient, scalar sound-speed, tensor-speed, and
tensor-friction gates for this explicit action.  It does not manufacture a
Planck/DESI joint posterior.  That separate empirical calculation requires
compiling the derived KGB functions into a Boltzmann solver and evaluating the
published CMB and large-scale-structure likelihoods.  No numerical likelihood
value is claimed until that executable run exists.
