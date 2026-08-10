# RFG-RXi: A Background-Preserving Foliation Completion

RFG-RXi is a new R-Universe action, distinct from the original regular RFG-R
reference action. It is motivated by the exact DAE result: the reference
branch has a finite-wavenumber zero of its curvature-constraint coefficient.
The extension is not fitted to LambdaCDM or to a CMB data product.

## Action

Let sigma_ij=K_ij-K gamma_ij/3. The cosmological action is

    S_Xi = S_RFG-R
     + M_Pl^2/2 int d^4x N sqrt(gamma) Xi [R3+sigma_ij sigma^ij],
    Xi=constant.

The operator is a lowest-spatial-derivative scalar of the preferred
foliation. Both R3 and sigma_ij vanish on flat FLRW. The RFG-R background
reconstruction and Friedmann equation remain exactly unchanged. The tensor
kinetic and gradient receive the same addition:

    Q_T = Q + Xi,    c_T^2=1.

The existing smooth Weyl switch multiplies the complete cosmological
Lagrangian, including this operator. In its local constant-zero region the
action and variations are exactly GR, so the stated local PPN matching remains
unchanged.

## Direct quadratic derivation

For one real Fourier mode, the independently symbolic ADM calculation gives

    Delta L_mode =
     Xi k^2 [6a^2 alpha zeta+3a^2 zeta^2+k^2 beta^2]/(3a).

After the same mode normalization used by the finite action and
s=k^2 beta/a^2,

    Delta L/a^3 = Xi [2(k^2/a^2) alpha zeta
                      +(k^2/a^2) zeta^2+s^2/3].

Thus the only finite matrices changed are

    Delta A_ss = 2Xi/3,
    Delta D_alpha,zeta = 2Xi k^2/a^2,
    Delta M_zeta,zeta = 2Xi k^2/a^2,

with exact rates dot(Delta D)=dot(Delta M)=-2H Delta D. The symbolic
identity is checked by validate_rfg_xi_completion.py.

Keeping the scalar traceless metric deformation until after variation gives
the full metric increments, in the same real-mode normalization,

```text
Delta E_lapse      =  2 Xi a k^2 zeta,
Delta E_shift      =  2 Xi beta k^4/(3a),
Delta E_trace      = -2 Xi a k^2 (alpha+zeta),
Delta E_traceless  = -2 Xi a k^4 (H beta+alpha+beta_dot+zeta)/3.
```

Consequently the action-defined traceless equation is

```text
dot(s) = -[3H + Qdot/(Q+Xi) + Q_X k^2/(3(Q+Xi)a^2)] s
         -(k^2/a^2)(alpha+zeta)
         +(Q_X/(Q+Xi))(k^2/a^2)(H alpha-dot(zeta))
         -Pi/(Q+Xi).
```

This is not the Einstein equation: the completion changes every scalar metric
Euler residual listed above. `validate_rfg_xi_metric_equations.py` derives
the four increments symbolically and verifies the executable shear equation
against the exact `Q -> Q+Xi` identity.

## Reference audit

Xi is a new dimensionless EFT coefficient. It is not inferred from any data in
this package. The unfitted normalization benchmark is Xi=1; Xi=2 is also
audited to show that the reported no-root result is not tied to a single
isolated benchmark. Both are constants independent of scale factor,
wavenumber, and observational data. Run:

    cd r_universe_completion
    PYTHONPATH=scripts python3 scripts/rfg_xi_completion.py
    PYTHONPATH=scripts python3 scripts/validate_rfg_xi_completion.py

The generated rfg_xi_completion_audit.csv contains both Xi=1 and Xi=2, each
on 49 scale factors in [1e-7,1] and 81 wavenumbers in k/H0 in [1e-4,1e6]. The
calculation returns:

    sampled curvature-constraint roots      = {1.0: 0, 2.0: 0}
    min Q_tensor                            = 1.736101654462e+00
    min |det A_(alpha,s)|                   = 5.991445560907e+00
    max transformed zeta kinetic residual   = 7.452e-15
    max reduced B antisymmetry residual      = 5.135e-10

The finite material kinetic block has inertia (4,0,0) at every audited point.
The final antisymmetry value is IEEE double-precision cancellation at the
smallest wavenumbers; the direct symbolic action has an exactly symmetric
mixing matrix.

## Scope

This removes the specific DAE singularity of the original RFG-R reference
action on the stated audit grid while preserving the R-Universe homogeneous
branch and luminal tensor speed. It is a mathematically explicit new theory,
not an empirical victory. It still requires a full photon--baryon--CDM--
neutrino initial-condition derivation, a converged hierarchy solver, gradient
and nonlinear stability tests, and a joint likelihood before any statement
about observations or comparison with LambdaCDM is possible.
