# Boltzmann implementation contract

The package supplies the complete background and Horndeski functions needed by
a linear Einstein--Boltzmann solver.  In the Bellini--Sawicki convention the
runtime functions are

\[
H(a)=H_0E(a),\qquad
\alpha_B(a),\qquad \alpha_K(a),\qquad
\alpha_M(a)=0,\qquad \alpha_T(a)=0,
\]

with standard minimally coupled baryons, CDM, photons, and neutrinos.  The
table `generated/ru_kgb_trajectory.csv` gives the exact values; the defining
functions are in `scripts/ru_kgb.py`, avoiding an interpolation ambiguity.

For a native run, the solver must evolve the coupled photon-baryon hierarchy,
massive/massless neutrinos, metric constraints, and the KGB scalar using these
functions.  It must recompute recombination and the drag scale, then evaluate
the selected CMB, BAO, RSD, and lensing likelihoods jointly with the same
parameter priors as LCDM.

The CMB script in this package is intentionally only a pre-recombination gate:
it computes the background sound horizon and the scalar-response size. It is
not called a CMB likelihood because it does not evolve `C_ell`.
