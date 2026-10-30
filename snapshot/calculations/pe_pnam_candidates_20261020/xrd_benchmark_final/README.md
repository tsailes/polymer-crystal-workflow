# Tashiro PE: position-only lattice consistency check

Compared 4 structures/cells with 45 published observed d spacings.
Target group: P n a m (Hall -P 2c 2n).
Distinct cell/group tables: 1. Q proximity threshold: 0.03 inverse Angstrom.

This output does not rank molecular packings or report a powder-profile R factor.
At fixed cell and target group, orientation/translation changes leave these allowed positions unchanged.

Each candidate folder contains allowed_positions.csv, observed_nearest_allowed.csv and result.json.
Signed residuals are predicted minus observed. Nearest-Q matches can reuse a predicted position.
target_symmetry_check verifies the supplied operations in the supplied origin/axis setting.
If that check fails, the allowed-position table is conditional on the requested target group and is not the candidate's verified extinction pattern.

Limitations:
- Only published, processed fiber-diffraction positions are compared; this is not a continuous powder profile.
- The target space group and cell are supplied priors, not inferred from these observations.
- Position matches contain no intensity, ADP, atomic-scattering-factor or refined-reference-coordinate information.
- At identical cell and target group, all packings have identical allowed lattice positions and position residuals.
- The reported Q tolerance is a descriptive threshold, not an experimental standard deviation or statistical confidence interval.
- Matching is observed-to-nearest-allowed-position, permits many-to-one assignments, and does not prove indexing.
- No penalty is assigned to unobserved predicted positions: the published list is not known to be complete.
- No full-pattern R factor, chi-square, probability, or packing ranking is produced.
- Extra structure-specific intensity cancellations and experimental detection thresholds are not modeled.
