# Friedlander Waveform — Formula, Units and Limits

## Positive-phase waveform

For a computed facade point, internal units are **seconds** for `tA` and
`t0`, **kPa·s** for reflected impulse `ir_alpha`, and **kPa** for reflected
pressure `Pr_alpha`. The CLI displays times in ms and impulse in both
kPa·ms and kPa·s; do not feed those displayed millisecond values back into
the solver without converting them.

```text
tau = (t - tA) / t0
P(t) = Pr_alpha * (1 - tau) * exp(-b * tau),  tA <= t <= tA + t0
P(t) = 0,                                      otherwise
```

The integrated positive-phase impulse is

```text
ir_alpha = Pr_alpha * t0 / b^2 * (b - 1 + exp(-b))
```

The `b` parameter is dimensionless. The solver attempts to match the
**reflected**, angle-dependent impulse, not the incident impulse.

## Solver and boundary flags

- `solve_friedlander_b(Pr_alpha, t0, ir_alpha)` uses Brent's method for one
  query, searching `b` in `[0.01, 50]`.
- `_solve_friedlander_b_batch()` uses bounded vectorized Newton iterations for
  a facade. Neither method silently claims a valid fit when the chart-derived
  triple has no root in the allowed bracket: it clips to a boundary and the
  computed result gets the `friedlander_b` flag. The scalar path also warns.
- In the Mach-front path, if the chart-derived duration cannot carry the
  reflected peak and impulse (`ir_alpha >= Pr_alpha*t0/2`), UFC's triangular
  fictitious duration is `t_rf = 2*ir_alpha/Pr_alpha` (Eq. 2-11). The code uses
  `t0 = 1.005*t_rf` so a finite positive `b` can be solved; the informational
  flag `t0_fictitious` records that adjustment. Do **not** describe the
  resulting finite-`b` waveform as an exact UFC triangle or as retaining the
  original chart-derived `t0`.
- Out-of-range chart inputs are clamped and flagged separately. `t0_fictitious`
  is informational; chart-limit and `friedlander_b` flags need attention when
  interpreting predicted curves.

## Distinct reference queries

For `W=200 kg`, `R=30 m`, `Hc=5 m`, normal incidence:

- `ufc-blast point` evaluates regular reflection and prints
  `Pr_alpha = 64.864 kPa`, `ir_alpha = 0.468400 kPa·s`,
  `tA = 54.598 ms`, `t0 = 19.818 ms`, `b = 1.0341`.
- `compute_facade(R=30, W=200, Hc=5, width=20, height=10, step=0.25)`
  applies the ground-evaluated Mach treatment below the triple-point
  boundary. Its 3,321-point grid has centre `(dx=0, dy=0)`
  `Pr_alpha ≈ 110.081271 kPa` and `ir_alpha ≈ 0.622702163 kPa·s`.
  This centre is **not** the regular point-command result.

These values are checks of the current method, not experimental blast data.
The model covers the positive phase only; it does not model the negative
(suction) phase or certify a structural component.
