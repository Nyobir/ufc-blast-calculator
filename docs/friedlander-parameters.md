# Friedlander Waveform — Formulas & Parameters

## The Friedlander Equation

Pressure time-history at a point on the facade:

```
P(t) = Pr_alpha * (1 - tau) * exp(-b * tau)     for tA <= t <= tA + t0
P(t) = 0                                         otherwise

where tau = (t - tA) / t0    (dimensionless time, 0 to 1)
```

### Parameters
| Symbol | Name | Source |
|--------|------|--------|
| Pr_alpha | Peak reflected pressure (kPa) | Cα(Ps0, α) x Ps0 from Figures 2-193 + 2-7 |
| tA | Arrival time (ms) | Figure 2-7 lookup by Z |
| t0 | Positive phase duration (ms) | Figure 2-7 lookup by Z |
| b | Decay coefficient (dimensionless) | Solved numerically from impulse |
| ir_alpha | Reflected impulse (kPa-ms) | Figure 2-194 lookup by (Ps0, α) |

## Impulse Integral (Analytical)

Integrating the Friedlander equation from 0 to t0:

```
ir = (Pr * t0) / b^2 * (b - 1 + exp(-b))
```

This is a **transcendental equation** in b — cannot be solved algebraically.

## b-Solver

### Single Point (Brent's method)
```python
def residual(b):
    return (Pr * t0) / b**2 * (b - 1 + exp(-b)) - ir

b = scipy.optimize.brentq(residual, 0.01, 50.0)
```

### Batch (Vectorized Newton's method)
Used by `compute_points_batch()` for 1000+ points simultaneously:
```python
b = np.full(N, 2.0)  # initial guess
for _ in range(100):
    eb = np.exp(-b)
    f  = (Pr * t0) / b**2 * (b - 1 + eb) - ir
    df = (Pr * t0) / b**3 * (-b + 2 - (2 + b) * eb)
    b  = b - f / df
```
Converges in ~10 iterations. **Important**: the impulse used is the **reflected** impulse ir_alpha (angle-dependent), not the incident impulse.

## Confirmed Test Values

From simulation validation categories:

| Category | Ps0 (kPa) | t0 (ms) | b | Shock velocity U (km/s) |
|----------|-----------|---------|---|------------------------|
| A (close) | 465 | 18.6 | 6.4 | 0.6659 |
| B (medium) | 270 | 14.8 | 4.8 | 0.5673 |
| C (far) | 165 | 14.5 | 4.2 | 0.5019 |

## Only Positive Phase

Currently only the positive phase is modeled. Negative phase Friedlander (suction) is not implemented.
