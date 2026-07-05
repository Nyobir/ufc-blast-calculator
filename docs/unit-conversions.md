# Unit Conversions — UFC Imperial to SI

All UFC 3-340-02 charts use imperial units. All internal storage and computation uses SI.

## Conversion Factors

| Parameter | Imperial (UFC) | SI (stored) | Factor | Notes |
|-----------|---------------|-------------|--------|-------|
| Pressure (Pr, Ps0) | psi | kPa | x 6.89476 | |
| Scaled distance (Z) | ft/lb^(1/3) | m/kg^(1/3) | x 0.396698 | |
| Scaled time (tA, t0) | ms/lb^(1/3) | ms/kg^(1/3) | x 1.3015 | **See trap below** |
| Scaled impulse (ir, is) | psi-ms/lb^(1/3) | kPa-ms/kg^(1/3) | x 8.974 | **See trap below** |

## The Denominator Trap

The scaling parameter W^(1/3) appears in the **denominator**: e.g., `tA/W^(1/3)`.

When converting lb^(1/3) to kg^(1/3) in the denominator:
```
1 lb = 0.453592 kg
1 lb^(1/3) = 0.7683 kg^(1/3)

tA [ms/lb^(1/3)] → tA [ms/kg^(1/3)]
Factor = 1/0.7683 = 1.3015   (NOT 0.7683!)
```

Getting this inverted makes all time and impulse values wrong by a factor of ~1.69x.

## Impulse Factor Derivation

```
impulse_factor = KPA_PER_PSI / Z_FACTOR
              = 6.89476 / 0.396698 * (0.396698)  ... wait, let's be precise:

psi-ms/lb^(1/3) → kPa-ms/kg^(1/3)
= 6.89476 [kPa/psi] × 1.3015 [lb^(1/3)/kg^(1/3)]
= 8.974
```

The Excel spreadsheet used 8.953 (close enough for digitization accuracy). The theoretically correct value is 8.9735.

## Unscaling to Physical Units

When computing final results from scaled table values:
```python
tA = tA_scaled * W**(1/3)       # ms (table stores ms/kg^(1/3))
t0 = t0_scaled * W**(1/3)       # ms
ir = ir_scaled * W**(1/3)       # kPa-ms (displayed in kPa-ms, NOT kPa-s)
```

## Key Constants in Code

```python
KPA_PER_PSI = 6.89476
Z_FACTOR = 0.396698               # ft/lb^(1/3) → m/kg^(1/3)
TIME_FACTOR = 1.3015              # ms/lb^(1/3) → ms/kg^(1/3)
IMPULSE_FACTOR = 8.974            # psi-ms/lb^(1/3) → kPa-ms/kg^(1/3)
HOB_THRESHOLD = 0.397             # m/kg^(1/3), free-air vs surface burst
```
