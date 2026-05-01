# Physics Regression Baselines

This directory holds pickled outputs from `src.protocal.Protocal.evolve()`
at fixed parameters. They define the frozen physics behavior the refactor
must preserve.

## Files

| File | Source |
|---|---|
| `qubit_static.pkl` | TransmonQubit basic properties |
| `ramsey_default.pkl` | `Protocal(type=1).evolve()` |
| `diff_echo_default.pkl` | `Protocal(type=2).evolve()` |
| `transient_default.pkl` | `Protocal(type=4).evolve()` |
| `cryoscope_default.pkl` | `Protocal(type=5).evolve()` |
| `predistortion_default.pkl` | P4 `PredistortionValidationWorkflow.run()` |
| `z_crosstalk_default.pkl` | P5 `ZCrosstalkWorkflow.run()` |

## Regenerating

Only regenerate when intentionally changing physics. Always justify the
change in the commit message.

```bash
python -m tests.regression.generate_baselines
```

Regression tests use `rtol=1e-6, atol=1e-9`, defined in
`tests/conftest.py`.
