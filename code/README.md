# Model Predictive Control — 02619

MATLAB experiments for the modified four-tank system. The assignment is in
[../documents/MPCexam2026.pdf](../documents/MPCexam2026.pdf).

## Run the experiments

From the `code/` directory, run these commands in MATLAB:

```matlab
run('problem02/simulate.m')
run('problem03/simulate.m')
run('problem04/simulate.m')
```

You can also change into a problem directory and run `simulate`. Each entry
script resolves its own paths, computes its dependencies, displays its figures,
and saves PNG (200 dpi), vector PDF, and MAT files in that problem's `results/`.
Figures use a light theme for readable report exports. Rerunning replaces files with the same names. Existing unrelated figures and
workspace variables are not cleared. Generated results are ignored by Git.

MATLAB with Control System Toolbox is required for problems 3–4 (`ss`, `tf`,
`pidtune`, `c2d`). Problem 2 needs only MATLAB. System Identification Toolbox is
optional: identification uses `tfest` when available and licensed, otherwise
continuous-time least-squares fitting with `fminsearch`.

To select the identification backend explicitly:

```matlab
addpath('common', 'problem04')
[identified, figures] = problem04Identification('least_squares'); % or 'toolbox'
[markov, markovFigures] = problem04Markov(identified);
```

## Shared code and interfaces

`common/` holds the physical parameters, routing, mass balance, observation and
output models, diffusion, equilibrium, linearization, integration, PI tuning,
identification, and Markov calculations. Problem functions hold experiment
configuration and plotting. Change physical defaults in `fourTankParameters`;
change operating points, solver settings, references, and noise in the relevant
experiment function. No class or positional parameter vector is required.

The parameter struct has fields `a` (four outlet areas), `A` (four tank areas),
`gamma` (two valve fractions), `g`, and `rho`. Units are centimetres, grams, and
seconds; pump and disturbance flows are cm³/s. States are the four tank masses.
Measured and controlled outputs are the levels of tanks 1 and 2; all four tank
levels are available separately.

- `fourTankDynamics(t, x, u, d, p)` is the shared drift for all three models.
- `fourTankOutput(x, p)` returns true controlled levels. `fourTankSensor(x, p, noise)`
  adds explicitly supplied measurement noise in centimetres.
- `fourTankSteadyState(u, d, p)` returns the operating point;
  `linearizeFourTank(op, p)` returns numeric Jacobians. Linear models use deviations
  around `op`, with both mass and level output matrices.
- `simulateFourTank(t, x0, u, d, p, options)` returns `t`, `x`, `h`, `y`, `z`, `u`,
  `d`, `p`, and `options`. Time is N×1; states/levels are 4×N; flows and lower-tank
  outputs are 2×N. Inputs and disturbances are held over each interval.
- Solver choices are `ode15s`, `ode45`, and `sde`. SDE options supply `sigma` (4×2)
  and Brownian increments `dW` (2×(N−1), already scaled by √Δt). Measurement noise
  is supplied as a 2×N array. Callers seed and draw randomness; the saved `seed`
  records that configuration. Constant diffusion makes the Milstein correction
  zero, so the update equals Euler–Maruyama with nonnegative mass projection.
- `markovParameters(G, Ts, count)` returns `H(:,:,1)=D` at k=0, then
  `H_k=C*A^(k−1)*B`. These are discrete impulse coefficients, without MATLAB's
  impulse-response 1/Ts scaling.

MAT files store a variable named `result`, including configuration, trajectories,
and computed models/controllers where relevant. Experiment collections use cell
arrays. Problem 3's comparison additionally stores absolute and deviation mass
states and both independently computed input histories; its reference and noise
settings are in grams, while `y` and `z` remain levels.

## Migration mapping and deliberate corrections

Sources are from `../../MPC/four-tank-model`; the source directory is untouched.

| Previous script(s) | Migrated experiment |
| --- | --- |
| `problem2_4_1.m`, `problem2_4.m`, `problem2_4_3.m` | `problem02Experiments`: deterministic, held stochastic disturbances, SDE |
| `problem3.m` | `problem03Analysis`: transfer functions and mass-coordinate PI tuning |
| `problem3_1.m` | `problem03Comparison`: independent noisy linear/nonlinear mass-feedback loops |
| `problem3_2.m` | `problem03LevelControl`: level-coordinate PI tuning and nonlinear control |
| `problem4_1.m` | `problem04Deterministic`: separate 10%, 25%, 50% pump steps |
| `problem4_2.m` | `problem04Stochastic`: the same steps at 0.5×, 1×, 2× noise |
| `problem4_4.m` | `problem04Identification`: normalized noisy 25% steps and fitted channels |
| `problem4_6.m` | `problem04Markov`: coefficients computed from those fitted channels |

The deterministic problem 2 experiment now has zero measurement noise.
Stochastic runs use explicit seeds; problem 4 noise comparisons reuse the same
Brownian and measurement realizations. Problem 3's comparison uses separate
controller states and feedback for each plant, and recomputes gains for its own
operating point instead of importing rounded gains from another experiment.
Problem 4's Markov calculation requires identified models, with no example
coefficient fallback or hidden workspace dependency.

Operating points are preserved, rather than standardized:

| Experiment | Pump flows | Mean disturbance flows |
| --- | --- | --- |
| Problem 2 and problem 4 | [300; 300] | [20; 20] |
| Problem 3 analysis | [100; 300] | [100; 100] |
| Problem 3 comparison | [300; 300] | [0; 0] |
| Problem 3 level control | [300; 300] | [100; 100] |

The original grids, step sizes, reference levels, and actuator bounds are
preserved. The comparison uses forward Euler with 0.5 s steps and pump bounds
[0,1000]; level control uses Tustin-discretized PI controllers, 2 s sampling,
`ode45`, and a lower pump bound of zero. Neither experiment introduces anti-windup.

## Coverage and validation

This is a migration of existing experiments. Problem 3 currently implements PI
control; the P and PID experiments are still missing. Problem 4 additionally uses
`problem04Report` to plot normalized responses for every step/noise combination
and measure fitted-model accuracy against separate deterministic experiments.
Its `normalized_steps.mat` contains the gains, response RMSEs, and Markov-tail
fractions used in the report. Later exam problems are not implemented here.
Old Python scaffolding, environments, generated images,
and the malformed scratch `model.m` are not imported.

For numerical and artifact checks, run:

```matlab
addpath('checks')
validateMigration
validateMigration(true) % also run all three entry scripts from root and locally
```

The extended check writes the usual ignored results and temporarily hides and
closes only the figures created during validation. Optional reference comparisons
use the old scripts when `../../MPC/four-tank-model` is available. Expected changes
from corrected feedback, gains, and noise settings are not required to match
legacy trajectories.
