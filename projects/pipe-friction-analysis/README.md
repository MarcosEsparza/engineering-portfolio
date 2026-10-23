# Pipe Friction Reanalysis in Python

A Python reanalysis of my Fall 2025 fluid friction lab (TecQuipment H408). I took the raw appendix data from the [original report](../../assets/docs/lab03-fluid-friction.pdf) and rebuilt the analysis as a small tested package. It uses one friction-factor convention throughout, the Colebrook equation instead of reading values off the Moody chart, and propagated measurement uncertainty.

**Tools:** Python (pandas, NumPy, SciPy, Matplotlib, pytest), GitHub Actions  
**Code:** [View on GitHub](https://github.com/MarcosEsparza/engineering-portfolio/tree/main/projects/pipe-friction-analysis)

[← Back to portfolio](../../)

---

## What the code does

1. Reads the raw readings (collection volume, time, and upstream/downstream piezometer levels) from [`data/h408_raw.csv`](data/h408_raw.csv).
2. Computes velocity, Reynolds number, head loss, and the **Darcy** friction factor from Darcy–Weisbach.
3. Compares each point with theory: **Blasius** for the smooth pipes and **Colebrook–White** (solved numerically with SciPy's `brentq`) for the smooth pipes and the rough pipe (ks/d = 0.036).
4. Propagates measurement uncertainty: ±1 mm per piezometer reading, ±0.5 s stopwatch timing, ±0.1 L collected volume. Since f ∝ h·t²/V², small head losses and short timings dominate.
5. Fits h ∝ uⁿ for each pipe as a physics check. Theory gives n ≈ 1.75 for smooth turbulent flow and n ≈ 2 for fully rough flow.
6. Unit tests check every relation against hand calculations from the report and known Moody chart values. GitHub Actions runs them on every push.

## Results

![Friction factor vs Reynolds number with uncertainty bars](results/friction_factor_vs_re.png)

![Head loss vs velocity on log-log axes](results/head_loss_scaling.png)

| Pipe | Re range | Measured f (avg) | Theory f (avg) | Measured / theory | Max uncertainty | Fitted n |
|---|---|---|---|---|---|---|
| 13.6 mm smooth | 9,800 – 16,900 | 0.051 | 0.028 | 1.8× | ±2% | 1.62 |
| 26.2 mm smooth | 4,200 – 5,900 | 0.142 | 0.037 | 3.9× | ±24% | 2.21 |
| 14 mm rough | 3,700 – 8,500 | 0.355 | 0.067 | 5.3× | ±8% | 1.99 |

## What I learned from redoing it

**13.6 mm smooth pipe.** The measured friction factor is about 1.8× the smooth-pipe prediction, and the offset holds steady across all five runs. Measurement uncertainty is only ±2%, so reading error can't explain it. A consistent offset like this points to something systematic, such as the effective length between tappings or an unaccounted-for minor loss. The original report reached the same conclusion (about 76% error).

**26.2 mm smooth pipe.** The head losses were only 6–13 mm of water, so each 1 mm piezometer reading causes ±12–24% uncertainty in f. The error bars make that visible, which the original tables didn't. Even so, the data sits about 4× above theory, so resolution alone doesn't explain the gap. With these flow rates, a better test would use a differential manometer or longer runs.

**14 mm rough pipe.** This was the most useful catch. The report found close agreement with the Moody chart, but the formula it used, f = hgd/(2Lu²), gives the **Fanning** friction factor. The Moody chart uses the **Darcy** factor, which is 4× larger. Comparing Darcy to Darcy, the measured values are about 5× the Colebrook prediction. The rest of the data still behaves well: the fitted exponent is n = 1.99, almost exactly what fully rough flow predicts. That points to a wrong geometric constant rather than bad readings. The recorded test length for this pipe is 0.2 m, while the other pipes use 0.912 m. With 0.912 m, the results land within about 17% of Colebrook. The next step would be to confirm the tapping spacing against the apparatus manual.

**Takeaway:** Stating the convention explicitly (Darcy vs. Fanning) and showing the uncertainty on every point turned two "unexplained error" results into specific, testable explanations.

## Run it yourself

```text
cd projects/pipe-friction-analysis
pip install -r requirements.txt
python analyze.py   # writes results/*.csv and results/*.png
pytest              # runs the unit tests
```

## Project layout

```text
pipe-friction-analysis/
├── analyze.py           # entry point: load data, compute, save tables and plots
├── data/h408_raw.csv    # raw readings transcribed from the lab report appendix
├── friction/
│   ├── core.py          # fluid relations: Re, Darcy–Weisbach, Blasius, Colebrook, uncertainty
│   └── analysis.py      # applies core.py to the dataset and summarizes each pipe
├── results/             # generated plots and CSV tables
└── tests/test_core.py   # pytest unit tests
```
