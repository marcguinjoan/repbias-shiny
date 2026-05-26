# REPBIAS Shiny App

**Unpacking Representation Bias: An Interactive Simulation Explorer**

Interactive Shiny application for the REPBIAS project (PID2021-128332NA-I00, PI: Marc Guinjoan), funded by the Spanish Ministry of Science and Innovation (MCIN/AEI/10.13039/501100011033/FEDER, UE).

Live app: https://mguinjoan.shinyapps.io/repbias/

---

## What this app does

The app lets researchers, students, and reviewers explore how three electoral institutions independently distort the ideological link between voters and legislators:

| Bias measure | Definition | Institution |
|---|---|---|
| **MRB** | IS − IP | Malapportionment |
| **SRB** | IR − IS | Seat-allocation formula + District Magnitude |
| **ThRB** | IR(threshold) − IR(base) | Electoral threshold |
| **TRB** | IR − IP | All institutions combined |

Where IP = ideology of population (population-weighted), IS = ideology of seats (seat-weighted), IR = ideology of representatives (seat-weighted after allocation).

Bias shares measure each institution's *relative weight*: `|MRB| / (|MRB| + |SRB| + |ThRB|)`.

---

## App structure

### Tab 1 — Explore Distributions
Draws from two pre-computed pools of 10,000+ simulations. Users filter by:
- Bias measure (MRB random/concentrated, SRB, ThRB)
- Polarization (Dalton index bands)
- Number of parties

Outputs: interactive scatter/density plot (plotly), marginal density strip, summary statistics table.

### Tab 2 — Scenario Designer
Runs **300 simulations on-the-fly** with user-defined institutional parameters:

| Parameter | Range | Notes |
|---|---|---|
| District Magnitude (DM) | 1–30 seats | All 100 districts identical |
| Malapportionment (MAL) | 0–0.45 | Samuels & Snyder (2001) index; deterministic construction |
| Malapportionment type | Random / Concentrated | Concentrated = one ideological bloc systematically over-represented |
| Electoral formula | D'Hondt, Sainte-Laguë, Mod-SL, Hare, Droop | Controls SRB |
| Electoral threshold | 0–20% | Applied uniformly per district |
| Number of parties | 2–6 or random | |
| Polarization | Random / Low / Medium / High | Controls Beta distribution shape for party ideologies |

Key design note on MAL: populations are constructed *deterministically* so that exactly the requested MAL value is achieved every run. With 100 equal-DM districts, seat shares are uniform (1/100), so setting pop_over = (1 − 2·MAL)·base and pop_under = (1 + 2·MAL)·base for 50 districts each yields exactly the target MAL.

Key formula: **MRB = MAL × (avg ideology of over-represented districts − avg ideology of under-represented districts)**. This means concentrated malapportionment can produce sizable MRB shares even at low MAL when SRB is small (e.g., high DM, many parties with balanced votes).

Outputs: distribution boxplot + jitter per institution, |TRB| density strip, summary statistics, collapsible scenario history table.

Buttons:
- **▶ Run**: generates new vote-share simulations + applies rules
- **↺ Refresh**: reuses same simulations, re-applies rules (fast — use when only changing DM, MAL, formula, threshold)
- **🎲 Randomize**: randomizes all parameters and runs immediately

---

## File structure

```
7_Shiny/
├── app_v2.R                    ← Main Shiny application
├── pool_mrb_generation.R       ← Generates pool_mrb_summary.rds (Tab 1, MRB/TRB)
├── pool_srb_generation.R       ← Generates pool_srb_districts.rds (Tab 1, SRB/ThRB)
├── pool_scenario_generation.R  ← Legacy pool script (no longer used by app)
├── run_batch.R                 ← Helper to run pool generation in batches
├── data/
│   ├── pool_mrb_summary.rds    ← NOT in git (regenerate: pool_mrb_generation.R)
│   └── pool_srb_districts.rds  ← NOT in git (regenerate: pool_srb_generation.R)
└── README.md
```

---

## Simulation design (shared across app and paper)

- **Territory**: 100 synthetic districts on a 10×10 grid
- **Ideology scale**: 0 (left) to 10 (right)
- **Party ideologies**: Beta(2.5, 2.5) × 10 — bell-shaped, centrist tendency
- **Vote shares**: Dirichlet(α = 1.5) — realistic district-level variation
- **Seat allocation**: D'Hondt by default; Tab 2 also offers Sainte-Laguë, Modified Sainte-Laguë (1.4), Hare (LR), Droop (LR)
- **Polarization (Dalton index)**: Σ vᵢ · |xᵢ − x̄|
- **ENEP (Laakso-Taagepera)**: 1 / Σ vᵢ²
- **MAL (Samuels & Snyder)**: 0.5 · Σ |sᵢ − pᵢ|

---

## Running locally

```r
# Install dependencies
install.packages(c("shiny", "tidyverse", "plotly", "gtools"))

# Generate data pools (one-time; takes several hours each)
source("pool_mrb_generation.R")   # produces data/pool_mrb_summary.rds
source("pool_srb_generation.R")   # produces data/pool_srb_districts.rds

# Launch app
shiny::runApp("app_v2.R")
```

The Scenario Designer tab works **without** the pool files (it computes on-the-fly). The Explore Distributions tab requires both `.rds` files.

---

## Deploying to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name="mguinjoan", token="...", secret="...")
rsconnect::deployApp(appDir = ".", appName = "repbias", appPrimaryDoc = "app_v2.R")
```

---

## Relationship to the paper

This app is the interactive companion to the paper *Unpacking Representation Bias: A Simulation Approach to Electoral Institutions*. The paper's static figures come from:

| Folder | Content |
|---|---|
| `4_Part 1/` | MRB analysis — 10,000 simulations |
| `5_Part 2A/` | SRB/ThRB by formula and DM |
| `5_Part 2B/` | Bias Share decomposition |
| `7_Shiny/` | This interactive app |

The LaTeX manuscript is at `C:\Users\1407091\Dropbox\Apps\Overleaf\Simulations malapportionment\main_v2.tex`.

---

## For AI assistants

If you are helping with this project, key things to know:

- All code and comments must be in **English**
- Save edits as `_v2` files; do not overwrite originals
- Use `detectCores() - 3` for parallelization in R
- R simulation files use `knitr::knit()` (no pandoc), not `rmarkdown::render()`
- The app file is `app_v2.R` (not `app.R`)
- Data pool files (`.rds`) are not in git — they must be generated locally
- The Scenario Designer computes on-the-fly; no pool file needed
- `MRB = IS − IP` (malapportionment effect), `SRB = IR − IS` (seat formula effect), `ThRB = IR(th) − IR(base)` (threshold effect)
- Bias share = |component| / (|MRB| + |SRB| + |ThRB|), computed per simulation then averaged
- High MRB share at low concentrated MAL is *correct* behavior when SRB is small (large DM, balanced multiparty votes)
