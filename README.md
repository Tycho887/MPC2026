# MPC 2026 report monorepo

The repository is imported directly into Overleaf. The main LaTeX document is
[mainReport.tex](mainReport.tex) at the repository root. Select it as the main
document in Overleaf and use the pdfLaTeX compiler.

| Path | Contents |
| --- | --- |
| `mainReport.tex`, `mcode.sty` | Report entry point and MATLAB listing style |
| `tex/`, `preamble/`, `ref/` | Chapters, report setup, and bibliography |
| `figures/`, `gfx/` | Report illustrations and template graphics |
| `tex/chapters/results/` | Committed figure snapshots used by the report |
| `code/` | MATLAB shared functions, problem experiments, and validation checks |
| `documents/` | Exam assignment PDF |

## Run MATLAB experiments

From the repository root in MATLAB:

```matlab
run('code/problem02/simulate.m')
run('code/problem03/simulate.m')
run('code/problem04/simulate.m')
```

Results are saved under `code/problemXX/results/` (ignored by Git). Existing
results were copied from the pre-merge workspace. See [code/README.md](code/README.md)
for model interfaces, toolbox requirements, and experiment coverage.

To validate:

```matlab
addpath('code/checks')
validateMigration
validateMigration(true) % also checks entry scripts and exported artifacts
```

## Build and update the report

Overleaf builds from the committed sources and figure snapshots; MATLAB is not
needed to compile the report. For a local build, run from the repository root:

```sh
latexmk -pdf mainReport.tex
```

To update simulation figures, copy the selected PDFs from
`code/problemXX/results/` into `tex/chapters/results/problemXX/`, review them,
and commit them before syncing to Overleaf. The report always uses these
snapshots so local and Overleaf builds use the same figures.

## Migration

The merged report remains the authoritative text. Earlier problem 2 and 3 TeX
drafts are preserved in `tex/drafts/pre-merge/` for reference and are not included
in the report build. MATLAB files, existing experiment results, and the assignment
were copied from `MPC-Exam-pre-merge` without importing its nested Git repository.
The MATLAB validation reference path was adjusted for `code/`.
`MPC` and `MPC-Exam-pre-merge` remain untouched; the monorepo does not require them
except for optional validation comparisons with the original lecture code.
