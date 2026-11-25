# MapBiomas Brazil | Nighttime Radiance

## About

This repository provides the Python implementation of MapBiomas Brazil Nighttime Radiance Collection 1, Version 1.0, released within MapBiomas Brazil Collection 11. The product covers Brazil annually from 2013 to 2025. General coordination is attributed to MapBiomas Brazil, technical coordination to PLAN-B, and institutional authorship to the MapBiomas Brazil Nighttime Radiance Team, as documented in the ATBD.

- [Algorithm Theoretical Basis Document (ATBD), August 2026](https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/09/ATBD-General-Collection-11-versao-ALAN_verFINAL.pdf)
- [MapBiomas Brazil data downloads](https://brasil.mapbiomas.org/en/downloads/)
- [Earth Engine source dataset](https://developers.google.com/earth-engine/datasets/catalog/NASA_VIIRS_002_VNP46A2)
- [NASA source product DOI](https://doi.org/10.5067/VIIRS/VNP46A2.002)

## Methodology

The workflow uses direct daily VIIRS radiance from NASA Black Marble VNP46A2 Collection 2. Each calendar year is processed independently. Internal QA flags and radiance thresholds screen observations; a median absolute deviation filter in natural-log space removes outliers. Quality-weighted lognormal estimation produces the expected annual radiance in linear physical units (nW cm⁻² sr⁻¹). No gap-filled radiance, auxiliary land-cover data or interannual smoothing is used in the annual estimator.

The distributed 250 m grid is obtained using nearest-neighbour allocation from nominal 500 m source information. This is not spatial super-resolution. The product is continuous radiance, not a land-cover classification or direct ground-level illuminance measurement.

## Repository structure

- `mapbiomas_viirs_collection1_final_1to1.py`: original production implementation, including the band dictionary and numerical reference test.
- `requirements.txt`: direct third-party dependency.
- `tests/test_self_test.py`: local test entry point.
- `RELEASE_NOTES.md`: Collection 1 release description.
- `.gitignore`: excludes local environments, sensitive configuration files and generated products.

## Requirements

- Python 3 (local numerical tests checked with Python 3.13.5).
- An Earth Engine-enabled Google Cloud project and an authorized Earth Engine account.
- A readable Earth Engine FeatureCollection defining the study boundary, and its version descriptor.
- A writable Earth Engine asset folder and/or Google Drive export destination.

Install the dependency in an isolated environment:

```bash
python -m venv .venv
```

Activate `.venv` before installing packages (Windows PowerShell: `.venv\Scripts\Activate.ps1`; macOS/Linux: `source .venv/bin/activate`). After activation, run `python -m pip install -r requirements.txt`. The exact production Earth Engine API version should be recorded by the team; this file does not assume an unverified version pin.

## Local verification

```bash
python mapbiomas_viirs_collection1_final_1to1.py --self-test
python -m unittest discover -s tests -v
```

The embedded tests do not require Earth Engine credentials. They check configuration, the 18-band schema and a numerical reference case; they do not validate server-side processing or spatial alignment.

## Execution

Replace `YOUR_GEE_PROJECT_ID`, `YOUR_PROJECT` and `BOUNDARY_VERSION` with your own authorized resources. Run `--help` to view all options. The following example builds and checks a single annual product without submitting export tasks; it requires Earth Engine access:

```bash
python mapbiomas_viirs_collection1_final_1to1.py --project YOUR_GEE_PROJECT_ID --roi-asset projects/YOUR_PROJECT/assets/BRAZIL_BOUNDARY --boundary-version BOUNDARY_VERSION --years 2025 --destination drive --dry-run
```

To submit a Drive export, remove `--dry-run` and optionally provide `--drive-folder`. For Asset export use `--destination asset --asset-folder projects/YOUR_PROJECT/assets/OUTPUT_FOLDER`; `--destination combined` exports to both. Use `--wait` to monitor tasks. Select the complete series with `--years 2013-2025`.

The public default grid is EPSG:5880 at 250 m. For exact release reproduction, supply the canonical CRS and six-value affine transform using `--crs` and `--crs-transform` instead of `--scale`. Do not alter the Collection 1 scientific defaults for reproduction of the released series. The canonical grid transform, boundary version and production dependency versions must be obtained from the release team.

## Outputs and reproducibility

Each annual image contains 18 co-registered bands. The main radiance layer is `expected_radiance`; the remaining bands describe log-space distribution, empirical quantiles, sampling support, quality and conditional uncertainty. The code defines their names, order, units and NoData rules. Drive exports are Cloud-Optimized GeoTIFFs; the script also supports Earth Engine assets. The GeoTIFF NoData value is -9999.

The public run manifest records source revision, parameters, input dataset, grid and verification information. An optional controlled log can contain internal identifiers and must not be committed. Use `--visual-check-status passed` only after a genuine manual product inspection. The conditional intervals are not a complete measure of uncertainty. Interpret radiance alongside observation-support bands such as `n_eff` and `valid_fraction`. See the ATBD and its supplementary validation report for the scientific assessment.

## Citation and data use

Please cite MapBiomas Brazil Nighttime Radiance Collection 1, Version 1.0 (2026), and the corresponding ATBD. Consult the [MapBiomas terms of use](https://brasil.mapbiomas.org/en/termos-de-uso/) for data attribution. The software licence and individual contributor credits must be confirmed by the project team before public release.
