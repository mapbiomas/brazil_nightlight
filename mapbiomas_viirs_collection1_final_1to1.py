#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MapBiomas Brazil – VIIRS Night-time Radiance Collection 1
=========================================================

FINAL repository-safe Python implementation of the annual VIIRS production
workflow documented in:

    MapBiomas Brazil Nighttime Radiance Collection 1
    Algorithm Theoretical Basis Document (ATBD), Version 1.0

Scientific scope
----------------
* Source: NASA Black Marble VNP46A2 Collection 2 in Google Earth Engine.
* Science band: DNB_BRDF_Corrected_NTL (DIRECT, not gap-filled).
* Period: annual products, one independent calendar-year estimate per output.
* No land-cover predictors, demographic variables, OBIA, trend model,
  neural reconstruction, second-sensor fusion or synthetic super-resolution.
* The source VIIRS information support remains nominally 500 m. A 250 m
  output lattice is used only for MapBiomas grid compatibility.

Central estimator
-----------------
For each retained daily radiance R_i > 0, with R0 = 1 nW cm^-2 sr^-1:

    x_i = ln(R_i / R0)

Candidate observations are screened using VNP46A2 internal QA. A robust
median/MAD filter identifies inliers:

    I_i = 1  if |x_i - median(x)| <= 3.5 * max(1.4826*MAD, 1e-6)

Accepted observation reliability is:

    q_i = q_cloud_confidence * q_cloud_mask_quality

where confident-clear = 1.00, probably-clear = 0.70, low cloud-mask quality
= 0.80, and medium/high cloud-mask quality = 1.00. Mandatory quality must
already equal zero and therefore does not receive a second graded weight.

A transparent Collection-1 QA-derived nominal log-space uncertainty is:

    u_i = 0.03 + 0.30 * (1 - q_i)

and the effective central-estimator weight is:

    w_i = I_i * q_i / (1 + u_i^2)

The weighted log mean, Kish effective sample size and finite-effective-sample
corrected weighted log variance are then used to estimate annual expected
radiance under a lognormal observation model:

    expected_radiance = R0 * exp(mu_w + 0.5 * s_w^2)

The principal output is therefore NOT stored in logarithmic units. The
logarithm is an internal estimation domain; expected_radiance is exported
back in linear physical radiance units (nW cm^-2 sr^-1).

Outputs
-------
Exactly 18 co-registered bands are exported, matching the Collection-1 ATBD:

  1. expected_radiance
  2. weighted_mean_log_radiance
  3. log_variance
  4. log_standard_deviation
  5. lognormal_cv
  6. median_radiance
  7. p10_radiance
  8. p90_radiance
  9. n_eff
 10. n_candidate
 11. n_inlier
 12. candidate_fraction
 13. valid_fraction
 14. outlier_fraction
 15. mean_reliability
 16. mean_retrieval_age_days
 17. ci95_low
 18. ci95_high

Export
------
The public workflow supports exactly the ATBD release modes: Earth Engine
Asset, Google Drive, or combined Asset + Drive export. Google Drive output is
a Cloud-Optimized GeoTIFF; Earth Engine adds the .tif extension to the file(s).
The default public grid is EPSG:5880 at 250 m. For an operational MapBiomas
release, provide the canonical affine transform with --crs-transform whenever
available so pixel size, origin and alignment are frozen explicitly.

Security
--------
No credentials, OAuth tokens, service-account keys, private project IDs,
asset IDs or local personal paths are embedded in this file.
"""

from __future__ import annotations

import argparse
import calendar
import datetime as dt
import hashlib
import json
import math
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

try:
    import ee
except ImportError:  # Allows credential-free/local --self-test without earthengine-api.
    ee = None


# =============================================================================
# Fixed Collection-1 definitions
# =============================================================================

DEFAULT_GEE_PROJECT = os.getenv("EARTHENGINE_PROJECT", "YOUR_GEE_PROJECT_ID")
DEFAULT_ROI_ASSET = "projects/YOUR_PROJECT/assets/BRAZIL_BOUNDARY"

DATASET = "NASA/VIIRS/002/VNP46A2"
SOURCE_DOI = "10.5067/VIIRS/VNP46A2.002"
RADIANCE_BAND = "DNB_BRDF_Corrected_NTL"

SOURCE_BANDS = (
    RADIANCE_BAND,
    "DNB_Lunar_Irradiance",
    "Latest_High_Quality_Retrieval",
    "Mandatory_Quality_Flag",
    "QF_Cloud_Mask",
    "Snow_Flag",
)

OUTPUT_BANDS = (
    "expected_radiance",
    "weighted_mean_log_radiance",
    "log_variance",
    "log_standard_deviation",
    "lognormal_cv",
    "median_radiance",
    "p10_radiance",
    "p90_radiance",
    "n_eff",
    "n_candidate",
    "n_inlier",
    "candidate_fraction",
    "valid_fraction",
    "outlier_fraction",
    "mean_reliability",
    "mean_retrieval_age_days",
    "ci95_low",
    "ci95_high",
)

assert len(OUTPUT_BANDS) == 18

SOFTWARE_VERSION = "1.0"

# Appendix A, serialized into annual image metadata and the public manifest.
# All bands are exported as Float32 in this implementation. Earth Engine assets
# retain masks; GeoTIFF exports encode masked pixels as NODATA (-9999).
BAND_DICTIONARY = (
    {"band": "expected_radiance", "unit": "nW cm^-2 sr^-1", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Quality-weighted lognormal annual expectation", "use": "Principal annual radiance layer"},
    {"band": "weighted_mean_log_radiance", "unit": "dimensionless", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Weighted mean of ln(R/R0)", "use": "Log-space location"},
    {"band": "log_variance", "unit": "dimensionless squared", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Effective-sample-corrected weighted log variance", "use": "Within-year log dispersion"},
    {"band": "log_standard_deviation", "unit": "dimensionless", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Square root of log variance", "use": "Log-space standard deviation"},
    {"band": "lognormal_cv", "unit": "dimensionless", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "sqrt(exp(log_variance) - 1)", "use": "Model-based relative dispersion"},
    {"band": "median_radiance", "unit": "nW cm^-2 sr^-1", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Median of retained daily radiances", "use": "Robust empirical location"},
    {"band": "p10_radiance", "unit": "nW cm^-2 sr^-1", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "10th percentile of retained daily radiances", "use": "Lower empirical distribution summary"},
    {"band": "p90_radiance", "unit": "nW cm^-2 sr^-1", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "90th percentile of retained daily radiances", "use": "Upper empirical distribution summary"},
    {"band": "n_eff", "unit": "effective observations", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Kish effective sample size after weighting", "use": "Weighted observation support"},
    {"band": "n_candidate", "unit": "observations", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Count after VIIRS QA and before MAD screening", "use": "Candidate support"},
    {"band": "n_inlier", "unit": "observations", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Count after MAD screening", "use": "Retained support"},
    {"band": "candidate_fraction", "unit": "proportion", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Candidate observations divided by calendar days", "use": "Pre-MAD temporal coverage"},
    {"band": "valid_fraction", "unit": "proportion", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Inlier observations divided by calendar days", "use": "Retained temporal coverage"},
    {"band": "outlier_fraction", "unit": "proportion", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "MAD-rejected candidates divided by candidate count", "use": "Robust-screening diagnostic"},
    {"band": "mean_reliability", "unit": "dimensionless", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Mean candidate reliability score", "use": "Average accepted-observation quality"},
    {"band": "mean_retrieval_age_days", "unit": "days", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Mean VNP46A2 retrieval-age value for candidates", "use": "Source provenance diagnostic"},
    {"band": "ci95_low", "unit": "nW cm^-2 sr^-1", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Lower conditional 95% bound", "use": "Diagnostic model interval"},
    {"band": "ci95_high", "unit": "nW cm^-2 sr^-1", "data_type": "float32", "nodata_rule": "masked in EE asset; -9999 in GeoTIFF", "description": "Upper conditional 95% bound", "use": "Diagnostic model interval"},
)

assert tuple(entry["band"] for entry in BAND_DICTIONARY) == OUTPUT_BANDS

R0 = 1.0  # nW cm^-2 sr^-1; makes ln(R/R0) dimensionless
NODATA = -9999.0
EPS = 1.0e-6


@dataclass(frozen=True)
class MethodParameters:
    """Collection-1 scientific parameters frozen into metadata/manifest."""

    # Hard candidate screening.
    min_radiance: float = 0.10
    lunar_max: float = 30.0
    min_cloud_mask_quality: int = 1
    max_cloud_confidence: int = 1
    mandatory_quality: int = 0

    # Candidate reliability q_i.
    probably_clear_weight: float = 0.70
    low_cloud_mask_quality_weight: float = 0.80

    # Robust log-space screening.
    mad_multiplier: float = 3.5
    mad_consistency_factor: float = 1.4826

    # Transparent QA-derived nominal log-space uncertainty used in w_i.
    uncertainty_floor_log: float = 0.03
    uncertainty_range_log: float = 0.30

    # Principal expected-radiance and CI release threshold.
    min_effective_observations: float = 2.0


# =============================================================================
# Reproducibility fingerprint and credential-free local tests
# =============================================================================


def source_code_revision() -> str:
    """Return a stable SHA-256 fingerprint of this exact source file."""
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def run_local_self_tests() -> dict:
    """
    Run ATBD numerical-reference and configuration tests without EE credentials.

    This is intentionally pure Python. It verifies the frozen 18-band schema,
    Collection-1 constants and an independently hard-coded numerical reference
    for the weighting / Kish n_eff / corrected log variance / lognormal mean / CI.
    """
    results = {}

    # Configuration reference tests.
    p = MethodParameters()
    results["configuration"] = {
        "output_band_count": len(OUTPUT_BANDS) == 18,
        "band_dictionary_order": tuple(d["band"] for d in BAND_DICTIONARY) == OUTPUT_BANDS,
        "direct_science_band": RADIANCE_BAND == "DNB_BRDF_Corrected_NTL",
        "mandatory_quality_zero": p.mandatory_quality == 0,
        "min_radiance_0_10": math.isclose(p.min_radiance, 0.10),
        "lunar_max_30": math.isclose(p.lunar_max, 30.0),
        "mad_3_5": math.isclose(p.mad_multiplier, 3.5),
        "mad_scale_1_4826": math.isclose(p.mad_consistency_factor, 1.4826),
        "min_neff_2": math.isclose(p.min_effective_observations, 2.0),
    }

    # Numerical reference: R=[1,2,4], q=[1.0,0.7,0.8], no MAD rejection.
    radiance = [1.0, 2.0, 4.0]
    reliability = [1.0, 0.7, 0.8]
    u = [p.uncertainty_floor_log + p.uncertainty_range_log * (1.0 - q) for q in reliability]
    w = [q / (1.0 + ui * ui) for q, ui in zip(reliability, u)]
    x = [math.log(r / R0) for r in radiance]
    sum_w = sum(w)
    sum_w2 = sum(v * v for v in w)
    n_eff = sum_w * sum_w / sum_w2
    mu = sum(wi * xi for wi, xi in zip(w, x)) / sum_w
    vb = sum(wi * (xi - mu) ** 2 for wi, xi in zip(w, x)) / sum_w
    s2 = vb * n_eff / (n_eff - 1.0)
    expected = R0 * math.exp(mu + 0.5 * s2)
    se = math.sqrt(s2 / n_eff)
    ci_low = R0 * math.exp(mu - 1.96 * se + 0.5 * s2)
    ci_high = R0 * math.exp(mu + 1.96 * se + 0.5 * s2)

    expected_reference = {
        "n_eff": 2.9294438339855664,
        "mu": 0.6357662626400151,
        "log_variance": 0.5217144665827740,
        "expected_radiance": 2.4513122561746314,
        "ci95_low": 1.0719525327414048,
        "ci95_high": 5.605595018190552,
    }
    actual = {
        "n_eff": n_eff,
        "mu": mu,
        "log_variance": s2,
        "expected_radiance": expected,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
    }
    numerical_pass = all(
        math.isclose(actual[k], expected_reference[k], rel_tol=1e-12, abs_tol=1e-12)
        for k in expected_reference
    )
    results["numerical_reference"] = {
        "passed": numerical_pass,
        "actual": actual,
        "reference": expected_reference,
    }

    config_pass = all(results["configuration"].values())
    results["passed"] = config_pass and numerical_pass
    if not results["passed"]:
        raise RuntimeError(f"Local Collection-1 self-tests failed: {results}")
    return results


def verify_source_schema() -> dict:
    """Verify required VNP46A2 bands against the live Earth Engine catalogue."""
    first = ee.Image(ee.ImageCollection(DATASET).first())
    available = list(first.bandNames().getInfo())
    missing = [band for band in SOURCE_BANDS if band not in available]
    if missing:
        raise RuntimeError(f"VNP46A2 source schema mismatch; missing bands: {missing}")
    return {
        "passed": True,
        "required_bands": list(SOURCE_BANDS),
        "missing_bands": [],
    }


# =============================================================================
# Earth Engine initialisation and ROI
# =============================================================================


def initialise_earth_engine(project: str, auth_mode: Optional[str] = None) -> None:
    """Initialise Earth Engine; authenticate interactively only when required."""
    if ee is None:
        raise RuntimeError(
            "earthengine-api is required for production/dry-run execution. "
            "Install it with: pip install earthengine-api"
        )
    if not project or project == "YOUR_GEE_PROJECT_ID":
        raise ValueError(
            "Provide a real Earth Engine Cloud project with --project or the "
            "EARTHENGINE_PROJECT environment variable."
        )

    try:
        ee.Initialize(project=project)
        ee.Number(1).getInfo()
        print(f"[AUTH] Earth Engine initialised with project: {project}")
        return
    except Exception as first_error:
        print(f"[AUTH] Existing Earth Engine credentials not usable: {first_error}")

    if auth_mode:
        ee.Authenticate(auth_mode=auth_mode)
    else:
        ee.Authenticate()

    ee.Initialize(project=project)
    ee.Number(1).getInfo()
    print(f"[AUTH] Earth Engine authenticated and initialised: {project}")



def load_roi_asset(asset_id: str) -> ee.Geometry:
    """Load the production polygon/FeatureCollection and return its geometry."""
    if not asset_id or "YOUR_PROJECT" in asset_id or "BRAZIL_BOUNDARY" in asset_id:
        raise ValueError(
            "Provide the production boundary with --roi-asset, e.g. "
            "projects/<project>/assets/<boundary>."
        )
    return ee.FeatureCollection(asset_id).geometry()


# =============================================================================
# Year handling
# =============================================================================


def parse_year_spec(spec: str) -> List[int]:
    """Parse a year expression such as '2013-2025,2027'."""
    if not spec or not spec.strip():
        raise ValueError("Year specification is empty.")

    years = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2:
                raise ValueError(f"Invalid year range: {token!r}")
            start, end = int(parts[0]), int(parts[1])
            if end < start:
                raise ValueError(f"Reversed year range is not allowed: {token!r}")
            years.update(range(start, end + 1))
        else:
            years.add(int(token))

    if not years:
        raise ValueError("No years were parsed.")
    return sorted(years)



def calendar_days(year: int) -> int:
    """Return the fixed calendar denominator used by ATBD support fractions."""
    return 366 if calendar.isleap(year) else 365



def dataset_date_bounds() -> Tuple[dt.date, dt.date]:
    """Read first/latest acquisition dates from the live Earth Engine catalogue."""
    collection = ee.ImageCollection(DATASET)
    min_ms = collection.aggregate_min("system:time_start").getInfo()
    max_ms = collection.aggregate_max("system:time_start").getInfo()
    if min_ms is None or max_ms is None:
        raise RuntimeError(f"Could not determine date bounds for {DATASET}.")

    first = dt.datetime.fromtimestamp(min_ms / 1000.0, tz=dt.timezone.utc).date()
    latest = dt.datetime.fromtimestamp(max_ms / 1000.0, tz=dt.timezone.utc).date()
    return first, latest



def source_image_count(year: int) -> int:
    """Count VNP46A2 source images in the annual temporal window for provenance."""
    start = ee.Date.fromYMD(year, 1, 1)
    end = ee.Date.fromYMD(year + 1, 1, 1)
    return int(ee.ImageCollection(DATASET).filterDate(start, end).size().getInfo())



def resolve_processing_years(
    requested_years: Optional[Sequence[int]],
    start_year: Optional[int],
    end_year: Optional[int],
    allow_incomplete: bool,
) -> Tuple[List[int], dict]:
    """
    Resolve production years against the live VNP46A2 catalogue.

    A publication year must be temporally closed in the source catalogue.
    The number of source images is recorded but is not forced to equal 365/366,
    because the ATBD defines 365/366 as the support-fraction denominator rather
    than claiming that every calendar day necessarily contributes a valid image.
    """
    first_date, latest_date = dataset_date_bounds()

    first_full_year = (
        first_date.year
        if first_date == dt.date(first_date.year, 1, 1)
        else first_date.year + 1
    )
    latest_closed_year = (
        latest_date.year
        if latest_date >= dt.date(latest_date.year, 12, 31)
        else latest_date.year - 1
    )

    if requested_years is not None:
        candidates = sorted(set(int(y) for y in requested_years))
    else:
        start = int(start_year) if start_year is not None else first_full_year
        end = int(end_year) if end_year is not None else latest_closed_year
        if end < start:
            raise ValueError(f"End year {end} is earlier than start year {start}.")
        candidates = list(range(start, end + 1))

    accepted: List[int] = []
    checks = []

    for year in candidates:
        n_source = source_image_count(year)
        closed = year <= latest_closed_year
        has_data = n_source > 0
        acceptable = has_data and (closed or allow_incomplete)

        checks.append(
            {
                "year": year,
                "source_image_count": n_source,
                "calendar_days": calendar_days(year),
                "calendar_year_closed": closed,
                "accepted": acceptable,
            }
        )

        if acceptable:
            accepted.append(year)
        else:
            print(
                f"[YEARS] Rejecting {year}: source_images={n_source}, "
                f"calendar_year_closed={closed}."
            )

    if not accepted:
        raise RuntimeError(
            "No acceptable processing years remain. Use --allow-incomplete-year "
            "only for explicitly labelled internal QC products."
        )

    return accepted, {
        "dataset_first_date": first_date.isoformat(),
        "dataset_latest_date": latest_date.isoformat(),
        "first_full_year": first_full_year,
        "latest_closed_year": latest_closed_year,
        "requested_candidates": candidates,
        "accepted_years": accepted,
        "allow_incomplete_year": allow_incomplete,
        "checks": checks,
    }


# =============================================================================
# VNP46A2 QA decoding and daily candidate construction
# =============================================================================


def bit_extract(image: ee.Image, from_bit: int, to_bit: int) -> ee.Image:
    """Extract an inclusive bit range from an integer QA band."""
    if from_bit < 0 or to_bit < from_bit:
        raise ValueError("Invalid bit range.")
    width = to_bit - from_bit + 1
    mask = (1 << width) - 1
    return image.rightShift(from_bit).bitwiseAnd(mask)



def preprocess_daily(image: ee.Image, p: MethodParameters) -> ee.Image:
    """
    Build the ATBD candidate observation for one VNP46A2 day.

    Returned candidate bands:
      ntl            direct DNB_BRDF_Corrected_NTL radiance
      logx           ln(ntl / R0)
      q              accepted-observation reliability
      u_log          nominal QA-derived log-space uncertainty
      base_w         q / (1 + u_log^2), before MAD inlier masking
      retrieval_age  Latest_High_Quality_Retrieval diagnostic
      obs            1 for each candidate observation
    """
    ntl = image.select(RADIANCE_BAND).toFloat().rename("ntl")
    lunar = image.select("DNB_Lunar_Irradiance").toFloat().rename("lunar")
    retrieval_age = (
        image.select("Latest_High_Quality_Retrieval")
        .toFloat()
        .rename("retrieval_age")
    )
    mqf = image.select("Mandatory_Quality_Flag").toInt().rename("mqf")
    qf = image.select("QF_Cloud_Mask").toInt().rename("qf")
    snow_flag = image.select("Snow_Flag").toInt().rename("snow_flag")

    day_night = bit_extract(qf, 0, 0)
    cloud_mask_quality = bit_extract(qf, 4, 5)
    cloud_confidence = bit_extract(qf, 6, 7)
    shadow = bit_extract(qf, 8, 8)
    cirrus = bit_extract(qf, 9, 9)
    snow_ice = bit_extract(qf, 10, 10)
    aurora = bit_extract(qf, 12, 12)
    lunar_eclipse = bit_extract(qf, 13, 13)

    candidate_mask = (
        day_night.eq(0)
        .And(cloud_mask_quality.gte(p.min_cloud_mask_quality))
        .And(cloud_confidence.lte(p.max_cloud_confidence))
        .And(shadow.eq(0))
        .And(cirrus.eq(0))
        .And(snow_ice.eq(0))
        .And(aurora.eq(0))
        .And(lunar_eclipse.eq(0))
        .And(snow_flag.eq(0))
        .And(mqf.eq(p.mandatory_quality))
        .And(lunar.lt(p.lunar_max))
        .And(ntl.gt(p.min_radiance))
    )

    # q_i: only graded factors that remain after mandatory-quality hard screening.
    cloud_conf_weight = (
        ee.Image(0.0)
        .where(cloud_confidence.eq(0), 1.0)
        .where(cloud_confidence.eq(1), p.probably_clear_weight)
    )

    cloud_mask_weight = (
        ee.Image(0.0)
        .where(cloud_mask_quality.eq(1), p.low_cloud_mask_quality_weight)
        .where(cloud_mask_quality.gte(2), 1.0)
    )

    q = cloud_conf_weight.multiply(cloud_mask_weight).rename("q")

    u_log = (
        ee.Image(p.uncertainty_floor_log)
        .add(
            ee.Image(p.uncertainty_range_log)
            .multiply(ee.Image(1.0).subtract(q))
        )
        .rename("u_log")
    )

    base_w = q.divide(ee.Image(1.0).add(u_log.pow(2))).rename("base_w")

    ntl_masked = ntl.updateMask(candidate_mask)
    logx = ntl_masked.divide(R0).log().rename("logx")
    obs = ee.Image(1.0).rename("obs").updateMask(candidate_mask)

    return (
        ee.Image.cat(
            ntl_masked,
            logx,
            q.updateMask(candidate_mask),
            u_log.updateMask(candidate_mask),
            base_w.updateMask(candidate_mask),
            retrieval_age.updateMask(candidate_mask),
            obs,
        )
        .copyProperties(image, ["system:time_start"])
    )


# =============================================================================
# Annual robust, quality-weighted lognormal Collection-1 estimator
# =============================================================================


def build_annual_product(year: int, roi: ee.Geometry, p: MethodParameters) -> ee.Image:
    """Build exactly one 18-band annual Collection-1 image."""
    start = ee.Date.fromYMD(year, 1, 1)
    end = ee.Date.fromYMD(year + 1, 1, 1)
    days = float(calendar_days(year))

    daily_raw = (
        ee.ImageCollection(DATASET)
        .select(list(SOURCE_BANDS))
        .filterBounds(roi)
        .filterDate(start, end)
    )

    daily = daily_raw.map(lambda img: preprocess_daily(ee.Image(img), p))

    # -------------------------------------------------------------------------
    # Candidate support before robust screening.
    # -------------------------------------------------------------------------
    n_candidate = daily.select("ntl").count().toFloat().rename("n_candidate")
    candidate_fraction = n_candidate.divide(days).rename("candidate_fraction")
    mean_reliability = daily.select("q").mean().rename("mean_reliability")
    mean_retrieval_age = (
        daily.select("retrieval_age").mean().rename("mean_retrieval_age_days")
    )

    # -------------------------------------------------------------------------
    # Median/MAD robust log-space inlier screening.
    # -------------------------------------------------------------------------
    log_median = daily.select("logx").median().rename("_log_median")

    def absolute_deviation(img: ee.Image) -> ee.Image:
        return (
            img.select("logx")
            .subtract(log_median)
            .abs()
            .rename("absdev")
            .copyProperties(img, ["system:time_start"])
        )

    mad_log = daily.map(absolute_deviation).median().rename("_mad_log")
    robust_scale = mad_log.multiply(p.mad_consistency_factor).max(EPS)
    inlier_limit = robust_scale.multiply(p.mad_multiplier)

    def retain_inliers(img: ee.Image) -> ee.Image:
        logx = img.select("logx")
        inlier = logx.subtract(log_median).abs().lte(inlier_limit)
        return (
            ee.Image.cat(
                img.select("ntl").updateMask(inlier),
                logx.updateMask(inlier),
                img.select("q").updateMask(inlier),
                img.select("u_log").updateMask(inlier),
                img.select("base_w").updateMask(inlier).rename("w"),
                ee.Image(1.0).rename("obs").updateMask(inlier),
            )
            .copyProperties(img, ["system:time_start"])
        )

    filtered = daily.map(retain_inliers)
    n_inlier = filtered.select("ntl").count().toFloat().rename("n_inlier")
    valid_fraction = n_inlier.divide(days).rename("valid_fraction")

    outlier_fraction = (
        n_candidate.subtract(n_inlier)
        .divide(n_candidate.max(1.0))
        .where(n_candidate.eq(0), 0.0)
        .rename("outlier_fraction")
    )

    # -------------------------------------------------------------------------
    # Weighted moments.
    # -------------------------------------------------------------------------
    sum_w = filtered.select("w").sum().rename("_sum_w")

    sum_w2 = filtered.map(
        lambda img: img.select("w").pow(2).rename("w2")
    ).sum().rename("_sum_w2")

    sum_wx = filtered.map(
        lambda img: img.select("w").multiply(img.select("logx")).rename("wx")
    ).sum().rename("_sum_wx")

    sum_wx2 = filtered.map(
        lambda img: (
            img.select("w")
            .multiply(img.select("logx").pow(2))
            .rename("wx2")
        )
    ).sum().rename("_sum_wx2")

    positive_weight = sum_w.gt(0.0).And(sum_w2.gt(0.0))
    sum_w_safe = sum_w.max(EPS)
    sum_w2_safe = sum_w2.max(EPS)

    n_eff = (
        sum_w.pow(2)
        .divide(sum_w2_safe)
        .where(positive_weight.Not(), 0.0)
        .rename("n_eff")
    )

    mu_w = (
        sum_wx.divide(sum_w_safe)
        .rename("weighted_mean_log_radiance")
        .updateMask(positive_weight)
    )

    biased_var = (
        sum_wx2.divide(sum_w_safe)
        .subtract(mu_w.pow(2))
        .max(0.0)
        .rename("_biased_log_variance")
    )

    variance_defined = positive_weight.And(n_eff.gt(1.0))
    correction = n_eff.divide(n_eff.subtract(1.0).max(EPS))
    log_variance = (
        biased_var.multiply(correction)
        .rename("log_variance")
        .updateMask(variance_defined)
    )

    log_sd = log_variance.sqrt().rename("log_standard_deviation")
    lognormal_cv = (
        log_variance.exp()
        .subtract(1.0)
        .max(0.0)
        .sqrt()
        .rename("lognormal_cv")
    )

    # -------------------------------------------------------------------------
    # Empirical inlier summaries.
    # -------------------------------------------------------------------------
    median_radiance = filtered.select("ntl").median().rename("median_radiance")

    percentiles = filtered.select("ntl").reduce(ee.Reducer.percentile([10, 90]))
    p10_radiance = percentiles.select("ntl_p10").rename("p10_radiance")
    p90_radiance = percentiles.select("ntl_p90").rename("p90_radiance")

    # -------------------------------------------------------------------------
    # Principal lognormal expectation and conditional interval.
    # -------------------------------------------------------------------------
    release_ok = variance_defined.And(n_eff.gte(p.min_effective_observations))

    expected_radiance = (
        mu_w.add(log_variance.multiply(0.5))
        .exp()
        .multiply(R0)
        .rename("expected_radiance")
        .updateMask(release_ok)
    )

    se_mu = log_variance.divide(n_eff.max(EPS)).sqrt().rename("_se_mu")

    ci95_low = (
        mu_w.subtract(se_mu.multiply(1.96))
        .add(log_variance.multiply(0.5))
        .exp()
        .multiply(R0)
        .rename("ci95_low")
        .updateMask(release_ok)
    )

    ci95_high = (
        mu_w.add(se_mu.multiply(1.96))
        .add(log_variance.multiply(0.5))
        .exp()
        .multiply(R0)
        .rename("ci95_high")
        .updateMask(release_ok)
    )

    # -------------------------------------------------------------------------
    # Exact 18-band ATBD schema and order.
    # -------------------------------------------------------------------------
    product = ee.Image.cat(
        expected_radiance,
        mu_w,
        log_variance,
        log_sd,
        lognormal_cv,
        median_radiance,
        p10_radiance,
        p90_radiance,
        n_eff.toFloat(),
        n_candidate,
        n_inlier,
        candidate_fraction.toFloat(),
        valid_fraction.toFloat(),
        outlier_fraction.toFloat(),
        mean_reliability.toFloat(),
        mean_retrieval_age.toFloat(),
        ci95_low,
        ci95_high,
    ).select(list(OUTPUT_BANDS)).clip(roi).toFloat()

    parameter_json = json.dumps(asdict(p), sort_keys=True, separators=(",", ":"))

    return product.set(
        {
            "product": "MapBiomas_Brazil_Nighttime_Radiance_Collection_1",
            "collection": "1",
            "version": "1.0",
            "year": year,
            "temporal_start": f"{year}-01-01",
            "temporal_end_exclusive": f"{year + 1}-01-01",
            "source_dataset": DATASET,
            "source_doi": SOURCE_DOI,
            "source_science_band": RADIANCE_BAND,
            "source_support": "VNP46A2 direct BRDF-corrected DNB radiance",
            "radiance_unit": "nW cm-2 sr-1",
            "log_reference_radiance": R0,
            "central_estimator": "quality_weighted_robust_lognormal_expectation",
            "quality_schema": (
                "night; cloud-mask-quality>=1; cloud-confidence<=1; no shadow; "
                "no cirrus; no snow/ice; no aurora; no lunar eclipse; "
                "Snow_Flag=0; Mandatory_Quality_Flag=0; lunar<30; radiance>0.10"
            ),
            "weight_model": "w=I*q/(1+u^2); u=0.03+0.30*(1-q)",
            "output_band_count": len(OUTPUT_BANDS),
            "output_band_order": ",".join(OUTPUT_BANDS),
            "output_grid_statement": (
                "250 m distributed grid is an alignment lattice; source VIIRS "
                "information support remains nominally 500 m"
            ),
            "resampling_statement": (
                "No bilinear/bicubic resampling is requested; Earth Engine "
                "default nearest-neighbour allocation is retained at export"
            ),
            "method_parameters_json": parameter_json,
            "system:time_start": start.millis(),
        }
    )


def attach_release_metadata(
    image: ee.Image,
    *,
    year: int,
    source_count: int,
    boundary_version: str,
    crs: str,
    scale: Optional[float],
    crs_transform: Optional[Sequence[float]],
    code_revision: str,
) -> ee.Image:
    """Attach Appendix A/C release metadata without exposing private identifiers."""
    if crs_transform is not None:
        transform_text = ",".join(f"{float(v):.15g}" for v in crs_transform)
        pixel_origin = f"x={float(crs_transform[2]):.15g},y={float(crs_transform[5]):.15g}"
        grid_definition = "explicit_affine_transform"
    else:
        transform_text = ""
        pixel_origin = "Earth Engine CRS/scale grid origin"
        grid_definition = "crs_plus_scale"

    return image.set(
        {
            "software_version": SOFTWARE_VERSION,
            "source_code_revision_sha256": code_revision,
            "source_image_count": int(source_count),
            "boundary_version": boundary_version,
            "output_crs": crs,
            "output_scale_m": float(scale) if scale is not None else -1.0,
            "output_crs_transform": transform_text,
            "pixel_origin": pixel_origin,
            "grid_definition": grid_definition,
            "resampling_rule": "nearest-neighbour/default Earth Engine allocation",
            "output_data_type": "Float32",
            "geotiff_nodata": NODATA,
            "band_dictionary_json": json.dumps(BAND_DICTIONARY, separators=(",", ":")),
            "reproducibility_schema": "ATBD_Appendix_A_and_C",
            "year": int(year),
        }
    )


# =============================================================================
# Export
# =============================================================================


def parse_crs_transform(text: Optional[str]) -> Optional[List[float]]:
    """Parse xScale,xShear,xOrigin,yShear,yScale,yOrigin."""
    if text is None:
        return None
    values = [float(v.strip()) for v in text.split(",") if v.strip()]
    if len(values) != 6:
        raise ValueError(
            "--crs-transform must contain exactly six comma-separated numbers."
        )
    return values



def _grid_kwargs(scale: Optional[float], crs_transform: Optional[Sequence[float]]) -> dict:
    if scale is not None and crs_transform is not None:
        raise ValueError("Use either --scale or --crs-transform, never both.")
    if scale is None and crs_transform is None:
        raise ValueError("Define the output grid with --scale or --crs-transform.")
    if crs_transform is not None:
        return {"crsTransform": list(crs_transform)}
    return {"scale": float(scale)}


def export_drive_geotiff(
    image: ee.Image,
    year: int,
    roi: ee.Geometry,
    drive_folder: str,
    crs: str,
    scale: Optional[float],
    crs_transform: Optional[Sequence[float]],
    file_dimensions: int,
) -> ee.batch.Task:
    """Create and start the ATBD Google Drive Cloud-Optimized GeoTIFF export."""
    description = f"MapBiomas_VIIRS_Nighttime_Radiance_C1_{year}"
    export_image = image.unmask(value=NODATA, sameFootprint=False)
    kwargs = dict(
        image=export_image,
        description=description,
        folder=drive_folder,
        fileNamePrefix=description,
        region=roi,
        crs=crs,
        maxPixels=10**13,
        shardSize=256,
        fileDimensions=file_dimensions,
        skipEmptyTiles=True,
        fileFormat="GeoTIFF",
        formatOptions={"cloudOptimized": True, "noData": NODATA},
    )
    kwargs.update(_grid_kwargs(scale, crs_transform))
    task = ee.batch.Export.image.toDrive(**kwargs)
    task.start()
    return task


def export_asset(
    image: ee.Image,
    year: int,
    roi: ee.Geometry,
    asset_folder: str,
    crs: str,
    scale: Optional[float],
    crs_transform: Optional[Sequence[float]],
) -> tuple[ee.batch.Task, str]:
    """Create and start the ATBD Earth Engine Asset export."""
    if not asset_folder:
        raise ValueError("--asset-folder is required for asset or combined export.")
    description = f"MapBiomas_VIIRS_Nighttime_Radiance_C1_{year}"
    asset_id = f"{asset_folder.rstrip('/')}/{description}"
    kwargs = dict(
        image=image,
        description=description,
        assetId=asset_id,
        region=roi,
        crs=crs,
        maxPixels=10**13,
    )
    kwargs.update(_grid_kwargs(scale, crs_transform))
    task = ee.batch.Export.image.toAsset(**kwargs)
    task.start()
    return task, asset_id


def submit_exports(
    image: ee.Image,
    year: int,
    roi: ee.Geometry,
    destination: str,
    drive_folder: str,
    asset_folder: Optional[str],
    crs: str,
    scale: Optional[float],
    crs_transform: Optional[Sequence[float]],
    file_dimensions: int,
) -> list[dict]:
    """Submit exactly the ATBD modes: drive, asset, or combined."""
    submitted = []
    if destination in {"drive", "combined"}:
        task = export_drive_geotiff(
            image, year, roi, drive_folder, crs, scale, crs_transform, file_dimensions
        )
        submitted.append({"kind": "drive", "task": task, "destination_identifier": drive_folder})
    if destination in {"asset", "combined"}:
        task, asset_id = export_asset(
            image, year, roi, asset_folder or "", crs, scale, crs_transform
        )
        submitted.append({"kind": "asset", "task": task, "destination_identifier": asset_id})
    if not submitted:
        raise ValueError(f"Unsupported destination: {destination}")
    return submitted


def wait_for_task(task: ee.batch.Task, poll_seconds: int = 30) -> dict:
    """Block until an Earth Engine export task finishes or fails."""
    while True:
        status = task.status()
        state = status.get("state", "UNKNOWN")
        if state in {"COMPLETED", "FAILED", "CANCELLED"}:
            if state != "COMPLETED":
                raise RuntimeError(
                    f"Earth Engine export {task.id} ended as {state}: "
                    f"{status.get('error_message', status)}"
                )
            return status
        print(f"[EXPORT] Task {task.id}: {state}")
        time.sleep(max(5, poll_seconds))


# =============================================================================
# Manifest and CLI
# =============================================================================


def write_manifest(path: Path, payload: dict) -> None:
    """Write a non-sensitive JSON provenance manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "MapBiomas Brazil Collection-1 annual VIIRS direct-radiance "
            "robust quality-weighted lognormal production."
        )
    )

    parser.add_argument("--project", default=DEFAULT_GEE_PROJECT)
    parser.add_argument(
        "--auth-mode",
        default=None,
        choices=["localhost", "notebook", "gcloud", "colab"],
    )
    parser.add_argument("--roi-asset", default=DEFAULT_ROI_ASSET)

    years = parser.add_argument_group("year selection")
    years.add_argument("--years", default=None, help="e.g. 2013-2025")
    years.add_argument("--start-year", type=int, default=2013)
    years.add_argument("--end-year", type=int, default=2025)
    years.add_argument(
        "--allow-incomplete-year",
        action="store_true",
        help="Internal QC only; publication products should use closed years.",
    )

    method = parser.add_argument_group("scientific parameters")
    method.add_argument("--min-radiance", type=float, default=0.10)
    method.add_argument("--lunar-max", type=float, default=30.0)
    method.add_argument("--probably-clear-weight", type=float, default=0.70)
    method.add_argument("--low-cloud-mask-quality-weight", type=float, default=0.80)
    method.add_argument("--mad-multiplier", type=float, default=3.5)
    method.add_argument("--uncertainty-floor-log", type=float, default=0.03)
    method.add_argument("--uncertainty-range-log", type=float, default=0.30)
    method.add_argument("--min-neff", type=float, default=2.0)

    export = parser.add_argument_group("ATBD export")
    export.add_argument("--destination", choices=["drive", "asset", "combined"], default="drive")
    export.add_argument("--drive-folder", default="MapBiomas_VIIRS_Collection1")
    export.add_argument(
        "--asset-folder",
        default=None,
        help="Earth Engine asset folder for --destination asset/combined.",
    )
    export.add_argument(
        "--crs",
        default="EPSG:5880",
        help="Default ATBD CRS; replace with canonical release CRS if required.",
    )

    grid = export.add_mutually_exclusive_group()
    grid.add_argument(
        "--scale",
        type=float,
        default=None,
        help="Output grid spacing in metres; defaults to 250 if no transform is supplied.",
    )
    grid.add_argument(
        "--crs-transform",
        default=None,
        help="Six-value affine transform: xScale,xShear,xOrigin,yShear,yScale,yOrigin.",
    )

    export.add_argument("--file-dimensions", type=int, default=8192)
    export.add_argument("--wait", action="store_true")
    export.add_argument(
        "--dry-run",
        action="store_true",
        help="Build and verify annual images/band schema without starting exports.",
    )
    export.add_argument("--poll-seconds", type=int, default=30)
    export.add_argument(
        "--manifest",
        default="viirs_collection1_run_manifest.json",
        help="Public non-sensitive JSON run manifest (no project/asset/task IDs).",
    )
    export.add_argument(
        "--controlled-log",
        default=None,
        help="Optional controlled JSON log; may contain project, asset and task identifiers.",
    )
    parser.add_argument(
        "--boundary-version",
        default=None,
        help="Required release descriptor/version of the boundary used for clipping.",
    )
    parser.add_argument(
        "--visual-check-status",
        choices=["not_run", "passed", "failed"],
        default="not_run",
        help="Record the manual visual product check required by Appendix C.",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run credential-free local numerical/configuration tests and exit.",
    )

    return parser



def validate_parameters(args: argparse.Namespace) -> None:
    if args.file_dimensions <= 0 or args.file_dimensions % 256 != 0:
        raise ValueError("--file-dimensions must be a positive multiple of 256.")
    if args.scale is not None and args.scale <= 0:
        raise ValueError("--scale must be positive.")
    if args.min_radiance <= 0:
        raise ValueError("--min-radiance must be > 0 before the logarithm.")
    if args.lunar_max <= 0:
        raise ValueError("--lunar-max must be positive.")
    if args.mad_multiplier <= 0:
        raise ValueError("--mad-multiplier must be positive.")
    if args.min_neff <= 1.0:
        raise ValueError("--min-neff must be > 1 for variance/CI definition.")
    if args.uncertainty_floor_log < 0 or args.uncertainty_range_log < 0:
        raise ValueError("Uncertainty parameters must be non-negative.")

    for name in (
        "probably_clear_weight",
        "low_cloud_mask_quality_weight",
    ):
        value = getattr(args, name)
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"--{name.replace('_', '-')} must be in [0, 1].")

    if not args.self_test:
        if not args.boundary_version or not args.boundary_version.strip():
            raise ValueError("--boundary-version is required for a reproducible production run.")
        if args.destination in {"asset", "combined"} and not args.asset_folder:
            raise ValueError("--asset-folder is required for --destination asset/combined.")



def main() -> int:
    args = build_parser().parse_args()
    validate_parameters(args)

    local_tests = run_local_self_tests()
    if args.self_test:
        print(json.dumps(local_tests, indent=2, ensure_ascii=False))
        print("[SELF-TEST] Collection-1 local numerical/configuration tests PASSED.")
        return 0

    initialise_earth_engine(args.project, args.auth_mode)
    source_schema_check = verify_source_schema()
    roi = load_roi_asset(args.roi_asset)

    explicit_years = parse_year_spec(args.years) if args.years else None
    selected_years, year_metadata = resolve_processing_years(
        requested_years=explicit_years,
        start_year=args.start_year,
        end_year=args.end_year,
        allow_incomplete=args.allow_incomplete_year,
    )

    parameters = MethodParameters(
        min_radiance=args.min_radiance,
        lunar_max=args.lunar_max,
        probably_clear_weight=args.probably_clear_weight,
        low_cloud_mask_quality_weight=args.low_cloud_mask_quality_weight,
        mad_multiplier=args.mad_multiplier,
        uncertainty_floor_log=args.uncertainty_floor_log,
        uncertainty_range_log=args.uncertainty_range_log,
        min_effective_observations=args.min_neff,
    )

    transform = parse_crs_transform(args.crs_transform)
    export_scale = None if transform is not None else (
        args.scale if args.scale is not None else 250.0
    )
    code_revision = source_code_revision()
    source_counts = {item["year"]: item["source_image_count"] for item in year_metadata["checks"]}

    spatial_alignment_check = {
        "status": "passed_explicit_affine_transform" if transform else "scale_only_configuration",
        "crs": args.crs,
        "scale_m": export_scale,
        "crs_transform": transform,
        "pixel_origin": (
            {"x": transform[2], "y": transform[5]} if transform else "Earth Engine CRS/scale grid origin"
        ),
        "resampling": "nearest-neighbour/default Earth Engine allocation",
    }

    print(f"[YEARS] Accepted years: {selected_years}")
    print(f"[SOURCE] {DATASET} :: {RADIANCE_BAND} (DIRECT radiance)")
    print("[SOURCE] Required VNP46A2 schema verified.")
    print(
        "[METHOD] QA -> ln(R/R0) -> robust MAD -> "
        "w=I*q/(1+u^2) -> Kish n_eff -> corrected log variance -> "
        "lognormal expected radiance."
    )
    print(
        "[OUTPUT] expected_radiance is LINEAR nW cm^-2 sr^-1 after "
        "lognormal back-transformation; it is not stored as log(radiance)."
    )
    print(f"[BANDS] Exactly {len(OUTPUT_BANDS)} bands: {', '.join(OUTPUT_BANDS)}")
    print(
        "[GRID] "
        + (f"crs={args.crs}, transform={transform}" if transform else f"crs={args.crs}, scale={export_scale} m")
    )
    print("[GRID] 250 m is an output lattice, not VIIRS super-resolution.")
    print(f"[EXPORT] ATBD destination mode: {args.destination}.")

    public_records = []
    controlled_records = []

    for year in selected_years:
        print(f"[PROCESS] Building Collection-1 annual product for {year}")
        product = build_annual_product(year, roi, parameters)
        product = attach_release_metadata(
            product,
            year=year,
            source_count=source_counts[year],
            boundary_version=args.boundary_version,
            crs=args.crs,
            scale=export_scale,
            crs_transform=transform,
            code_revision=code_revision,
        )

        # Server-side output schema verification before any export.
        actual_bands = product.bandNames().getInfo()
        if actual_bands != list(OUTPUT_BANDS):
            raise RuntimeError(
                "Unexpected output schema. "
                f"Expected {list(OUTPUT_BANDS)}, got {actual_bands}."
            )

        public_record = {
            "year": year,
            "source_image_count": source_counts[year],
            "output_name_prefix": f"MapBiomas_VIIRS_Nighttime_Radiance_C1_{year}",
            "schema_check": "passed",
            "destination_mode": args.destination,
            "export_submitted": False,
            "export_status": "not_submitted" if args.dry_run else "pending_submission",
        }

        if args.dry_run:
            public_record["mode"] = "dry_run"
            public_records.append(public_record)
            print(f"[DRY-RUN] {year}: image and exact 18-band schema verified; no export submitted.")
            continue

        submissions = submit_exports(
            image=product,
            year=year,
            roi=roi,
            destination=args.destination,
            drive_folder=args.drive_folder,
            asset_folder=args.asset_folder,
            crs=args.crs,
            scale=export_scale,
            crs_transform=transform,
            file_dimensions=args.file_dimensions,
        )

        public_record["export_submitted"] = True
        public_record["export_status"] = "submitted"
        public_record["submitted_outputs"] = [s["kind"] for s in submissions]

        controlled_record = {
            "year": year,
            "project": args.project,
            "roi_asset": args.roi_asset,
            "boundary_version": args.boundary_version,
            "destination_mode": args.destination,
            "submissions": [],
        }

        for submission in submissions:
            task = submission["task"]
            state = task.status().get("state", "UNKNOWN")
            controlled_record["submissions"].append(
                {
                    "kind": submission["kind"],
                    "destination_identifier": submission["destination_identifier"],
                    "task_id": task.id,
                    "state_at_submission": state,
                }
            )
            print(f"[EXPORT] {year} {submission['kind']}: task_id={task.id}, state={state}")

        if args.wait:
            final_states = []
            for submission, controlled_submission in zip(submissions, controlled_record["submissions"]):
                final_status = wait_for_task(submission["task"], args.poll_seconds)
                final_state = final_status.get("state", "UNKNOWN")
                controlled_submission["final_state"] = final_state
                final_states.append(final_state)
            public_record["export_status"] = (
                "completed" if all(s == "COMPLETED" for s in final_states) else "failed"
            )

        public_records.append(public_record)
        controlled_records.append(controlled_record)

    verification = {
        "local_configuration_test": "passed",
        "local_numerical_reference_test": "passed",
        "source_schema_check": "passed",
        "annual_output_schema_checks": "passed",
        "spatial_alignment_check": spatial_alignment_check,
        "visual_product_check": args.visual_check_status,
    }

    public_manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "software_version": SOFTWARE_VERSION,
        "source_code_revision_sha256": code_revision,
        "script": Path(__file__).name,
        "product": "MapBiomas Brazil Nighttime Radiance Collection 1",
        "collection": 1,
        "version": "1.0",
        "source_dataset": DATASET,
        "source_doi": SOURCE_DOI,
        "source_science_band": RADIANCE_BAND,
        "source_bands": list(SOURCE_BANDS),
        "year_resolution": year_metadata,
        "boundary_version": args.boundary_version,
        "method_parameters": asdict(parameters),
        "methodology": {
            "input_boundary": "VIIRS-only; boundary used for filtering/clipping/export only",
            "hard_mask": (
                "night; cloud-mask-quality>=1; cloud-confidence<=1; no shadow; "
                "no cirrus; no snow/ice; no aurora; no lunar eclipse; "
                "Snow_Flag=0; Mandatory_Quality_Flag=0; lunar<30; direct radiance>0.10"
            ),
            "log_transform": "x=ln(R/R0), R0=1 nW cm^-2 sr^-1",
            "robust_filter": "|x-median(x)| <= 3.5*max(1.4826*MAD,1e-6)",
            "reliability": "q=cloud_confidence_weight*cloud_mask_quality_weight",
            "nominal_log_uncertainty": "u=0.03+0.30*(1-q)",
            "effective_weight": "w=I*q/(1+u^2)",
            "effective_sample_size": "n_eff=(sum(w)^2)/sum(w^2)",
            "central_estimator": "R0*exp(mu_w+0.5*s_w^2)",
            "conditional_ci": "R0*exp(mu_w +/- 1.96*sqrt(s_w^2/n_eff) + 0.5*s_w^2)",
            "support_fraction_denominator": "365 or 366 calendar days",
        },
        "output": {
            "band_dictionary": list(BAND_DICTIONARY),
            "band_count": len(OUTPUT_BANDS),
            "principal_band": "expected_radiance",
            "principal_band_units": "nW cm^-2 sr^-1",
            "principal_band_storage_domain": "linear physical radiance",
            "internal_model_domain": "natural-log radiance ratio",
            "drive_format": "Cloud-Optimized GeoTIFF",
            "earth_engine_asset_supported": True,
            "destination_mode": args.destination,
            "nodata_geotiff": NODATA,
            "crs": args.crs,
            "scale_m": export_scale,
            "crs_transform": transform,
            "resampling": "nearest-neighbour/default Earth Engine allocation",
            "spatial_support_note": "250 m lattice does not create spatial information beyond nominal 500 m VNP46A2 support",
        },
        "verification": verification,
        "annual_exports": public_records,
        "security": {
            "credentials_embedded": False,
            "public_manifest_omits_project_id": True,
            "public_manifest_omits_asset_ids": True,
            "public_manifest_omits_task_ids": True,
            "public_manifest_omits_local_paths": True,
        },
        "scope_exclusions": [
            "gap-filled radiance in primary estimator",
            "land cover",
            "population",
            "roads/buildings",
            "OBIA",
            "trend modelling",
            "ecological attribution",
            "second-sensor fusion",
            "neural reconstruction",
            "synthetic super-resolution",
        ],
    }

    write_manifest(Path(args.manifest), public_manifest)
    print(f"[DONE] Public manifest written: {Path(args.manifest).name}")

    if args.controlled_log:
        controlled_log = {
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "software_version": SOFTWARE_VERSION,
            "source_code_revision_sha256": code_revision,
            "project": args.project,
            "roi_asset": args.roi_asset,
            "boundary_version": args.boundary_version,
            "destination_mode": args.destination,
            "drive_folder": args.drive_folder if args.destination in {"drive", "combined"} else None,
            "asset_folder": args.asset_folder if args.destination in {"asset", "combined"} else None,
            "annual_exports": controlled_records,
            "verification": verification,
        }
        write_manifest(Path(args.controlled_log), controlled_log)
        print(f"[DONE] Controlled production log written: {Path(args.controlled_log).name}")

    if args.dry_run:
        print("[DONE] Dry-run completed successfully; no export tasks were submitted.")
    else:
        print("[DONE] ATBD export task(s) submitted successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
