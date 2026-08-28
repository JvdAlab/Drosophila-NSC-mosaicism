# Script to compute per-regulon differential-activity metrics for the pySCENIC pipeline
# Post-DVC stage: consumes the DVC-stage outputs, not DVC-tracked.

# %% Import libraries
import argparse
import ast
import logging
import re
import warnings
import yaml
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from pathlib import Path

warnings.filterwarnings("ignore")


# %% Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


# %% Analysis constants
CONDITION_COL = "condition"
REF = "Control"
ALT = "ND75-KD"
SAFE_CONDITION_NAMES = {REF: "Control", ALT: "ND75_KD"}
LDH_GENE = "Ldh"


# %% Set up config
def load_config(config_path):
    """Load pipeline configuration from a YAML file."""
    logger.info(f"Loading configuration from: {config_path}")
    if not Path(config_path).exists():
        raise FileNotFoundError(f"Config file not found at {config_path}")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# %% I/O helpers
def extract_tf_name(regulon_col: str) -> str:
    """Strip SCENIC's trailing (+)/(-) suffix, preserving gene-name parentheses."""
    return re.sub(r"\([+-]\)$", "", regulon_col).strip()


def load_ctx(path: Path) -> pd.DataFrame:
    """Load the cisTarget ctx-pruning CSV with its 3-row compound header."""
    df = pd.read_csv(path, header=[0, 1], index_col=[0, 1])
    df.columns = [
        b if a.startswith("Unnamed") else f"{a}.{b}" if b else a for a, b in df.columns
    ]
    df = df.rename(
        columns={
            "Enrichment.AUC": "AUC",
            "Enrichment.NES": "NES",
            "Enrichment.MotifSimilarityQvalue": "MotifSimilarityQvalue",
            "Enrichment.OrthologousIdentity": "OrthologousIdentity",
            "Enrichment.Annotation": "Annotation",
            "Enrichment.Context": "Context",
            "Enrichment.TargetGenes": "TargetGenes",
            "Enrichment.RankAtMax": "RankAtMax",
        }
    )
    df.index.names = ["TF", "MotifID"]
    return df.reset_index()


def parse_target_genes(raw) -> set:
    """Parse a ctx `TargetGenes` cell into a set of gene names."""
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return set()
    if isinstance(raw, str):
        try:
            parsed = ast.literal_eval(raw)
            if isinstance(parsed, (list, tuple)):
                return {x[0] if isinstance(x, (list, tuple)) else x for x in parsed}
            if isinstance(parsed, (set, frozenset)):
                return set(parsed)
        except Exception:
            pass
        return set(re.findall(r"[A-Za-z][A-Za-z0-9_\-\.]+", raw))
    return set()


# %% Statistics helpers
def rank_biserial_r(u_stat: float, n1: int, n2: int) -> float:
    """r > 0 -> more active in ND75-KD (alt); r < 0 -> more active in Control (ref)."""
    return (2.0 * u_stat) / (n1 * n2) - 1.0


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Pooled-SD Cohen's d (a=ND75-KD, b=Control). Positive -> higher in ND75-KD."""
    na, nb = len(a), len(b)
    pooled = np.sqrt(
        ((na - 1) * a.std(ddof=1) ** 2 + (nb - 1) * b.std(ddof=1) ** 2) / (na + nb - 2)
    )
    return (a.mean() - b.mean()) / pooled if pooled > 0 else 0.0


# %% Regulon metrics
def compute_regulon_metrics(auc_mtx, ctx_df, adata, rescue_map: dict) -> pd.DataFrame:
    """Per-regulon DRA metrics (Mann-Whitney U, rank-biserial r, Cohen's d, BH-FDR),
    AUC summaries, NES/motif counts, and TF-expression-vs-AUC Spearman (rescue-mapped)."""
    cells_ref = adata.obs_names[adata.obs[CONDITION_COL] == REF].tolist()
    cells_alt = adata.obs_names[adata.obs[CONDITION_COL] == ALT].tolist()
    common = [c for c in auc_mtx.index if c in adata.obs_names]
    auc_sub = auc_mtx.loc[common]
    cells_ref = [c for c in cells_ref if c in common]
    cells_alt = [c for c in cells_alt if c in common]

    records = []
    for reg_col in auc_sub.columns:
        auc_vals = auc_sub[reg_col].values
        auc_ref = auc_sub.loc[cells_ref, reg_col].values
        auc_alt = auc_sub.loc[cells_alt, reg_col].values

        u, p_mw = mannwhitneyu(auc_alt, auc_ref, alternative="two-sided")
        rb_r = rank_biserial_r(u, len(auc_alt), len(auc_ref))
        d = cohens_d(auc_alt, auc_ref)

        glob_median = np.median(auc_vals)
        frac_active_ref = float(np.mean(auc_ref > glob_median))
        frac_active_alt = float(np.mean(auc_alt > glob_median))

        tf_name = extract_tf_name(reg_col)
        if tf_name in adata.var_names:
            lookup_name = tf_name
        else:
            rescued = rescue_map.get(tf_name)
            lookup_name = rescued if rescued in adata.var_names else None

        if lookup_name is not None:
            expr = adata[common, lookup_name].X
            expr = np.asarray(expr.todense()).flatten() if hasattr(expr, "todense") else np.asarray(expr).flatten()
            rho, p_rho = spearmanr(expr, auc_sub[reg_col].values)
        else:
            rho, p_rho = np.nan, np.nan

        ctx_row = ctx_df[ctx_df["TF"] == tf_name]
        nes = ctx_row["NES"].max() if len(ctx_row) else np.nan
        n_motifs = len(ctx_row)

        if "TargetGenes" in ctx_df.columns and len(ctx_row):
            all_targets = set()
            for raw in ctx_row["TargetGenes"]:
                all_targets |= parse_target_genes(raw)
            n_targets = len(all_targets)
        else:
            n_targets = np.nan

        records.append(
            dict(
                regulon=reg_col,
                TF=tf_name,
                mw_pval=p_mw,
                rank_biserial_r=rb_r,
                cohens_d=d,
                mean_auc=float(np.mean(auc_vals)),
                median_auc=float(np.median(auc_vals)),
                frac_active_ref=frac_active_ref,
                frac_active_alt=frac_active_alt,
                delta_frac_active=frac_active_alt - frac_active_ref,
                NES=nes,
                n_motifs=n_motifs,
                tf_auc_spearman=rho,
                tf_auc_spearman_p=p_rho,
                n_targets=n_targets,
            )
        )

    df = pd.DataFrame(records)
    mask = df["mw_pval"].notna()
    _, fdr, _, _ = multipletests(df.loc[mask, "mw_pval"], method="fdr_bh")
    df.loc[mask, "mw_fdr"] = fdr
    return df


def correlate_with_gene(auc_mtx, adata, gene, tf_list) -> pd.DataFrame:
    """Spearman(regulon AUC, `gene` expression) per regulon: overall + per condition + BH-FDR."""
    common = [c for c in auc_mtx.index if c in adata.obs_names]
    expr = adata[common, gene].X
    expr = np.asarray(expr.todense()).flatten() if hasattr(expr, "todense") else np.asarray(expr).flatten()

    tf_to_col = {extract_tf_name(c): c for c in auc_mtx.columns}

    records = []
    for tf in tf_list:
        col = tf_to_col.get(tf)
        if col is None:
            continue

        auc_vals = auc_mtx.loc[common, col].values
        rho, pval = spearmanr(auc_vals, expr)
        row = {"TF": tf, "regulon": col, "rho_overall": rho, "pval_overall": pval}

        obs = adata[common].obs[CONDITION_COL].values
        for cond in [REF, ALT]:
            mask = obs == cond
            r, p = (
                spearmanr(auc_vals[mask], expr[mask]) if mask.sum() > 10 else (np.nan, np.nan)
            )
            safe_cond = SAFE_CONDITION_NAMES[cond]
            row[f"rho_{safe_cond}"] = r
            row[f"pval_{safe_cond}"] = p

        records.append(row)

    df = pd.DataFrame(records)
    mask = df["pval_overall"].notna()
    _, fdr, _, _ = multipletests(df.loc[mask, "pval_overall"], method="fdr_bh")
    df.loc[mask, "fdr_overall"] = fdr
    df["delta_rho"] = df[f"rho_{SAFE_CONDITION_NAMES[ALT]}"] - df[f"rho_{SAFE_CONDITION_NAMES[REF]}"]
    return df.sort_values("rho_overall", key=abs, ascending=False).reset_index(drop=True)


# %% Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pyscenic_config.yaml")
    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)

    p = cfg["paths"]
    run_cfg = cfg["run"]

    # Resolve paths
    OUT_DIR = Path(run_cfg["output_dir"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = OUT_DIR / "04_compute_regulon_metrics.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info(f"Logging to {log_file}")

    FILES = {
        "H5AD": Path(p["work_dir"]) / p["h5ad"],
        "RESCUE_MAP": Path(p["res_dir"]) / p["gene_map"],
        "CTX_RESULTS": OUT_DIR / "ctx_pruning_results.csv",
        "AUC_RESULTS": OUT_DIR / "auc_matrix.csv",
        "METRICS_OUT": OUT_DIR / "regulon_metrics.csv",
        "LDH_CORRELATION_OUT": OUT_DIR / "regulon_Ldh_correlation.csv",
        "ANNOTATED_OUT": OUT_DIR / "top_regulons_annotated.csv",
    }

    logger.info("Loading AnnData, AUC matrix, ctx-pruning results, and rescue map...")
    adata = ad.read_h5ad(FILES["H5AD"])

    observed_conditions = set(adata.obs[CONDITION_COL].dropna().unique())
    if observed_conditions != {REF, ALT}:
        raise ValueError(
            f"Expected condition values {sorted({REF, ALT})}, "
            f"found {sorted(observed_conditions)}. Update the AnnData metadata first."
        )

    auc_mtx = pd.read_csv(FILES["AUC_RESULTS"], index_col=0)
    ctx_df = load_ctx(FILES["CTX_RESULTS"])
    rescue_df = pd.read_csv(FILES["RESCUE_MAP"])
    rescue_map = dict(zip(rescue_df["SCENIC_DB_Symbol"], rescue_df["Original_Symbol"]))
    logger.info(
        f"{adata.n_obs} cells x {adata.n_vars} genes; "
        f"{adata.obs[CONDITION_COL].value_counts().to_dict()}"
    )

    if LDH_GENE not in adata.var_names:
        raise ValueError(f"'{LDH_GENE}' not found in adata.var_names.")

    logger.info("Computing per-regulon DRA metrics...")
    metrics = compute_regulon_metrics(auc_mtx, ctx_df, adata, rescue_map)
    metrics.to_csv(FILES["METRICS_OUT"], index=False)
    logger.info(f"Wrote {FILES['METRICS_OUT']} ({len(metrics)} regulons)")

    logger.info(f"Computing Spearman correlation with {LDH_GENE}...")
    ldh_cor = correlate_with_gene(auc_mtx, adata, LDH_GENE, metrics["TF"].tolist())
    ldh_cor.to_csv(FILES["LDH_CORRELATION_OUT"], index=False)
    logger.info(f"Wrote {FILES['LDH_CORRELATION_OUT']} ({len(ldh_cor)} TFs)")

    logger.info("Building annotated top-regulons table...")
    top = metrics.reindex(
        metrics["rank_biserial_r"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)
    top["direction"] = np.where(top["rank_biserial_r"] > 0, "↑ ND75-KD", "↓ ND75-KD")

    ldh_rename = {
        "rho_overall": "ldh_rho_overall",
        "fdr_overall": "ldh_fdr",
        "delta_rho": "ldh_delta_rho",
    }
    for col in ldh_cor.columns:
        if col.startswith("rho_") and col != "rho_overall":
            ldh_rename[col] = f"ldh_{col}"

    top = top.merge(
        ldh_cor[list(ldh_rename.keys()) + ["TF"]].rename(columns=ldh_rename),
        on="TF",
        how="left",
    )
    top.to_csv(FILES["ANNOTATED_OUT"], index=False)
    logger.info(f"Wrote {FILES['ANNOTATED_OUT']} ({len(top)} regulons)")

    n_rescued = metrics["TF"].isin(rescue_map.keys()).sum()
    n_recovered = (
        metrics["TF"].isin(rescue_map.keys()) & metrics["tf_auc_spearman"].notna()
    ).sum()
    logger.info(
        f"Rescue map applied: {n_recovered}/{n_rescued} rescue-mappable TFs "
        f"recovered a tf_auc_spearman value."
    )
    logger.info("04_compute_regulon_metrics complete.")
