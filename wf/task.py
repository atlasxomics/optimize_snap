import glob
import itertools
import logging
import os
import snapatac2 as snap
import subprocess

from typing import List

from latch import message
from latch.resources.tasks import large_task
from latch.types import LatchDir

import wf.plotting as pl
import wf.preprocessing as pp
import wf.utils as utils


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


@large_task
def snap_task(
    runs: List[utils.Run],
    genome: utils.Genome,
    project_name: str,
    tile_size: List[int] = 5000,
    n_features: List[int] = 25000,
    resolution: List[float] = 1.0,
    varfeat_iters: List[int] = 1,
    leiden_iters: int = -1,
    min_cluster_size: int = 20,
    min_tss: float = 2.0,
    min_frags: int = 10,
) -> LatchDir:

    samples = [run.run_id for run in runs]
    conditions = list({run.condition for run in runs})

    # Get channels for specifying plot point size, use max for now...
    channels = max({utils.get_channels(run) for run in runs})

    # Set 'groups' list for differential analysis
    groups = ["cluster"]
    if len(samples) > 1:
        groups.append("sample")
    if len(conditions) > 1:
        groups.append("condition")
    logging.info(f"Comparing features amoung groups {groups}.")

    qc_metrics = ["n_fragment", "log10_frags", "tsse"]

    genome = genome.value  # Convert to str

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    if min_frags == 0:
        logging.warrning("Minimum fragments set to 0.")
        message(
            typ="warning",
            data={"title": "min_frags", "body": "Minimum fragments set to 0."}
        )

    # Build sets of parameters --
    sets = list(
        itertools.product(tile_size, n_features, resolution, varfeat_iters)
    )
    logging.info(f"Iterating through paramter sets {sets}...")

    # Create AnnData objects --------------------------------------------------
    logging.info("Creating AnnData objects...")
    adatas = pp.make_anndatas(runs, genome, min_frags=min_frags)
    adatas = pp.filter_adatas(adatas, min_tss=min_tss)

    # Iterate through parameter sets ------------------------------------------
    count = 1
    for set in sets:
        ts, vf, cr, vi = set
        logging.info(
            f"Set {count}: tile size {ts}, variable features {vf}, clustering \
            resolution {cr}, variable feature iterations {vi}"
        )
        cr_str = str(cr).replace(".", "-")
        set_dir = f"{out_dir}/set{count}_ts{ts}-vf{vf}-cr{cr_str}-vi{vi}"
        os.makedirs(set_dir, exist_ok=True)

        logging.info(f"Adding tile matrix to objects with tile size {ts}...")
        snap.pp.add_tile_matrix(adatas, bin_size=ts)

        if len(samples) > 1:
            logging.info("Combining objects...")
            adata = pp.combine_anndata(adatas, samples, filename="combined")
        else:
            adata = adatas[0]

        logging.info(
            f"Selecting features with {vf} features and {vi} clustering \
            iteration(s)"
        )
        snap.pp.select_features(adata, n_features=vf, max_iter=vi)

        logging.info(
            f"Performing dimensionality reduction with resolution {cr}..."
        )
        adata = pp.add_clusters(adata, cr, leiden_iters, min_cluster_size)

        adata = pp.add_spatial(adata)  # Add spatial coordinates to tixels

        adata.write(f"{set_dir}/combined.h5ad")

        # bedgraphs --
        for group in groups:
            snap.ex.export_coverage(
                adata, groupby=group, suffix=f"{group}.bedgraph.zst"
            )

        # Plotting --
        pl.plot_umaps(adata, groups, "umap.pdf")
        pl.plot_spatial(
            adata,
            samples,
            "cluster",
            "spatial_dim.pdf",
            pt_size=utils.pt_sizes[channels]["dim"]
        )

        pdfs = glob.glob("*.pdf")
        bgs = glob.glob("*.zst")
        subprocess.run(["mv"] + pdfs + bgs + set_dir)

        count += 1

    figures_dir = f"{out_dir}/figures"

    pl.plot_spatial_qc(
        adata,
        samples,
        qc_metrics,
        f"{figures_dir}/spatial_qc.pdf",
        pt_size=utils.pt_sizes[channels]["qc"]
    )
    snap.pl.frag_size_distr(
        adata, interactive=False, out_file=f"{figures_dir}/fragment_size.pdf"
    )
    snap.pl.tsse(
        adata, interactive=False, out_file=f"{figures_dir}/tss_frags.pdf"
    )

    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}")
