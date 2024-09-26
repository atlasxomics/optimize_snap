from typing import List, Optional

from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

from wf.task import opt_task
from wf.utils import Run, Genome


metadata = LatchMetadata(
    display_name="optimize_snap",
    author=LatchAuthor(
        name="James McGann",
        email="jamesm@atlasxomics.com",
        github="github.com/atlasxomics"
    ),
    repository="https://github.com/atlasxomics/optimize_snap",
    license="MIT",
    parameters={
        "runs": LatchParameter(
            display_name="runs",
            description="List of runs to be analyzed; each run must contain a \
                run_id and fragments.tsv file; optional: condition, tissue \
                position file for filtering on/off tissue. Note that multiple \
                Conditions must be separted by '_' (i.e., Female-control).",
            batch_table_column=True
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Name of output directory in snap_opt_outs/",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome for runs.",
            batch_table_column=True,
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the \
                tile matrix.",
            batch_table_column=True,
        ),
        "n_features": LatchParameter(
            display_name="number of features",
            description="Number of features to be selected as 'most \
                accessible' in tile matrix.",
            batch_table_column=True
        ),
        "resolution": LatchParameter(
            display_name="clustering resolution",
            description="Clustering resolution for Leiden algorithm; higher \
                values result in more clusters.",
            batch_table_column=True
        ),
        "varfeat_iters": LatchParameter(
            display_name="variable features iterations",
            description="Iterations performed when selecting variable \
                features for tile matrix. 'If greater than 1, this function \
                will perform iterative clustering and feature selection based \
                on variable features found using previous clustering results. \
                This is similar to the procedure implemented in ArchR... \
                Default value is 1, which means no iterative clustering is \
                 performed.'- SnapATAC2 docs",
            batch_table_column=True
        ),
        "min_cluster_size": LatchParameter(
            display_name="minimum cells per cluster",
            description="Minimum number of cells in a cluster.",
            batch_table_column=True,
            hidden=True
        ),
        "min_tss": LatchParameter(
            display_name="minimum TSS",
            description="The minimum numeric transcription start site (TSS) \
                enrichment score required for a cell to pass filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="The minimum number of mapped fragments required  \
                per cell to pass filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "pt_size": LatchParameter(
            display_name="Override point size",
            description="Point size for spatial plot of clustering. \
                Recommendations: 50x:75, 96x:5, 220:5.",
            batch_table_column=True,
            hidden=True
        ),
        "qc_pt_size": LatchParameter(
            display_name="Override spatial QC point size",
            description="Point size for spatial plot of clustering. \
                Recommendations: 50x:25, 96x:1, 220:0.5.",
            batch_table_column=True,
            hidden=True
        ),
    }
)


@workflow(metadata)
def opt_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    tile_size: List[int] = [5000],
    n_features: List[int] = [25000],
    resolution: List[float] = [1.0],
    varfeat_iters: List[int] = [1],
    min_cluster_size: int = 20,
    min_tss: float = 2.0,
    min_frags: int = 10,
    pt_size: Optional[int] = None,
    qc_pt_size: Optional[int] = None,
) -> LatchDir:
    """Determine optimal input parameters for spatial epigenetic analysis.
    * Variable Feature Iterations
        * https://github.com/kaizhang/SnapATAC2/issues/109
    * Number of Features
        * https://github.com/kaizhang/SnapATAC2/issues/111
    """

    results = opt_task(
        runs=runs,
        genome=genome,
        project_name=project_name,
        tile_size=tile_size,
        n_features=n_features,
        resolution=resolution,
        varfeat_iters=varfeat_iters,
        min_cluster_size=min_cluster_size,
        min_tss=min_tss,
        min_frags=min_frags,
        pt_size=pt_size,
        qc_pt_size=qc_pt_size
    )

    return results
