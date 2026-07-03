"""nilhmm: HMM-based introgression calling for NIL populations"""

from .core import (
    call_introgressions, introgression_hmm,
    call_introgressions_counts, introgression_hmm_counts,
)
from .io import read_vcf, read_vcf_counts, write_results
from .grid_search import optimize_parameters

__version__ = "0.1.0"
__all__ = [
    "call_introgressions", "introgression_hmm",
    "call_introgressions_counts", "introgression_hmm_counts",
    "read_vcf", "read_vcf_counts", "write_results", "optimize_parameters",
]
