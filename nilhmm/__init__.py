"""nilhmm: HMM-based introgression calling for NIL populations"""

from .core import call_introgressions, introgression_hmm
from .io import read_vcf, write_results
from .grid_search import optimize_parameters

__version__ = "0.1.0"
__all__ = ["call_introgressions", "introgression_hmm", "read_vcf", "write_results", "optimize_parameters"]
