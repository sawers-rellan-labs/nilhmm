#!/bin/bash
cd "/Users/fvrodriguez/Library/CloudStorage/GoogleDrive-frodrig4@ncsu.edu/My Drive/repos/nilhmm"
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md to reflect current project structure

- Remove references to deleted files (bzea_vcf_introgression_caller.py, File_S14.Step_2_nNIL_introgression_calls_from_chip_data.py)
- Update to reference current files (scripts/call_bzea_introgressions.py, nilhmm/core.py)
- Reorganize document structure to match actual package layout
- Add coverage-specific parameter defaults and package usage examples
- Maintain all technical HMM implementation details"
