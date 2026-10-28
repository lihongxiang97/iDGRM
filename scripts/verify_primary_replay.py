"""Replay saved primary evidence; compare all labels, not biological accuracy."""
import argparse
import json
from pathlib import Path

import pandas as pd

from idgrm import IDGRMConfig, __version__
from idgrm.classifier import classify_primary_pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    comparisons = []
    for dataset in ('pear_hapA', 'arabidopsis', 'pear_hapB'):
        root = Path(args.analysis) / dataset / 'primary'
        previous = pd.read_csv(root / 'primary_classifications.tsv', sep='\t')
        evidence = pd.read_csv(root / 'tissue_evidence.tsv', sep='\t')
        metadata = json.loads((root / 'primary_run_metadata.json').read_text(encoding='utf-8'))
        config = IDGRMConfig(**{k: v for k, v in metadata['configuration'].items()
                               if k in IDGRMConfig.__dataclass_fields__})
        current = classify_primary_pairs(evidence, previous, config)
        comparison = previous[['pair_id', 'primary_class']].merge(
            current[['pair_id', 'primary_class']], on='pair_id', validate='one_to_one', suffixes=('_previous', '_current'))
        changed = int(comparison.primary_class_previous.ne(comparison.primary_class_current).sum())
        assert len(comparison) == len(previous) and changed == 0
        current.to_csv(output / f'{dataset}_primary.tsv', sep='\t', index=False)
        comparisons.append(dict(dataset=dataset, pairs=len(comparison), changed=changed))
        print(comparisons[-1], flush=True)
    (output / 'replay_verification.json').write_text(json.dumps(dict(
        software_version=__version__, purpose='Version regression, not accuracy validation',
        comparisons=comparisons), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
