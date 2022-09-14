from amis import (
    get_model_name,
    reconstruct_multi_models,
)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Recommend substitutions to a wildtype sequence'
    )
    parser.add_argument('sequence', type=str,
                        help='Wildtype sequence')
    parser.add_argument(
        '--model-names',
        type=str,
        default=[ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ],
        nargs='+',
        help='Type of language model (e.g., esm1b, esm1v1)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=None
        help='alpha stringency parameter'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()

    models = [
        get_model_name(model_name) for model_name in args.model_names
    ]

    mutations_models = reconstruct_multi_models(args.sequence, models, alpha=args.alpha)
    for k, v in sorted(mutations_models.items(), key=lambda item: -item[1]):
        mut_str = f'{k[1]}{k[0] + 1}{k[2]}'
        print(f'{mut_str}\t{v}')
