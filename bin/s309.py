from utils import *
from amis import (
    get_model,
    encode,
    decode,
    deep_mutational_scan,
    compare,
    evolve,
    interpolate,
    extrapolate,
    reconstruct,
    reconstruct_prose,
)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='S309 analysis')
    parser.add_argument('--namespace', type=str, default='s309',
                        help='Model namespace')
    parser.add_argument('--model-name', type=str, default='esm1b',
                        help='Type of language model (e.g., esm1b, esm-msa)')
    parser.add_argument('--lnorm', action='store_true', default=False,
                        help='Softmax length normalization')
    args = parser.parse_args()
    return args


def load_data():
    fname = 'data/s309/s309.fa'

    seqs, name_seqs = {}, {}
    for record in SeqIO.parse(fname, 'fasta'):
        seq = record.seq
        name = record.description
        assert(name not in name_seqs)
        
        seqs[seq] = {
            'name': name,
        }
        name_seqs[name] = seq

    return seqs, name_seqs

if __name__ == '__main__':
    args = parse_args()

    model = get_model(args)

    seqs, name_seqs = load_data()

    signal_peptide = str(SeqIO.read('data/constant/signal_peptide.fa', 'fasta').seq)
    heavy_seq = str(SeqIO.read('data/constant/igg_heavy.fa', 'fasta').seq)
    light_seq = str(SeqIO.read('data/constant/igg_light_kappa.fa', 'fasta').seq)

    vh = str(name_seqs['S309_VH'])
    vl = str(name_seqs['S309_VL'])
    
    lsig = len(signal_peptide)
    lhvy = len(heavy_seq)
    llgt = len(light_seq)

    new = reconstruct(vh, model, decode_kwargs={ 'exclude': 'unnatural' })
    compare(vh, new, namespace='S309 VH')
    print('')

    new = reconstruct(vl, model, decode_kwargs={ 'exclude': 'unnatural' })
    compare(vl, new, namespace='S309 VL')
