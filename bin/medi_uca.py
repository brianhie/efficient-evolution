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
    parser = argparse.ArgumentParser(description='MEDI UCA analysis')
    parser.add_argument('--namespace', type=str, default='MEDI UCA',
                        help='Model namespace')
    parser.add_argument('--model-name', type=str, default='esm1b',
                        help='Type of language model (e.g., esm1b, esm-msa)')
    args = parser.parse_args()
    return args


def load_data():
    fname = 'data/medi8852/uca_amis.fa'

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

    vh = (
        'QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARGGHITIFGVNIDAFDIWGQGTMVTVSS'
    )
    vl = (
        'DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK'
    )

    new = reconstruct(vh, model, decode_kwargs={ 'exclude': 'unnatural' })
    compare(vh, new, namespace='MEDI_UCA VH')
    print('')

    new = reconstruct(vl, model, decode_kwargs={ 'exclude': 'unnatural' })
    compare(vl, new, namespace='MEDI_UCA VK')
