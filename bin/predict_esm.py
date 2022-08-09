# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# Modified by Brian Hie (brianhie@stanford.edu).

import argparse
import pathlib

import torch

from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import itertools
from typing import List, Tuple
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns

from utils import deep_mutational_scan


def label_row(row, sequence, token_probs, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"

    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)

    # add 1 for BOS
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()


def compute_pppl(row, sequence, model, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"

    # modify the sequence
    sequence = sequence[:idx] + mt + sequence[(idx + 1) :]

    # encode the sequence
    data = [
        ("protein1", sequence),
    ]

    batch_converter = alphabet.get_batch_converter()

    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)

    # compute probabilities at each position
    log_probs = []
    for i in range(1, len(sequence) - 1):
        batch_tokens_masked = batch_tokens.clone()
        batch_tokens_masked[0, i] = alphabet.mask_idx
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens_masked.cuda())["logits"], dim=-1)
        log_probs.append(token_probs[0, i, alphabet.get_idx(sequence[i])].item())  # vocab size
    return sum(log_probs)


def predict_esm(
        sequence: str,
        model_locations: List[str],
        scoring_strategy: str = 'wt-marginals',
        mutation_col: str = 'mutant',
        offset_idx: int = 0,
        nogpu: bool = False,
        verbose: int = 0,
):
    # Conduct the deep mutational scan.
    data = [
        f'{wt}{pos + offset_idx}{mt}'
        for pos, wt, mt in deep_mutational_scan(sequence)
    ]
    df = pd.DataFrame(data, columns=[ mutation_col ])

    # inference for each model
    for model_location in model_locations:
        model, alphabet = pretrained.load_model_and_alphabet(model_location)
        model.eval()
        if torch.cuda.is_available() and not nogpu:
            model = model.cuda()
            if verbose:
                print("Transferred model to GPU")

        batch_converter = alphabet.get_batch_converter()

        data = [
            ("protein1", sequence),
        ]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)

        if scoring_strategy == "wt-marginals":
            with torch.no_grad():
                token_probs = torch.log_softmax(
                    model(batch_tokens.cuda())["logits"], dim=-1
                )
            df[model_location] = df.apply(
                lambda row: label_row(
                    row[mutation_col],
                    sequence,
                    token_probs,
                    alphabet,
                    offset_idx,
                ),
                axis=1,
            )
        elif scoring_strategy == "masked-marginals":
            all_token_probs = []
            if verbose:
                iterator = tqdm(range(batch_tokens.size(1)))
            else:
                iterator = range(batch_tokens.size(1))
            for i in iterator:
                batch_tokens_masked = batch_tokens.clone()
                batch_tokens_masked[0, i] = alphabet.mask_idx
                with torch.no_grad():
                    token_probs = torch.log_softmax(
                        model(batch_tokens_masked.cuda())["logits"], dim=-1
                    )
                all_token_probs.append(token_probs[:, i])  # vocab size
            token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
            df[model_location] = df.apply(
                lambda row: label_row(
                    row[mutation_col],
                    sequence,
                    token_probs,
                    alphabet,
                    offset_idx,
                ),
                axis=1,
            )
        elif scoring_strategy == "pseudo-ppl":
            tqdm.pandas()
            df[model_location] = df.progress_apply(
                lambda row: compute_pppl(
                    row[mutation_col], sequence, model, alphabet, offset_idx
                ),
                axis=1,
            )

    return df


if __name__ == "__main__":
    seqs_abs = {
        'medi_vh': 'QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARSGHITVFGVNVDAFDMWGQGTMVTVSS',
        'uca_vh': 'QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARGGHITIFGVNIDAFDIWGQGTMVTVSS',
        'mab114_vh': 'EVQLVESGGGLIQPGGSLRLSCAASGFALRMYDMHWVRQTIDKRLEWVSAVGPSGDTYYADSVKGRFAVSRENAKNSLSLQMNSLTAGDTAIYYCVRSDRGVAGLFDSWGQGILVTVSS',
        'mU_vh': 'EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYDMHWVRQATGKGLEWVSAIGTAGDTYYPGSVKGRFTISRENAKNSLYLQMNSLRAGDTAVYYCVRSDRGVAGLFDSWGQGTLVTVSS',
        's309_vh': 'QVQLVQSGAEVKKPGASVKVSCKASGYPFTSYGISWVRQAPGQGLEWMGWISTYNGNTNYAQKFQGRVTMTTDTSTTTGYMELRRLRSDDTAVYYCARDYTRGAWFGESLIGGFDNWGQGTLVTVSS',
        'r7_vh': 'QVQLVESGGGVVQPGRSLRLSCAASGFTFSNYAMYWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRTEDTAVYYCASGSDYGDYLLVYWGQGTLVTVSS',
        'c143_vh': 'EVQLVESGGGLVQPGGSLRLSCAASGFSVSTKYMTWVRQAPGKGLEWVSVLYSGGSDYYADSVKGRFTISRDNSKNALYLQMNSLRVEDTGVYYCARDSSEVRDHPGHPGRSVGAFDIWGQGTMVTVSS',
        
        'medi_vl': 'DIQMTQSPSSLSASVGDRVTITCRTSQSLSSYTHWYQQKPGKAPKLLIYAASSRGSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK',
        'uca_vl': 'DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK',
        'mab114_vl': 'DIQMTQSPSSLSASVGDRITITCRASQAFDNYVAWYQQRPGKVPKLLISAASALHAGVPSRFSGSGSGTHFTLTISSLQPEDVATYYCQNYNSAPLTFGGGTKVEIK',
        'mU_vl': 'DIQMTQSPSSLSASVGDRVTITCRASQGISNYLAWYQQKPGKVPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQKYNSAPLTFGGGTKVEIK',
        's309_vl': 'EIVLTQSPGTLSLSPGERATLSCRASQTVSSTSLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQHDTSLTFGGGTKVEIK',
        'r7_vl': 'QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQSEDEADYYCNSLTSISTWVFGGGTKLTVL',
        'c143_vl': 'QSALTQPASVSGSPGQSITISCTGTSNDVGSYTLVSWYQQYPGKAPKLLIFEGTKRSSGISNRFSGSKSGNTASLTISGLQGEDEADYYCCSYAGASTFVFGGGTKLTVL',
    }

    model_locations = [
        'esm1b_t33_650M_UR50S',
        'esm1v_t33_650M_UR90S_1',
        'esm1v_t33_650M_UR90S_2',
        'esm1v_t33_650M_UR90S_3',
        'esm1v_t33_650M_UR90S_4',
        'esm1v_t33_650M_UR90S_5',
    ]

    data = []
    for seq_name in seqs_abs:
        for model_location in model_locations:
            df_wt_marginals = predict_esm(
                seqs_abs[seq_name], [ model_location ],
                scoring_strategy='wt-marginals',
            )
            df_masked_marginals = predict_esm(
                seqs_abs[seq_name], [ model_location ],
                scoring_strategy='masked-marginals',
            )
            corr, _ = ss.spearmanr(
                df_wt_marginals[model_location],
                df_masked_marginals[model_location],
            )
            print(f'{seq_name}\t{model_location}\t{corr}')
            plt.figure()
            plt.scatter(
                df_wt_marginals[model_location],
                df_masked_marginals[model_location],
            )
            plt.xlabel('WT marginals')
            plt.ylabel('Masked marginals')
            plt.title(f'{seq_name}, {model_location}, Spearman r = {corr:.03}')
            plt.savefig(
                f'figures/esm_marginal_scatter/{seq_name}_{model_location}.png',
                dpi=500
            )
            plt.close()

            data.append([
                seq_name,
                model_location,
                corr
            ])

    df = pd.DataFrame(data, columns=[ 'name', 'model', 'corr' ])
    plt.figure()
    sns.barplot(data=df, x='name', y='corr', hue='model')
    plt.ylim([ .8, 1.01 ])
    plt.savefig('figures/esm_marginal_summary.svg')
    plt.close()
