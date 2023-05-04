import numpy as np
import scipy.special
import sys
import torch

def err_model(name):
    raise ValueError('Model {} not supported'.format(name))

def get_model(args):
    if args.model_name == 'esm1b':
        from fb_model import FBModel
        model = FBModel(
            'esm1b_t33_650M_UR50S',
            repr_layer=[-1],
        )
    elif args.model_name.startswith('esm1v'):
        from fb_model import FBModel
        model = FBModel(
            'esm1v_t33_650M_UR90S_' + args.model_name[-1],
            repr_layer=[-1],
        )
    elif args.model_name == 'esm-msa':
        from fb_model import FBModel
        model = FBModel(
            'esm_msa1_t12_100M_UR50S',
            repr_layer=[-1],
        )
    elif args.model_name == 'prose':
        from prose_model import ProseModel
        model = ProseModel()
    else:
        err_model(args.model_name)

    return model

def get_model_name(name):
    if name == 'esm1b':
        from fb_model import FBModel
        model = FBModel(
            'esm1b_t33_650M_UR50S',
            repr_layer=[-1],
        )
    elif name.startswith('esm1v'):
        from fb_model import FBModel
        model = FBModel(
            'esm1v_t33_650M_UR90S_' + name[-1],
            repr_layer=[-1],
        )
    elif name == 'esm-msa':
        from fb_model import FBModel
        model = FBModel(
            'esm_msa1_t12_100M_UR50S',
            repr_layer=[-1],
        )
    elif name == 'prose':
        from prose_model import ProseModel
        model = ProseModel()
    else:
        err_model(name)

    return model

def encode(seq, model):
    return model.encode(seq)

def decode(embedding, model, exclude=set()):
    if exclude == 'unnatural':
        exclude = set([
            'B', 'J', 'O', 'U', 'X', 'Z', '-', '.',
        ])
    
    logits = model.decode(embedding)

    if 'esm_msa1_t' in model.name_:
        logits = logits[0]
        embedding = embedding[0]

    assert(logits.shape[0] == embedding.shape[0])
    assert(logits.shape[1] == len(model.alphabet_.all_toks))

    valid_idx = [
        idx for idx, tok in enumerate(model.alphabet_.all_toks)
    ]
    logits = logits[:, valid_idx]

    argmax = ''.join([
        model.alphabet_.all_toks[valid_idx[tok_idx]]
        if ('<' not in model.alphabet_.all_toks[valid_idx[tok_idx]] and
            model.alphabet_.all_toks[valid_idx[tok_idx]] not in exclude) else '.'
        for tok_idx in np.argmax(logits, 1)
    ])

    return argmax

def deep_mutational_scan(seq, model):
    if model.name_ == 'prose':
        from prose.alphabets import Uniprot21
        from prose.utils import pack_sequences
        alphabet = Uniprot21()
        x = [ torch.from_numpy(alphabet.encode(seq.encode())).long() ]
        x, _ = pack_sequences(x)
        logits = model.model_(x).data.cpu().detach().numpy()

        for i in range(logits.shape[0]):
            for j in range(logits.shape[1]):
                pos, mt = i + 1, alphabet.decode(j).decode('utf-8')
                val = logits[i, j]
                print(f'{pos}\t{mt}\t{val}')
    else:
        logits = model.decode(model.encode(seq))

        for i in range(logits.shape[0]):
            for j in range(logits.shape[1]):
                pos, mt = i + 1, model.alphabet_.all_toks[j]
                val = logits[i, j]
                print(f'{pos}\t{mt}\t{val}')

def reconstruct(seq, model, encode_kwargs={}, decode_kwargs={}):
    if model.name_ == 'prose':
        return reconstruct_prose(seq, model)
    
    return decode(
        encode(seq, model, **encode_kwargs),
        model, **decode_kwargs
    )

def soft_reconstruct(seq, model, alpha=1., offset=1):
    if model.name_ == 'prose':
        raise NotImplementedError('Does not support prose reconstruction')

    exclude = set([
        'B', 'J', 'O', 'U', 'X', 'Z', '-', '.',
    ])

    logits = model.predict_sequence_prob(seq)
    probs = scipy.special.softmax(logits, axis=1)

    mutations = []
    
    for i in range(probs.shape[0] - 1):
        if i == 0:
            continue
        pos = i - offset
        wt_j = model.alphabet_.tok_to_idx[seq[pos]]
        wt_prob = probs[i, wt_j]
        for j in range(probs.shape[1]):
            mt = model.alphabet_.all_toks[j]
            if mt in exclude or '<' in mt:
                continue
            if j == wt_j:
                continue
            mt_prob = probs[i, j]
            if mt_prob > alpha * wt_prob:
                mutations.append((pos, seq[pos], mt))

    return mutations

def reconstruct_prose(seq, model):
    from prose.alphabets import Uniprot21
    from prose.utils import pack_sequences
    
    alphabet = Uniprot21()
    x = [ torch.from_numpy(alphabet.encode(seq.encode())).long() ]
    x, _ = pack_sequences(x)
    
    logits = model.model_(x).data.cpu().detach().numpy()

    return ''.join([
        alphabet.decode(np.argmax(logits[i])).decode('utf-8')
        for i in range(logits.shape[0])
    ])

def compare(seq_old, seq_new, start=0, end=None, namespace=None):
    #print(f'Old: {seq_old}')
    #print(f'New: {seq_new}')
    if namespace is not None:
        sys.stdout.write(f'{namespace} mutations: ')
    for idx, (ch_old, ch_new) in enumerate(
            zip(seq_old, seq_new)
    ):
        if idx < start:
            continue
        if end is not None and idx >= end:
            continue
        if ch_new == '.':
            continue
        if ch_old != ch_new:
            sys.stdout.write(f'{ch_old}{idx - start + 1}{ch_new}, ')
    sys.stdout.write('\n')

def diff(seq_old, seq_new, start=0, end=None):
    different_muts = []
    for idx, (ch_old, ch_new) in enumerate(
            zip(seq_old, seq_new)
    ):
        if idx < start:
            continue
        if end is not None and idx >= end:
            continue
        if ch_new == '.':
            continue
        if ch_old != ch_new:
            different_muts.append((idx, ch_old, ch_new))
    return different_muts
    
def reconstruct_multi_models(
        wt_seq,
        model_names=[
            'esm1b',
            'esm1v1',
            'esm1v2',
            'esm1v3',
            'esm1v4',
            'esm1v5',
        ],
        alpha=None,
        return_names=False,
):
    mutations_models, mutations_model_names = {}, {}
    for model_name in model_names:
        model = get_model_name(model_name)
        if alpha is None:
            wt_new = reconstruct(
                wt_seq, model, decode_kwargs={ 'exclude': 'unnatural' }
            )
            mutations_model = diff(wt_seq, wt_new)
        else:
            mutations_model = soft_reconstruct(
                wt_seq, model, alpha=alpha,
            )
        for mutation in mutations_model:
            if mutation not in mutations_models:
                mutations_models[mutation] = 0
                mutations_model_names[mutation] = []
            mutations_models[mutation] += 1
            mutations_model_names[mutation].append(model.name_)
        del model

    if return_names:
        return mutations_models, mutations_model_names

    return mutations_models

def evolve(seq_uca, model, n_generations=1):
    print(f'Gen 0: {seq_uca}')
    seq_curr = seq_uca
    for gen in range(n_generations):
        seq_ml = reconstruct(seq_curr, model)
        compare(seq_curr, seq_ml)
        if seq_ml == seq_curr:
            print(f'Converged at generation {gen}')
            break
        print(f'Gen {gen + 1}: {seq_ml}')
        seq_curr = seq_ml

def interpolate(baseline, target, model, n_steps=15,
                start=0, end=None):
    for alpha in np.linspace(0., 1., n_steps):
        midpoint = (alpha * encode(baseline, model)) + \
                   ((1 - alpha) * encode(target, model))
        new_seq = decode(midpoint, model)
        compare(target, new_seq, start, end)

def extrapolate(baseline, target, model, n_steps=15,
                start=0, end=None):
    for alpha in np.linspace(1., 5., n_steps):
        delta = (encode(target, model) - encode(baseline, model))
        new_seq = decode((delta * alpha) + encode(target, model), model)
        compare(target, new_seq, start, end)
