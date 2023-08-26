import torch
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained
import numpy as np
import warnings

class FBModel(object):
    def __init__(self, name, repr_layer=[-1]):
        self.name_ = name
        self.repr_layer_ = repr_layer

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            model, alphabet = pretrained.load_model_and_alphabet(name)

        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
        self.model_ = model
        self.alphabet_ = alphabet
        self.unk_idx_ = alphabet.tok_to_idx['<unk>']
        self.offset_ = 1

        assert(all(
            -(model.num_layers + 1) <= i <= model.num_layers
            for i in [ -1 ]
        ))
        self.repr_layers_ = [
            (i + model.num_layers + 1) % (model.num_layers + 1)
            for i in [ -1 ]
        ]

    def predict_sequence_prob(self, seq):
        output = []

        batch_size = 1022 # ESM max sequence context.

        for batchi in range(((len(seq) - 1) // batch_size) + 1):
            start = batchi * batch_size
            end = (batchi + 1) * batch_size
    
            seqs = [ seq[start:end] ]
            labels = [ 'seq0' ]
    
            dataset = FastaBatchedDataset(labels, seqs)
            batches = dataset.get_batch_indices(batch_size, extra_toks_per_seq=1)
            data_loader = torch.utils.data.DataLoader(
                dataset, collate_fn=self.alphabet_.get_batch_converter(),
                batch_sampler=batches
            )
    
            with torch.no_grad():
                for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                    if torch.cuda.is_available():
                        toks = toks.to(device='cuda', non_blocking=True)
                    out = self.model_(
                        toks,
                        repr_layers=self.repr_layers_,
                        return_contacts=False
                    )
                    logits = out['logits'].to(device='cpu').numpy()
    
            output.append(logits[0])
    
        for i, x in enumerate(output):
            if len(output) > 1 and i > 0:
                output[i] = output[i][1:]
            if len(output) > 1 and i < len(output) - 1:
                output[i] = output[i][:-1]
    
        return np.vstack(output)

    def encode(self, seq_encode):
        seqs = [ seq_encode ]
        labels_full = [ 'seq0' ]
    
        embedded_seqs = {}
    
        max_len = max([ len(seq) for seq in seqs ])
        window_size = 1022
        n_windows = ((max_len - 1) // window_size) + 1
        for window_i in range(n_windows):
            start = window_i * window_size
            end = (window_i + 1) * window_size
    
            seqs_window = [ seq[start:end] for seq in seqs ]
    
            dataset = FastaBatchedDataset(labels_full, seqs_window)
            batches = dataset.get_batch_indices(4096, extra_toks_per_seq=1)
            data_loader = torch.utils.data.DataLoader(
                dataset, collate_fn=self.alphabet_.get_batch_converter(),
                batch_sampler=batches
            )
    
            with torch.no_grad():
                for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                    if isinstance(labels[0], list):
                        labels = labels[0]
                    if torch.cuda.is_available():
                        toks = toks.to(device='cuda', non_blocking=True)
                    out = self.model_(
                        toks,
                        repr_layers=self.repr_layers_,
                        return_contacts=False
                    )
                    representations = {
                        layer: t.to(device='cpu')
                        for layer, t in out['representations'].items()
                    }
    
                    for i, label in enumerate(labels):
                        seq_idx = int(label[3:])
                        seq = seqs[seq_idx]
                        assert(len(representations.items()) == 1)
                        for _, t in representations.items():
                            if 'esm_msa1_t' in self.name_:
                                representation = t[i, :, 1:len(strs[0][i]) + 1]
                            else:
                                representation = t[i, 1:len(strs[i]) + 1]
                        if seq not in embedded_seqs:
                            embedded_seqs[seq] = []
                        embedded_seqs[seq].append({
                            f'embedding{window_i}': representation.numpy()
                        })
    
        for seq in embedded_seqs:
            embeddings = []
            for window_i in range(n_windows):
                for meta in embedded_seqs[seq]:
                    if f'embedding{window_i}' in meta:
                        embedding_i = meta[f'embedding{window_i}']
                        if n_windows > 1 and window_i > 0:
                            embedding_i = embedding_i[1:]
                        if n_windows > 1 and window_i < n_windows - 1:
                            embedding_i = embedding_i[:-1]
                        embeddings.append(embedding_i)
                        break
            embedding = np.vstack(embeddings)
            embedded_seqs[seq] = [ { 'embedding': embedding } ]
    
        return embedded_seqs[seq_encode][0]['embedding']

    def decode(self, embedding):
        with torch.no_grad():
            embedding = torch.from_numpy(embedding)
            if torch.cuda.is_available():
                embedding = embedding.to(device='cuda', non_blocking=True)
            logits = self.model_.lm_head(embedding)
            logits = logits.to(device='cpu').numpy()
        return logits
