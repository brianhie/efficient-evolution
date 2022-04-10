from anndata import AnnData
import evolocity
import pandas as pd
from scipy.sparse import coo_matrix

def load_nearest(fname, cutoff=None):
    seqs = []
    with open(fname) as f:
        for line in f:
            similarity, seqs_str = line.rstrip().split('\t')
            similarity = float(similarity)
            if cutoff is not None and similarity < cutoff:
                continue
            seqs_str = seqs_str.replace("'", "").strip('[]')
            seqs += seqs_str.split()
    return seqs

if __name__ == '__main__':
    seqs_abs = {
        'medi_vh': 'QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARSGHITVFGVNVDAFDMWGQGTMVTVSS',
        'medi_uca_vh': 'QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARGGHITIFGVNIDAFDIWGQGTMVTVSS',
        'mab114_vh': 'EVQLVESGGGLIQPGGSLRLSCAASGFALRMYDMHWVRQTIDKRLEWVSAVGPSGDTYYADSVKGRFAVSRENAKNSLSLQMNSLTAGDTAIYYCVRSDRGVAGLFDSWGQGILVTVSS',
        'mab114_uca_vh': 'EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYDMHWVRQATGKGLEWVSAIGTAGDTYYPGSVKGRFTISRENAKNSLYLQMNSLRAGDTAVYYCVRSDRGVAGLFDSWGQGTLVTVSS',
        's309_vh': 'QVQLVQSGAEVKKPGASVKVSCKASGYPFTSYGISWVRQAPGQGLEWMGWISTYNGNTNYAQKFQGRVTMTTDTSTTTGYMELRRLRSDDTAVYYCARDYTRGAWFGESLIGGFDNWGQGTLVTVSS',
        'regn10987_vh': 'QVQLVESGGGVVQPGRSLRLSCAASGFTFSNYAMYWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRTEDTAVYYCASGSDYGDYLLVYWGQGTLVTVSS',
        'c143_vh': 'EVQLVESGGGLVQPGGSLRLSCAASGFSVSTKYMTWVRQAPGKGLEWVSVLYSGGSDYYADSVKGRFTISRDNSKNALYLQMNSLRVEDTGVYYCARDSSEVRDHPGHPGRSVGAFDIWGQGTMVTVSS',
        
        'medi_vl': 'DIQMTQSPSSLSASVGDRVTITCRTSQSLSSYTHWYQQKPGKAPKLLIYAASSRGSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK',
        'medi_uca_vl': 'DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK',
        'mab114_vl': 'DIQMTQSPSSLSASVGDRITITCRASQAFDNYVAWYQQRPGKVPKLLISAASALHAGVPSRFSGSGSGTHFTLTISSLQPEDVATYYCQNYNSAPLTFGGGTKVEIK',
        'mab114_uca_vl': 'DIQMTQSPSSLSASVGDRVTITCRASQGISNYLAWYQQKPGKVPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQKYNSAPLTFGGGTKVEIK',
        's309_vl': 'EIVLTQSPGTLSLSPGERATLSCRASQTVSSTSLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQHDTSLTFGGGTKVEIK',
        'regn10987_vl': 'QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQSEDEADYYCNSLTSISTWVFGGGTKLTVL',
        'c143_vl': 'QSALTQPASVSGSPGQSITISCTGTSNDVGSYTLVSWYQQYPGKAPKLLIFEGTKRSSGISNRFSGSKSGNTASLTISGLQGEDEADYYCCSYAGASTFVFGGGTKLTVL',
    }

    for name in seqs_abs:
        wt_seq = seqs_abs[name]

        nearest_seqs = load_nearest(
            f'target/uniref/uniref90_{name}.txt',
        )

        adata = AnnData(
            X=coo_matrix((len(nearest_seqs) + 1, 1)),
            obs={ 'seq': [ wt_seq ] + nearest_seqs },
        )

        evolocity.tl.onehot_msa(adata, reference=0)

        X_onehot = adata.obsm['X_onehot']
        vocabulary = adata.uns['onehot_vocabulary']

        X_counts = X_onehot.reshape((
            len(nearest_seqs) + 1,
            len(wt_seq),
            len(vocabulary),
        )).sum(0)

        data = []
        for i in range(X_counts.shape[0]):
            wt = wt_seq[i]
            n_aligned = sum([
                counts for idx, counts in enumerate(X_counts[i])
                if vocabulary[idx] != '-'
            ])
            for j in range(X_counts.shape[1]):
                mt = vocabulary[j]
                counts = X_counts[i, j]
                frac = counts / len(nearest_seqs)
                frac_aligned = counts / n_aligned
                data.append([ i + 1, wt, mt, counts, frac, frac_aligned ])
                
        df = pd.DataFrame(data, columns=[
            'pos',
            'wt',
            'mt',
            'counts',
            'fraction',
            'fraction_aligned',
        ])

        df.to_csv(f'target/uniref/uniref90_{name}_msa.txt', sep='\t')
