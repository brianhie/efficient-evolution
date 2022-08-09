from anndata import AnnData
import evolocity
import pandas as pd
from scipy.sparse import coo_matrix

def load_nearest(gene):
    imgt_vdj_fname = 'data/imgt/vdj_information.tsv'
    seqs = []
    with open(imgt_vdj_fname) as f:
        for line in f:
            fields = line.rstrip().split()
            if fields[1] == gene:
                seqs.append(fields[3])
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

    seqs_vdj = {
        'medi_vh': [ 'IGHV6-1', 'IGHD3-3', 'IGHJ3' ],
        'medi_uca_vh': [ 'IGHV6-1', 'IGHD3-3', 'IGHJ3' ],
        'mab114_vh': [ 'IGHV3-13', 'IGHD6-19', 'IGHJ4' ],
        'mab114_uca_vh': [ 'IGHV3-13', 'IGHD6-19', 'IGHJ4' ],
        's309_vh': [ 'IGHV1-18', '', 'IGHJ4' ],
        'regn10987_vh': [ 'IGHV3-30', '', 'IGHJ4' ],
        'c143_vh': [ 'IGHV3-66', '', 'IGHJ3' ],

        'medi_vl': [ 'IGKV1-39', 'IGKJ1' ],
        'medi_uca_vl': [ 'IGKV1-39', 'IGKJ1' ],
        'mab114_vl': [ 'IGKV1-27', 'IGKJ4' ],
        'mab114_uca_vl': [ 'IGKV1-27', 'IGKJ4' ],
        's309_vl': [ 'IGKV3-20', 'IGKJ4' ],
        'regn10987_vl': [ 'IGLV2-14', 'IGLJ3', ],
        'c143_vl': [ 'IGLV2-23', 'IGLJ3' ],
    }

    for name in seqs_abs:
        print(f'Considering {name}...')
        
        wt_seq = seqs_abs[name]

        data = []
        for allele in seqs_vdj[name]:
            if not allele:
                continue
            
            nearest_seqs = load_nearest(allele)
            print(f'Allele: {allele}, depth {len(nearest_seqs)}')
            
            adata = AnnData(
                X=coo_matrix((len(nearest_seqs) + 1, 1)), # Empty placeholder.
                obs={ 'seq': [ wt_seq ] + nearest_seqs },
            )
            
            evolocity.tl.onehot_msa(adata, reference=0)
            
            X_onehot = adata.obsm['X_onehot']
            vocabulary = adata.uns['onehot_vocabulary']
            tok_to_idx = { vocabulary[idx]: idx for idx in vocabulary }
            X_counts = X_onehot.reshape((
                len(nearest_seqs) + 1,
                len(wt_seq),
                len(vocabulary),
            )).sum(0)

            for i in range(X_counts.shape[0]):
                wt = wt_seq[i]
                n_aligned = sum([
                    counts for idx, counts in enumerate(X_counts[i])
                    if vocabulary[idx] != '-'
                ])
                if n_aligned == 0:
                    wt_frac = -1
                else:
                    wt_frac = X_counts[i, tok_to_idx[wt]] / n_aligned
                for j in range(X_counts.shape[1]):
                    mt = vocabulary[j]
                    counts = X_counts[i, j]
                    frac = counts / len(nearest_seqs)
                    ratio = frac / wt_frac
                    data.append([ i + 1, wt, mt, allele, counts, frac, ratio ])
                
        df = pd.DataFrame(data, columns=[
            'pos',
            'wt',
            'mt',
            'allele',
            'counts',
            'fraction',
            'likelihood_ratio',
        ])

        df.to_csv(f'target/imgt/imgt_{name}_msa.txt', sep='\t')
