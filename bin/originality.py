from utils import *
from amis import reconstruct_multi_models

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

def load_ratios(name, database='abysis'):
    if database == 'abysis':
        fname = f'target/{database}/{database}_counts_{name}.txt'
    elif database == 'uniref':
        if name.startswith('uca_'):
            name = 'medi_' + name
        if name.startswith('mU_'):
            name = 'mab114_uca_' + name.split('_')[-1]
        if name.startswith('r7_'):
            name = 'regn10987_' + name.split('_')[-1]
        fname = f'target/{database}/{database}90_{name}_msa.txt'
    else:
        raise ValueError(f'{name} not supported')
    
    pos_max = {}
    with open(fname) as f:
        f.readline()
        for line in f:
            fields = line.rstrip().split()
            pos, mt, frac = int(fields[1]), fields[3], float(fields[-2])
            if (pos not in pos_max) or \
               (pos in pos_max and frac > pos_max[pos]):
                pos_max[pos] = frac
                    
    mutation_ratios = {}
    with open(fname) as f:
        f.readline()
        for line in f:
            fields = line.rstrip().split()
            pos, wt, mt = fields[1:4]
            pos = int(pos)
            mutation_ratios[(pos, wt, mt)] = max(
                float(fields[-2]) / pos_max[pos],
                float(fields[-1])
            )
    return mutation_ratios

if __name__ == '__main__':
    model_names = [ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ]

    for name in seqs_abs:
        seq = seqs_abs[name]
        print(f'\n{name}')
        mutations_models = reconstruct_multi_models(seq, model_names, alpha=0.5)

        abysis_likelihood_ratios = load_ratios(name, 'abysis')
        uniref_likelihood_ratios = load_ratios(name, 'uniref')
        
        for (pos, wt, mt), n_models in sorted(mutations_models.items(), key=lambda item: -item[1]):
            pos += 1
            abysis_ratio = abysis_likelihood_ratios[(pos, wt, mt)]
            uniref_ratio = uniref_likelihood_ratios[(pos, wt, mt)]
            fields = [ pos, wt, mt, n_models, abysis_ratio, uniref_ratio ]
            print('\t'.join([ str(field) for field in fields ]))
