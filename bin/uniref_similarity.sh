python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARSGHITVFGVNVDAFDMWGQGTMVTVSS \
       > target/uniref/uniref90_medi_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       DIQMTQSPSSLSASVGDRVTITCRTSQSLSSYTHWYQQKPGKAPKLLIYAASSRGSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK \
       > target/uniref/uniref90_medi_vl.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARGGHITIFGVNIDAFDIWGQGTMVTVSS \
       > target/uniref/uniref90_medi_uca_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRTFGQGTKVEIK \
       > target/uniref/uniref90_medi_uca_vl.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       EVQLVESGGGLIQPGGSLRLSCAASGFALRMYDMHWVRQTIDKRLEWVSAVGPSGDTYYADSVKGRFAVSRENAKNSLSLQMNSLTAGDTAIYYCVRSDRGVAGLFDSWGQGILVTVSS \
       > target/uniref/uniref90_mab114_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       DIQMTQSPSSLSASVGDRITITCRASQAFDNYVAWYQQRPGKVPKLLISAASALHAGVPSRFSGSGSGTHFTLTISSLQPEDVATYYCQNYNSAPLTFGGGTKVEIK \
       > target/uniref/uniref90_mab114_vl.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYDMHWVRQATGKGLEWVSAIGTAGDTYYPGSVKGRFTISRENAKNSLYLQMNSLRAGDTAVYYCVRSDRGVAGLFDSWGQGTLVTVSS \
       > target/uniref/uniref90_mab114_uca_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       DIQMTQSPSSLSASVGDRVTITCRASQGISNYLAWYQQKPGKVPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQKYNSAPLTFGGGTKVEIK \
       > target/uniref/uniref90_mab114_uca_vl.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       QVQLVQSGAEVKKPGASVKVSCKASGYPFTSYGISWVRQAPGQGLEWMGWISTYNGNTNYAQKFQGRVTMTTDTSTTTGYMELRRLRSDDTAVYYCARDYTRGAWFGESLIGGFDNWGQGTLVTVSS \
       > target/uniref/uniref90_s309_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       EIVLTQSPGTLSLSPGERATLSCRASQTVSSTSLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQHDTSLTFGGGTKVEIK \
       > target/uniref/uniref90_s309_vl.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       QVQLVESGGGVVQPGRSLRLSCAASGFTFSNYAMYWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRTEDTAVYYCASGSDYGDYLLVYWGQGTLVTVSS \
       > target/uniref/uniref90_regn10987_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQSEDEADYYCNSLTSISTWVFGGGTKLTVL \
       > target/uniref/uniref90_regn10987_vl.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       EVQLVESGGGLVQPGGSLRLSCAASGFSVSTKYMTWVRQAPGKGLEWVSVLYSGGSDYYADSVKGRFTISRDNSKNALYLQMNSLRVEDTGVYYCARDSSEVRDHPGHPGRSVGAFDIWGQGTMVTVSS \
       > target/uniref/uniref90_c143_vh.txt &

python bin/uniref_similarity.py \
       data/uniref/uniref202003/uniref_2020_03.fasta \
       QSALTQPASVSGSPGQSITISCTGTSNDVGSYTLVSWYQQYPGKAPKLLIFEGTKRSSGISNRFSGSKSGNTASLTISGLQGEDEADYYCCSYAGASTFVFGGGTKLTVL \
       > target/uniref/uniref90_c143_vl.txt &

wait
