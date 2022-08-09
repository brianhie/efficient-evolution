# Efficient evolution from general protein language models

Scripts for running the analysis described in the paper ["Efficient evolution of human antibodies from general protein language models and sequence information alone"](https://www.biorxiv.org/content/10.1101/2022.04.10.487811v1).

### Running the model

To evaluate the model on a new sequence, clone this repository and run
```bash
python bin/recommend.py [sequence]
```
where `[sequence]` is the wildtype protein sequence you want to evolve. The script will output a list of substitutions and the number of recommending language models.

### Paper analysis scripts

To reproduce the analysis in the paper, first download and extract data with the commands:
```bash
wget https://zenodo.org/record/6968342/files/data.tar.gz
tar xvf data.tar.gz
```

To acquire mutations to a given antibody, run the command
```bash
bash bin/eval_models.sh [antibody_name]
```
where `[antibody_name]` is one of `medi8852`, `medi_uca`, `mab114`, `mab114_uca`, `s309`, `regn10987`, or `c143`.

DMS experiments can be run with the command
```bash
bash bin/dms.sh
```
