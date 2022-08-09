# Efficient evolution from protein language models

Scripts for reproducing results in the paper ["Efficient evolution of human antibodies from general protein language models and sequence information alone"](https://www.biorxiv.org/content/10.1101/2022.04.10.487811v1).

To download data, run the commands:
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
