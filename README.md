# Efficient evolution from general protein language models

Scripts for running the analysis described in the paper ["Efficient evolution of human antibodies from general protein language models"](https://www.nature.com/articles/s41587-023-01763-2).

### Running the model

To evaluate the model on a new sequence, clone this repository and run
```bash
python bin/recommend.py [sequence]
```
where `[sequence]` is the wildtype protein sequence you want to evolve. The script will output a list of substitutions and the number of recommending language models.

To recommend mutations to antibody variable domain sequences, we have simply run the above script separately on the heavy and light chain sequences.

We have also made a [Google Colab](https://colab.research.google.com/drive/18QLOmi5yNb1i9wztAzv981Wgk2E4IP4q?usp=sharing) notebook available. However, this notebook requires a full download and installation of the language models for each run. When making many predictions, we recommend the local installation above, as this will allow you to cache and reuse the models.

### Paper analysis scripts

To reproduce the analysis in the paper, first download and extract data with the commands: 

```bash
git clone https://github.com/brianhie/efficient-evolution.git
cd efficient-evolution

wget https://zenodo.org/record/6968342/files/data.tar.gz
tar xvf data.tar.gz
```

### Install dependancy

Using conda package manager is easiest one.

```bash
conda create -f enviroments.yml
conda activate evo
(evo)python bin/recommend.py [your_sequence]
```

> When your first run. It tried download esm model and it takes time. 


### To acquire mutations to a given antibody.

run the command

```bash
(evo)bash bin/eval_models.sh [antibody_name]
```
where `[antibody_name]` is one of `medi8852`, `medi_uca`, `mab114`, `mab114_uca`, `s309`, `regn10987`, or `c143`.

DMS experiments can be run with the command
```bash
(evo)bash bin/dms.sh
```
