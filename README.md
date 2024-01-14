# BeProf
 A comprehensive computational benchmark for evaluating deep learning-based protein function prediction approaches

## Overview
The evaluation for Fmax, Smin, AUPR, IC_AUPR and DP_AUPR. The details can be obtained from our article.

## Usage

You can directly run our code on your results as follows:
```
python evaluation.py --predict pred_file --output_path out_path --true true_file --background bk_file --go go_file --metrics metric_idx

arguments:
    predict: your own predicted result file
    output_path: the output directory
    true: the ground truth of test proteins
    background: the information of known proteins
    go: relations between GO terms
    metrics: the idx of metrics need to be calculated
```
Example
```
python evaluation.py --predict './Example/predict_result.pkl' --output_path './result' --true './Example/true_test.pkl' --background './Example/all_protein_information.pkl' --go './Example/go.obo' --metrics '0,1,2,3,4'
```

### Details of parameters and data format

- `pred_file`：your own predicted result file (`.pkl` is required). It is a dictionary as follows：

```
- Protein
    - tag (bp/mf/cc)
        - GO term id -> predicted probability

Example:
{
    '5MP1_HUMAN': {
        'bp': {'GO:0008150': 0.2505, 
               'GO:0000003': 0.2223, 
               'GO:0009653': 0.2223, ...},
        'cc': {'GO:0005575': 1.0, 
               'GO:0005737': 1.0, 
               'GO:0032991': 0.0281, ...},
        'mf': {'GO:0003674': 0.2505,
               'GO:0003676': 0.0281,
               'GO:1901363': 0.0281, ...}}
}
```

- `true_file`: the ground truth of test proteins (`.pkl` is required). It is a dictionary as follows：
```
- Protein
    - tag (all_bp/all_mf/all_cc) -> GO term ids (set format)

Example:
{
    '5MP1_HUMAN': {
        'all_bp': {'GO:0051171', 'GO:0010608', 'GO:0006446', 'GO:0060255', ...}, 
        'all_cc': {'GO:0005575', 'GO:0110165', 'GO:0005622', 'GO:0005737'}, 
        'all_mf': set() }
}
```

- `bk_file`: the background of all known protein, which is used to calculate the Information Content (`IC`) and `Smin` (`.pkl` is required). The data format is the same as `test_data_separate_file`.

- `go_file`: the relations between GO terms and can be downloaded directly from https://current.geneontology.org/ontology/go.obo (`.obo` is required).

- `metric_idx`: The metrics need to be calculated.
```
0 stands for F_max
1 stands for Smin
2 stands for AUPR
3 stands for ICAUPR
4 stands for DPAUPR
```
