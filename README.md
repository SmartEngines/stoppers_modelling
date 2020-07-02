# Fast stopping rule modelling experiments

This repository contains experimental data and code for the fast approximate modelling of the combination result for implementing a stopping rule for video stream text recognition.

To calculate the estimation values with the evaluated methods and measure execution time, use the `precalculations.ipynb` notebook (function `precalculate(field_types, method, parallel_processes)`). To only run the analysis and construct the plots (using the precalculated values), run `comparison.ipynb`. For our experiments the code was executed used Python 3.7.4 running under Jupyter notebook version 6.0.1.

Description of each source code file and data directory follows.

## Code overview

1\. `metrics.py` - python module containing implementations for basic functions required for the experiments;

2\. `combination.py` - python module containing implementation of the text string combination algorithm, described in [[this article](https://mmp.susu.ru/article/en/578)], used for computing the expected distance estimation using the [[base stopping method](https://arxiv.org/abs/1910.04107)];

3\. `combination_with_estimation.py` - python module containing a modified implementation of the text string combination algorithm, used for computing the expected distance estimation using Method A and Method B;

4\. `treap.py` - python module containing implementation of a balanced binary search tree - [[treap with random priorities](http://doi.acm.org/10.1145/277651.277660)];

5\. `precalculations.ipynb` - Jupyter notebook with code for precalculations of the error level (distance from the combined result to the correct text field value), estimation value, and timing;

6\. `comparison.ipynb` - Jupyter notebook with code for comparing the stopping methods and constructing plots (Figures 3, 4, and 5).

7\. `prepare_ic15_yvt.ipynb` - Jupyter notebook with code for preparing IC15-Train and YVT datasets.


## Overview of data directories

1\. `data_<DATASET>/` - directory with text field recognition results, used as a dataset for the experiments. Document fields were taken from [[MIDV-500 dataset](http://www.computeroptics.smr.ru/eng/KO/Annot/KO43-5/430515e.html)] and [[MIDV-2019 dataset](https://arxiv.org/abs/1910.04009)], recognized using text field recognition subsystem of [[Smart IDReader](https://doi.org/10.1109/ICDAR.2017.347)]. Fields are grouped into four field types (`data_<DATASET>/date` - numeric dates, `data_<DATASET>/docnum` - document numbers, `data_<DATASET>/latin` - Latin name components, `data_<DATASET>/mrz` - MRZ lines). Each field clip is stored in a Pickle format and has the filename format `<FIELDTYPE>_<CLIPID>_<FIELDNAME>.pkl`. Arbitrary text objects were taken from [[IC15 Text in Videos dataset](https://doi.org/10.1109/ICDAR.2015.7333942)] and [[YouTube Video Text dataset](https://doi.org/10.1109/WACV.2014.6836024)], and recognized using the pretrained @clovaai model [[available here](https://github.com/clovaai/deep-text-recognition-benchmark)]. Arbitrary text objects have field type `none`.
  The Pickle data format for each field clip is the following:
  ```
  {
    "clip_id": "<CLIPID>",
    "ideal": "<CORRECT_FIELD_VALUE>",
    "field_name": "<FIELDNAME>",
    "field_type": "<FIELDTYPE(date/docnum/latin/mrz)>",
    "clip": [ // list of per-frame recognition results
      [ // per-frame recognition result - list of character classification results
        { // character classification result
          "A": 0.9435345435, // "<CLASS_LABEL>": <MEMBERSHIP_ESTIMATION>
          "B": 0.0224234234,
          ...
        },
        ...
      ],
      ...
    ]
  }
  ```
2\. `precalc_base_<DATASET>/` - directory with precalculated values for analysis. [[Base stopping method](https://arxiv.org/abs/1910.04107)] is used. Data files are generated using the corresponding functions in `precalculations.ipynb`. Each processed field clip is stored in JSON format and has the filename format `<FIELDTYPE>_<CLIPID>_<FIELDNAME>_base_precalc.json`. 
  The JSON data format for each field clip is the following:
  ```
  [ // list of stages (one per each processed frame)
    [ // data for i-th stage
      0.02934239847,  // Normalized Generalized Levenshtein Distance from combined result to the correct value (error level)
      0.00001239872,  // Sum of modelled distances (see expected distance estimation (eq. 2))
      0.00000023423,  // Time required to obtain combined result (seconds)
      0.00693821907   // Time required to calculate the estimation value (seconds)
    ],
    ...
  ]
  ```

3\. `precalc_summation_<DATASET>/` - directory with precalculated values for analysis. Method A is used. Data files are generated using the corresponding functions in `precalculations.ipynb`. The format is the same as `precalc_base_<DATASET>/`.

4\. `precalc_treap_<DATASET>/` - directory with precalculated values for analysis. Method B is used. Data files are generated using the corresponding functions in `precalculations.ipynb`. The format is the same as `precalc_base_<DATASET>/`.