# Kasumi benchmark

R scripts used to analyse publicly available data, benchmark Kasumi against baselines and related methods and generate plots for the Kasumi paper.

For package requirements see [utils.R](utils.R). igraph version 1.5.1 is required to reproduce the exact results from the paper. Aternatively, a Dockerfile is provided for convenience.

Requires a development version of mistyR (>= 1.99.10) that can be obtained from [jtanevski/mistyR](https://github.com/jtanevski/mistyR).
The Kasumi R package can be obtained from [jtanevski/kasumi](https://github.com/jtanevski/kasumi).

Streamlined examples of running Kasumi on the CTCL CODEX dataset and the PDAC SMI dataset to reproduce the results shown in the paper can found in the examples folder.

To reproduce the analysis of the [DCIS dataset](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8792442/), run the following Python code in the `Image Data/Segmetation_Outlines_and_Labels_Mendeley/` folder in order to generate the required csv files

```Python
import glob
import re
from skimage import io
from skimage.measure import regionprops_table
import pandas

for f in glob.glob("*_labels.tiff"):
    mask = io.imread(f)
    rp = regionprops_table(mask,  properties=('label', 'centroid'))
    pandas.DataFrame(rp).to_csv(re.findall("\d+", f)[0] + ".csv", index = False)
```