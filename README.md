# ggslicer
Prepare and visualize volumetric imaging data using `ggplot2`. 

`ggslicer` provides a variety of helper functions for plotting brain slices in R using `ggplot2`, built around [`RMINC`](https://github.com/Mouse-Imaging-Centre/RMINC). It thus requires `RMINC` to function. Given the dependency on `RMINC`, it only works with MINC files for now. The `minc-toolkit` provides `nii2mnc` and `mnc2nii` that efficiently and correctly converts between MINC and NIFTI file formats. 

A complete re-write of `ggslicer` built around the ITK framework is in progress, and is nearing completion, and will be made available as an R package, as opposed the the current set of functions that need to be sourced. This will allow for file format-independent reading of data and plotting. Regardless of this, consider using MINC due to its efficient and flexible storage of images and metadata (it is built on top of HDF5), handy features (minchistory), and well-developed ecosystem (`minc-toolkit`, `minc-stuffs`). 

## Usage

First source the files:

```
source("plotting_functions/plotting_functions.R")
source("plotting_functions/plotting_functions_labels.R")
```

Then, use any of the functions available in the global environment.

Four main functions allow for efficient wrangling of MINC data (stored as files on disk) into a data frame / tibble:

- `prepare_anatomy()`: read a minc file and return a data frame
- `prepare_masked_anatomy()`: read a minc file and mask file and return a data frame with mask values included
- `prepare_contours()`: read a minc file associated with a template/stats image and return a set of points that form contours at specified intensity levels
- `prepare_label_contours()`: read a minc file associated with an atlas/label image and return a set of points that form contours around specified regions

The output of each of these are a list, with one of the elements being the data frame to be plotted.
