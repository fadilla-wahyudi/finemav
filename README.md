# *FineMAV*
This stand-alone program implements the Fine-Mapping of Adaptive Variation (*FineMAV*) statistic for detection of positively selected variants.

The *FineMAV* score of the derived allele for each SNP is calculated by multiplying three metrics:

- Derived allele purity -- measure of population differentiation that ranges from 0 to 1, where a score of 1 indicates that all the derived alelles fall into a single population. 
- Derived allele frequency 
- The Combined Annotation-Dependent Depletion PHRED-scaled C-score (CADD_PHRED) -- a measure of functionality

For more information about how this program works, please visit our [wiki](https://github.com/fadilla-wahyudi/finemav/wiki/1.-Home).

## Installation

### GitHub
First, clone this repository.
```
git clone https://github.com/fadilla-wahyudi/finemav/
```
Give the binary file executable permissions.
```
cd finemav
chmod +x finemav
```
To access the file globally, move it to one of the directories in your PATH. To check the PATH variable, type in `echo $PATH`.
Update this


### GUI version
GUI versions are available for Linux and MacOS:

- [Linux](https://drive.google.com/file/d/1xBhQGpUhVd02kyIuevVIuqac4zJ_13Tm/view?usp=sharing)
- [MacOS](https://drive.google.com/file/d/1hHp1SFps89_pFRPPPGCnQHTXuzsEXu5t/view?usp=sharing)


## Citation
If you utilise this program, please cite the following papers together:
>Szpak, M., Mezzavilla, M., Ayub, Q. et al. *FineMAV*: prioritizing candidate genetic variants driving local adaptations in human populations. Genome Biol **19**, 5 (2018). https://doi.org/10.1186/s13059-017-1380-2


