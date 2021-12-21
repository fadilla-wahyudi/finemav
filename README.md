# *FineMAV*
This stand-alone program implements the Fine-Mapping of Adaptive Variation (*FineMAV*) statistic for detection of positively selected variants.

The *FineMAV* score of the derived allele for each SNP is calculated by multiplying three metrics:

- Derived allele purity -- measure of population differentiation that ranges from 0 to 1, where a score of 1 indicates that all the derived alelles fall into a single population. 
- Derived allele frequency 
- The Combined Annotation-Dependent Depletion PHRED-scaled C-score (CADD_PHRED) -- a measure of functionality

For more information about this program, including installing and running the software on Linux, please visit our [wiki](https://github.com/fadilla-wahyudi/finemav/wiki).

## News

### GUI version
![#89f015](https://via.placeholder.com/15/89f015/000000?text=+) 4th May 2021 : We recently released the GUI version of *FineMAV* for Linux. For more information on how to install and launch the GUI version of the software, please visit our [wiki](https://github.com/fadilla-wahyudi/finemav/wiki/5.-GUI-version-(Linux)).

- [Linux](https://drive.google.com/file/d/1ch8d5rGR9S-_bPrkTNGflKaBGMFaerwk/view?usp=sharing)


## Citation
If you utilise this program, please cite the following papers together:
> Wahyudi, F., Aghakhanian, F., Rahman, S. et al. Prioritising positively selected variants in whole-genome sequencing data using *FineMAV*. BMC Bioinformatics **22**, 604 (2021). https://doi.org/10.1186/s12859-021-04506-9
>Szpak, M., Mezzavilla, M., Ayub, Q. et al. *FineMAV*: prioritizing candidate genetic variants driving local adaptations in human populations. Genome Biol **19**, 5 (2018). https://doi.org/10.1186/s13059-017-1380-2
