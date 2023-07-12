# NanoMito
Assembling mitochondrial genomes from native sequencing Nanopore data

Required dependencies:
- NanoFilt https://github.com/wdecoster/nanofilt (De Coster et al., 2018)
- Minimap2 https://github.com/lh3/minimap2 (Li, 2018)
- Samtools https://github.com/samtools/samtools (Li et al., 2009)
- Flye assembler https://github.com/fenderglass/Flye (Kolmogorov et al., 2019)

Optional dependencies: 
- Guppy_basecaller https://community.nanoporetech.com/downloads
- Canu assembler https://github.com/marbl/canu (Koren et al., 2017)

Adjusting Flye

Flye was developed for long read data and uses a minimum read overlap of 1000 bases. When assembling mitochondrial genomes, especially from degraded samples, shorter reads might still be valuable. The 1000 bases minimum can be circumvented by adjusting the "main.py" file in flye's installation folder.
Find the following piece of code: 
```sh
parser.add_argument("-m", "--min-overlap", dest="min_overlap", metavar="int",
                       type=lambda v: check_int_range(v, 1000, 10000),
                       default=None, help="minimum overlap between reads [auto]")
```
Change `(v, 1000, 10000)` to `(v, 100, 10000)` and your minimum read length overlap is set to 100 bases. Use really not recommended when you do have proper long reads. Short “long reads” reads might skew coverage locally ruining your assembly. Flye was probably hard coded like that for a reason!

Be sure required dependencies are installed. The script is not meant to be run as is, adjust file locations to your specific data locations. Between consecutive rounds of assembly, draft assemblies need to be curated manually and set as input for the next round.

Quality check
After successful mitochondrial assembly the quality can be checked by annotating the genome on Mitos webserver, http://mitos.bioinf.uni-leipzig.de (Bernt et al., 2013). The annotation will show missing and duplicate elements of the mitochondrial genome.
