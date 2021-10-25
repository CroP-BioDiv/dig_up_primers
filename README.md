# Dig up primers

Workflow for finding SSRs primers on one or more assemblies of related species.


## Instalation

Workflow is implemented as python3 script [dig_up_primers.py](dig_up_primers.py). Script uses external software:

* [MISA](https://webblast.ipk-gatersleben.de/misa/),
* [RepeatMasker](https://www.repeatmasker.org/),
* [Primer3](https://primer3.org/),
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).


## Strategy

Strategy is, a standard one:

* locate SSRs,
* design PCR primers around SSRs,
* filter primers by wanted characteristics,
* check primers that amplify exactly once in each assembly.


### SSR location

SSR location is done with perl script MISA. Reported SSRs are additionally filtered by removing SSRs that are

* complex repeat,
* too long,
* to be close to other SSR or sequnece end.


### Primer design

PCR primer design is done with Primer3 program. For each SSR, region of sequence sorrounding it is extracted, and Primer3 is run on that data.
Additional Primer3 restriction parameters are set to filter regions without promissing primers.


### Low complexity check

Each region around SSR is screened for low complexity with RepeatMasker program. SSR which regions are found to be of low complexity (except SSR part) are removed from further processing.


### Amplification

In silico aplification is done with blastn program (BLAST). Assemblies are searched for both primer sequences, located near one to the other, in the right direction.
If any assembly doesn't contain primer location, or contains more than one primer location, than that primer is discarded.


## Usage

Command accept lot of arguments, which are mostly passed to external programs. Arguments have specified defaults, set to resonable values. Check help:
```
python3 dig_up_primers.py -h
```


The most simple usage to specify project's directory (preferably new one) and assemblies to work with:
```
python3 dig_up_primers.py -w <project_directory> -a <path_to_fasta_file> [-a <path_to_fasta_file>]*
```

Script creates project directory if needed and stores values of command's arguments in it. Subsequent scipt calls in project's directory will continue calculation (if stopped).


## Projects structure

Project directory is structures into subdirectories and files covering each step of a process. Subdirectories contain directory for each assembly file.


## Citation

Beier S, Thiel T, Münch T, Scholz U, Mascher M (2017) MISA-web: a web server for microsatellite prediction. Bioinformatics 33 2583–2585. [dx.doi.org/10.1093/bioinformatics/btx198](http://dx.doi.org/10.1093/bioinformatics/btx198)

Thiel T, Michalek W, Varshney R, Graner A (2003) Exploiting EST databases for the development and characterization of gene-derived SSR-markers in barley (_Hordeum vulgare_ L.). Theoretical and Applied Genetics 106 (3): 411-422

Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG. Primer3--new capabilities and interfaces.
Nucleic Acids Res. 2012 Aug 1;40(15):e115.

Koressaar T and Remm M. Enhancements and modifications of primer design program Primer3. Bioinformatics 2007;23(10):1289-1291.

Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 <http://www.repeatmasker.org>.

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.