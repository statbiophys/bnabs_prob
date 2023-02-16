# bnabs_prob

##### Cosimo Lupo &#169; 2018-2023

##### License

Free use of the present code is granted under the terms of the [GNU General Public License version 3](https://www.gnu.org/licenses/quick-guide-gplv3.html) (GPLv3).

##### Contacts

For any issue, question or bug, please write [us](mailto:cosimo.lupo89@gmail.com) an email.

---

### Summary

The present code, companion of the manuscript *"Determining probabilities of HIV-1 bNAb development in healthy and chronically infected individuals"*, is conceptually divided in two sections:

1) `ig`, where we process B-Cell Receptor (BCR) repertoires from either healthy or chronically infected patientes, grouped into three cohorts: "healthy\_control", "hiv1", and "hcv". Sequences are firstly annotated through igBlast, then sorted and grouped by cohorts. Secondly, IGoR infers their recombination and evolution statistics, finally producing cohort-specific models.

2) `bnabs`, where we annotate bnabs sequences, analyze their features and put them in relation with their probability of generation and evolution, according to the cohort-specific models inferred above.

Test datasets of 5'000 BCR heavy- and light-chain IgG sequences for two healthy donors are provided under `ig/datasets/`, allowing to directly run the present code. Under `bnabs/sequences`, instead, we provided heavy- and light-chain sequences for the bnabs analyzed in the manuscript, together with a summary of their features relevant for the present analysis (e.g. their neutralization breadth). Some IGoR models, inferred on the whole cohort of healthy patients, are also provided under the folder `templates/igor_models/inferred`, again for testing purposes.

All the scripts, written in Python3, can hence be run with the attached test datasets. However, they rely on a local installation of **igBlast** and **IGoR** softwares. Expand the sections below for installation and configuration details. V(D)J templates for igBlast and standard IGoR models are also shipped with this code, under the folder `templates`.

Finally, though most of the Python packages used in the scripts are quite common (e.g. `numpy` or `pandas`), some others are less frequent and could not be alredy present in standard Python3 distributions. It can be the case for the following packages:

- [biopython](https://pypi.org/project/biopython/)
- [PyYAML](https://pypi.org/project/PyYAML/)
- [ATrieGC](https://pypi.org/project/atriegc/)

that can be pip-installed as usual.

---

### igBlast installation

[igBlast](https://www.ncbi.nlm.nih.gov/igblast/index.cgi) is a powerful and versatile Ig annotation software. We used the following releases:
	
- blastn: [2.9.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/)
- igblastn: [1.13.0](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.13.0/)
	
with V,D,J templates extracted from [IMGT](https://www.imgt.org), formatted and attached to this release (under `templates`).

Please go through the following of this section for the instructions on how to install igBlast and produce the templates in the desired format.

<details><summary>Installation details</summary>
<p>

##### Installation of the Blast command line tools

Go to the webpage:
[https://www.ncbi.nlm.nih.gov/books/NBK279671/](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
and follow the instructions. It is needed for eg building the database afterwards. A better guide can be found at: [https://ncbi.github.io/igblast/](https://ncbi.github.io/igblast/) (recommended).

Unwrap the tar.gz file through the command:

```
tar xvzf file_name
```

##### Installation of igBlast

Go to the webpage:
[https://ncbi.github.io/igblast/cook/How-to-set-up.html](https://ncbi.github.io/igblast/cook/How-to-set-up.html)
for the main instructions.

Unwrap the `.tar.gz` file through the command:

```
tar xvzf file_name
```

Change the permissions to directories and files downloaded, through the command:

```
chmod -R u+rw *
```

##### Setting paths

Paths can be saved into `.bashrc` file through commands like:

- `export PATH=$PATH:$HOME/igBlast/ncbi-blast-2.9.0+/bin:$HOME/igBlast/ncbi-igblast-1.13.0/bin`

- `export BLASTDB=$HOME/igBlast/blastdb`

- `export IGDATA=$HOME/igBlast`

so to be able to use igBlast commands (e.g. `igblastn`) independently from the current working directory. Otherwise, the full path of such commands (e.g. `$HOME/igBlast/ncbi-igblast-1.13.0/bin/igblastn`) has to be used.

##### Making the database

IG databases can be downloaded from IMGT at:
[http://www.imgt.org/vquest/refseqh.html](http://www.imgt.org/vquest/refseqh.html), at the section 'IG "V-REGION", "D-REGION", "J-REGION", "C-GENE exon" sets'. Ungapped germline genes should be downloaded from the column 'F+ORF+all P'.
Each page should be opened, sequences copied and then pasted into a new file in the germline tree.

- Another huge database from IMGT can be downloaded at:
[http://www.imgt.org/download/GENE-DB/](http://www.imgt.org/download/GENE-DB/)
but it includes several species in the same files, and for each of them, also several other genes that do not appear in the above database.
Again, be careful to download ungapped genes from the 'F+ORF+all P' section.

- Sequences too short can cause problems when creating the database; a work-around is to shorten comments preceeding the sequence.
Also, one could also keep just the name of the gene and nothing else.
To do that on the IMGT database, use the following command: `awk '{if(gsub(">",">")==1){split($0,a,"|"); print ">"a[2]}else{print $1}}' FILE_IN > FILE_OUT`. This issue should have been fixed in the newer version of igBlast.

- Otherwise (and more easily), run the following command contained into the Blast script `edit_imgt_file.pl`: `./edit_imgt_file.pl imgt_file > my_seq_file` or, even better, rely on the following script (based on the one above) that both filters the IMGT nomenclature and pastes together different lines of the same sequence, putting a space between them: `./filterAndSort_IMGT_templates.sh`, each time choosing in the header of the script the kind of database you want to build.

- Database can be made through the command `makeblastdb -in database_file.fasta -parse_seqids -dbtype nucl` for each of the three gene types (V,D,J).

</p>
</details>

---

### IGoR installation

After a first sequence annotation and a quality filtering through igBlast, BCR sequences are then ready to be analyzed through [IGoR](https://github.com/statbiophys/IGoR) software, which allows to infer V(D)J recombination related processes from sequencing data, as well as hyper-mutation statistics.

The underlying methodology and some biologically relevant results are described in the following [paper](https://www.nature.com/articles/s41467-018-02832-w):

- Quentin Marcou, Thierry Mora, Aleksandra M. Walczak. *"High-throughput immune repertoire analysis with IGoR"*. Nature Communications 9, 561 (2018). 

We relied on a local installation of the [1.4.1](https://github.com/statbiophys/IGoR/releases/tag/1.4.1) release.

Please go through the following of this section for the instructions on how to install and use IGoR.

<details><summary>Installation details</summary>
<p>

The complete documentation of IGoR can be found [here](https://statbiophys.github.io/IGoR/), with all the installation details, a list of known issues and suggested solutions, and an exhaustive usage guide.

The desired release of IGoR can be downloaded from [GitHub](https://github.com/statbiophys/IGoR/releases).

##### On Linux platforms

Once in the IGoR root directory, three key commands have to be launched, one after the other, to install IGoR on Linux platforms:

```
./configure
make
make install
```

For a user-level installation, if e.g. the user has no administrator privilegese, the flag `--prefix` has to appended to the `./ configure` command. For example, in order to install IGoR under the user home directory:
`./configure --prefix=$HOME`

Then, `make` and `make install` steps can be executed smoothly as before.

At the end of the installation procedure, IGoR’s executable will appear under the `igor_src` folder.

In case of multiple IGoR versions, it is recommended to explicitly use the full path to the executable of the desired version, e.g. `$HOME/igor_1.4.1/igor_src/igor`. Otherwise, if correctly exported during the `make install` step, the `igor` command will be accessible from any location without any full path specification.

##### On MacOS platforms

The installation steps previously mentioned are based on OpenMP-compatible compilers, as e.g. GNU `gcc`. Unfortunately, the default Apple compiler, though still callable through the `gcc` command, does not belong to this class of compilers and hence will cause a [fatal error](https://statbiophys.github.io/IGoR/#troubleshoots).

The suggested workaround is to install [Homebrew](https://brew.sh) and then the GNU gcc compiler through it (gcc 7.x versions are recommended, as their compatibility with the IGoR release 1.4.1 used here is granted):

```
brew install gcc@7
```

The freshly installed compiler will not overwrite the default Apple compiler, so `gcc` command will still refer to the latter. To this aim, when launching the `./configure` command for installing IGoR, a further flag has to be used, so to explicitly tell which compiler has to be used during the installation. Referring to releases 7.x of GNU gcc:

```
./configure CC="gcc-7" CXX="g++-7"
```

`make` and `make install` command can then be executed smoothly.

</p>
</details>

---

### The `ig` section

There are two main scripts, `launch_annotation.py` and `launch_igor.py`, respectively for sequence annotation through igBlast and IGoR evaluation and model inference.

The analysis run by each script can be customized through a dedicated `.yaml` config file. Among the possible accessible parameters, please notice (and modify, if needed) the full path of both igBlast and IGoR executables.

The **annotation** step can be invoked as:

```
python3 launch_annotation.py
```

with settings fully customizable from the `config_annotation.yaml` file. The user can modify the path of input data, choose the type of chain (heavy, kappa or lambda) and the cohort such data come from, and also if a certain step of the script has to be run or not (e.g., the user can choose to run only the core of igBlast annotation, leaving the parsing of igBlast intermediate output for later).

Its output consists in:

- a `.igBlast_statistics` file, csv-formatted, with annotation results for each sequence (e.g., best-scoring V template, number and position of point-mutations, and so on);
- a set of `.csv` and `.fasta` files, where annotated sequences are sorted, grouped and filtered according to customizable criterions (e.g., include sequences with indels or not, divide in-frame sequences from out-of-frame, and so on).

Both kinds of files are produced separately for each donor, and at the cohort level, by grouping together sequences from different donors belonging to the same cohort under examination.

The **IGoR evaluation/inference** step, launched through the command:

```
python3 launch_igor.py
```

is again fully customizable from the `config_igor.yaml` file (e.g., the user can choose to run only the alignment step, leaving the inference or the final evaluation for later).

The output consists into a set of folders (`aligns`, `evaluate`, `inference`, `output`) containing IGoR results, as described in its [documentation](https://statbiophys.github.io/IGoR/). In particular, the two files under the `inference` folder, `final_parms.txt` and `final_marginals.txt`, are the ones containing the inferred model to be used later for bnabs evaluation. To this aim, they will also be copied automatically under the `template/igor_models/inferred` folder.

Also, a summary file is produced by combining for each sequence the results from igBlast annotation and IGoR evaluation, potentially useful for a deeper analysis at the single-sequence level, or also to extract summary statistics at the donor or cohort level (e.g., the distribution of the fraction of point-mutated positions).

The "hiv1" cohort can be further stratified by antiretroviral therapy (ART) treatment:

- "hiv1\_ART\_OFF"
- "hiv1\_ART\_ON"

and by serum neutralization breadth:

- "hiv1\_non\_neutralizers"
- "hiv1\_weak\_neutralizers"
- "hiv1\_intermediate\_neutralizers"
- "hiv1\_top\_neutralizers"

Related IGoR models can be inferred on these sub-cohorts, by modifying the `cohort` parameter in the `config_igor.yaml ` file.

The auxiliary script `funcs_lineages.py` can be used to infer lineages in 
annotated datasets. This file also contains RAxML commands used to 
reconstruct phylogenies and ancestral states in largest lineages, and 
functions to analyze phylogenies to quantify skewedness.


---

### The `bnabs` section

In this section, for each bnab are reported the heavy- and light-chain sequences, together with some key features (e.g. their neutralization potency), taken from the literature.

Though already annotated, they can be re-annotated through the command:

```
python3 launch_annotation.py
```

in order to extract annotation results in the desired format and to prepare `.csv` and `.fasta` files for the following IGoR evaluation step. Again, through the `config_annotation.yaml` file, the user can customize the path of the input file, choose the type of chain to be analyzed, and so on.

Output files are exactly analogous with those of the `ig` annotation step.

Finally, by means of IGoR, it is possible to evaluate recombination and evolution probabilities of bnabs. Once chosen the desired model in the `config_igor.yaml` file (among the default IGoR ones or those inferred in the `ig` step on the three cohorts), the command:

```
python3 launch_igor.py
```

allows to get the aforementioned probabilities for each bnab. The output is of the same kind of that obtained in the `ig` step (apart from the `inference` folder, since here bnabs are just *evaluated* according to some IGoR models, and not used for further model inference).

Finally, as for the `ig` step, a summary file is produced, combining for each bnab the results from igBlast annotation and from IGoR evaluation according to a certain model, recombination and hyper-mutation probabilities, and neutralization properties. It's this set of files that is eventually used for the final bnab analysis, i.e. the assessment of the correlation between their probability of being generated and developed, and their neutralization properties.
