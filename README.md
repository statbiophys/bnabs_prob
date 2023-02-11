# bnabs_prob

The present code, companion of the manuscript *``Determining probabilities of HIV-1 bNAb development in healthy and chronically infected individuals''*, is conceptually divided in two sections:

1) `ig`, where we process B-Cell Receptor (BCR) repertoires from both healthy and chronically infected patientes. Sequences are firstly annotated through igBlast, then sorted and grouped by cohorts. Secondly, IGoR infers their recombination statistics, finally producing cohort-specific recombination models.

2) `bnabs`, where we analyze bnabs features and put them in relation with their probability of generation and evolution according to the cohort-specific models inferred above.

Test datasets of 5'000 BCR heavy- and light-chain sequences for an healthy donor are provided under `ig/healthy_IgG`, allowing to directly run the present code. Under `bnabs/sequences`, instead, we provided heavy- and light-chain sequences for the bnabs analyzed in the manuscript, together with a summary of their features relevant for the present analysis (e.g. their neutralization breadth). Some IGoR models, inferred on the whole cohort of healthy patients, are also provided under the folder `templates/igor_models`, again for testing purposes.

All the scripts, written in Python3, can hence be run with the attached test datasets. However, they rely on a local installation of igBlast and IGoR softwares. See following sections for installation and configuration details. V(D)J templates for igBlast and standard IGoR models are also provided here, under the folder `templates`.

-

### igBlast installation

[igBlast](https://www.ncbi.nlm.nih.gov/igblast/index.cgi) is a powerful and versatile Ig annotation software. We used the following releases:

- blastn: [2.9.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/)
- igblastn: [1.13.0](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.13.0/)

with V,D,J templates extracted from [IMGT](https://www.imgt.org), formatted and attached to this release.

Please go through the following of this section for the instructions on how to install igBlast and produce the templates in the desired format.

##### Installation of the Blast command line tools

Go to the webpage:
[https://www.ncbi.nlm.nih.gov/books/NBK279671/](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
and follow the instructions. It is needed for eg building the database afterwards. A better guide can be found at: [https://ncbi.github.io/igblast/](https://ncbi.github.io/igblast/) (recommended).

Unwrap the tar.gz file through the command: `tar xvzf file_name`

##### Installation of igBlast

Go to the webpage:
[https://ncbi.github.io/igblast/cook/How-to-set-up.html](https://ncbi.github.io/igblast/cook/How-to-set-up.html)
for the main instructions.

Unwrap the `.tar.gz` file through the command: `tar xvzf file_name`

Change the permissions to directories and files downloaded, through the command: `chmod -R u+rw *`

##### Setting paths

Paths can be saved into `.bashrc` file through commands like:

- `export PATH=$PATH:$HOME/igBlast/ncbi-blast-2.9.0+/bin:$HOME/igBlast/ncbi-igblast-1.13.0/bin`

- `export BLASTDB=$HOME/igBlast/blastdb`

- `export IGDATA=$HOME/igBlast`

##### Making the database

IG databases can be downloaded from IMGT at:
[http://www.imgt.org/vquest/refseqh.html](http://www.imgt.org/vquest/refseqh.html), at the section 'IG "V-REGION", "D-REGION", "J-REGION", "C-GENE exon" sets'. Ungapped germline genes should be downloaded from the column 'F+ORF+all P'.
Each page should be opened, sequences copied and then pasted into a new file in the germline tree.

- Another huge database from IMGT can be downloaded at:
[http://www.imgt.org/download/GENE-DB/](http://www.imgt.org/download/GENE-DB/)
but it includes several species in the same files, and for each of them, also several 'strange' genes that do not appear in the above database.
Again, be careful to download ungapped genes from the 'F+ORF+all P' section.

- Sequences too short can cause problems when creating the database; a work-around is to shorten comments preceeding the sequence.
Also, one could also keep just the name of the gene and nothing else.
To do that on the IMGT database, use the following command:
`awk '{if(gsub(">",">")==1){split($0,a,"|"); print ">"a[2]}else{print $1}}' FILE_IN > FILE_OUT`
This issue should have been fixed in the newer version of igBlast.

- Otherwise (and more easily), run the following command contained into the Blast script "edit_imgt_file.pl":
`./edit_imgt_file.pl imgt_file > my_seq_file`
or, even better, rely on the following script (based on the one above) that both filters the IMGT nomenclature and pastes together different lines of the same sequence, putting a space between them:
`./filterAndSort_IMGT_templates.sh`
each time choosing in the header of the script the kind of database you want to build.

- Database can be made through the command
`makeblastdb -in database_file.fasta -parse_seqids -dbtype nucl`
for each of the three gene types (V,D,J).

-

### IGoR installation

After a first sequence annotation and a quality filtering through igBlast, the dataset is now ready to be analyzed through [IGoR](https://github.com/statbiophys/IGoR) software, which allows to infer V(D)J recombination related processes from sequencing data.

The underlying methodology and some biologically relevant results are described in the following [paper](https://www.nature.com/articles/s41467-018-02832-w):

- Quentin Marcou, Thierry Mora, Aleksandra M. Walczak. *''High-throughput immune repertoire analysis with IGoR''*. Nature Communications 9, 561 (2018). 

We relied on a local installation of the [1.4.1](https://github.com/statbiophys/IGoR/releases/tag/1.4.1) release. Please refer to its [documentation](https://statbiophys.github.io/IGoR/) for a complete installation and usage guide.

-

### The `ig` section

There are two main scripts, `launch_annotation.py` and `launch_igor.py`, respectively for sequences annotation through igBlast and IGoR evaluation and model inference. The analysis run by each script can be customized through a dedicated `.yaml` config file.

The **annotation** step can be invoked as:

`python3 launch_annotation.py`

with settings fully customizable from the `config_annotation.yaml` file. The user can modify the path of input data, choose the type of chain (heavy, kappa or lambda) and the cohort such data come from, and also if a certain step of the script has to be run or not (e.g., the user can choose to run only the core of igBlast annotation, leaving the parsing of igBlast intermediate output for later).

Its output consists in:

- a `.igBlast_statistics` file, csv-formatted, with annotation results for each sequence (e.g., best-scoring V template, number and position of point-mutations, and so on);
- a set of `.csv` and `.fasta` files, where annotated sequences are sorted, grouped and filtered according to customizable criterions (e.g., include sequences with indels or not, divide in-frame sequences from out-of-frame, and so on).

Both kinds of files are produced separately for each donor, and at the cohort level, by grouping together sequences from different donors belonging to the same cohort under examination.

The **IGoR evaluation/inference** step, launched through the command:

`python3 launch_igor.py`

is again fully customizable from the `config_igor.yaml` file (e.g., the user can choose to run only the alignment step, leaving the inference or the final evaluation for later).

The output consists into a set of folders (`aligns`, `evaluate`, `inference`, `output`) containing IGoR results, as described in its [documentation](https://statbiophys.github.io/IGoR/). In particular, the two files under the `inference` folder, `final_parms.txt` and `final_marginals.txt`, are the ones containing the inferred model to be used later for bnabs evaluation.

Also, a summary file is produced by combining for each sequence the results from igBlast annotation and IGoR evaluation, potentially useful for a deeper analysis at the single-sequence level, or also to extract summary statistics at the donor or cohort level (e.g., the distribution of the fraction of point-mutated positions).

-

### The `bnabs` section

In this section, for each bnab are reported the heavy- and light-chain sequences, together with some key features (e.g. their neutralization potency), taken from the literature.

Though already annotated, they can be re-annotated through the command:

`python3 launch_annotation.py`

in order to extract annotation results in the desired format and to prepare `.csv` and `.fasta` files for the following IGoR evaluation step. Again, through the `config_annotation.yaml` file, the user can customize the path of the input file, choose the type of chain to be analyzed, and so on.

Output files are exactly analogous with those of the `ig` annotation step.