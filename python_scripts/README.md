# Instruction for Script_for_search_metagenomes.py

This script can be used to search and download the list of metagenomes available in the NCBI database. You can get acquainted with the parameters passed to the script via the-h flag:

## Using -h flag

```
$ python3 Script_for_search_metagenomes.py -h
usage: Script_for_search_metagenomes.py [-h] --output OUTPUT

Search all metagenomes in the NCBI database and return the result table.

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        Output file name
```
## Example

``` $ python3 Script_for_search_metagenomes.py -o file_name```

# Instruction for Download_assembly.py

This script can be used to download the selected metagenome from the NCBI database. You must submit a csv table to the input after the script runs "Script_for_search_metagenomes.py". You can get acquainted with the parameters passed to the script via the-h flag

## Using -h flag
```
$ python3 Download_assembly.py -h
usage: Download_assembly.py [-h] --term TERM --input INPUT --output_dir
                            OUTPUT_DIR

Downloading metagenomes by term.

optional arguments:
  -h, --help            show this help message and exit
  --term TERM, -t TERM  Term, string value
  --input INPUT, -i INPUT
                        Input .csv file
  --output_dir OUTPUT_DIR, -d OUTPUT_DIR
                        The path to the directory where you want to
                        save the files, without "/"
```
## Example
```
$ python3 Download_assembly.py -t metagenome_name -i list_metagenomes.csv -d /path/to/dir
```

# Instruction for Trees_builder.py

Use Python for this script


# Instruction Remove_gap.py

This script is necessary to remove gaps in "fasta" file. You can get acquainted with the parameters passed to the script via the-h flag

## Using -h flag

```
$ python3 Remove_gap.py -h
usage: Remove_gap.py [-h] --input INPUT --output OUTPUT

Remove gap in "fasta" file.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input file
  --output OUTPUT, -o OUTPUT
                        Output file
```
## Example
```
python3 Remove_gap.py -i input_file.fasta -o output_file.fasta
