# MIRROR

Workflow of the **MIRROR (Mimicking Inverted Repeats to Recruit ADAR via Engineered Oligoribonucleotide)**.

<iframe src="./img/MIRROR.png" width="600" height="400"></iframe>

## Python packages and softwares

The script was tested on a centos machine (CentOS Linux release 7.8.2003).

|        python version        | 3.8   |
| :--------------------------: | ----- |
|            pandas            | 1.5.1 |
|          Biopython           | 1.79  |
| RNA python apis of ViennaRNA | 4.0.2 |
|          RNAhybrid           | 2.1.2 |

## Usage

The main script, **Alu_mimic.py**, generates gRNAs with different structural features extracted from the substrates. The **substrates** directory contains the identified inverted Alu pairs used for gRNA generation.

You can use *python Alu_mimic.py -h* to view the help information. Before running the script, ensure the necessary packages and software are installed, and replace the RNAhybrid path in the script with the path to your own RNAhybrid installation.

The output includes a CSV file containing all the generated gRNAs and basical information, and a log file showing structural information predcited by RNAhybrid. You can use *less -S* to view it on a linux machine.

<img src="./img/log.png" alt="workflow" width="300" height="200">

 

