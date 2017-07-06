# AssemblyChallenge
Short-read assembler using de Bruijn graphs.

The assembly process is divided into 3 phases:

Phase I - Auto Cutoff Detection

Build a de Bruijn graph from the (uncorrected) reads. Simplify the graph using tip removal and bubble collapse. Get
the k-mer coverage distribution from the graph and fit a Gauss-Poisson mixture model to identify a suitable
coverage cutoff.

Phase II - Read Error Correction

Error-correct the input reads using the cutoff threshold obtained in Phase I. The read error-correction consists
of two runs of spectral alignment followed by one run of the bestReplacement algorithm which tries to find and
correct the most likely error every read.

Phase III - Final De Bruijn Graph Construction

Re-build the de Bruijn graph on the error corrected reads. Simplify the graph using tip and bubble removal. Finally
apply a coverage cutoff.

# Usage
    java LaunchAssembly [reads.fasta] [output.fasta] [size] [k-mer]

where

*   `[reads.fasta]` - input fasta file

*   `[output.fasta]` - output file (fasta)

*   `[size]` - k-mer size (odd integer)

*   `[k-mer]` - (optional) print incoming and outgoing k-mers of [k-mer]

# Compilation
From project directory

    mkdir bin

    javac -classpath src/ src/LaunchAssembly.java -d bin

# JDK
Java SE Development Kit 8, Update 131 (JDK 8u131).




