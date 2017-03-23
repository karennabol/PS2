## Gene Regulatory Network Generator in Julia (JuGRN)

### Introduction ###
JuGRN is a code generation system that transforms simple descriptions of the connectivity of gene regulatory networks into model code written in the [Julia](http://julialang.org) programming language.

### Installation and Requirements
You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/JuGRN-Generator.git

or

	$ git clone https://github.com/varnerlab/JuGRN-Generator.git

To execute a code generation job, Julia must be installed on your machine along with the Julia packages ``ArgParse`` and ``JSON``.
Julia can be downloaded/installed on any platform.
The required [Julia](http://julialang.org) packages can be installed by executing the commands:

	julia> Pkg.add("ArgParse")

and

	julia> Pkg.add("JSON")

in the Julia REPL.  

### How do I generate model code? ###
To generate a GRN model, issue the command ``make_julia_model.jl`` from the command line (outside of the REPL):

	$ julia make_julia_model.jl -m <input path> -o <output path> -s <host type>

The ``make_julia_model.jl`` command takes four command line arguments:

Argument | Required | Default | Description
--- | --- | --- | ---
-m | Yes	| none | Path to model input file (your \*.net file)
-o | No	| current directory | Path where files are written
-s | No	| bacterial | Host type (bacterial \| mammalian)

### Format for the GRN model input file ###
JuGRN-Generator transforms structured flat files into GRN model code. JuGRN-Generator takes flat files of the form:

~~~
// ----------------------------------------------------------------------- //
// JuGRN interactions -
//
// Record:
// actor {activate(*),induce(*) | inhibit(*),repress(*)} (target,...)
// ---------------------------------------------------------------------- //

// three gene memory network -
gene_1 induces (gene_2,gene_3)
gene_2 activates gene_3
gene_3 activates gene_2

~~~

The model specification file (by default given the filename `Network.net`) defines the biology of the model that gets generated.
# PS1
