
# Ribbon (Ribosomal Bring On)

version 1.2

Karel Estrada && Enrique Merino

# INTRODUCTION
Ribbons extracts all copies of mRNA ribosomal 16s from bacteria or archaea genomes, using conserved motifs in the 16S mRNA ribosomal to detect them in the genomes and evaluate the true regions that contain a copy of this.

The main purpose of Ribbon is to extract complete copies of the 16s, including its 3'prime, where we will find the signature of the anti-Shine-Dalgarno sequence.

# REQUIREMENTS:

-Perl v5.8 or above

-Perl Modules: 

    Getopt::Std

    List::Util

The best way to install a perl module would be CPAN

-MEME Suite (tested on v5.0.1) software (http://meme-suite.org/doc/download.html)


Depending on which shell you have, please add these environmental variables:

In bash shell, could be .bashrc or .bash_profile add:

    export MEMEribbon=/path/where/you/installed/ribbon_1.2/db

    PATH=$PATH:/path/where/you/installed/ribbon_1.2/bin


Please, make sure that everything is installed and running properly.
Once you defined the environmental variables you can run the script "ribbon" to chech if everything it's okay. You can also run the example in the "example" folder, like this:


    $ cd example

    $ ribbon ecoli.fna ecoli


if everything worked well, the md5sum firm of generated file "ecoli.rib_complete_tail.fna" must be equal to the md5sum firm in the "md5sum_ecoli.rib_complete_tail.txt" file.


Any questions:

Karel Estrada
karel[at]ibt.unam.mx / kjestradag@gmail.com
