# psycloned_nemo

This repo simply contains the PSyclone-processed NEMO source files for the
ORCA1_SI3_GO8 configuration on the Met Office dev_r10037_GPU branch.
PSyclone has been used to add OpenACC directives to the code base, intended
to be used with the 'managed memory' option of the PGI compiler.

Assuming you have cloned/downloaded the code to the root NEMO directory 
(containing the `makenemo` script) and it is in a directory named
`psycloned_nemo` then to build NEMO do:

    ./makenemo -n ORCA1_SI3_GO8 -r SPITZ12 -m YOUR_ARCH_FILE -e ${PWD}/psycloned_nemo del_key 'key_iomput key_mpp_mpi'

Note that YOUR_ARCH_FILE must contain appropriate flags for OpenACC
compilation with managed memory.
