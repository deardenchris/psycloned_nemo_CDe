# psycloned_nemo

This repo simply contains the PSyclone-processed NEMO source files for the
ORCA1_GO8 configuration (no MPI, no sea-ice) on the Met Office NEMO_4.0_mirror_SI3_GPU  branch.
PSyclone has been used to add OpenACC directives to the code base, intended
to be used with the 'managed memory' option of the PGI compiler.

The code was processed using PSyclone master branch (v2.0.0), with manual tweaks added to improve GPU performance. 

Assuming you have cloned/downloaded the code to the root NEMO directory 
(containing the `makenemo` script) and it is in a directory named
`psycloned_nemo_CDe` then to build NEMO do:

    ./makenemo -n MO_GO8_GPU_ocean_only -r SPITZ12 -m YOUR_ARCH_FILE -e ${PWD}/psycloned_nemo_CDe -j 4 add_key "IEEE_IS_NAN=ieee_is_nan key_nosignedzero" del_key "key_iomput key_mpp_mpi key_si3"

Note that YOUR_ARCH_FILE must contain appropriate flags for OpenACC
compilation with managed memory.
