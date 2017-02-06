#!/bin/bash

#
# run number code
#  XZAAAYNN
#    X == 1 zexp, 2 def
#    Z == 1 for anitnu, 2 for nu
#  AAA == A for target (dominant)
#    Y == 1 for mixed material, == 0 for pure
#   NN == run number
#

# this is silliness, python not required given all the run number structure we ended up with...

python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/default_2_12_2/ gntp.21\(.*\)ghep.root default_nubar_qe_like_scint.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/default_2_12_2/ gntp.220120\(.*\)ghep.root default_nu_qe_like_carbon.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/default_2_12_2/ gntp.220121\(.*\)ghep.root default_nu_qe_like_scint.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/default_2_12_2/ gntp.220560\(.*\)ghep.root default_nu_qe_like_iron.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/default_2_12_2/ gntp.222070\(.*\)ghep.root default_nu_qe_like_lead.txt

python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/z_exp_2_12_2/ gntp.11\(.*\)ghep.root zexp_nubar_qe_like_scint.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/z_exp_2_12_2/ gntp.120120\(.*\)ghep.root zexp_nu_qe_like_carbon.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/z_exp_2_12_2/ gntp.120121\(.*\)ghep.root zexp_nu_qe_like_scint.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/z_exp_2_12_2/ gntp.120560\(.*\)ghep.root zexp_nu_qe_like_iron.txt
python walk_for_file_list.py /pnfs/minerva/persistent/users/perdue/GENIE_files/qe_like/z_exp_2_12_2/ gntp.122070\(.*\)ghep.root zexp_nu_qe_like_lead.txt
