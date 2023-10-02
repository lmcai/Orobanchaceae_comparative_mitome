sed "s/xxx/$1/g" codeml.ctl >codeml_run.ctl
./paml-4.10.7/bin/codeml codeml_run.ctl
mv mlc $1.mlc
