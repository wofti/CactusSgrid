#!/bin/bash

read -p "Type the name of the project: " XXX

mkdir $XXX
cd $XXX

git clone giter@quark.physics.fau.edu:sgrid
cp ../MyConfig sgrid/
cd sgrid/
make git_clone
cd ../

echo ''
echo 'Sgrid downladed'
echo ''

cp ../dns.* ./
curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/ET_2022_05/GetComponents
chmod a+x GetComponents
./GetComponents dns.th
cp ../ll_comp.sh Cactus/
cd ..
mv $XXX ../../../

echo ''
echo 'Cactus downloaded'
echo 'All done'
echo 'Things to do now'
echo ''
echo '1. Go to the Sgrid directory and try to compile it by typing "make". If failed, you must properly set MyConfog and try again.'
echo ''
echo '2. Set the right path to the Sgrid library in your dns.cfg file. It probably supposed to be'
echo 'LIBDIRS = path_to_this_project/sgrid/lib, and then the library itself:'
echo 'LIBS = sgrid'
echo ''
echo '3. Go to Cactus directory and compile it using ll_comp or by typing:'
echo './simfactory/bin/sim setup-silent'
echo './simfactory/bin/sim build -j12 --thornlist ../dns.th --optionlist ../dns.cfg'
echo ''
echo '4. Set the right path to the Sgrid-generated initial data directory in dns.par. and build the skeleton of the run by typing:'
echo './simfactory/bin/sim create-run name_of_the_run --parfile ../dns.par'
echo ''
echo '5. See Cactus/arrangements/CactusSgrid/DNSdata/README for details'
echo ''
