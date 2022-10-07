# DNS_preparer

The bash script which downloads sgrid and cactus and all their components necessary to compile cactus with DNSdata thorn.

One should simply execute "prepare_DNSdata.sh". The other files will be copied into the right places for further edition.

--------------------------------- Explonation of the "other files" ---------------------------------


- MyConfig    	- it is a file required by Sgrid to be compiled with a proper options. You are setting here the compiler, libraries, etc.
- dns.th      	- it is a thorns list required by Cactus (and by GetComponents script to download all necessary components to compile Cactus.
- dns.cfg		- it is a configuration file which contains a compilation options. In particular, the Sgrid library path is set in it.
- dns.par		- it is an example of Cactus input parameter file. There is a place in it to set the parameters related to DNSdata thorn. In particular, you must set the path to Sgrid-generated initial data in it.
- ll_comp.sh	- it is just a sbath submission script which you may use to compile Cactus. Of course, it is a template to be modified before using. 


--------------------------------- What to do to make it working ---------------------------------

There are a couple of things to do to make it working. Please be aware that further steps may depend on the machine you are working at and your private modules, libraries, compilers, etc.

1. Go to the Sgrid directory and try to compile it by typing "make". If failed, you must properly set MyConfog and try again.

2. Set the right path to the Sgrid library in your dns.cfg file. It probably supposed to be
"LIBDIRS = path_to_this_project/sgrid/lib", and then the library itself:
"LIBS = sgrid"

3. Go to Cactus directory and compile it using ll_comp.sh or by typing:
"./simfactory/bin/sim setup-silent"
"./simfactory/bin/sim build -j12 --thornlist ../dns.th --optionlist ../dns.cfg"

4. Set the right path to the Sgrid-generated initial data directory in dns.par file: 
"DNSdata::sgrid_datadir = "path_to_ID_directory" ", 
and then build the skeleton of the run by typing:
"./simfactory/bin/sim create-run name_of_the_run --parfile ../dns.par"

5. See Cactus/arrangements/CactusSgrid/DNSdata/README for details related to the properties of the DNSdata thorn.

