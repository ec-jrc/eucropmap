universe       = docker
docker_image   = jeoreg.cidsn.jrc.it:5000/jeodpp-htcondor/d5-eucropmap_py3:1.1 
executable     = ./start_process_2levels.sh 
transfer_input_files = /eos/jeodpp/data/projects/REFOCUS/cropclassif/start_process.sh,/eos/jeodpp/data/projects/REFOCUS/cropclassif/s1_classify_2levels.py
arguments      = /eos/jeodpp/data/projects/REFOCUS/cropclassif/list_rasters_eu_stratum_all2.lst $(ClusterID) $(ProcId) /eos/jeodpp/data/projects/REFOCUS/data/S1_GS/S1_classif_v7/ /eos/jeodpp/data/projects/REFOCUS/classification/result_v7/models/mask/RFmodel_LUCAS_[2]_level_1_all-polygons_janv-jul2018_15122020best /eos/jeodpp/data/projects/REFOCUS/classification/result_v7/models/crops/RFmodel_LUCAS_[2]_level_2_all-polygons_janv-jul2018_15122020best-smallgridter

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
#job_machine_attrs = Machine  
#job_machine_attrs_history_length = 5           
#requirements = target.machine =!= MachineAttrMachine1 && target.machine =!= MachineAttrMachine2

request_memory = 200GB
request_cpus   = 1
output         = /eos/jeodpp/htcondor/processing_logs/REFOCUS/cropclassif.$(ClusterId).$(ProcId).out 
error          = /eos/jeodpp/htcondor/processing_logs/REFOCUS/cropclassif.$(ClusterId).$(ProcId).err 
log            = /eos/jeodpp/htcondor/processing_logs/REFOCUS/cropclassif.$(ClusterId).$(ProcId).log
batch_name = "EUcropmap_S2"
queue infile from /eos/jeodpp/data/projects/REFOCUS/cropclassif/list_rasters_eu_stratum_all2.lst 


