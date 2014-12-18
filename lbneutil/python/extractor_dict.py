#!/usr/bin/env python
import sys, getopt
import os
import subprocess
from subprocess import Popen, PIPE
import time
import samweb_cli
#from samweb_client.utility import fileEnstoreChecksum
import ast
import project_utilities, root_metadata
import json

def getmetadata(inputfile):
        # Set up the experiment name for samweb Python API
	samweb = samweb_cli.SAMWebClient(experiment=project_utilities.get_experiment())

	# Extract metadata into a pipe.
	local = project_utilities.path_to_local(inputfile)
	if local != '':
		proc = subprocess.Popen(["sam_metadata_dumper", local], stdout=subprocess.PIPE)
	else:
		url = project_utilities.path_to_url(inputfile)
		proc = subprocess.Popen(["sam_metadata_dumper", url], stdout=subprocess.PIPE)
	if local != '' and local != inputfile:
		os.remove(local)

	mdtext=''
	for line in proc.stdout.readlines():
		if line[-3:-1] != ' ,':
			mdtext = mdtext + line.replace(", ,",",")
	mdtop = json.JSONDecoder().decode(mdtext)
	if len(mdtop.keys()) == 0:
		print 'No top-level key in extracted metadata.'
		sys.exit(1)
	file_name = mdtop.keys()[0]
	mdart = mdtop[file_name]

	# define an empty python dictionary which will hold sam metadata.
	# Some fields can be copied directly from art metadata to sam metadata.
	# Other fields require conversion.
	md = {}

	# Loop over art metadata.
	for mdkey in mdart.keys():
		mdval = mdart[mdkey]

		# Skip some art-specific fields.

		if mdkey == 'file_format_version':
			pass
		elif mdkey == 'file_format_era':
			pass

		# Ignore primary run_type field (if any).
		# Instead, get run_type from runs field.

		elif mdkey == 'run_type':
			pass

		# Ignore data_stream for now.

		elif mdkey == 'data_stream':
			pass

		# Ignore process_name for now.

		elif mdkey == 'process_name':
			pass

		# Application family/name/version.

		elif mdkey == 'applicationFamily':
			if not md.has_key('application'):
				md['application'] = {}
			md['application']['family'] = mdval
		elif mdkey == 'StageName':
			if not md.has_key('application'):
				md['application'] = {}
			md['application']['name'] = mdval
		elif mdkey == 'applicationVersion':
			if not md.has_key('application'):
				md['application'] = {}
			md['application']['version'] = mdval

		# Parents.

		elif mdkey == 'parents':
			mdparents = []
			for parent in mdval:
				parent_dict = {'file_name': parent}
				mdparents.append(parent_dict)
			md['parents'] = mdparents

		# Other fields where the key or value requires minor conversion.

		elif mdkey == 'first_event':
			md[mdkey] = mdval[2]
		elif mdkey == 'last_event':
			md[mdkey] = mdval[2]
		elif mdkey == 'lbneMCGenerators':
			md['lbne_MC.generators']  = mdval
		elif mdkey == 'lbneMCOscillationP':
			md['lbne_MC.oscillationP']  = mdval
		elif mdkey == 'lbneMCTriggerListVersion':
			md['lbne_MC.trigger-list-version']  = mdval
		elif mdkey == 'lbneMCBeamEnergy':
			md['lbne_MC.beam_energy']  = mdval
		elif mdkey == 'lbneMCBeamFluxID':
			md['lbne_MC.beam_flux_ID']  = mdval
		elif mdkey == 'lbneMCName':
			md['lbne_MC.name']  = mdval
		elif mdkey == 'lbneMCDetectorType':
			md['lbne_MC.detector_type']  = mdval
		elif mdkey == 'lbneMCNeutrinoFlavors':
			md['lbne_MC.neutrino_flavors']  = mdval
		elif mdkey == 'lbneMCMassHierarchy':
			md['lbne_MC.mass_hierarchy']  = mdval
		elif mdkey == 'lbneMCMiscellaneous':
			md['lbne_MC.miscellaneous']  = mdval
		elif mdkey == 'lbneMCGeometryVersion':
			md['lbne_MC.geometry_version']  = mdval
		elif mdkey == 'lbneMCOverlay':
			md['lbne_MC.overlay']  = mdval
		elif mdkey == 'lbneDataRunMode':
			md['lbne_data.run_mode']  = mdval
		elif mdkey == 'lbneDataDetectorType':
			md['lbne_data.detector_type']  = mdval
		elif mdkey == 'lbneDataName':
			md['lbne_data.name']  = mdval

		# For all other keys, copy art metadata directly to sam metadata.
		# This works for run-tuple (run, subrun, runtype) and time stamps.

		else:
			md[mdkey] = mdart[mdkey]

	# Get the other meta data field parameters						
	md['file_name'] =  inputfile.split("/")[-1]
	md['file_size'] =  os.path.getsize(inputfile)
	md['crc'] = root_metadata.fileEnstoreChecksum(inputfile)
	return md

if __name__ == "__main__":
	md = getmetadata(str(sys.argv[1]))
	#print md	
	mdtext = json.dumps(md, indent=2, sort_keys=True)
	print mdtext
	sys.exit(0)	
