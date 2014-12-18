#!/usr/bin/env python
#----------------------------------------------------------------------
#
# Name: project_utilities.py
#
# Purpose: A python module containing various experiment-specific
#          python utility functions.
#
# Created: 28-Oct-2013  H. Greenlee
#
#----------------------------------------------------------------------

import os

# Don't fail (on import) if samweb is not available.

try:
    import samweb_cli
except ImportError:
    pass

def get_dropbox(filename):
    
    # Get metadata.
    
    md = {}
    exp = 'lbne'
    if os.environ.has_key('SAM_EXPERIMENT'):
        exp = os.environ['SAM_EXPERIMENT']
    samweb = samweb_cli.SAMWebClient(experiment=exp)
    try:
        md = samweb.getMetadata(filenameorid=filename)
    except:
        pass

    # Extract the metadata fields that we need.
    
    file_type = ''

    if md.has_key('file_type'):
        file_type = md['file_type']

    if not file_type:
        raise RuntimeError, 'Missing or invalid metadata for file %s.' % filename

    # Construct dropbox path.

    path = '/lbne/data/lbnepro/dropbox/%s' % file_type
    return path

# Return fcl configuration for experiment-specific sam metadata.

def get_sam_metadata(project, stage):
    result = 'services.user.FileCatalogMetadataLBNE: {\n'
    for key in project.parameters:
        result = result + '  %s: "%s"\n' % (key, project.parameters[key])
    for key in stage.parameters:
        result = result + '  %s: "%s"\n' % (key, stage.parameters[key])
    result = result + '  StageName: "%s"\n' % stage.name
    result = result + '}\n'
    return result

# Function to return path to the setup_lbne.sh script

def get_setup_script_path():

    OASIS_DIR="/cvmfs/oasis.opensciencegrid.org/lbne/products/"
    FERMIAPP_DIR="/grid/fermiapp/lbne/software/"

    if os.path.isfile(FERMIAPP_DIR+"setup_lbne.sh"):
        setup_script = FERMIAPP_DIR+"setup_lbne.sh"
    elif os.path.isfile(OASIS_DIR+"setup_lbne.sh"):
        setup_script = OASIS_DIR+"setup_lbne.sh"
    else:
        raise RuntimeError, "Could not find setup script at "+FERMIAPP_DIR+" or "+OASIS_DIR

    return setup_script

# Construct dimension string for project, stage.

def dimensions(project, stage):

    dim = 'file_type %s' % project.file_type
    dim = dim + ' and data_tier %s' % stage.data_tier
    for key in project.parameters:
        if key == 'MCName':
            dim = dim + ' and lbne_MC.name %s' % project.parameters[key]
        if key == 'DataName':
            dim = dim + ' and lbne_data.name %s' % project.parameters[key]
    dim = dim + ' and version %s' % project.release_tag
    dim = dim + ' and application %s' % stage.name
    dim = dim + ' and availability: anylocation'
    return dim
