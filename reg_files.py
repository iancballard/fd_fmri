import os
import sys
import time
import shutil
import os.path as op
from textwrap import dedent
import argparse

import matplotlib as mpl
mpl.use("Agg")

import nipype
from nipype import Node, SelectFiles, DataSink, IdentityInterface

import lyman
import lyman.workflows as wf
from lyman import tools

project = lyman.gather_project_info()

# Set roots of output storage
data_dir = project["data_dir"]
exp_name = 'ser_8mm'
exp = lyman.gather_experiment_info('ser_8mm', None)
subj_source = tools.make_subject_source(['fd_104'])
analysis_dir = op.join(project["analysis_dir"], exp_name)
working_dir = op.join(project["working_dir"], exp_name)
nipype.config.set("execution", "crashdump_dir", project["crash_dir"])

# Is this a model or timeseries registration?
regtype =  "model"
space = 'mni'
smoothing = 'unsmoothed'

# Are we registering across experiments?
cross_exp = False

subject_id = 'fd_104'

# Retrieve the right workflow function for registration
# Get the workflow function dynamically based on the space
warp_method = project["normalization"]
flow_name = "%s_%s_reg" % (space, regtype)
reg, reg_input, reg_output = wf.create_reg_workflow(flow_name,
													space,
													regtype,
													warp_method,
													False,
													cross_exp)

# Define a smoothing info node here. Use an iterable so that running
# with/without smoothing doesn't clobber working directory files
# for the other kind of execution
smooth_source = Node(IdentityInterface(fields=["smoothing"]),
					 iterables=("smoothing", [smoothing]),
					 name="smooth_source")

# Set up the registration inputs and templates
reg_templates = dict(
	masks="{subject_id}/preproc/run_*/functional_mask.nii.gz",
	means="{subject_id}/preproc/run_*/mean_func.nii.gz",
					 )

if regtype == "model":
	# First-level model summary statistic images
	reg_base = "{subject_id}/mvpa/"
	reg_templates.update(dict(
		search=op.join(reg_base, "searchlight_state_radius8.nii.gz")))

# Native anatomy to group anatomy affine matrix and warpfield
if space == "mni":
	aff_ext = "mat" if warp_method == "fsl" else "txt"
	reg_templates["warpfield"] = op.join(data_dir, "{subject_id}",
										 "normalization/warpfield.nii.gz")
	reg_templates["affine"] = op.join(data_dir, "{subject_id}",
									  "normalization/affine." + aff_ext)

# Rigid (6dof) functional-to-anatomical matrices
rigid_stem = "{subject_id}/preproc/run_*/func2anat_"
if warp_method == "ants" and space == "mni":
	reg_templates["rigids"] = rigid_stem + "tkreg.dat"
else:
	reg_templates["rigids"] = rigid_stem + "flirt.mat"
reg_lists = reg_templates.keys()

# Rigid matrix from anatomy to target experiment space
targ_analysis_dir = op.join(project["analysis_dir"], 'ser_8mm')
reg_templates["first_rigid"] = op.join(targ_analysis_dir,
									   "{subject_id}", "preproc",
									   "run_1", "func2anat_flirt.mat")

# Define the registration data source node
reg_source = Node(SelectFiles(reg_templates,
							  force_lists=reg_lists,
							  base_directory=analysis_dir),
				  "reg_source")

# Registration inputnode
reg_inwrap = tools.InputWrapper(reg, subj_source,
								reg_source, reg_input)
reg_inwrap.connect_inputs()

# The source node also needs to know about the smoothing on this run
reg.connect(smooth_source, "smoothing", reg_source, "smoothing")

# Set up the registration output and datasink
reg_sink = Node(DataSink(base_directory=analysis_dir), "reg_sink")

reg_outwrap = tools.OutputWrapper(reg, subj_source,
								reg_sink, reg_output)
reg_outwrap.set_subject_container()
reg_outwrap.sink_outputs("reg.%s" % space)

# Reg has some additional substitutions to strip out iterables
# and rename the timeseries file
reg_subs = [("_smoothing_", "")]
reg_outwrap.add_regexp_substitutions(reg_subs)

reg.base_dir = working_dir

reg.run()
