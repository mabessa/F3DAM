%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Generate postprocessing file                                           %
%                                                                         %
%  This function is called by the code: STEP1_DataDriven_main_code.m      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is part of the data-driven framework and follows the      %
% LICENSE AGREEMENT included in te HEADING of the file:                   %
%                                                                         %
%        STEP1_DataDriven_main_code.m                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function generate_postproc_file(DoE_point,DiscInput_point,UserSettings,jDoE,iDiscInp,jDoE_dir,postproc_dir,abq_inp_file_noExt,CPUid)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
Lx = UserSettings.analyses.RVE_dimensions(1);
Ly = UserSettings.analyses.RVE_dimensions(2);

nmat_phases = UserSettings.analyses.Number_Mat_Phases;


% alpha=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Python File for Post-processing Abaqus simulations               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Started creation of Python file');
%
postproc_scriptfile_name = strcat(postproc_dir,'/script_RVE_postproc_CPU',num2str(CPUid),'.py');
fid = fopen(postproc_scriptfile_name,'wt');
%
fprintf(fid,'#=====================================================================#\n');
fprintf(fid,'#\n');
fprintf(fid,'# Created by M.A. Bessa on %s\n',datestr(now));
fprintf(fid,'#=====================================================================#\n');
fprintf(fid,'from odbAccess import *\n');
fprintf(fid,'import os\n');
fprintf(fid,'import numpy\n');
fprintf(fid,'import collections\n');
fprintf(fid,'from copy import deepcopy\n');
fprintf(fid,'try:\n');
fprintf(fid,'    import cPickle as pickle  # Improve speed\n');
fprintf(fid,'except ValueError:\n');
fprintf(fid,'    import pickle\n');
fprintf(fid,'\n');
fprintf(fid,'def dict_merge(a, b):\n');
fprintf(fid,'    #recursively merges dict''s. not just simple a[''key''] = b[''key''], if\n');
fprintf(fid,'    #both a and b have a key who''s value is a dict then dict_merge is called\n');
fprintf(fid,'    #on both values and the result stored in the returned dictionary.\n');
fprintf(fid,'    if not isinstance(b, dict):\n');
fprintf(fid,'        return b\n');
fprintf(fid,'    result = deepcopy(a)\n');
% fprintf(fid,'    result = a\n');
fprintf(fid,'    for k, v in b.iteritems():\n');
fprintf(fid,'        if k in result and isinstance(result[k], dict):\n');
fprintf(fid,'                result[k] = dict_merge(result[k], v)\n');
fprintf(fid,'        else:\n');
fprintf(fid,'            result[k] = deepcopy(v)\n');
% fprintf(fid,'            result[k] = v\n');
fprintf(fid,'    return result\n');
fprintf(fid,'\n');
fprintf(fid,'#\n');
fprintf(fid,'# This output point corresponds to the following DoE point and a DiscInput point:\n');
fprintf(fid,'stringDoE = ''DoE''+str(%i)\n',jDoE);
fprintf(fid,'stringDiscInput = ''DiscInput''+str(%i)\n',iDiscInp);
fprintf(fid,'#\n');
fprintf(fid,'# Set the work directory for ABAQUS (where the odb file is):\n');
fprintf(fid,'RVEdir=''%s''\n',jDoE_dir);
fprintf(fid,'\n');
% fprintf(fid,'BCdir= RVEdir + ''/BC%i''\n',iBC);
% fprintf(fid,'#\n');
fprintf(fid,'# Set directory where you want the post-processing file to be generated\n');
fprintf(fid,'postproc_dir=''%s''\n',postproc_dir);
fprintf(fid,'#\n');
fprintf(fid,'# Name of the file job\n');
fprintf(fid,'jobname=''%s''\n',abq_inp_file_noExt);
fprintf(fid,'#\n');
for jmat_phase=1:nmat_phases
    %
    phase_name = strcat('Phase_',num2str(jmat_phase));
    %
    fprintf(fid,'# Define the name of the %s part\n',phase_name);
    fprintf(fid,'%s_part=''FinalRVE-1.%s''\n',phase_name,phase_name);
    fprintf(fid,'#\n');
end
%
fprintf(fid,'#\n');
fprintf(fid,'odbfile=jobname+''.odb'' # Define name of this .odb file\n');
fprintf(fid,'#\n');
fprintf(fid,'os.chdir(RVEdir)\n');
%
if UserSettings.gen_options.run_simul_option == 1 % Then the postprocessing file is associated to the CPU number
    fprintf(fid,'file_postproc_path = postproc_dir+''/''+''RVEs_postprocessing_variables_CPU%s.p''\n',num2str(CPUid));
else % Then the postprocessing file generated previously is unique and independent of CPU number
    fprintf(fid,'file_postproc_path = postproc_dir+''/''+''RVEs_postprocessing_variables.p''\n');
end
%
fprintf(fid,'# Flag saying post-processing file exists:\n');
fprintf(fid,'try:\n');
fprintf(fid,'    # Try to load a previous post-processing file with all information:\n');
fprintf(fid,'    readfile_postproc = open(file_postproc_path, ''r'')\n');
fprintf(fid,'    RVEs_data_old = pickle.load(readfile_postproc)\n');
fprintf(fid,'    readfile_postproc.close()\n');
fprintf(fid,'    postproc_exists = 1\n');
fprintf(fid,'except Exception, e:\n');
fprintf(fid,'    postproc_exists = 0 # Flag saying that there is no post-processing file saved previously\n');
fprintf(fid,'\n');
fprintf(fid,'#\n');
%
fprintf(fid,'try:\n');
fprintf(fid,'    # Try to open this RVE odb file\n');
fprintf(fid,'    RVEodb = openOdb(path=odbfile) # Open .odb file\n');
fprintf(fid,'except Exception, e:\n');
fprintf(fid,'    print >> sys.stderr, ''does not exist''\n');
fprintf(fid,'    print >> sys.stderr, ''Exception: %%s'' %% str(e)\n');
fprintf(fid,'    sys.exit(1) # Exit the code because there is nothing left to do!\n');
fprintf(fid,'\n');
fprintf(fid,'#\n');
fprintf(fid,'# Determine the number of steps in the output database.\n');
fprintf(fid,'mySteps = RVEodb.steps\n');
fprintf(fid,'numSteps = len(mySteps)\n');
fprintf(fid,'#\n');
for jmat_phase=1:nmat_phases
    %
    phase_name = strcat('Phase_',num2str(jmat_phase));
    %
    fprintf(fid,'%s_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''%s'']\n',phase_name,upper(phase_name));
    fprintf(fid,'%s_elSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].elementSets[''%s'']\n',phase_name,upper(phase_name));
end
%
fprintf(fid,'botEdge_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''BOTEDGE'']\n');
fprintf(fid,'botEdge_elSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].elementSets[''BOTEDGE'']\n');
fprintf(fid,'rightEdge_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''RIGHTEDGE'']\n');
fprintf(fid,'rightEdge_elSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].elementSets[''RIGHTEDGE'']\n');
fprintf(fid,'topEdge_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''TOPEDGE'']\n');
fprintf(fid,'topEdge_elSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].elementSets[''TOPEDGE'']\n');
fprintf(fid,'leftEdge_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''LEFTEDGE'']\n');
fprintf(fid,'leftEdge_elSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].elementSets[''LEFTEDGE'']\n');
%
fprintf(fid,'entireRVE_nSet = RVEodb.rootAssembly.nodeSets['' ALL NODES'']\n');
fprintf(fid,'entireRVE_elSet = RVEodb.rootAssembly.elementSets['' ALL ELEMENTS'']\n');
if UserSettings.gen_options.pbcs_option == 1 % Only if using PBCs
    fprintf(fid,'dummy1_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''%s'']\n',upper(UserSettings.analyses.dummy1));
    fprintf(fid,'dummy2_nSet = RVEodb.rootAssembly.instances[''FINALRVE-1''].nodeSets[''%s'']\n',upper(UserSettings.analyses.dummy2));
end
fprintf(fid,'Lx = %6.5E # RVE dimension along x\n',Lx);
fprintf(fid,'Ly = %6.5E # RVE dimension along y\n',Ly);
%
fprintf(fid,'#\n');
% fprintf(fid,'RVEframe = RVEodb.steps[mySteps.keys()[0]].frames[0] # Undeformed config.\n');
%
fprintf(fid,'# For each step, obtain the following:\n');
fprintf(fid,'#     1) The step key.\n');
fprintf(fid,'#     2) The number of frames in the step.\n');
fprintf(fid,'#     3) The increment number of the last frame in the step.\n');
fprintf(fid,'#\n');
fprintf(fid,'totalNumFrames = 0\n');
fprintf(fid,'for iStep in range(numSteps):\n');
fprintf(fid,'    stepKey = mySteps.keys()[iStep]\n');
fprintf(fid,'    step = mySteps[stepKey]\n');
fprintf(fid,'    numFrames = len(step.frames)\n');
fprintf(fid,'    totalNumFrames = totalNumFrames + numFrames\n');
fprintf(fid,'    #\n');
fprintf(fid,'\n');
fprintf(fid,'# Preallocate quantities for speed\n');
fprintf(fid,'RVEframe = RVEodb.steps[mySteps.keys()[0]].frames[0] # Undeformed config.\n');
%
if UserSettings.gen_options.analysis_option ~= 2
    %% NOT a steady state dynamics simulation
    %
    % Compute volume and volume at integration points for the ENTIRE RVE
    fprintf(fid,'# Extract volume at integration point in ENTIRE RVE:\n');
    fprintf(fid,'ivolField = RVEframe.fieldOutputs[''IVOL'']\n');
    fprintf(fid,'ivolSubField = ivolField.getSubset(region=entireRVE_elSet, position=INTEGRATION_POINT)\n');
    if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
        for jmat_phase=1:nmat_phases
            %
            phase_name = strcat('Phase_',num2str(jmat_phase));
            %
            fprintf(fid,'ivolSubField_%s = ivolField.getSubset(region=%s_elSet, position=INTEGRATION_POINT)\n',phase_name,phase_name);
        end
    end
    fprintf(fid,'#\n');
    fprintf(fid,'ivol=numpy.zeros(( len(ivolSubField.values) ))\n');
    fprintf(fid,'tot_vol = 0.0\n');
    if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
        for jmat_phase=1:nmat_phases
            %
            phase_name = strcat('Phase_',num2str(jmat_phase));
            %
            fprintf(fid,'ivol_%s=numpy.zeros(( len(ivolSubField_%s.values) ))\n',phase_name,phase_name);
            fprintf(fid,'tot_vol_%s = 0.0\n',phase_name);
        end
    end
    fprintf(fid,'for i in range(0, len(ivolSubField.values)):\n');
    fprintf(fid,'    ivol[i] = ivolSubField.values[i].data # Volume for i-th integration point\n');
    fprintf(fid,'    tot_vol = tot_vol + ivol[i] # total volume\n');
    fprintf(fid,'\n');
    if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
        for jmat_phase=1:nmat_phases
            %
            phase_name = strcat('Phase_',num2str(jmat_phase));
            %
            fprintf(fid,'for i in range(0, len(ivolSubField_%s.values)):\n',phase_name);
            fprintf(fid,'    ivol_%s[i] = ivolSubField_%s.values[i].data # Volume for i-th integration point\n',phase_name,phase_name);
            fprintf(fid,'    tot_vol_%s = tot_vol_%s + ivol_%s[i] # total volume\n',phase_name,phase_name,phase_name);
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'# finished computing volume at integration points and total volume\n');
    fprintf(fid,'#\n');
    %
    % fprintf(fid,'evolField = RVEframe.fieldOutputs[''EVOL'']\n');
    % fprintf(fid,'evolSubField = evolField.getSubset(region=entireRVE_elSet, position=WHOLE_ELEMENT)\n');
    % fprintf(fid,'#\n');
    % fprintf(fid,'evol=numpy.zeros(( len(evolSubField.values) ))\n');
    % fprintf(fid,'tot_vol = 0.0\n');
    % fprintf(fid,'for i in range(0, len(evolSubField.values)):\n');
    % fprintf(fid,'    evol[i] = evolSubField.values[i].data # Volume for i-th element\n');
    % fprintf(fid,'    tot_vol = tot_vol + evol[i] # total volume\n');
    % fprintf(fid,'\n');
    % fprintf(fid,'# finished computing element volume\n');
    % fprintf(fid,'#\n');
    %
    %
    % First consider quantities at the INTEGRATION POINTS:
    for ivar = 1:length(UserSettings.postproc.vars_int_pts)
        fprintf(fid,'%s_Field = RVEframe.fieldOutputs[''%s'']\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'%s_SubField = %s_Field.getSubset(region=entireRVE_elSet, position=INTEGRATION_POINT)\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'#\n');
        fprintf(fid,'if isinstance(%s_SubField.values[0].data,float):\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'    # Then variable is a scalar\n');
        fprintf(fid,'    av_%s = numpy.zeros(( totalNumFrames ))\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'else:\n');
        fprintf(fid,'    # Variable is an array\n');
        fprintf(fid,'    av_%s = numpy.zeros(( totalNumFrames,len(%s_SubField.values[0].data) ))\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'\n');
        %
        if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
            for jmat_phase=1:nmat_phases
                %
                phase_name = strcat('Phase_',num2str(jmat_phase));
                %
                fprintf(fid,'%s_SubField_%s = %s_Field.getSubset(region=%s_elSet, position=INTEGRATION_POINT)\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'#\n');
                fprintf(fid,'if isinstance(%s_SubField_%s.values[0].data,float):\n',UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'    # Then variable is a scalar\n');
                fprintf(fid,'    av_%s_%s = numpy.zeros(( totalNumFrames ))\n',UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'else:\n');
                fprintf(fid,'    # Variable is an array\n');
                fprintf(fid,'    av_%s_%s = numpy.zeros(( totalNumFrames,len(%s_SubField_%s.values[0].data) ))\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'\n');
            end
        end
    end
    % Then consider quantities at the WHOLE ELEMENT:
    for ivar = 1:length(UserSettings.postproc.vars_whole_els)
        fprintf(fid,'%s_Field = RVEframe.fieldOutputs[''%s'']\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'%s_SubField = %s_Field.getSubset(region=entireRVE_elSet, position=WHOLE_ELEMENT)\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'#\n');
        fprintf(fid,'if isinstance(%s_SubField.values[0].data,float):\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'    # Then variable is a scalar\n');
        fprintf(fid,'    av_%s = numpy.zeros(( totalNumFrames ))\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'else:\n');
        fprintf(fid,'    # Variable is an array\n');
        fprintf(fid,'    av_%s = numpy.zeros(( totalNumFrames,len(%s_SubField.values[0].data) ))\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'\n');
        %
        if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
            for jmat_phase=1:nmat_phases
                %
                phase_name = strcat('Phase_',num2str(jmat_phase));
                %
                fprintf(fid,'%s_SubField_%s = %s_Field.getSubset(region=%s_elSet, position=WHOLE_ELEMENT)\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'#\n');
                fprintf(fid,'if isinstance(%s_SubField_%s.values[0].data,float):\n',UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'    # Then variable is a scalar\n');
                fprintf(fid,'    av_%s_%s = numpy.zeros(( totalNumFrames ))\n',UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'else:\n');
                fprintf(fid,'    # Variable is an array\n');
                fprintf(fid,'    av_%s_%s = numpy.zeros(( totalNumFrames,len(%s_SubField_%s.values[0].data) ))\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'\n');
            end
        end
    end
    % Then consider quantities at the NODES:
    for ivar = 1:length(UserSettings.postproc.vars_nodes)
        fprintf(fid,'%s_Field = RVEframe.fieldOutputs[''%s'']\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'%s_SubField = %s_Field.getSubset(region=entireRVE_nSet, position=NODAL)\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'#\n');
        fprintf(fid,'if isinstance(%s_SubField.values[0].data,float):\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'    # Then variable is a scalar\n');
        fprintf(fid,'    av_%s = numpy.zeros(( totalNumFrames ))\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'else:\n');
        fprintf(fid,'    # Variable is an array\n');
        fprintf(fid,'    av_%s = numpy.zeros(( totalNumFrames,len(%s_SubField.values[0].data) ))\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'\n');
        %
        if UserSettings.postproc.local_averages == 1 % In case we want to have the separate average behavior at each material phase
            for jmat_phase=1:nmat_phases
                %
                phase_name = strcat('Phase_',num2str(jmat_phase));
                %
                fprintf(fid,'%s_SubField_%s = %s_Field.getSubset(region=%s_nSet, position=NODAL)\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'#\n');
                fprintf(fid,'if isinstance(%s_SubField_%s.values[0].data,float):\n',UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'    # Then variable is a scalar\n');
                fprintf(fid,'    av_%s_%s = numpy.zeros(( totalNumFrames ))\n',UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'else:\n');
                fprintf(fid,'    # Variable is an array\n');
                fprintf(fid,'    av_%s_%s = numpy.zeros(( totalNumFrames,len(%s_SubField_%s.values[0].data) ))\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'\n');
            end
        end
    end
    %
    %
else
    %% This is a STEADY STATE DYNAMICS simulation, so IVOL and EVOL are not available to compute averages...
    %
    fprintf(fid,'tot_vol = Lx * Ly\n'); % Just manually compute the total area (volume in 3D)
    for jmat_phase=1:nmat_phases
        %
        phase_name = strcat('Phase_',num2str(jmat_phase));
        %
        fprintf(fid,'tot_vol_%s = ''Not available''\n',phase_name); % This simulation type doesn't output local volumes
    end
    %
end

%% Compute average strain and stress measures from Displacements and Reaction Forces (faster)

% And finally consider quantities that lead to nominal stress P and deformation gradient F (if using PBCs):
if UserSettings.gen_options.pbcs_option == 1 % Then we are using Periodic Boundary Conditions (PBCs)
    postproc_variables_pbcs = {'U', 'RF'};
    for ivar = 1:length(postproc_variables_pbcs)
        fprintf(fid,'%s_Field = RVEframe.fieldOutputs[''%s'']\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'%s_dummy1_SubField = %s_Field.getSubset(region=dummy1_nSet, position=NODAL)\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'%s_dummy2_SubField = %s_Field.getSubset(region=dummy2_nSet, position=NODAL)\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'#\n');
        fprintf(fid,'if isinstance(%s_dummy1_SubField.values[0].data,float):\n',postproc_variables_pbcs{ivar});
        fprintf(fid,'    # Then variable is a scalar\n');
        fprintf(fid,'    %s_dummy1 = numpy.zeros(( totalNumFrames ))\n',postproc_variables_pbcs{ivar});
        fprintf(fid,'    %s_dummy2 = numpy.zeros(( totalNumFrames ))\n',postproc_variables_pbcs{ivar});
        if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
            fprintf(fid,'    %s_dummy1_imag = numpy.zeros(( totalNumFrames ))\n',postproc_variables_pbcs{ivar});
            fprintf(fid,'    %s_dummy2_imag = numpy.zeros(( totalNumFrames ))\n',postproc_variables_pbcs{ivar});
        end
        fprintf(fid,'else:\n');
        fprintf(fid,'    # Variable is an array\n');
        fprintf(fid,'    %s_dummy1 = numpy.zeros(( totalNumFrames,len(%s_dummy1_SubField.values[0].data) ))\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'    %s_dummy2 = numpy.zeros(( totalNumFrames,len(%s_dummy2_SubField.values[0].data) ))\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
            fprintf(fid,'    %s_dummy1_imag = numpy.zeros(( totalNumFrames,len(%s_dummy1_SubField.values[0].data) ))\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
            fprintf(fid,'    %s_dummy2_imag = numpy.zeros(( totalNumFrames,len(%s_dummy2_SubField.values[0].data) ))\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        end
        fprintf(fid,'\n');
    end
    %
else % NOT using PBCs
    % So we will need the reaction force on the left edge (RF_leftEdge)
    fprintf(fid,'RF_leftEdge_Field = RVEframe.fieldOutputs[''RF'']\n');
    fprintf(fid,'RF_leftEdge_SubField = RF_leftEdge_Field.getSubset(region=leftEdge_nSet, position=NODAL)\n');
    fprintf(fid,'#\n');
    fprintf(fid,'if isinstance(RF_leftEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'    # Then variable is a scalar\n');
    fprintf(fid,'    sum_RF_leftEdge = numpy.zeros(( totalNumFrames ))\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'    sum_imagRF_leftEdge = numpy.zeros(( totalNumFrames ))\n');
    end
    fprintf(fid,'else:\n');
    fprintf(fid,'    # Variable is an array\n');
    fprintf(fid,'    sum_RF_leftEdge = numpy.zeros(( totalNumFrames,len(RF_leftEdge_SubField.values[0].data) ))\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'    sum_imagRF_leftEdge = numpy.zeros(( totalNumFrames,len(RF_leftEdge_SubField.values[0].conjugateData) ))\n');
    end
    fprintf(fid,'\n');
    % And RF on bottom edge (RF_botEdge)
    fprintf(fid,'RF_botEdge_Field = RVEframe.fieldOutputs[''RF'']\n');
    fprintf(fid,'RF_botEdge_SubField = RF_botEdge_Field.getSubset(region=botEdge_nSet, position=NODAL)\n');
    fprintf(fid,'#\n');
    fprintf(fid,'if isinstance(RF_botEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'    # Then variable is a scalar\n');
    fprintf(fid,'    sum_RF_botEdge = numpy.zeros(( totalNumFrames ))\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'    sum_imagRF_botEdge = numpy.zeros(( totalNumFrames ))\n');
    end
    fprintf(fid,'else:\n');
    fprintf(fid,'    # Variable is an array\n');
    fprintf(fid,'    sum_RF_botEdge = numpy.zeros(( totalNumFrames,len(RF_botEdge_SubField.values[0].data) ))\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'    sum_imagRF_botEdge = numpy.zeros(( totalNumFrames,len(RF_botEdge_SubField.values[0].conjugateData) ))\n');
    end
    fprintf(fid,'\n');
    %
    % And the displacement on the right edge (U_rightEdge)
    fprintf(fid,'U_rightEdge_Field = RVEframe.fieldOutputs[''U'']\n');
    fprintf(fid,'U_rightEdge_SubField = U_rightEdge_Field.getSubset(region=rightEdge_nSet, position=NODAL)\n');
    fprintf(fid,'#\n');
    fprintf(fid,'if isinstance(U_rightEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'    # Then variable is a scalar\n');
    fprintf(fid,'    av_U_rightEdge = numpy.zeros(( totalNumFrames ))\n');
    fprintf(fid,'else:\n');
    fprintf(fid,'    # Variable is an array\n');
    fprintf(fid,'    av_U_rightEdge = numpy.zeros(( totalNumFrames,len(U_rightEdge_SubField.values[0].data) ))\n');
    fprintf(fid,'\n');
    %
    % And the displacement U on the top edge (U_topEdge)
    fprintf(fid,'U_topEdge_Field = RVEframe.fieldOutputs[''U'']\n');
    fprintf(fid,'U_topEdge_SubField = U_topEdge_Field.getSubset(region=topEdge_nSet, position=NODAL)\n');
    fprintf(fid,'#\n');
    fprintf(fid,'if isinstance(U_topEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'    # Then variable is a scalar\n');
    fprintf(fid,'    av_U_topEdge = numpy.zeros(( totalNumFrames ))\n');
    fprintf(fid,'else:\n');
    fprintf(fid,'    # Variable is an array\n');
    fprintf(fid,'    av_U_topEdge = numpy.zeros(( totalNumFrames,len(U_topEdge_SubField.values[0].data) ))\n');
    fprintf(fid,'\n');
end
%
fprintf(fid,'# Loop over Steps and Frames to compute average quantities in RVE\n');
fprintf(fid,'eye = numpy.identity(2)\n');
fprintf(fid,'defGrad = numpy.zeros(( totalNumFrames,2,2 ))\n');
fprintf(fid,'nomP = numpy.zeros(( totalNumFrames,2,2 ))\n');
if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
    fprintf(fid,'nomP_imag = numpy.zeros(( totalNumFrames,2,2 ))\n');
    fprintf(fid,'frequency = numpy.zeros( totalNumFrames )\n');
    fprintf(fid,'YoungsE1_real = numpy.zeros( totalNumFrames )\n');
    fprintf(fid,'YoungsE1_imag = numpy.zeros( totalNumFrames )\n');
end
fprintf(fid,'jacobian = numpy.zeros(totalNumFrames)\n');
% If using Green strain we can output it as well as the PK2 stress: 
if UserSettings.gen_options.strain_option == 1 % Then using we are using Green-strain:
    fprintf(fid,'Green_strain = numpy.zeros(( totalNumFrames,2,2 ))\n');
    fprintf(fid,'PK2 = numpy.zeros(( totalNumFrames,2,2 ))\n');
end
fprintf(fid,'stepTotalTime = numpy.zeros(numSteps)\n');
fprintf(fid,'previousFrame = 0\n');
fprintf(fid,'numFrames = 0\n');
%
fprintf(fid,'for iStep in range(numSteps):\n');
fprintf(fid,'    previousFrame = previousFrame + numFrames\n');
fprintf(fid,'    stepKey = mySteps.keys()[iStep]\n');
fprintf(fid,'    step = mySteps[stepKey]\n');
fprintf(fid,'    stepTotalTime[iStep] = step.timePeriod\n');
fprintf(fid,'    numFrames = len(step.frames)\n');
fprintf(fid,'    #\n');
fprintf(fid,'    for iFrame_step in range (0,numFrames):\n');
fprintf(fid,'        iFrame = previousFrame+iFrame_step\n'); % iFrame_step is the frame at this step
fprintf(fid,'        RVEframe = RVEodb.steps[stepKey].frames[iFrame]\n');
if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
    fprintf(fid,'        frequency[iFrame] = RVEframe.frequency\n');
end
fprintf(fid,'        #\n');
%
if UserSettings.gen_options.analysis_option ~= 2
    %% NOT a steady state dynamics simulation
    % First consider quantities at the INTEGRATION POINTS:
    for ivar = 1:length(UserSettings.postproc.vars_int_pts)
        fprintf(fid,'        # Variable: %s\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'        %s_Field = RVEframe.fieldOutputs[''%s'']\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'        %s_SubField = %s_Field.getSubset(region=entireRVE_elSet, position=INTEGRATION_POINT)\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        if isinstance(%s_SubField.values[0].data,float):\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'            # Then variable is a scalar:\n');
        fprintf(fid,'            # Loop over every element to compute average\n');
        fprintf(fid,'            for i in range(0, len(%s_SubField.values)):\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'                av_%s[iFrame] = av_%s[iFrame] + %s_SubField.values[i].data*ivol[i]\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'            #\n');
        fprintf(fid,'            av_%s[iFrame] = av_%s[iFrame]/tot_vol\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'        else:\n');
        fprintf(fid,'            # Variable is an array:\n');
        fprintf(fid,'            # Loop over every element to compute average\n');
        fprintf(fid,'            for j in range(0, len(%s_SubField.values[0].data)):\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'                for i in range(0, len(%s_SubField.values)):\n',UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'                    av_%s[iFrame][j] = av_%s[iFrame][j] + %s_SubField.values[i].data[j]*ivol[i]\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'                #\n');
        fprintf(fid,'                av_%s[iFrame][j] = av_%s[iFrame][j]/tot_vol\n',UserSettings.postproc.vars_int_pts{ivar},UserSettings.postproc.vars_int_pts{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        # Finished computing average for this variable!\n');
        %
        if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
            for jmat_phase=1:nmat_phases
                %
                phase_name = strcat('Phase_',num2str(jmat_phase));
                %
                fprintf(fid,'        # Variable: %s\n',UserSettings.postproc.vars_int_pts{ivar});
                fprintf(fid,'        %s_SubField_%s = %s_Field.getSubset(region=%s_elSet, position=INTEGRATION_POINT)\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'        #\n');
                fprintf(fid,'        if isinstance(%s_SubField_%s.values[0].data,float):\n',UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'            # Then variable is a scalar:\n');
                fprintf(fid,'            # Loop over every element to compute average\n');
                fprintf(fid,'            for i in range(0, len(%s_SubField_%s.values)):\n',UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'                av_%s_%s[iFrame] = av_%s_%s[iFrame] + %s_SubField_%s.values[i].data*ivol_%s[i]\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name,phase_name);
                fprintf(fid,'            #\n');
                fprintf(fid,'            av_%s_%s[iFrame] = av_%s_%s[iFrame]/tot_vol_%s\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name,phase_name);
                fprintf(fid,'        else:\n');
                fprintf(fid,'            # Variable is an array:\n');
                fprintf(fid,'            # Loop over every element to compute average\n');
                fprintf(fid,'            for j in range(0, len(%s_SubField_%s.values[0].data)):\n',UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'                for i in range(0, len(%s_SubField_%s.values)):\n',UserSettings.postproc.vars_int_pts{ivar},phase_name);
                fprintf(fid,'                    av_%s_%s[iFrame][j] = av_%s_%s[iFrame][j] + %s_SubField_%s.values[i].data[j]*ivol_%s[i]\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name,phase_name);
                fprintf(fid,'                #\n');
                fprintf(fid,'                av_%s_%s[iFrame][j] = av_%s_%s[iFrame][j]/tot_vol_%s\n',UserSettings.postproc.vars_int_pts{ivar},phase_name,UserSettings.postproc.vars_int_pts{ivar},phase_name,phase_name);
                fprintf(fid,'        #\n');
                fprintf(fid,'        # Finished computing average over the %s for this variable!\n',phase_name);
            end
        end
    end
    % Then consider quantities at the WHOLE ELEMENT:
    for ivar = 1:length(UserSettings.postproc.vars_whole_els)
        fprintf(fid,'        # Variable: %s\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'        %s_Field = RVEframe.fieldOutputs[''%s'']\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'        %s_SubField = %s_Field.getSubset(region=entireRVE_elSet, position=WHOLE_ELEMENT)\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        if isinstance(%s_SubField.values[0].data,float):\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'            # Then variable is a scalar:\n');
        fprintf(fid,'            # Loop over every element to compute average\n');
        fprintf(fid,'            for i in range(0, len(%s_SubField.values)):\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'                av_%s[iFrame] = av_%s[iFrame] + %s_SubField.values[i].data\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'            #\n');
        fprintf(fid,'            av_%s[iFrame] = av_%s[iFrame]/tot_vol\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'        else:\n');
        fprintf(fid,'            # Variable is an array:\n');
        fprintf(fid,'            # Loop over every element to compute average\n');
        fprintf(fid,'            for j in range(0, len(%s_SubField.values[0].data)):\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'                for i in range(0, len(%s_SubField.values)):\n',UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'                    av_%s[iFrame][j] = av_%s[iFrame][j] + %s_SubField.values[i].data[j]\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'                #\n');
        fprintf(fid,'                av_%s[iFrame][j] = av_%s[iFrame][j]/tot_vol\n',UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar},UserSettings.postproc.vars_whole_els{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        # Finished computing average for this variable!\n');
        %
        if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
            for jmat_phase=1:nmat_phases
                %
                phase_name = strcat('Phase_',num2str(jmat_phase));
                %
                fprintf(fid,'        %s_SubField_%s = %s_Field.getSubset(region=%s_elSet, position=WHOLE_ELEMENT)\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'        #\n');
                fprintf(fid,'        if isinstance(%s_SubField_%s.values[0].data,float):\n',UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'            # Then variable is a scalar:\n');
                fprintf(fid,'            # Loop over every element to compute average\n');
                fprintf(fid,'            for i in range(0, len(%s_SubField_%s.values)):\n',UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'                av_%s_%s[iFrame] = av_%s_%s[iFrame] + %s_SubField_%s.values[i].data\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'            #\n');
                fprintf(fid,'            av_%s_%s[iFrame] = av_%s_%s[iFrame]/tot_vol_%s\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name,phase_name);
                fprintf(fid,'        else:\n');
                fprintf(fid,'            # Variable is an array:\n');
                fprintf(fid,'            # Loop over every element to compute average\n');
                fprintf(fid,'            for j in range(0, len(%s_SubField_%s.values[0].data)):\n',UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'                for i in range(0, len(%s_SubField_%s.values)):\n',UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'                    av_%s_%s[iFrame][j] = av_%s_%s[iFrame][j] + %s_SubField_%s.values[i].data[j]\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name);
                fprintf(fid,'                #\n');
                fprintf(fid,'                av_%s_%s[iFrame][j] = av_%s_%s[iFrame][j]/tot_vol_%s\n',UserSettings.postproc.vars_whole_els{ivar},phase_name,UserSettings.postproc.vars_whole_els{ivar},phase_name,phase_name);
                fprintf(fid,'        #\n');
                fprintf(fid,'        # Finished computing average over the %s for this variable!\n',phase_name);
            end
        end
    end
    % Then consider quantities at the NODES:
    for ivar = 1:length(UserSettings.postproc.vars_nodes)
        fprintf(fid,'        # Variable: %s\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'        %s_Field = RVEframe.fieldOutputs[''%s'']\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'        %s_SubField = %s_Field.getSubset(region=entireRVE_nSet, position=NODAL)\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        if isinstance(%s_SubField.values[0].data,float):\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'            # Then variable is a scalar:\n');
        fprintf(fid,'            # Loop over every element to compute average\n');
        fprintf(fid,'            for i in range(0, len(%s_SubField.values)):\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'                av_%s[iFrame] = av_%s[iFrame] + %s_SubField.values[i].data\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'            #\n');
        fprintf(fid,'            av_%s[iFrame] = av_%s[iFrame]/len(%s_SubField.values)\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'        else:\n');
        fprintf(fid,'            # Variable is an array:\n');
        fprintf(fid,'            # Loop over every element to compute average\n');
        fprintf(fid,'            for j in range(0, len(%s_SubField.values[0].data)):\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'                for i in range(0, len(%s_SubField.values)):\n',UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'                    av_%s[iFrame][j] = av_%s[iFrame][j] + %s_SubField.values[i].data[j]\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'                #\n');
        fprintf(fid,'                av_%s[iFrame][j] = av_%s[iFrame][j]/len(%s_SubField.values)\n',UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar},UserSettings.postproc.vars_nodes{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        # Finished computing average for this variable!\n');
        %
        if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
            for jmat_phase=1:nmat_phases
                %
                phase_name = strcat('Phase_',num2str(jmat_phase));
                %
                fprintf(fid,'        %s_SubField_%s = %s_Field.getSubset(region=%s_nSet, position=NODAL)\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'        #\n');
                fprintf(fid,'        if isinstance(%s_SubField_%s.values[0].data,float):\n',UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'            # Then variable is a scalar:\n');
                fprintf(fid,'            # Loop over every element to compute average\n');
                fprintf(fid,'            for i in range(0, len(%s_SubField_%s.values)):\n',UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'                av_%s_%s[iFrame] = av_%s_%s[iFrame] + %s_SubField_%s.values[i].data\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'            #\n');
                fprintf(fid,'            av_%s_%s[iFrame] = av_%s_%s[iFrame]/len(%s_SubField_%s.values)\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'        else:\n');
                fprintf(fid,'            # Variable is an array:\n');
                fprintf(fid,'            # Loop over every element to compute average\n');
                fprintf(fid,'            for j in range(0, len(%s_SubField_%s.values[0].data)):\n',UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'                for i in range(0, len(%s_SubField_%s.values)):\n',UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'                    av_%s_%s[iFrame][j] = av_%s_%s[iFrame][j] + %s_SubField_%s.values[i].data[j]\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'                #\n');
                fprintf(fid,'                av_%s_%s[iFrame][j] = av_%s_%s[iFrame][j]/len(%s_SubField_%s.values)\n',UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name,UserSettings.postproc.vars_nodes{ivar},phase_name);
                fprintf(fid,'        #\n');
                fprintf(fid,'        # Finished computing average over the %s for this variable!\n',phase_name);
            end
        end
    end
    %
end % end of NOT a steady state simulation
%
%
% And finally consider quantities that lead to nominal stress P and deformation gradient F (if using PBCs):
if UserSettings.gen_options.pbcs_option == 1 % Then we are using Periodic Boundary Conditions (PBCs)
    for ivar = 1:length(postproc_variables_pbcs)
        fprintf(fid,'        # Variable: %s\n',postproc_variables_pbcs{ivar});
        fprintf(fid,'        %s_Field = RVEframe.fieldOutputs[''%s'']\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'        %s_dummy1_SubField = %s_Field.getSubset(region=dummy1_nSet, position=NODAL)\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'        %s_dummy2_SubField = %s_Field.getSubset(region=dummy2_nSet, position=NODAL)\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'        #\n');
        fprintf(fid,'        if isinstance(%s_dummy1_SubField.values[0].data,float):\n',postproc_variables_pbcs{ivar});
        fprintf(fid,'            # Then variable is a scalar:\n');
        fprintf(fid,'            %s_dummy1[iFrame] = %s_dummy1_SubField.values[0].data\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'            %s_dummy2[iFrame] = %s_dummy2_SubField.values[0].data\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
            fprintf(fid,'            %s_dummy1_imag[iFrame] = %s_dummy1_SubField.values[0].conjugateData\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
            fprintf(fid,'            %s_dummy2_imag[iFrame] = %s_dummy2_SubField.values[0].conjugateData\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        end
        fprintf(fid,'            #\n');
        fprintf(fid,'        else:\n');
        fprintf(fid,'            # Variable is an array:\n');
        fprintf(fid,'            for j in range(0, len(%s_dummy1_SubField.values[0].data)):\n',postproc_variables_pbcs{ivar});
        fprintf(fid,'                %s_dummy1[iFrame][j] = %s_dummy1_SubField.values[0].data[j]\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        fprintf(fid,'                %s_dummy2[iFrame][j] = %s_dummy2_SubField.values[0].data[j]\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
            fprintf(fid,'                %s_dummy1_imag[iFrame][j] = %s_dummy1_SubField.values[0].conjugateData[j]\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
            fprintf(fid,'                %s_dummy2_imag[iFrame][j] = %s_dummy2_SubField.values[0].conjugateData[j]\n',postproc_variables_pbcs{ivar},postproc_variables_pbcs{ivar});
        end
        fprintf(fid,'                #\n');
        fprintf(fid,'        #\n');
        fprintf(fid,'        # Finished saving this variable at the dummy nodes!\n');
    end
    % Now compute the deformation gradient:
    fprintf(fid,'        # Now compute the deformation gradient, jacobian and nominal stress:\n');
    fprintf(fid,'        for j in range(0,2):\n');
    fprintf(fid,'            defGrad[iFrame][0][j] = U_dummy1[iFrame][j] + eye[0][j]\n');
    fprintf(fid,'            defGrad[iFrame][1][j] = U_dummy2[iFrame][j] + eye[1][j]\n');
    % Important note: The stress in 2D is Force/length. However, I confirmed that we need to
    %divide the Reaction force at the dummy node by the Length associated to that dummy node
    %because we provided a displacement as the STRAIN value (not a real displacement), i.e.
    %the constraint equations take the input displacement we assign to the dummy node and
    %multiply it by the Length of the RVE along that direction.
    % In conclusion: RF_dummyi_j/L_i/L_j => RF_dummyi_j/area
    % Additional note: in 3D it is the same thing (divide by volume)
    fprintf(fid,'            nomP[iFrame][0][j] = RF_dummy1[iFrame][j]/tot_vol\n');
    fprintf(fid,'            nomP[iFrame][1][j] = RF_dummy2[iFrame][j]/tot_vol\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'            # Since this is a steady state dynamics simulation, compute the imaginary parts:\n');
        fprintf(fid,'            nomP_imag[iFrame][0][j] = RF_dummy1_imag[iFrame][j]/tot_vol\n');
        fprintf(fid,'            nomP_imag[iFrame][1][j] = RF_dummy2_imag[iFrame][j]/tot_vol\n');
        fprintf(fid,'        # Compute the Young''s moduli in the 1 direction: \n');
        fprintf(fid,'        if U_dummy1[iFrame][0] == 0.0:\n');
        fprintf(fid,'            YoungsE1_real[iFrame] = 0.0 \n');
        fprintf(fid,'            YoungsE1_imag[iFrame] = 0.0 \n');
        fprintf(fid,'        else:\n');
        fprintf(fid,'            YoungsE1_real[iFrame] = nomP[iFrame][0][0] / U_dummy1[iFrame][0] \n');
        fprintf(fid,'            YoungsE1_imag[iFrame] = nomP_imag[iFrame][0][0] / U_dummy1[iFrame][0] \n');
    end
    fprintf(fid,'        jacobian[iFrame] = numpy.linalg.det(defGrad[iFrame][:][:])\n');
    % If using Green-strain:
    if UserSettings.gen_options.strain_option == 1 % Then using we are using Green-strain:
        fprintf(fid,'        Green_strain[iFrame,:,:] = 0.5*( numpy.dot(numpy.transpose(defGrad[iFrame,:,:]),defGrad[iFrame,:,:])-numpy.identity(2) )\n');
        fprintf(fid,'        PK2[iFrame,:,:] = numpy.dot( nomP[iFrame,:,:],numpy.linalg.inv(numpy.transpose(defGrad[iFrame,:,:])) )\n');
    end
    %
    %
else % NOT using PBCs
    % So we will need the reaction force on the left edge (RF_leftEdge)
    fprintf(fid,'        RF_leftEdge_Field = RVEframe.fieldOutputs[''RF'']\n');
    fprintf(fid,'        RF_leftEdge_SubField = RF_leftEdge_Field.getSubset(region=leftEdge_nSet, position=NODAL)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        if isinstance(RF_leftEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'            # Then variable is a scalar:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for i in range(0, len(RF_leftEdge_SubField.values)):\n');
    fprintf(fid,'                sum_RF_leftEdge[iFrame] = sum_RF_leftEdge[iFrame] + RF_leftEdge_SubField.values[i].data\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'                sum_imagRF_leftEdge[iFrame] = sum_imagRF_leftEdge[iFrame] + RF_leftEdge_SubField.values[i].conjugateData\n');
    end
    fprintf(fid,'            #\n');
    %     fprintf(fid,'            sum_RF_leftEdge[iFrame] = sum_RF_leftEdge[iFrame]/len(RF_leftEdge_SubField.values)\n');
    fprintf(fid,'        else:\n');
    fprintf(fid,'            # Variable is an array:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for j in range(0, len(RF_leftEdge_SubField.values[0].data)):\n');
    fprintf(fid,'                for i in range(0, len(RF_leftEdge_SubField.values)):\n');
    fprintf(fid,'                    sum_RF_leftEdge[iFrame][j] = sum_RF_leftEdge[iFrame][j] + RF_leftEdge_SubField.values[i].data[j]\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'                    sum_imagRF_leftEdge[iFrame][j] = sum_imagRF_leftEdge[iFrame][j] + RF_leftEdge_SubField.values[i].conjugateData[j]\n');
    end
    fprintf(fid,'                #\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        # Finished computing sum for this variable!\n');
    % And RF on bottom edge (RF_leftEdge)
    fprintf(fid,'        RF_botEdge_Field = RVEframe.fieldOutputs[''RF'']\n');
    fprintf(fid,'        RF_botEdge_SubField = RF_botEdge_Field.getSubset(region=botEdge_nSet, position=NODAL)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        if isinstance(RF_botEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'            # Then variable is a scalar:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for i in range(0, len(RF_botEdge_SubField.values)):\n');
    fprintf(fid,'                sum_RF_botEdge[iFrame] = sum_RF_botEdge[iFrame] + RF_botEdge_SubField.values[i].data\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'                sum_imagRF_botEdge[iFrame] = sum_imagRF_botEdge[iFrame] + RF_botEdge_SubField.values[i].data\n');
    end
    fprintf(fid,'            #\n');
    %     fprintf(fid,'            sum_RF_botEdge[iFrame] = sum_RF_botEdge[iFrame]/len(RF_botEdge_SubField.values)\n');
    fprintf(fid,'        else:\n');
    fprintf(fid,'            # Variable is an array:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for j in range(0, len(RF_botEdge_SubField.values[0].data)):\n');
    fprintf(fid,'                for i in range(0, len(RF_botEdge_SubField.values)):\n');
    fprintf(fid,'                    sum_RF_botEdge[iFrame][j] = sum_RF_botEdge[iFrame][j] + RF_botEdge_SubField.values[i].data[j]\n');
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'                    sum_imagRF_botEdge[iFrame][j] = sum_imagRF_botEdge[iFrame][j] + RF_botEdge_SubField.values[i].conjugateData[j]\n');
    end
    fprintf(fid,'                #\n');
    %     fprintf(fid,'                sum_RF_botEdge[iFrame][j] = sum_RF_botEdge[iFrame][j]/len(RF_botEdge_SubField.values)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        # Finished computing sum for this variable!\n');
    %
    % And the displacement U on the right edge (U_rightEdge)
    fprintf(fid,'        U_rightEdge_Field = RVEframe.fieldOutputs[''U'']\n');
    fprintf(fid,'        U_rightEdge_SubField = U_rightEdge_Field.getSubset(region=rightEdge_nSet, position=NODAL)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        if isinstance(U_rightEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'            # Then variable is a scalar:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for i in range(0, len(U_rightEdge_SubField.values)):\n');
    fprintf(fid,'                av_U_rightEdge[iFrame] = av_U_rightEdge[iFrame] + U_rightEdge_SubField.values[i].data\n');
    fprintf(fid,'            #\n');
    fprintf(fid,'            av_U_rightEdge[iFrame] = av_U_rightEdge[iFrame]/len(U_rightEdge_SubField.values)\n');
    fprintf(fid,'        else:\n');
    fprintf(fid,'            # Variable is an array:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for j in range(0, len(U_rightEdge_SubField.values[0].data)):\n');
    fprintf(fid,'                for i in range(0, len(U_rightEdge_SubField.values)):\n');
    fprintf(fid,'                    av_U_rightEdge[iFrame][j] = av_U_rightEdge[iFrame][j] + U_rightEdge_SubField.values[i].data[j]\n');
    fprintf(fid,'                #\n');
    fprintf(fid,'                av_U_rightEdge[iFrame][j] = av_U_rightEdge[iFrame][j]/len(U_rightEdge_SubField.values)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        # Finished computing average for this variable!\n');
    %
    % And the displacement on the top edge (U_topEdge)
    fprintf(fid,'        U_topEdge_Field = RVEframe.fieldOutputs[''U'']\n');
    fprintf(fid,'        U_topEdge_SubField = U_topEdge_Field.getSubset(region=topEdge_nSet, position=NODAL)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        if isinstance(U_topEdge_SubField.values[0].data,float):\n');
    fprintf(fid,'            # Then variable is a scalar:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for i in range(0, len(U_topEdge_SubField.values)):\n');
    fprintf(fid,'                av_U_topEdge[iFrame] = av_U_topEdge[iFrame] + U_topEdge_SubField.values[i].data\n');
    fprintf(fid,'            #\n');
    fprintf(fid,'            av_U_topEdge[iFrame] = av_U_topEdge[iFrame]/len(U_topEdge_SubField.values)\n');
    fprintf(fid,'        else:\n');
    fprintf(fid,'            # Variable is an array:\n');
    fprintf(fid,'            # Loop over every element to compute average\n');
    fprintf(fid,'            for j in range(0, len(U_topEdge_SubField.values[0].data)):\n');
    fprintf(fid,'                for i in range(0, len(U_topEdge_SubField.values)):\n');
    fprintf(fid,'                    av_U_topEdge[iFrame][j] = av_U_topEdge[iFrame][j] + U_topEdge_SubField.values[i].data[j]\n');
    fprintf(fid,'                #\n');
    fprintf(fid,'                av_U_topEdge[iFrame][j] = av_U_topEdge[iFrame][j]/len(U_topEdge_SubField.values)\n');
    fprintf(fid,'        #\n');
    fprintf(fid,'        # Finished computing average for this variable!\n');
    %
    %
    
    % Now compute the deformation gradient:
    fprintf(fid,'        # Now compute the principal average strains (applied to the RVE):\n');
    fprintf(fid,'        epsilon1 = av_U_rightEdge[iFrame][0] / Lx \n');
    fprintf(fid,'        epsilon2 = av_U_topEdge[iFrame][1] / Ly \n');
    fprintf(fid,'        # Compute average strains (applied to the RVE) in original frame:\n');
    fprintf(fid,'        alpha = RVEs_data_old[stringDiscInput][stringDoE][''alpha_NoPBCs''] \n');
    fprintf(fid,'        epsilon = numpy.zeros(( 2,2 ))\n');
    fprintf(fid,'        epsilon[0][0] = epsilon1*numpy.power(cos(-alpha),2) + epsilon2*numpy.power(sin(-alpha),2) \n'); % epsilon(1,1)
    fprintf(fid,'        epsilon[1][1] = epsilon1*numpy.power(sin(-alpha),2) + epsilon2*numpy.power(cos(-alpha),2) \n'); % epsilon(2,2)
    fprintf(fid,'        epsilon[0][1] = (epsilon2-epsilon1)*sin(-alpha)*cos(-alpha) \n'); % epsilon(1,2)
    fprintf(fid,'        epsilon[1][0] = epsilon[0][1] \n'); % epsilon(2,1)
    fprintf(fid,'        # Now compute the deformation gradient, jacobian and nominal stress:\n');
    
    fprintf(fid,'        for j in range(0,2):\n');
    fprintf(fid,'            for k in range(0,2):\n');
    fprintf(fid,'                defGrad[iFrame][j][k] = epsilon[j][k] + eye[j][k]\n');
    
    %     fprintf(fid,'        defGrad[iFrame][:][:] = epsilon[:][:] + eye[0:1][0:1]\n');
    fprintf(fid,'        #\n');
    % Now compute nominal stress in the original frame:
    fprintf(fid,'        # Now compute the average nominal stresses in the rotated frame:\n');
    fprintf(fid,'        nomP1 = -sum_RF_leftEdge[iFrame][0] / Ly \n');
    fprintf(fid,'        nomP2 = -sum_RF_botEdge[iFrame][1] / Lx \n');
    fprintf(fid,'        nomP[iFrame][0][0] = nomP1*numpy.power(cos(-alpha),2) + nomP2*numpy.power(sin(-alpha),2) \n'); % epsilon(1,1)
    fprintf(fid,'        nomP[iFrame][1][1] = nomP1*numpy.power(sin(-alpha),2) + nomP2*numpy.power(cos(-alpha),2) \n'); % epsilon(2,2)
    fprintf(fid,'        nomP[iFrame][0][1] = (nomP2-nomP1)*sin(-alpha)*cos(-alpha) \n'); % epsilon(1,2)
    fprintf(fid,'        nomP[iFrame][1][0] = nomP[iFrame][0][1] \n'); % epsilon(2,1)
    %
    if UserSettings.gen_options.analysis_option == 2 % This is a steady state dynamics simulation
        fprintf(fid,'        # Since this is a steady state dynamics simulation, compute the imaginary parts:\n');
        fprintf(fid,'        nomP1_imag = -sum_imagRF_leftEdge[iFrame][0] / Ly \n');
        fprintf(fid,'        nomP2_imag = -sum_imagRF_botEdge[iFrame][1] / Lx \n');
        fprintf(fid,'        nomP_imag[iFrame][0][0] = nomP1_imag*numpy.power(cos(-alpha),2) + nomP2_imag*numpy.power(sin(-alpha),2) \n'); % epsilon(1,1)
        fprintf(fid,'        nomP_imag[iFrame][1][1] = nomP1_imag*numpy.power(sin(-alpha),2) + nomP2_imag*numpy.power(cos(-alpha),2) \n'); % epsilon(2,2)
        fprintf(fid,'        nomP_imag[iFrame][0][1] = (nomP2_imag-nomP1_imag)*sin(-alpha)*cos(-alpha) \n'); % epsilon(1,2)
        fprintf(fid,'        nomP_imag[iFrame][1][0] = nomP_imag[iFrame][0][1] \n'); % epsilon(2,1)
        fprintf(fid,'        # Compute the Young''s moduli in the 1 direction: \n');
        fprintf(fid,'        YoungsE1_real[iFrame] = nomP[iFrame][0][0] / epsilon[0][0] \n');
        fprintf(fid,'        YoungsE1_imag[iFrame] = nomP_imag[iFrame][0][0] / epsilon[0][0] \n');
    end
    %     sigma(1,1) = sigma1*cos(-alpha)^2 + sigma2*sin(-alpha)^2;
    %     sigma(2,2) = sigma1*sin(-alpha)^2 + sigma2*cos(-alpha)^2;
    %     sigma(1,2) = (sigma2 - sigma1)*sin(-alpha)*cos(-alpha);
    %     sigma(2,1) = sigma(1,2);
    fprintf(fid,'        # Finally compute the Jacobian:\n');
    fprintf(fid,'        jacobian[iFrame] = numpy.linalg.det(defGrad[iFrame][:][:])\n');
    %
    % If using Green-strain:
    if UserSettings.gen_options.strain_option == 1 % Then using we are using Green-strain:
        fprintf(fid,'        Green_strain[iFrame,:,:] = 0.5*( numpy.dot(numpy.transpose(defGrad[iFrame,:,:]),defGrad[iFrame,:,:])-numpy.identity(2) )\n');
        fprintf(fid,'        PK2[iFrame,:,:] = numpy.dot( nomP[iFrame,:,:],numpy.linalg.inv(numpy.transpose(defGrad[iFrame,:,:])) )\n');
    end
    %
end
%
fprintf(fid,'\n');
fprintf(fid,'#\n');
%
fprintf(fid,'# Save all variables to a single structured variable with all the data\n');
dictionary_string = '';
% First consider quantities at the INTEGRATION POINTS:
for ivar = 1:length(UserSettings.postproc.vars_int_pts)
    if ivar == 1
        add_string = [ '''',UserSettings.postproc.vars_int_pts{ivar},''':av_',UserSettings.postproc.vars_int_pts{ivar} ];
    else
        add_string = [ ', ','''',UserSettings.postproc.vars_int_pts{ivar},''':av_',UserSettings.postproc.vars_int_pts{ivar} ];
    end
    dictionary_string = strcat(dictionary_string,add_string);
end
% Then consider quantities at the WHOLE ELEMENT:
for ivar = 1:length(UserSettings.postproc.vars_whole_els)
    if ivar == 1 && isempty(UserSettings.postproc.vars_int_pts) == 1
        add_string = [ '''',UserSettings.postproc.vars_whole_els{ivar},''':av_',UserSettings.postproc.vars_whole_els{ivar} ];
    else
        add_string = [ ', ','''',UserSettings.postproc.vars_whole_els{ivar},''':av_',UserSettings.postproc.vars_whole_els{ivar} ];
    end
    dictionary_string = strcat(dictionary_string,add_string);
end
% Then consider quantities at the NODES:
for ivar = 1:length(UserSettings.postproc.vars_nodes)
    if ivar == 1 && isempty(UserSettings.postproc.vars_int_pts) == 1 ...
            && isempty(UserSettings.postproc.vars_whole_els) == 1
        add_string = [ '''',UserSettings.postproc.vars_nodes{ivar},''':av_',UserSettings.postproc.vars_nodes{ivar} ];
    else
        add_string = [ ', ','''',UserSettings.postproc.vars_nodes{ivar},''':av_',UserSettings.postproc.vars_nodes{ivar} ];
    end
    dictionary_string = strcat(dictionary_string,add_string);
end
%
% Do the same thing for the local averages (if asked):
if UserSettings.postproc.local_averages == 1 % In case we want separate the average behavior at each material phase
    for jmat_phase=1:nmat_phases
        %
        phase_name = strcat('Phase_',num2str(jmat_phase));
        %
        % First consider quantities at the INTEGRATION POINTS:
        for ivar = 1:length(UserSettings.postproc.vars_int_pts)
            add_string = [ ', ','''',UserSettings.postproc.vars_int_pts{ivar},'_',phase_name,''':av_',UserSettings.postproc.vars_int_pts{ivar},'_',phase_name ];
            dictionary_string = strcat(dictionary_string,add_string);
        end
        % Then consider quantities at the WHOLE ELEMENT:
        for ivar = 1:length(UserSettings.postproc.vars_whole_els)
            add_string = [ ', ','''',UserSettings.postproc.vars_whole_els{ivar},'_',phase_name,''':av_',UserSettings.postproc.vars_whole_els{ivar},'_',phase_name ];
            dictionary_string = strcat(dictionary_string,add_string);
        end
        % Then consider quantities at the NODES:
        for ivar = 1:length(UserSettings.postproc.vars_nodes)
            add_string = [ ', ','''',UserSettings.postproc.vars_nodes{ivar},'_',phase_name,''':av_',UserSettings.postproc.vars_nodes{ivar},'_',phase_name ];
            dictionary_string = strcat(dictionary_string,add_string);
        end
    end
end % no local averages requested
%
% Now take care of the variables that are always outputted, independent of user choice:
permanent_vars = '''step_total_time'':stepTotalTime, ''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol';
%
% If using Green-strain:
if UserSettings.gen_options.strain_option == 1 % Then using we are using Green-strain:
    permanent_vars = strcat(permanent_vars,', ''Green_strain'':Green_strain, ''PK2'':PK2');
end
%
%
if UserSettings.postproc.local_averages == 1
    for jmat_phase=1:nmat_phases
        %
        phase_name = strcat('Phase_',num2str(jmat_phase));
        %
        permanent_vars = strcat(permanent_vars,', ''total_vol_',phase_name,''':tot_vol_',phase_name);
%         if interphase_material == 1 % interphase exists:
%             permanent_vars = strcat(permanent_vars,', ''total_vol_IP'':tot_vol_IP');
%         end
    end
end
%
if UserSettings.gen_options.analysis_option == 2 % This is a steady state simulation
    permanent_vars = strcat(permanent_vars,', ''P_imag'':nomP_imag, ''YoungsE1_real'':YoungsE1_real, ''YoungsE1_imag'':YoungsE1_imag, ''frequency'':frequency');
end
%
if (length(UserSettings.postproc.vars_int_pts)+length(UserSettings.postproc.vars_whole_els)+length(UserSettings.postproc.vars_nodes)) > 0
    % If user asked for some post-processing variables then:
    fprintf(fid,'RVE_variables = {%s, %s}\n',dictionary_string,permanent_vars);
    %     if UserSettings.gen_options.pbcs_option == 1
    %         if UserSettings.postproc.local_averages == 1
    %             if interphase_material ~= 0
    %                 fprintf(fid,'RVE_variables = {%s, ''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix, ''total_vol_particle'':tot_vol_particle, ''total_vol_IP'':tot_vol_IP, ''step_total_time'':stepTotalTime}\n',dictionary_string);
    %             else % no interphase:
    %                 fprintf(fid,'RVE_variables = {%s, ''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix, ''total_vol_particle'':tot_vol_particle, ''step_total_time'':stepTotalTime}\n',dictionary_string);
    %             end
    %         else % no local averages
    %             fprintf(fid,'RVE_variables = {%s, ''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol, ''step_total_time'':stepTotalTime}\n',dictionary_string);
    %         end
    %     else % no pbcs
    %         if UserSettings.postproc.local_averages == 1
    %             if interphase_material ~= 0
    %                 fprintf(fid,'RVE_variables = {%s, ''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix, ''total_vol_particle'':tot_vol_particle, ''total_vol_IP'':tot_vol_IP, ''step_total_time'':stepTotalTime}\n',dictionary_string);
    %             else % no interphase:
    %                 fprintf(fid,'RVE_variables = {%s, ''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix,  ''total_vol_particle'':tot_vol_particle, ''step_total_time'':stepTotalTime}\n',dictionary_string);
    %             end
    %         else % no local averages
    %             fprintf(fid,'RVE_variables = {%s, ''total_vol'':tot_vol, ''step_total_time'':stepTotalTime}\n',dictionary_string);
    %         end
    %     end
else % No post-processing variables were asked by the user!
    fprintf(fid,'RVE_variables = {%s}\n',permanent_vars);
    %     if UserSettings.gen_options.pbcs_option == 1
    %         if UserSettings.postproc.local_averages == 1
    %             if interphase_material ~= 0
    %                 fprintf(fid,'RVE_variables = {''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix, ''total_vol_particle'':tot_vol_particle, ''total_vol_IP'':tot_vol_IP, ''step_total_time'':stepTotalTime}\n');
    %             else % no interphase:
    %                 fprintf(fid,'RVE_variables = {''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix,  ''total_vol_particle'':tot_vol_particle, ''step_total_time'':stepTotalTime}\n');
    %             end
    %         else % no local averages
    %             fprintf(fid,'RVE_variables = {''P'':nomP, ''F'':defGrad, ''J'':jacobian, ''total_vol'':tot_vol, ''step_total_time'':stepTotalTime}\n');
    %         end
    %     else % no pbcs
    %         if UserSettings.postproc.local_averages == 1
    %             if interphase_material ~= 0
    %                 fprintf(fid,'RVE_variables = {''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix, ''total_vol_particle'':tot_vol_particle, ''total_vol_IP'':tot_vol_IP, ''step_total_time'':stepTotalTime}\n');
    %             else % no interphase:
    %                 fprintf(fid,'RVE_variables = {''total_vol'':tot_vol, ''total_vol_matrix'':tot_vol_matrix,  ''total_vol_particle'':tot_vol_particle, ''step_total_time'':stepTotalTime}\n');
    %             end
    %         else % no local averages
    %             fprintf(fid,'RVE_variables = {''total_vol'':tot_vol, ''step_total_time'':stepTotalTime}\n');
    %         end
    %     end
end
fprintf(fid,'#\n');
fprintf(fid,'RVEs_data_update = {stringDiscInput : {stringDoE: RVE_variables} }\n');
fprintf(fid,'if postproc_exists == 1:\n');
% fprintf(fid,'    RVEs_data = merge_dicts(RVEs_data_old, RVEs_data_update)\n');
fprintf(fid,'    RVEs_data = dict_merge(RVEs_data_old, RVEs_data_update)\n');
% fprintf(fid,'    RVEs_data = RVEs_data_old\n');
% fprintf(fid,'    RVEs_data.update(RVEs_data_update)\n');
fprintf(fid,'else:\n');
fprintf(fid,'   RVEs_data = RVEs_data_update\n');
fprintf(fid,'\n');
fprintf(fid,'# Save post-processing information to pickle file:\n');
fprintf(fid,'writefile_postproc = open(file_postproc_path, ''w'')\n');
fprintf(fid,'pickle.dump(RVEs_data, writefile_postproc)\n');
fprintf(fid,'writefile_postproc.close()\n');
fprintf(fid,'\n');
fprintf(fid,'# End of file.\n');
%
% NOTE: To access the list "E" of the RVE1 and BC1: [x['E'] for x in RVEs_data['RVE1']['BC1']]
%
%
fclose(fid);
%
% disp(' ');
% disp('Creation of Post-Processing Python file COMPLETED');
% disp('Elapsed Time [min]: ');
% disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Python Script in ABAQUS CAE                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp(['Started Abaqus odb post-processing for DoE point ',int2str(jDoE),' of Discrete Input point ',int2str(iDiscInp)]);
disp(' ');
system([UserSettings.analyses.abaqus_path,' cae noGUI=',postproc_scriptfile_name]);
% disp(' ');
% disp('RVEs_postprocessing_variables.p file UPDATED for this RVE and B.C.');
% disp('Elapsed Time [min]: ');
% disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
