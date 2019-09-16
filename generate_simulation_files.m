%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Generate the input files for the FEM simulation                        %
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
function generate_simulation_files(DoE_point,DiscInput_point,UserSettings,jDoE,iDiscInp,jDoE_dir,postproc_dir,abq_inp_file_noExt,run_file,CPUid)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Lx = UserSettings.analyses.RVE_dimensions(1);
Ly = UserSettings.analyses.RVE_dimensions(2);

nmat_phases = UserSettings.analyses.Number_Mat_Phases;

% Vector with the IDs of the "nmaterials" material phases
mat_phases_IDs = [DiscInput_point{1}];
% 
% % Cell with the strings (names) of the files with the properties for each
% % material phase.
% mat_phases_input_files = {DiscInput_points{2}};
%
applied_strain = zeros(1,3);
applied_strain(1) = UserSettings.analyses.FixedInput_values{1};
applied_strain(2) = UserSettings.analyses.FixedInput_values{2};
applied_strain(3) = UserSettings.analyses.FixedInput_values{3};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert Boundary Conditions

%         alpha = zeros(max(selectBCs_jRVE),1); % transformation angle used for NON periodic boundary conditions
%         for iBC=selectBCs_jRVE
epsilon = zeros(2,2); % initialize small strain for this iBC
% epsilon1 = 0.0; epsilon2 = 0.0; % principal strains used for NON periodic boundary conditions
% Calculate the deformation gradient to apply the Boundary Conditions
% according to the strain measure considered in the microstructure file:
if UserSettings.gen_options.strain_option == 0 % Strain provided is equal to the small strain:
    epsilon(1,1) = applied_strain(1);
    epsilon(2,2) = applied_strain(2);
    epsilon(1,2) = applied_strain(3);
    epsilon(2,1) = epsilon(1,2);
elseif UserSettings.gen_options.strain_option == 1 % Strain provided is Green-Lagrange strain
    E11 = applied_strain(1);
    E22 = applied_strain(2);
    E12 = applied_strain(3);
%     Sin_theta = input_BC(iBC,3); % In case the third component is the sin of angle theta:
%     E12 = 1/2*sqrt(2*E11+1)*sqrt(2*E22+1)*Sin_theta;
    %
    EGreen = [E11,E12; E12,E22];
    % Deformation gradient
    F = sqrtm(2*EGreen+eye(2));
    epsilon = 1/2*(F+F')-eye(2); % small strain to be applied to the RVE
end
%
% If not using periodic boundary conditions, then the easiest way is to
% convert to principal strains space (BE CAREFUL IF THE RVE IS ANISOTROPIC!!!!)
if UserSettings.gen_options.pbcs_option ~= 1 % NOT using PBCs
    % In 2D:
%                 I1 = trace(epsilon);
%                 I2 = epsilon(1,1)*epsilon(2,2) - epsilon(1,2)^2;
%                 I3 = 0.0;
    %
%                 phi = 1/3*acos((2*I1^3-9*I1*I2+27*I3)/(2*(I1^2-3*I2)^(3/2)));
%                 sqrt_term = sqrt(((epsilon(1,1)-epsilon(2,2))/2+epsilon(1,2)^2)^2);
%                 epsilon1 = (epsilon(1,1)+epsilon(1,2))/2 + sqrt_term;
%                 epsilon2 = (epsilon(1,1)+epsilon(1,2))/2 - sqrt_term;
%                 epsilon2 = I1/3+2/3*sqrt(I1^2-3*I2)*cos(phi+2*pi()/3);
%                 epsilon2 = I1/3+2/3*sqrt(I1^2-3*I2)*cos(phi+4*pi()/3);
    %
    % Compute the transformation angle of the principal strains
    if epsilon(1,2) == 0
        alpha = 0;
        %
        epsilon1 = epsilon(1,1);
        epsilon2 = epsilon(2,2);
    else
        alpha = 1/2*atan( epsilon(1,2)/((epsilon(1,1)-epsilon(2,2))/2) );
        %
        epsilon1 = epsilon(1,1)*cos(alpha)^2 + epsilon(2,2)*sin(alpha)^2 + 2*epsilon(1,2)*sin(alpha)*cos(alpha);
        epsilon2 = epsilon(1,1)*sin(alpha)^2 + epsilon(2,2)*cos(alpha)^2 - 2*epsilon(1,2)*sin(alpha)*cos(alpha);
        epsilon_shear = (epsilon(2,2)-epsilon(1,1))*sin(alpha)*cos(alpha) +epsilon(1,2)*(cos(alpha)^2-sin(alpha)^2);
        if abs(epsilon_shear) > 1e-5
            disp('Check rotation to principal strains frame... Shear is not zero! Hit any key to continue:'); pause
        end
    end
    %
    % Rotating to the principal angle leads to zero shear and
    % the principal strains:
    %
    % Note: if using non-periodic boundary conditions, the RVE is subjected
    % to the principal strains (no shear strain) and then the homogenized
    % stress is transformed back to the original frame (see generate_postproc_file.m)
    %
    % In 3D:
%                 I1 = trace(epsilon);
%                 I2 = epsilon(1,1)*epsilon(2,2)+epsilon(2,2)*epsilon(3,3)+epsilon(1,1)*epsilon(3,3) -...
%                     epsilon(1,2)^2-epsilon(2,3)^2-epsilon(3,1)^2;
%                 I3 = det(epsilon);
%                 %
%                 phi = 1/3*acos((2*I1^3-9*I1*I2+27*I3)/(2*(I1^2-3*I2)^(3/2)));
%                 epsilon1 = I1/3+2/3*sqrt(I1^2-3*I2)*cos(phi);
%                 epsilon2 = I1/3+2/3*sqrt(I1^2-3*I2)*cos(phi+2*pi()/3);
%                 epsilon3 = I1/3+2/3*sqrt(I1^2-3*I2)*cos(phi+4*pi()/3);
%
    % Just save the angle ALPHA for later use in postprocessing the results
    %
%     disp(' ');
%     disp('Started creation of Python file');
    %
    postproc_scriptfile_name = strcat(postproc_dir,'/script_write_alpha_NoPBCs_CPU',num2str(CPUid),'.py');
    fid = fopen(postproc_scriptfile_name,'wt');
    %
    fprintf(fid,'#=====================================================================#\n');
    fprintf(fid,'#\n');
    fprintf(fid,'# Created by M.A. Bessa on %s\n',datestr(now));
    fprintf(fid,'#=====================================================================#\n');
%     fprintf(fid,'from odbAccess import *\n');
%     fprintf(fid,'import os\n');
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
%     fprintf(fid,'# Set the work directory for ABAQUS (where the odb file is):\n');
%     fprintf(fid,'RVEdir=''%s''\n',jDoE_dir);
%     fprintf(fid,'\n');
    % fprintf(fid,'BCdir= RVEdir + ''/BC%i''\n',iBC);
    % fprintf(fid,'#\n');
    fprintf(fid,'# Set directory where you want the post-processing file to be generated\n');
    fprintf(fid,'postproc_dir=''%s''\n',postproc_dir);
    fprintf(fid,'#\n');
%     fprintf(fid,'# Name of the file job\n');
%     fprintf(fid,'jobname=''%s''\n',abq_inp_file_noExt);
%     fprintf(fid,'#\n');
%     for jmat_phase=1:nmat_phases
%         %
%         phase_name = strcat('Phase_',num2str(jmat_phase));
%         %
%         fprintf(fid,'# Define the name of the %s part\n',phase_name);
%         fprintf(fid,'%s_part=''FinalRVE-1.%s''\n',phase_name,phase_name);
%         fprintf(fid,'#\n');
%     end
%     %
%     fprintf(fid,'#\n');
%     fprintf(fid,'odbfile=jobname+''.odb'' # Define name of this .odb file\n');
%     fprintf(fid,'#\n');
%     fprintf(fid,'os.chdir(RVEdir)\n');
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
    fprintf(fid,'alpha = %6.5E \n',alpha);
    alpha_var = '''alpha_NoPBCs'':alpha';
    fprintf(fid,'RVE_variables = {%s}\n',alpha_var);
    fprintf(fid,'#\n');
    fprintf(fid,'stringDoE = ''DoE''+str(%i)\n',jDoE);
    fprintf(fid,'stringDiscInput = ''DiscInput''+str(%i)\n',iDiscInp);
    % fprintf(fid,'stringBC = ''BC''+str(%i)\n',iBC);
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
    % NOTE: RVEs_data['DiscInp1']['DoE1']['alpha_NoPBCs']
    %
    %
    fclose(fid);
    % Run this little python code:
    system(['python ',postproc_dir,'/script_write_alpha_NoPBCs_CPU',num2str(CPUid),'.py']);
end


% Write the simulations input files:
%
disp(' '); disp(' ');
disp('Started writing the main input file');
%
disp(' ');
string=strcat('Analysis of Input point ',num2str(iDiscInp,'%i'),' for DoE point ',num2str(jDoE,'%i'));
disp(string);
disp(' ');
%
%% Write lines that are common to every load case:
%
abq_inp_file_fullpath = strcat(jDoE_dir,'/',abq_inp_file_noExt,'.inp');
fid = fopen(abq_inp_file_fullpath,'wt');
%
% Heading
%
fprintf(fid,'** ***************************************************************************\n');
fprintf(fid,'**   Created by M.A. Bessa on %s\n',datestr(now));
fprintf(fid,'**\n');
fprintf(fid,'**   %s :\n',string);
fprintf(fid,'**   Applied deformation: %2.4f \n',[epsilon(1,1),epsilon(1,2),epsilon(2,2)]);
if UserSettings.gen_options.analysis_option == 1
    fprintf(fid,'** INPUT file for IMPLICIT STATIC analysis\n');
elseif UserSettings.gen_options.analysis_option == 2
    fprintf(fid,'** INPUT file for IMPLICIT STEADY STATE analysis\n');
elseif UserSettings.gen_options.analysis_option == 3
    fprintf(fid,'** INPUT file for EXPLICIT analysis\n');
end
fprintf(fid,'**\n');
fprintf(fid,'** ***************************************************************************\n');
fprintf(fid,'**\n');
%
% Including other files
%
if UserSettings.gen_options.pbcs_option == 1 % Periodic Boundary Conditions
    fprintf(fid,'** Include file with RVE mesh and Periodic Boundary Conditions:\n');
    fprintf(fid,'*INCLUDE, INPUT=include_mesh_FinalRVE.inp\n');
else
    fprintf(fid,'** Include file with RVE mesh (WITHOUT Periodic Boundary Conditions):\n');
    fprintf(fid,'*INCLUDE, INPUT=include_mesh_FinalRVE_noPBCs.inp\n');
end
%
% Section controls
fprintf(fid,'** Section controls:\n');
for jmat_phase=1:nmat_phases
    %
    phase_name = strcat('Phase_',num2str(jmat_phase));
    %
    %
    if UserSettings.gen_options.analysis_option == 1 || 2 % If running implicit
        % Matrix section controls
    %     fprintf(fid,'*SECTION CONTROLS, NAME=ctrls_Matrix, Hourglass=STIFFNESS\n');
        fprintf(fid,'*SECTION CONTROLS, NAME=ctrls_%s, Hourglass=ENHANCED\n',phase_name);
        fprintf(fid,'1., 1., 1.\n');
    else % if running explicit
        fprintf(fid,'*SECTION CONTROLS, NAME=ctrls_%s, Hourglass=RELAX STIFFNESS\n',phase_name);
        fprintf(fid,'1., 1., 1.\n');
    end
    %
end
%
% Amplitude Curve
if UserSettings.gen_options.analysis_option ~= 2 % Not in Implicit Steady State Dynamics
    fprintf(fid,'** Amplitude for the load curve:\n');
    fprintf(fid,'*Amplitude, name=AMP-1, DEFINITION=SMOOTH STEP\n');
    fprintf(fid,'0.,   0.,     %6.5E,      1.0\n',UserSettings.analyses.total_time);
end
%
% Include material properties
for jmat_phase=1:nmat_phases
    %
    material_name = strcat('Mat_Phase_',num2str(jmat_phase));
    material_type = mat_phases_IDs(jmat_phase);
    fprintf(fid,'** Include file with properties of Phase %i:\n',jmat_phase);
    if material_type == 1
        propsfile_name = strcat('include_props_',material_name,'_iso_elas.inp');
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_Matrix_iso_elas.inp\n');
    elseif material_type == 2
        propsfile_name = strcat('include_props_',material_name,'_ortho_elas.inp');
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_Matrix_ortho_elas.inp\n');
    elseif material_type == 3
        propsfile_name = strcat('include_props_',material_name,'_viscoelas.inp');
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_Matrix_viscoelas.inp\n');
    elseif material_type == 4
        propsfile_name = strcat('include_props_',material_name,'_parab_plas.inp');
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_Matrix_parab_plas.inp\n');
    elseif material_type == 5
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_cohesive.inp\n');
    elseif material_type == 6
        propsfile_name = strcat('include_props_',material_name,'_hyperelas_arruda.inp');
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_Matrix_hyperelas_arruda.inp\n');
    elseif material_type == 7
        propsfile_name = strcat('include_props_',material_name,'_hyperelas_neo.inp');
%         fprintf(fid,'*INCLUDE, INPUT=../../include_props_Matrix_hyperelas_neo.inp\n');
    else
        error('Selected material model does not exist');
    end
    %
    % Include material file:
    fprintf(fid,'*INCLUDE, INPUT=%s\n',propsfile_name);
end
%
% fprintf(fid,'** Include file with Particle properties:\n');
% if particle_material == 1
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_Particle_iso_elas.inp\n');
% elseif particle_material == 2
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_Particle_ortho_elas.inp\n');
% elseif particle_material == 3
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_Particle_viscoelas.inp\n');
% elseif particle_material == 4
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_Particle_parab_plas.inp\n');
% elseif particle_material == 6
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_Particle_hyperelas_arruda.inp\n');
% elseif particle_material == 7
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_Particle_hyperelas_neo.inp\n');
% else
%     disp('Selected Particle model does not exist');
% end
% %
% fprintf(fid,'** Include file with IP properties:\n');
% if interphase_material == 0
%     % No interphase exists
% elseif interphase_material == 1
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_iso_elas.inp\n');
% elseif interphase_material == 2
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_ortho_elas.inp\n');
% elseif interphase_material == 3
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_viscoelas.inp\n');
% elseif interphase_material == 4
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_parab_plas.inp\n');
% elseif interphase_material == 5
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_cohesive.inp\n');
% elseif interphase_material == 6
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_hyperelas_arruda.inp\n');
% elseif interphase_material == 7
%     fprintf(fid,'*INCLUDE, INPUT=../../include_props_IP_hyperelas_neo.inp\n');
% else
%     disp('Selected IP model does not exist');
% end
% %
fprintf(fid,'**\n');
%
fprintf(fid,'**Initial Conditions\n');
fprintf(fid,'*INITIAL CONDITIONS, TYPE=TEMPERATURE\n');
for jmat_phase=1:nmat_phases
    %
    phase_name = strcat('Phase_',num2str(jmat_phase));
    fprintf(fid,'FinalRVE-1.%s, %3g\n',phase_name,UserSettings.analyses.T0_global);
    % fprintf(fid,'FinalRVE-1.Particle, %3g\n',UserSettings.analyses.T0_global);
    % if interphase_material ~=0
    %     fprintf(fid,'FinalRVE-1.IPLayer, %3g\n',UserSettings.analyses.T0_global);
    % end
end
%
% STEP
%
if UserSettings.gen_options.analysis_option == 1 % Running implicit static analysis
    fprintf(fid,'** IMPLICIT STATIC Load Step 1 --------------------------------------\n');
    fprintf(fid,'*Step, name=Step-1, nlgeom=Yes\n');
    fprintf(fid,'*Static, Stabilize\n');
    fprintf(fid,' %6.5E,%6.5E,%6.5E,%6.5E\n',UserSettings.analyses.init_time_inc,UserSettings.analyses.total_time,UserSettings.analyses.min_time_inc,UserSettings.analyses.max_time_inc);
    % If you want to run dynamic:
%     fprintf(fid,'*STEP, INC=1000000, UNSYMM=YES\n');
%     fprintf(fid,'*DYNAMIC, HAFTOL=1E-3\n');
%     fprintf(fid,' 0.0001,1,1E-15,0.0001\n');
elseif UserSettings.gen_options.analysis_option == 2 % Running implicit steady state analysis
    fprintf(fid,'** IMPLICIT STEADY STATE Load Step 1 --------------------------------------\n');
    fprintf(fid,'*Step, name=Step-1, nlgeom=NO, perturbation\n');
    fprintf(fid,'*Steady State Dynamics, direct, friction damping=NO\n');
    fprintf(fid,' %6.5E, %6.5E, %i, %2.4f\n',UserSettings.analyses.f_initial,UserSettings.analyses.f_final,UserSettings.analyses.NoF,UserSettings.analyses.Bias);
elseif UserSettings.gen_options.analysis_option == 3 % Running explicit dynamic analysis
    fprintf(fid,'** EXPLICIT Load Step 1 --------------------------------------\n');
    fprintf(fid,'*Step, name=Step-1, nlgeom=YES\n');
    fprintf(fid,'*Dynamic, Explicit\n');
    fprintf(fid,', 1.,, 0.0001\n');
    fprintf(fid,'*Bulk Viscosity\n');
    fprintf(fid,'0.06, 1.2\n');
    fprintf(fid,'** Mass Scaling: Semi-Automatic\n');
    fprintf(fid,'**               Whole Model\n');
    fprintf(fid,'*Variable Mass Scaling, dt=1e-05, type=below min, frequency=1\n');
end
%
if UserSettings.gen_options.pbcs_option == 1
    fprintf(fid,'** FOR PERIODIC BOUNDARY CONDITIONS\n');
    fprintf(fid,'** Note: Node %s is "dummy1", i.e. Epsilon_1i\n',UserSettings.analyses.dummy1);
    fprintf(fid,'** Note: Node %s is "dummy2", i.e. Epsilon_2j\n',UserSettings.analyses.dummy2);
else
    fprintf(fid,'** Note: NOT USING PERIODIC BOUNDARY CONDITIONS\n');
end
%
%% Set the Periodic Boundary Conditions
%
%
if UserSettings.gen_options.pbcs_option == 1
    fprintf(fid,'*BOUNDARY,OP=NEW\n');
    fprintf(fid,'** Note: Epsilon_11 = %6.5E\n',epsilon(1,1));
    fprintf(fid,' FinalRVE-1.%s, 1, 1, %6.5E\n',UserSettings.analyses.dummy1,epsilon(1,1));
    if UserSettings.gen_options.analysis_option ~= 2 % Not a steady state dynamics simulation
        fprintf(fid,'** Note: Epsilon_22 = %6.5E\n',epsilon(2,2));
        fprintf(fid,' FinalRVE-1.%s, 2, 2, %6.5E\n',UserSettings.analyses.dummy2,epsilon(2,2));
        fprintf(fid,'** Note: Epsilon_12 = Epsilon_21 = %6.5E\n',epsilon(1,2));
        fprintf(fid,' FinalRVE-1.%s, 2, 2, %6.5E\n',UserSettings.analyses.dummy1,epsilon(1,2));
        fprintf(fid,' FinalRVE-1.%s, 1, 1, %6.5E\n',UserSettings.analyses.dummy2,epsilon(2,1));
        fprintf(fid,'** Note: Nodes restricting rigid body motion\n');
        fprintf(fid,'** FinalRVE-1.xSupport, 1, 1\n');
        fprintf(fid,'** FinalRVE-1.ySupport, 2, 2\n');
    end
else % No PBCs
%     fprintf(fid,'*BOUNDARY,OP=NEW\n');
%     fprintf(fid,'** Note: Nodes restricting rigid body motion\n');
%     fprintf(fid,'** FinalRVE-1.xSupport, 1, 1\n');
%     fprintf(fid,'** FinalRVE-1.ySupport, 2, 2\n');
%         fprintf(fid,'FinalRVE-1.LeftEdge, 1, 1\n');
%         fprintf(fid,'FinalRVE-1.VertexLB, 1, 2\n');
    if UserSettings.gen_options.analysis_option == 3
        fprintf(fid,'*BOUNDARY,OP=NEW, Amplitude=AMP-1\n'); % In case it is explicit use Amplitude curve
    else
        fprintf(fid,'*BOUNDARY,OP=NEW\n'); % If implicit no need for amplitude curve
    end
    % 11 direction
    fprintf(fid,'FinalRVE-1.RightEdge, 1, 1, %6.5E\n',epsilon1*Lx);
    fprintf(fid,'FinalRVE-1.LeftEdge, 1, 1, %6.5E\n',0.0);
    % 22 direction
    if UserSettings.gen_options.analysis_option ~= 2 % Not a steady state dynamics simulation
        fprintf(fid,'FinalRVE-1.TopEdge, 2, 2, %6.5E\n',epsilon2*Ly);
    end
    fprintf(fid,'FinalRVE-1.BotEdge, 2, 2, %6.5E\n',0.0);
end
%
%
%% Write remaining lines that are common to every load case
if UserSettings.gen_options.analysis_option ~= 2
fprintf(fid,'*OUTPUT, FIELD, TIME INTERVAL=%6.5E\n',UserSettings.analyses.output_time_interval);
else
    fprintf(fid,'*OUTPUT, FIELD\n');
end
fprintf(fid,'*ELEMENT OUTPUT\n');
fprintf(fid,' S, E, LE, ELEN, ELEDEN, ENER, IVOL, EVOL,\n');
% fprintf(fid,' S, E, EE, ELEN, IVOL, EVOL, ELEDEN, ENER, ELEN, LE, PE, PEEQ, PEMAG\n');
for jmat_phase=1:nmat_phases
    %
    phase_name = strcat('Phase_',num2str(jmat_phase));
    material_type = mat_phases_IDs(jmat_phase);
    if material_type == 4 % VUMAT
        fprintf(fid,'*ELEMENT OUTPUT, ELSET=FinalRVE-1.%s\n',phase_name);
        fprintf(fid,' SDV\n');
    elseif material_type == 5 % Cohesive element
        fprintf(fid,'*ELEMENT OUTPUT, ELSET=FinalRVE-1.%s\n',phase_name);
        fprintf(fid,' SDEG\n');
    end
end
%
fprintf(fid,'*NODE OUTPUT\n');
fprintf(fid,' RF, U\n');
if UserSettings.gen_options.analysis_option ~= 2
    fprintf(fid,'*OUTPUT, HISTORY, TIME INTERVAL=%6.5E\n',UserSettings.analyses.output_time_interval);
else
    fprintf(fid,'*OUTPUT, HISTORY\n');
end
fprintf(fid,'*ENERGY OUTPUT\n');
fprintf(fid,' ALLIE, ALLAE, ALLCD, ALLWK, ALLKE, ALLPD, ALLSE, ETOTAL\n');
%if UserSettings.gen_options.analysis_option == 1 % Running implicit static
%    fprintf(fid,'*CONTROLS,PARAMETERS=FIELD\n');
%    fprintf(fid,' 0.01, , , , 0.05\n');
%    fprintf(fid,'*CONTROLS,PARAMETERS=TIME INCREMENTATION\n');
%    fprintf(fid,' 10, 10, 4, 500, 10, 4, 30, 50, 10, 3\n');
%    fprintf(fid,' .25, .5, .75, .85\n');
%end
%fprintf(fid,'*CONTROLS,PARAMETERS=FIELD,FIELD=DISPLACEMENT\n');
%fprintf(fid,' 0.10, 1., , , 0.10\n');
fprintf(fid,'*END STEP\n');
%
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create run file for computer cluster ARES                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_file_fullpath = strcat(jDoE_dir,'/',run_file); % run file for ARES cluster
fid = fopen(run_file_fullpath,'wt');
%
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'#PBS -V\n');
fprintf(fid,'#PBS -N Samples\n');
fprintf(fid,'#PBS -q batch\n');
fprintf(fid,'#PBS -l nodes=%i:ppn=%i\n',UserSettings.analyses.nnodes,UserSettings.analyses.ncpus);
fprintf(fid,'#PBS -l walltime=999:99:99\n');
fprintf(fid,'\n');
fprintf(fid,'echo Working directory is $PBS_O_WORKDIR\n');
fprintf(fid,'cd $PBS_O_WORKDIR\n');
fprintf(fid,'echo Running on host `hostname`\n');
fprintf(fid,'echo Time is `date`\n');
fprintf(fid,'echo Directory is `pwd`\n');
fprintf(fid,'echo The following processors are allocated to this job:\n');
fprintf(fid,'echo `cat $PBS_NODEFILE`\n');
fprintf(fid,'NP=`wc -l < $PBS_NODEFILE`\n');
fprintf(fid,'\n');
fprintf(fid,'abq6142 double job=%s cpus=$NP interactive\n',abq_inp_file_noExt);
fprintf(fid,'\n');
fprintf(fid,'# These environment variables are available in every batch job\n');
fprintf(fid,'#\n');
fprintf(fid,'# $PBS_ENVIRONMENT set to PBS_BATCH to indicate that the job is a batch job; otherwise,\n');
fprintf(fid,'#                  set to PBS_INTERACTIVE to indicate that the job is a PBS interactive job\n');
fprintf(fid,'# $PBS_JOBID       the job identifier assigned to the job by the batch system\n');
fprintf(fid,'# $PBS_JOBNAME     the job name supplied by the user\n');
fprintf(fid,'# $PBS_NODEFILE    the name of the file that contains the list of nodes assigned to the job\n');
fprintf(fid,'# $PBS_QUEUE       the name of the queue from which the job is executed\n');
fprintf(fid,'# $PBS_O_HOME      value of the HOME variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_LANG      value of the LANG variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_LOGNAME   value of the LOGNAME variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_PATH      value of the PATH variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_MAIL      value of the MAIL variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_SHELL     value of the SHELL variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_TZ        value of the TZ variable in the environment in which qsub was executed\n');
fprintf(fid,'# $PBS_O_HOST      the name of the host upon which the qsub command is running\n');
fprintf(fid,'# $PBS_O_QUEUE     the name of the original queue to which the job was submitted\n');
fprintf(fid,'# $PBS_O_WORKDIR   the absolute path of the current working directory of the qsub command\n');
fprintf(fid,'#\n');
fprintf(fid,'# End of example PBS script\n');
%
fclose(fid);
%
% end
%
%     tmp=pwd;
%     cd(dir_fea_name);
%     copyfile('../../../dir_umat_final_thick.f','dir_umat.f');
%     copyfile('../../../yield_c.inp');
%     copyfile('../../../yield_t.inp');
%     copyfile('../../../run.pbs');
%     copyfile('../../../abaqus_v6_to_copy.env','abaqus_v6.env');
%
% unix(strcat(abaqus_path,' j=run ask_delete=OFF int cpus=4'));    cd(tmp);
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Run simulations automatically?                                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if run_simul_option == 1
%     cd(dir_fea_name);
% %     system(string);
%     % Run simulations in local desktop
%     string = strcat('abq6132 job=',abq_inp_file_noExt,' cpus=',num2str(UserSettings.analyses.ncpus),' double interactive ask_delete=OFF');
%     system(string);
% elseif run_simul_option == 2
%     % Run simulations in ARES cluster
%     string = strcat('qsub ',run_file);
%     system(string);
% end
%
%
% disp(' ');
% disp('run.inp file created');
% disp('Elapsed Time [min]: ');
% disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
