%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%            Matlab Routine for data-driven discovery of material         %
%               properties and constitutive behavior                      %
%                                                                         %
%                March 2017 - version 1.0                                 %
%                                                                         %
%                Code developed by:                                       %
%                Miguel A. Bessa - M.A.Bessa@tudelft.nl                   %
%              	                                                          %
% >---------------------------------------------------------------------< %
% > LICENSE:                                                            < %
% >                                                                     < %
% > The data-driven code is proprietary software provided under the     < %
% > following license:                                                  < %
% >                                                                     < %
% >    1) Non-commercial use (e.g. academic research):                  < %
% >      |                                                              < %
% >      | 1.1) Permission to use: The user is allowed to use this      < %
% >      |                         software free of charge and edit the < %
% >      |                         code for non-commercial purposes.    < %
% >      |                          However, the user cannot distribute < %
% >      |                          this software (main functions and   < %
% >      |                          auxiliary functions) to anyone else.< %
% >      |                          Any new user needs to contact       < %
% >      |                          Miguel A. Bessa directly to         < %
% >      |                          acquire the code and have permission< %
% >      |                          to use it.                          < %
% >      |  1.2) Attribution: Any work resulting from using part or     < %
% >      |                    the entire software (e.g., research       < %
% >      |                    publications, presentations, posters, etc)< %
% >      |                    has to acknowledge the usage of the       < %
% >      |                    software and cite the appropriate         < %
% >      |                    reference [1]:                            < %
% >      |                                                              < %
% >      |  [1] Bessa, M. A., et al. "A framework for data-driven       < %
% >      |      analysis of materials under uncertainty: Countering     < %
% >      |      the curse of dimensionality." Computer Methods in       < %
% >      |      Applied Mechanics and Engineering 320 (2017): 633-667.  < %
% >                                                                     < %
% >    2) Any other type of usage (for commercial purposes) is by       < %
% >       default prohibited. For this type of usage, please            < %
% >       contact Miguel A. Bessa.                                      < %
% >---------------------------------------------------------------------< %
%                                                                         %
% NOTE: If you received a copy of this code from a someone else, please   %
%       contact Miguel A. Bessa to get permission to use it.              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clc; close all; clear all; fclose('all');
%
tic;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Option Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Essential options:
%
UserSettings.gen_options.matlab_cpus = 2; % Number of "MATLAB workers" (CPUs used to run MATLAB code)
    %
UserSettings.gen_options.selectDoEs = 1:2; % DoE points that will be analyzed.
    % -> selectDoEs = vector with RVE numbers
    % Example 1 - Analyze 120 consecutive DoEs: selectDoEs = 1:1:120
    % Example 2 - Analyze DoEs 34, 53 and 111: selectDoEs = [34 , 53 , 111];
    % Example 3 - Analyze only DoE 27: selectDoEs = 27
    % Note: the code will create separate folders for each DoE being analyzed
    %
UserSettings.gen_options.selectDiscInputs = 1:2; % Discrete Input points that will be analyzed.
   % This parameter can be used the same way as the selectDOEs.
   %
UserSettings.gen_options.analysis_option = 1; % Type of analysis to run:
   % -> 1-Implicit, static
   % -> 2-Implicit, steady state dynamics
   % -> 3-Explicit
   %
UserSettings.gen_options.pbcs_option = 1; % Use Periodic Boundary Conditions? 0-NO, 1-YES
   % Important Note: in EXPLICIT using periodic boundary conditions is EXTREMELY slow
   % (different implementation needed)
   %
UserSettings.gen_options.strain_option = 1; % Choose strain measure to apply Boundary Conditions:
   % -> 0: microstructure_input_file provides SMALL
   %       STRAIN
   % -> 1: microstructure_input_file provides
   %       GREEN-LAGRANGE STRAIN
   %
UserSettings.gen_options.mesh_option = 1; % Generate mesh and auxiliary input files;
   % -> 0-NO
   % -> 1-YES
   %
UserSettings.gen_options.abq_inputs_option = 1; % Generates FEA inputs with B.C.s and materials:
   % -> 0-NO
   % -> 1-YES
   %
UserSettings.gen_options.run_simul_option = 1; % Do you want to run the simulations automatically?
   % -> 0-NO
   % -> 1-YES, but running in local desktop
   % -> 2-YES, but running in ARES cluster
   %
UserSettings.gen_options.postproc_option = 1; % Do RVE postprocessing?
   % -> 0-NO
   % -> 1-YES
   %
UserSettings.gen_options.convert_p2mat_option = 1; % Do you want to convert the python post-processing
   %file to a .mat file that MATLAB can read?
   % -> 0-NO
   % -> 1-YES
   %
UserSettings.gen_options.delodb_option = 0; % Delete RVE odb file?
   % -> 0-NO
   % -> 1-YES
   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for defining the Design of Experiments (Folder: 1_DoEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the Design of Experiments (DoE) with the continuous design variables:
%
% Example: parameters that define changes in geometry, material properties, etc.
%

% UserSettings.microstructure_input_file = 'INFO_10000RVEs_5Samples'; % Name of the file with all the sample microstructures (without ".mat")
%


% Name of the DoE file (to create or existing already) [without ".mat"]:
UserSettings.DoE.file_name = 'DOE_Bas_example';
UserSettings.DoE.size = 1000; % Number of DoE points
UserSettings.DoE.vars        = {'C1' , 'C2'}; % Input design variables names
% UserSettings.DoE_vars        = {'theta' , 'h'}; % Input design variables names
% UserSettings.DoE_vars        = {'theta' , 'h' , 't'}; % Input design variables names
                                                %this is just for reference and to
                                                %confirm the dimension of the input
                                                %design space: dim = length(vars)
UserSettings.DoE.LowerBounds = [0.05 , 0.05]; % (1 x dim vector) The lower bounds of each DoE variables.
UserSettings.DoE.UpperBounds = [0.3 , 0.3]; %(1 x dim vector) The upper bounds of each DoE variable.
% UserSettings.DoE_LowerBounds = [10.0  , 2.0e-3]; % (1 x dim vector) The lower bounds of each DoE variables.
% UserSettings.DoE_UpperBounds = [315.0 , 6.0e-3]; %(1 x dim vector) The upper bounds of each DoE variable.
% UserSettings.DoE_LowerBounds = [10.0  , 2.0e-3 , 50.0e-6]; % (1 x dim vector) The lower bounds of each DoE variables.
% UserSettings.DoE_UpperBounds = [315.0 , 8.0e-3 , 500.0e-6]; %(1 x dim vector) The upper bounds of each DoE variable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for defining the Discrete Input variables (folder: 2_Input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define DISCRETE inputs (set points with discrete properties that change)
%
% Examples: materials for the different phases, Loading conditions, etc. 
% Define the Discrete Inputs that you would like to vary discretely (Examples:
% materials, loading directions, etc.)
UserSettings.DiscreteInput.file_name = 'INPUT_Bas_example'; % Name of the file with all the discrete INPUTS (without ".mat")
% UserSettings.DiscreteInput.size = 3; % Number of Discrete Input points (sets of fixed values)
%
% The input file will have all the discrete inputs.
%If USER wants to define different material types for each phase it can start with
%the Material ID (material type) for each phase present in the RVE followed by the name of
%the file containing that material's properties.
% The material type can be one of the following:
    % -> 1-ISOTROPIC ELASTICITY
    % -> 2-ORTHOTROPIC ELASTICITY (Not very useful in 2D)
    % -> 3-VISCOELASTICITY
    % -> 4-PLASTICITY (VUMAT only for Explicit!)
    % -> 5-COHESIVE BEHAVIOR (Does not make sense for the bulk phases!)
    % -> 6-ARRUDA-BOYCE HYPERELASTICITY
    % -> 7-NEO-HOOKEAN HYPERELASTICITY
UserSettings.DiscreteInput.vars = {'Phase1_Mat1_ID', 'Phase1_Mat1_filename'}; % Discrete input variables

UserSettings.DiscreteInput.points = {};
% Include input properties for first input point:
UserSettings.DiscreteInput.points(1,:) = { 6 , 'Rubber_hyperelasticity_arruda' };
% Include input properties for second input point:
UserSettings.DiscreteInput.points(2,:) = { 7 , 'Rubber_hyperelasticity_neohookean' };
% -> Note: EACH ROW OF UserSettings.DiscreteInput.points REPRESENTES A INPUT POINT (columns correspond to each variable) 
%



% Remaining Fixed properties can be grouped into a single CELL variable for convenience: 
% Example: USER may want to define the material properties for all the
% phases only once if different materials are not being considered.
%
% Fixed Inputs:
UserSettings.analyses.FixedInput_vars = {'Epsilon_11' , 'Epsilon_22' , 'Epsilon_12'};
UserSettings.analyses.FixedInput_values = { 1.0 , -0.1 , 0.2 };
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options relative to the (FEM) analyses (folder: 3_Analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define ABAQUS executable path:
%abaqus_path = '/opt/abaqus/Commands/abq6132';
UserSettings.analyses.abaqus_path = 'abaqus';
%
UserSettings.analyses.folder = 'Bas_example_2'; % Name of the folder where all the files for all the simulations will be written
%
% Number of phases in RVE (THIS IS IMPORTANT! If this number is wrong, then
%the material properties of each phase will not be read accordingly)
UserSettings.analyses.Number_Mat_Phases = 1;
% RVE dimentions (in X direction and Y direction):
UserSettings.analyses.RVE_dimensions = [3.5,3.5];
% RVE center:
UserSettings.analyses.RVE_center = [3.5/2, 3.5/2];
% Approximate mesh size (this value should be sufficiently small for appropriate meshing):
UserSettings.analyses.mesh_size = 0.02;
% Mesh refinement factor:
UserSettings.analyses.mesh_refinement_factor = 1.25;
% Number of iterations to find appropriate mesh (scalar):
UserSettings.analyses.mesh_iter = 3;
%
% Number of CPUs to run FEM analyses:
UserSettings.analyses.ncpus = 1; % Number of CPUs to use in ARES cluster OR in LOCAL DESKTOP
UserSettings.analyses.nnodes = 1; % Number of nodes in ARES cluster only
% -> NOTE THAT THIS IS DIFFERENT FROM THE NUMBER OF MATLAB WORKERS, since
%    each Matlab worker (CPU) will be calling simulations with ncpus*nnodes (if
%    in cluster) or ncpus (if local desktop)
%
if UserSettings.gen_options.analysis_option == 1 % Implicit static
    UserSettings.analyses.init_time_inc = 0.02; % Initial time increment for simulation
    UserSettings.analyses.total_time = 1.0; % Total time of simulation
    UserSettings.analyses.min_time_inc = 1E-05; % Minimum time increment allowed
    UserSettings.analyses.max_time_inc = 0.02; % Maximum time increment allowed
    UserSettings.analyses.output_time_interval = 0.1; % Time interval to output to odb file
elseif UserSettings.gen_options.analysis_option == 2 % Steady state dynamics
    UserSettings.analyses.f_initial = 1e-7;   % initial frequency
    UserSettings.analyses.f_final = 1e7;   % final frequency
    UserSettings.analyses.NoF = 200;     % number of frequency points
    UserSettings.analyses.Bias = 1;     % bias parameter
end
%
% Global temperature of analysis:
UserSettings.analyses.T0_global = 120;
%
% Name of Node set for Dummy node related to Epsilon_1i
UserSettings.analyses.dummy1 = 'LR_point';
%
% Name Node set for Dummy node related to Epsilon_2i
UserSettings.analyses.dummy2 = 'TB_point';
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options relative to the Postprocessing of the analyses (folder: 4_Postprocessing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% What are the post-processing variables that you want to average over the RVE?
UserSettings.postproc.vars_int_pts = {}; % Variables that exist at INTEGRATION POINTS! (It can be an empty cell)
UserSettings.postproc.vars_whole_els = {'ELSE'}; % Variables that exist at the WHOLE ELEMENT! (It can be an empty cell)
UserSettings.postproc.vars_nodes = {}; % Variables that exist at the WHOLE ELEMENT! (It can be an empty cell)
%
UserSettings.postproc.local_averages = 0; % if local_averages = 1 it will also calculate averages in each material phase
%
%
%% END OF USER INPUTS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of code (in principle, USER does not need to change anything  %
% in this file from this point forward)                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save folder where this file is located
pwd;
main_dir = pwd;
%
% Save path to microstructure generation files
DoE_dir = strcat(cd,'/1_DoEs');
if exist(DoE_dir,'dir') == 0, mkdir(DoE_dir); end
% Save path to Input files for material properties
input_dir = strcat(cd,'/2_Inputs');
if exist(input_dir,'dir') == 0, mkdir(input_dir); end
% Create a folder for this RVE, if it doesn't already exist:
analyses_dir = strcat(cd,'/3_Analyses/',UserSettings.analyses.folder);
if exist(analyses_dir,'dir') == 0, mkdir(analyses_dir); end
% Create a folder for postprocessing, if it doesn't already exist:
postproc_dir = strcat(cd,'/4_Postprocessing/',UserSettings.analyses.folder);
if exist(postproc_dir,'dir') == 0, mkdir(postproc_dir); end
%
% Create file to indicate if there were ERRORS in generating the materials (ERROR_materials):
% Delete old ERROR_materials if it exists:
errorfile_name = strcat(analyses_dir,'/ERROR_FILE');
if exist(errorfile_name,'file') == 2, delete(errorfile_name); end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create or LOAD the DoE file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Design of Experiments (DoE):
if exist([DoE_dir,'/',UserSettings.DoE.file_name,'.mat'],'file') == 0
    %
    DoE_dim = length(UserSettings.DoE.vars);
    DoE.points = generate_DoE(DoE_dim, UserSettings.DoE.size, UserSettings.DoE.LowerBounds, UserSettings.DoE.UpperBounds);
    DoE.size = UserSettings.DoE.size;
    DoE.vars = UserSettings.DoE.vars;
    DoE.LowerBounds = UserSettings.DoE.LowerBounds;
    DoE.UpperBounds = UserSettings.DoE.UpperBounds;
    %
    % Save the RVEs_data structure to a MATLAB format:
    savefile = strcat(DoE_dir,'/',UserSettings.DoE.file_name,'.mat');
    save(savefile,'DoE');
else
    disp(' ');
    disp('DoE file already exists! File NOT updated');
    %
%     errorfile_name = strcat(analyses_dir,'/ERROR_FILE'); % Write ERROR file
%     fid = fopen(errorfile_name,'a+');
%     fprintf(fid,'DoE file already exists! File NOT updated\n');
%     fclose(fid);
    loadfile = strcat(DoE_dir,'/',UserSettings.DoE.file_name,'.mat');
    load(loadfile,'DoE');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create or LOAD the Input file (with discrete input properties)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([input_dir,'/',UserSettings.DiscreteInput.file_name,'.mat'],'file') == 0
    %
%     Input.size = UserSettings.DiscreteInput.size;
    DiscInput.vars = UserSettings.DiscreteInput.vars;
    DiscInput.points = UserSettings.DiscreteInput.points;
    %
    % Save the RVEs_data structure to a MATLAB format:
    savefile = strcat(input_dir,'/',UserSettings.DiscreteInput.file_name,'.mat');
    save(savefile,'DiscInput');
else
    disp(' ');
    disp('INPUT file already exists! File NOT updated');
    %
    loadfile = strcat(input_dir,'/',UserSettings.DiscreteInput.file_name,'.mat');
    load(loadfile,'DiscInput');
end
%
%
% Delete postprocessing files that may have been used previously with
% different CPUs:
if UserSettings.gen_options.postproc_option == 1
    % Delete postprocessing files for different CPUs (these were just
    % temporary)
    cd(postproc_dir);
%             string = ['rm ',abq_inp_file{jDoE},'*.odb'];
    string = ['rm ',strcat('STRUCTURES_postprocessing_variables_CPU'),'*'];

    system(string);
    %
    cd(main_dir);
end
%
% Restart MATLAB parallel pool:
delete(gcp('nocreate'));
parpool('local', UserSettings.gen_options.matlab_cpus); % USE WORKERS
%
string = 'rm tempname.mat';
system(string);
%
% Create file to indicate if there were ERRORS in generating the materials (ERROR_materials):
% Delete old ERROR_materials if it exists:
errorfile_name = strcat(analyses_dir,'/ERROR_materials');
if exist(errorfile_name,'file') == 2, delete(errorfile_name); end
% %
%
for iDiscInp=UserSettings.gen_options.selectDiscInputs
    %
    iDiscInput_dir = strcat(analyses_dir,'/DiscInput_point',int2str(iDiscInp));
    % Create a folder for this set of Fixed Variables, if it doesn't already exist:
    if exist(iDiscInput_dir,'dir') == 0, mkdir(iDiscInput_dir); end
    %
    parfor jDoE = UserSettings.gen_options.selectDoEs
%     for jDoE = UserSettings.gen_options.selectDoEs
        %
        t = getCurrentTask();
        CPUid(jDoE) = t.ID;
%         CPUid(jDoE) = 1;
        %
        % Create a folder for this DoE, if it doesn't already exist:
        jDoE_dir = strcat(iDiscInput_dir,'/DoE_point',int2str(jDoE));
        if exist(jDoE_dir,'dir') == 0, mkdir(jDoE_dir); end
        %
        % Delete old ERROR_FILE if it exists:
        errorfile_name = strcat(jDoE_dir,'/ERROR_FILE');
        if exist(errorfile_name,'file') == 2, delete(errorfile_name); end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generates Files with Mesh, PBCs and Material Properties             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if UserSettings.gen_options.mesh_option == 1
            % Before executing anything, delete previous .rpy and .rec files of ABAQUS
            string = 'rm abaqus.rpy* abaqus*.rec';
            system(string);
            %
            % Call function: generate_mesh
            generate_mesh(DoE.points(jDoE,:),DiscInput.points(iDiscInp,:),UserSettings,jDoE,iDiscInp,jDoE_dir);
        end
        %
        % Call function that reads the material properties
        generate_materials(DoE.points(jDoE,:),DiscInput.points(iDiscInp,:),UserSettings,input_dir,jDoE_dir);
        %
        %
        abq_inp_file_noExt = {};
        run_file = {};
        %
        if UserSettings.gen_options.analysis_option == 1
            abq_inp_file_noExt{jDoE} = strcat('DoE',int2str(jDoE),'_implicit_static'); % file without extension
            run_file{jDoE} = strcat('run_DoE',int2str(jDoE),'_implicit_static.sh'); % run file for ARES cluster
        elseif UserSettings.gen_options.analysis_option == 2
            abq_inp_file_noExt{jDoE} = strcat('DoE',int2str(jDoE),'_implicit_steadystate'); % file without extension
            run_file{jDoE} = strcat('run_DoE',int2str(jDoE),'_implicit_steadystate.sh'); % run file for ARES cluster
        elseif UserSettings.gen_options.analysis_option == 3
            abq_inp_file_noExt{jDoE} = strcat('DoE',int2str(jDoE),'_explicit'); % file without extension
            run_file{jDoE} = strcat('run_DoE',int2str(jDoE),'_explicit.sh'); % run file for ARES cluster
        else % Implicit analysis
            disp(strcat('Analysis option selected is unavailable'));
            errorfile_name = strcat(DoE_point_dir,'/ERROR_FILE');
            fid = fopen(errorfile_name,'a+');
            fprintf(fid,'\nAnalysis option selected is unavailable');
            fclose(fid);
            error('Analysis option unavailable');
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate FEA files for each simulation of the RVE and B.C.s         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if UserSettings.gen_options.abq_inputs_option == 1
            % Before executing anything, delete previous .rpy and .rec files of ABAQUS
            string = 'rm abaqus.rpy* abaqus*.rec';
            system(string);
            % Call function: generate_simulation_files
            generate_simulation_files(DoE.points(jDoE,:),DiscInput.points(iDiscInp,:),UserSettings,jDoE,iDiscInp,jDoE_dir,postproc_dir,abq_inp_file_noExt{jDoE},run_file{jDoE},CPUid(jDoE)); % Applied Boundary Condition number iBC
            %
        end
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run simulation automatically?                                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if UserSettings.gen_options.run_simul_option == 1
            cd(jDoE_dir);
            %     system(string);
            % Run simulations in local desktop
            string = [UserSettings.analyses.abaqus_path,' job=',abq_inp_file_noExt{jDoE},...
                ' cpus=',num2str(UserSettings.analyses.ncpus),' double interactive ask_delete=OFF'];
            system(string);
        elseif UserSettings.gen_options.run_simul_option == 2
            cd(jDoE_dir);
%                     string = ['chmod +x ',runfile_name{iBC}];
%                     system(string);
            % Run simulations in ARES cluster
            string = ['qsub ',run_file{jDoE}];
            system(string);
        end
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate Postprocessing file for all RVEs                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        cd(main_dir); % Return to main folder
        %
        if UserSettings.gen_options.postproc_option == 1
            % Before executing anything, delete previous .rpy and .rec files of ABAQUS
            string = 'rm abaqus.rpy* abaqus*.rec';
            system(string);
            % Call function: generate_postproc_file
%             for iBC=selectBCs_jRVE
                generate_postproc_file(DoE.points(jDoE,:),DiscInput.points(iDiscInp,:),UserSettings,jDoE,iDiscInp,jDoE_dir,postproc_dir,abq_inp_file_noExt{jDoE},CPUid(jDoE)); % Applied Boundary Condition number iBC
%             end
            %
        end
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete ODB files for this RVE?                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if UserSettings.gen_options.delodb_option == 1
            % Delete odb for each BC
            cd(jDoE_dir);
    %             string = ['rm ',abq_inp_file{jDoE},'*.odb'];
            string = ['rm ',strcat('DoE',int2str(jDoE)),'*'];

            system(string);
            %
        end
        %
        cd(main_dir); % Return to main folder
        %
    end
end
%
%
% Since we are using parallel computing, we have to merge the Pickle
% postprocessing files from each CPU:
if UserSettings.gen_options.postproc_option == 1
    generate_unique_postproc_file(UserSettings,postproc_dir);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the python post-processing file RVEs_postprocessing_variables.p %
% to a MATLAB post-processing file RVEs_postprocessing_variables.m        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if UserSettings.gen_options.convert_p2mat_option == 1
    % Read Python homogenized variables from the pickle file generated
    file_postproc = strcat(postproc_dir,'/RVEs_postprocessing_variables.p');
    file_load = loadpickle(file_postproc);
    RVEs_data = file_load.RVEs_data;
    clear file_load;
    % The structured variable RVEs_data constains all the postprocessing
    %information from all the RVEs analyzed.
%     if pbcs_option == 1
%         DiscInput_strings = fieldnames(RVEs_data);
%         nDiscInputs = numel(DiscInput_strings);
%         for iDiscInp= 1:nDiscInputs
%             DoE_strings = fieldnames(RVEs_data.(DiscInput_strings{iDiscInp}));
%             nDoEs = numel(DoE_strings);
%             for iDoE = 1:nDoEs
%                 Green_strain_NEW = zeros(size(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).F));
%                 PK2_NEW = zeros(size(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).P));
%                 for iFrame = 1:size( RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).F, 1)
%                     F_local = reshape(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).F(iFrame,:,:),size(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).F,2),size(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).F,3)); % Deformation gradient for this frame
%                     P_local = reshape(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).P(iFrame,:,:),size(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).P,2),size(RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).P,3)); % nominal stress for this frame
%                     % Compute Green-Lagrange strain:
%                     Green_strain_NEW(iFrame,:,:) = 1/2*(F_local'*F_local-eye(size(F_local,1)));
%                     % Compute 2nd P-K stress:
%                     PK2_NEW(iFrame,:,:) = P_local*inv(F_local');
%                 end
%                 Green_strain_NEW-RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).Green_strain
%                 PK2_NEW - RVEs_data.(DiscInput_strings{iDiscInp}).(DoE_strings{iDoE}).PK2
%             end
%         end
%     end
    %
%end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to identify some simulations that had problems or simulations that  %
    % did not run at all.                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Code to check the simulation results
    % This code determines the simulations that need to be run again.
    %
%     RVEs_data_file = strcat(postproc_dir,'/RVEs_postprocessing_variables.mat');
%     load(RVEs_data_file); % Load variable with microstructures
    %
    tolerance = 1e-3; % Tolerance used to check total simulation time as well as if the RVEs have the correct initial volume
    %
    %%
    icount = 0;
    DiscInput_strings = fieldnames(RVEs_data);
    nDiscInput = numel(DiscInput_strings);
    %
    total_number_of_expected_DoEs = DoE.size;
    %
    %
    for iDiscInput = 1:nDiscInput
        DoE_strings = fieldnames(RVEs_data.(DiscInput_strings{iDiscInput}));
        nDoEs = numel(DoE_strings);
        unfinished_simulations.(DiscInput_strings{iDiscInput}) = [];
        missing_simulations.(DiscInput_strings{iDiscInput}) = zeros(total_number_of_expected_DoEs,3);
        missing_simulations.(DiscInput_strings{iDiscInput})(1:total_number_of_expected_DoEs,1) = 1:total_number_of_expected_DoEs;
        for jDoE=1:nDoEs
            % Extract just the number of the RVE (not the entire string);
            Str = DoE_strings(jDoE);
            Str = Str{1};
            Key = 'DoE';
            Index = strfind(Str, Key);
            jDoE_integer = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
            %
            % Eliminate this RVE from the list of RVEs that were not meshed
            missing_simulations.(DiscInput_strings{iDiscInput})( find(missing_simulations.(DiscInput_strings{iDiscInput})(:,1)==jDoE_integer), : ) = [];
            %
            if UserSettings.gen_options.analysis_option ~= 2 % Not a steady state dynamics simulation
                total_vol = RVEs_data.(DiscInput_strings{iDiscInput}).(DoE_strings{jDoE}).total_vol;
                if abs(UserSettings.analyses.total_time-1.0) > tolerance
                % Comment line above and uncomment next if the RVEs do not have voids (checks volume consistency)
%                 if abs(UserSettings.analyses.total_time-1.0) > tolerance || abs(total_vol-prod(UserSettings.analyses.RVE_dimensions)) > tolerance
                    icount = icount + 1;
                    unfinished_simulations.(DiscInput_strings{iDiscInput})(icount,1) = jDoE_integer;
                    % Save total time and total volume
                    unfinished_simulations.(DiscInput_strings{iDiscInput})(icount,2) = UserSettings.analyses.total_time;
                    unfinished_simulations.(DiscInput_strings{iDiscInput})(icount,3) = total_vol;
                end
            end
            %
        end
        %
        % Gather all the simulations that may have problems
        Bad_DoEs.(DiscInput_strings{iDiscInput}) =  [missing_simulations.(DiscInput_strings{iDiscInput}); unfinished_simulations.(DiscInput_strings{iDiscInput})];
        % Sort the elements:
        Bad_DoEs.(DiscInput_strings{iDiscInput}) = sortrows(Bad_DoEs.(DiscInput_strings{iDiscInput}),1);
    end
    %
    % Save the RVEs_data structure to a MATLAB format:
    savefile = strcat(postproc_dir,'/RVEs_postprocessing_variables.mat');
    save(savefile,'RVEs_data','Bad_DoEs');
    %
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of MAIN Routine                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
