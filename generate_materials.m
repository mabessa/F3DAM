%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Generate material properties                                           %
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
function generate_materials(DoE_point,DiscInput_point,UserSettings,input_dir,jDoE_dir)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Material Properties From input Text File                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmat_phases = UserSettings.analyses.Number_Mat_Phases;

% Vector with the IDs of the "nmaterials" material phases
mat_phases_IDs = [DiscInput_point{1}];

% Cell with the strings (names) of the files with the properties for each
% material phase.
mat_phases_input_files = {DiscInput_point{2}};

analyses_dir = jDoE_dir;
%
for jmat_phase=1:nmat_phases
    %% Choose the material:
    input_file = strcat(input_dir,'/',mat_phases_input_files{jmat_phase});
    material_type = mat_phases_IDs(jmat_phase);
    material_name = strcat('Mat_Phase_',num2str(jmat_phase));
    %
    if material_type == 0
        material_type = 999; % No material model for the matrix is not accepted. So an error will occur.
    elseif material_type == 5
        material_type = 999; % Cohesive model for the matrix does not make sense. So an error will occur.
    end
    %
    %
    %% Write the include file for this material:
    %
    if material_type == 1 % Material is Isotropic Elastic
        %
        % Read input file
        [Dens,E,niu,...
            alpha,T0] = read_iso_elas(input_file,analyses_dir);
        %
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_iso_elas.inp');
        fid = fopen(propsfile_name,'w');
        %
        fprintf(fid,'**\n');
        fprintf(fid,'** %s material properties:\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*ELASTIC\n');
        fprintf(fid,' %6g, %6.4f, 0.\n',E,niu);
        fprintf(fid,'*EXPANSION, ZERO=%4g\n',T0);
        fprintf(fid,' %6.5E, 0.\n',alpha);
        fprintf(fid,' %6.5E, 200.\n',alpha);
        fprintf(fid,'**\n');
        %
        fclose(fid);
        %
        %
    elseif material_type == 2 % Orthotropic Elasticity (No damage implemented)
        %
        [Dens,E1,E2,niu12,niu21,...
            G12,G23,niu23,alpha1,...
            alpha2,T0,GcT,Xt,Xc] = read_ortho_elas(input_file,analyses_dir);
        %
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_ortho_elas.inp');
        fid = fopen(propsfile_name,'w');
        %
        fprintf(fid,'**\n');
        fprintf(fid,'** %s material properties:\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*ELASTIC, TYPE=ENGINEERING CONSTANTS\n');
        fprintf(fid,' %7g, %7g, %7g, %6.4f, %6.4f, %6.4f, %7g, %7g,\n',E1,E2,E2,...
            niu12,niu12,niu23,G12,G12);
        fprintf(fid,' %7g , 0.\n',G23);
        fprintf(fid,'*EXPANSION, TYPE=ORTHO, ZERO=%4g\n',T0);
        fprintf(fid,' %6.5E, %6.5E, %6.5E, 0.\n',alpha1,alpha2,alpha2);
        fprintf(fid,' %6.5E, %6.5E, %6.5E, 200.\n',alpha1,alpha2,alpha2);
        fprintf(fid,'**\n');
        %
        fclose(fid);
        %
        %
    elseif material_type == 3 % Viscoelasticity
        %
        [Dens,E,niu,...
            alpha,T0,...
            freq,Ep,Epp] = read_viscoelas(input_file,analyses_dir);
        %
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_viscoelas.inp');
        fid = fopen(propsfile_name,'w');
        %
        freq = fliplr(freq); % Reverse the order of elements
        Ep  = fliplr(Ep); % Reverse the order of elements
        Epp = fliplr(Epp); % Reverse the order of elements
        TanDelta = Epp./Ep;
        [m , location] = max(TanDelta); % maximum TanDelta and index location
        E0 = max(Ep); % Max Ep
        K0 = E0/(3*(1-2*niu)); % Now transform Young's Modulus into Shear Modulus assuming Bulk Modulus is constant
        E_inft = min(Ep); % Min Ep
        niu_inft = (3*K0-E_inft)/(6*K0); % long time Young's Modulus
        %
        Ecomp = complex(Ep,Epp);
        Gcomp = 3*K0*Ecomp./(9*K0-Ecomp);
        Gp = real(Gcomp);
        Gpp = imag(Gcomp);
        G_inft = min(Gp);
        %
        % Write material properties include file:
        fprintf(fid,'**\n');
        fprintf(fid,'** %s material properties:\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*ELASTIC, moduli=LONG TERM\n');
        fprintf(fid,' %6g, %6.6f\n',E_inft/(10^12),niu_inft);
        fprintf(fid,'*EXPANSION, ZERO=%4g\n',T0);
        fprintf(fid,' %6.5E, 0.\n',alpha);
        fprintf(fid,' %6.5E, 200.\n',alpha);
        %
        if jmat_phase ~=3 % Not the interphase
            %
            fprintf(fid,'*VISCOELASTIC, frequency=TABULAR\n');
            for i = 1:length(Gpp)
                fprintf(fid,' %6.5E, %6.5E, 0., 0., %6.5E\n',Gpp(i)/G_inft,1-Gp(i)/G_inft,freq(i)./(2*pi()));
            end
            fprintf(fid,'**\n');
            %
        else % Interphase with a chosen shifting and broadening
            % Do shifting
            for i =1:length(freq)
                f(i) = log10(freq(i))-R1;
            end
            % % Do broadening
            for i=1:length(f)
                fb(i) = ( f(i)-f(location) )*B1+f(location); % Broadening
            end
            %
            for i=1:length(fb)
                newfreq(i) = 10^fb(i);
            end
            %
            fprintf(fid,'*VISCOELASTIC, frequency=TABULAR\n');
            for i = 1:length(Gpp)
                fprintf(fid,' %6.5E, %6.5E, 0., 0., %6.5E\n',Gpp(i)/G_inft,1-Gp(i)/G_inft,newfreq(i)/(2*pi()));
            end
            fprintf(fid,'**\n');
        end
        %
        fclose(fid);
        %
        %
    elseif material_type == 4 % Plasticity using Paraboloidal model
        %
        % Read material input:
        [E,niu,alpha,T0,niup,...
            GcT,Xt,Xc,Dens,...
            Tref,EPRref,...
            alpha_c,alpha_t,...
            beta_c,beta_t,...
            NCC,YCC,NCT,YCT] = read_plas_paraboloid(input_file,analyses_dir);
        %
        % Write material properties to include in input file
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_parab_plas.inp');
        fid = fopen(propsfile_name,'w');
        %
        fprintf(fid,'**\n');
        % Material roperties
        fprintf(fid,'** %s material properties:\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*DEPVAR\n');
        fprintf(fid,' 22\n');
        fprintf(fid,'*USER MATERIAL, CONSTANTS=13\n');
        fprintf(fid,'** E, NIU, ALPHA, NIUp, XT, XC, Gc, EPRref\n');
        fprintf(fid,' %6g, %6.4f, %6.5e, %6.4f, %6g, %6g, %6.4f, %6.4f\n',...
            E,niu,alpha,niup,Xt,Xc,GcT,EPRref);
        fprintf(fid,'** Tref, ALPHA_C, ALPHA_T, BETA_C, BETA_T\n');
        fprintf(fid,' %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n',...
            Tref,alpha_c,alpha_t,beta_c,beta_t);
        fprintf(fid,'**\n');
        %
        fclose(fid);
        %
        % Write Hardening Curves to a file to be included in the USER MATERIAL
        %
        hardfile_name = strcat(analyses_dir,'/TWO_hardening_laws_for_vumat.f');
        fid = fopen(hardfile_name,'w');
        %
        fprintf(fid,'C Hardening Curves for VUMAT\n');
        fprintf(fid,'C ==================================================================\n');
        fprintf(fid,'C\n');
        fprintf(fid,'      PARAMETER (NCC=%d,NCT=%d)\n',round(NCC),round(NCT));
        fprintf(fid,'      DIMENSION YCC(NCC,2),YCT(NCT,2)\n');
        %
        % Write Compression hardening curve
        fprintf(fid,'C\n');
        fprintf(fid,'C     Compression hardening curve\n');
        fprintf(fid,'      data (YCC(i,1), i=1,NCC) /\n');
        for i=1:1:(NCC-1)
            fprintf(fid,'     *  %6.9f,\n',YCC(i,1));
        end
        fprintf(fid,'     *  %6.9f /\n',YCC(NCC,1));
        %
        fprintf(fid,'      data (YCC(i,2), i=1,NCC) /\n');
        for i=1:1:(NCC-1)
            fprintf(fid,'     *  %6.9f,\n',YCC(i,2));
        end
        fprintf(fid,'     *  %6.9f /\n',YCC(NCC,2));
        %
        % Write Tension hardening curve
        fprintf(fid,'C\n');
        fprintf(fid,'C     Tension hardening curve\n');
        fprintf(fid,'      data (YCT(i,1), i=1,NCT) /\n');
        for i=1:1:(NCT-1)
            fprintf(fid,'     *  %6.9f,\n',YCT(i,1));
        end
        fprintf(fid,'     *  %6.9f /\n',YCT(NCT,1));
        %
        fprintf(fid,'      data (YCT(i,2), i=1,NCT) /\n');
        for i=1:1:(NCT-1)
            fprintf(fid,'     *  %6.9f,\n',YCT(i,2));
        end
        fprintf(fid,'     *  %6.9f /\n',YCT(NCT,2));
        %
        fprintf(fid,'C\n');
        %
        fclose(fid);
        %
    elseif material_type == 5 % Cohesive law
        %
        % Read material input:
        [Dens,Enn,Ess,Ett,...
         sig_nn,sig_ss,sig_tt,...
         Gc_nn,Gc_ss,Gc_tt,...
         d_stab] = read_cohesive(input_file,analyses_dir);
        %
        % Write material properties to include in input file
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_cohesive.inp');
        fid = fopen(propsfile_name,'w');
        %
        % Interface properties (cohesive elements between fibres and matrix)
        fprintf(fid,'** %s material\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*ELASTIC, TYPE=TRACTION\n');
        fprintf(fid,' %6.5e, %6.5e, %6.5e\n',Enn,Ess,Ett);
        fprintf(fid,'*DAMAGE INITIATION, CRITERION=MAXS\n');
        fprintf(fid,' %6.5e, %6.5e, %6.5e\n',sig_nn, sig_ss, sig_tt);
        fprintf(fid,'*DAMAGE EVOLUTION, TYPE=ENERGY, MIXED MODE BEHAVIOR=BK, MODE MIX RATIO=ENERGY, POWER=1.45, SOFTENING=LINEAR\n');
        fprintf(fid,' %6.5e, %6.5e, %6.5e\n',Gc_nn,Gc_ss,Gc_tt);
        fprintf(fid,'*DAMAGE STABILIZATION\n');
        fprintf(fid,' %6.5e\n',d_stab);
        %
        %
    elseif material_type == 6 % Material is Hyperelastic (Arruda-Boyce model)
        %
        % Read input file
        [Dens,mu,lambda,...
         D,alpha,T0] = read_hyperelas_arruda(input_file,analyses_dir);
        %
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_hyperelas_arruda.inp');
        fid = fopen(propsfile_name,'w');
        %
        fprintf(fid,'**\n');
        fprintf(fid,'** %s material properties:\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*HYPERELASTIC, arruda-boyce\n');
        fprintf(fid,' %6g, %6.4f, %6.4f\n',mu,lambda,D);
        fprintf(fid,'*EXPANSION, ZERO=%4g\n',T0);
        fprintf(fid,' %6.5E, 0.\n',alpha);
        fprintf(fid,' %6.5E, 200.\n',alpha);
        fprintf(fid,'**\n');
        %
        fclose(fid);
        %
        %
    elseif material_type == 7 % Material is Hyperelastic (Neo-Hookean model)
        %
        % Read input file
        [Dens,C10,D1,...
         alpha,T0] = read_hyperelas_neo(input_file,analyses_dir);
        %
        propsfile_name = strcat(analyses_dir,'/include_props_',material_name,'_hyperelas_neo.inp');
        fid = fopen(propsfile_name,'w');
        %
        fprintf(fid,'**\n');
        fprintf(fid,'** %s material properties:\n',material_name);
        fprintf(fid,'*MATERIAL, NAME=%s\n',material_name);
        fprintf(fid,'*DENSITY\n');
        fprintf(fid,' %6.5E\n',Dens);
        fprintf(fid,'*HYPERELASTIC, neo hooke\n');
        fprintf(fid,' %6g, %6.4f\n',C10,D1);
        fprintf(fid,'*EXPANSION, ZERO=%4g\n',T0);
        fprintf(fid,' %6.5E, 0.\n',alpha);
        fprintf(fid,' %6.5E, 200.\n',alpha);
        fprintf(fid,'**\n');
        %
        fclose(fid);
        %
        %
    else % NO MATERIAL MODEL
        disp('Selected %s material model does not exist',material_name);
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nSelected %s material model does not exist\n',material_name);
        fclose(fid);
        %     return;
    end
    %
end % end loop for the material phases
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Miscellaneous Information                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% disp(' ');
% disp('Generation of include material .inp files COMPLETED');
% disp('Elapsed Time [min]: ');
% disp(toc/60);
% status_mesh = 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Read isotropic elasticity properties for this material
function [Dens,E,niu,...
          alpha,T0] = read_iso_elas(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end % Skip the first 13 lines
    Dens = str2num(fgetl(fid)); % Density - Dens
    fgetl(fid);
    E = str2num(fgetl(fid)); % Young's modulus
    fgetl(fid);
    niu = str2num(fgetl(fid)); % Poisson ratio
    fgetl(fid); 
    alpha = str2num(fgetl(fid)); % Coefficient of thermal expansion
    fgetl(fid);
    T0 = str2num(fgetl(fid)); % Initial temperature (also the stress free temperature, where there is no thermal expansion) - T0
    %
    fclose(fid);
end
%% Read ORTHOTROPIC elasticity + DAMAGE properties for this material
function [Dens,E1,E2,niu12,niu21,...
          G12,G23,niu23,alpha1,...
          alpha2,T0,GcT,Xt,Xc] = read_ortho_elas(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end
    Dens = str2num(fgetl(fid));
    fgetl(fid);
    E1 = str2num(fgetl(fid));
    E2 = str2num(fgetl(fid));
    fgetl(fid);
    niu12 = str2num(fgetl(fid));
    fgetl(fid);
    niu21 = niu12*E2/E1;
    G12 = str2num(fgetl(fid));
    G23 = str2num(fgetl(fid));
    fgetl(fid);
    niu23 = E2/(2*G23)-1;
    alpha1 = str2num(fgetl(fid));
    alpha2 = str2num(fgetl(fid));
    fgetl(fid);
    T0 = str2num(fgetl(fid)); % Stress free temperature, where there is no thermal expansion - T0
    fgetl(fid);
    GcT = str2num(fgetl(fid));
    fgetl(fid);
    Xt = str2num(fgetl(fid));
    Xc = str2num(fgetl(fid));
    %
    fclose(fid);
end
%% Read material with VISCOELASTICITY law
function [Dens,E,niu,...
          alpha,T0,...
          freq,Ep,Epp] = read_viscoelas(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end % Skip the first 13 lines
    Dens = str2num(fgetl(fid)); % Density - Dens
    fgetl(fid);
    E = str2num(fgetl(fid)); % Young's modulus
    fgetl(fid);
    niu = str2num(fgetl(fid)); % Poisson ratio
    fgetl(fid); 
    alpha = str2num(fgetl(fid)); % Coefficient of thermal expansion
    fgetl(fid);
    T0 = str2num(fgetl(fid)); % Initial temperature (also the stress free temperature, where there is no thermal expansion) - T0
    %
    for i=1:1:7, fgetl(fid); end % Skip the next 7 lines
    stat = 1; ncount = 0;
    while stat == 1
        tline=fgetl(fid);
        if tline ==-1
            break;
        else
            [data , stat] = str2num(tline);
            ncount = ncount + 1;
            freq(ncount) = data(1,1); % frequency
            Ep(ncount) = data(1,2); % E' (E prime)
            Epp(ncount) = data(1,3); % E'' (E double prime)
        end
    end
    %
    fclose(fid);
end
%% Read material with PARABOLOID PLASTICITY LAW
function [E,niu,alpha,T0,niup,...
          GcT,Xt,Xc,Dens,...
          Tref,EPRref,...
          alpha_c,alpha_t,...
          beta_c,beta_t,...
          NCC,YCC,NCT,YCT] = read_plas_paraboloid(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end % Skip the first 13 lines
    Dens = str2num(fgetl(fid)); % Density - Dens
    fgetl(fid);
    E = str2num(fgetl(fid)); % Young's modulus
    fgetl(fid);
    niu = str2num(fgetl(fid)); % Poisson ratio
    fgetl(fid); 
    alpha = str2num(fgetl(fid)); % Coefficient of thermal expansion
    fgetl(fid);
    T0 = str2num(fgetl(fid)); % Stress free temperature, where there is no thermal expansion - T0
    fgetl(fid);
    niup = str2num(fgetl(fid)); % Plastic Poisson Coefficient - niup
    fgetl(fid);
    GcT = str2num(fgetl(fid)); % Energy Release Rate - GcT
    fgetl(fid);
    Xt = str2num(fgetl(fid)); % Tensile Strength - Xt
    Xc = str2num(fgetl(fid)); % Compressive Strength - Xc
    fgetl(fid);
    Tref = str2num(fgetl(fid)); % Reference temperature for plasticity laws - Tref
    fgetl(fid);
    EPRref = str2num(fgetl(fid)); % Reference equivalent plastic strain - EPRref
    fgetl(fid);
    alpha_c = str2num(fgetl(fid)); % Parameter for strain-rate dependency in compression - ALPHA_C
    fgetl(fid);
    alpha_t = str2num(fgetl(fid)); % Parameter for strain-rate dependency in tension - ALPHA_T
    fgetl(fid);
    beta_c = str2num(fgetl(fid)); % Parameter for strain-rate dependency in compression - BETA_C
    fgetl(fid);
    beta_t = str2num(fgetl(fid)); % Parameter for strain-rate dependency in compression - BETA_T
    %
    % Read Hardening Curves From Text File
    %
    for i=1:1:13, fgetl(fid); end % Skip 13 lines
    % Compression Hardening curve
    NCC = str2num(fgetl(fid)); fgetl(fid);
    YCC = zeros(NCC,2);
    % Get equivalent plastic strain:
    for i=1:1:NCC
        YCC(i,1)=str2num(fgetl(fid));
    end
    fgetl(fid); fgetl(fid);
    % Get stresses:
    for i=1:1:NCC
        YCC(i,2)=str2num(fgetl(fid));
    end
    %
    %
    % Tension Hardening Curve
    for i=1:1:4, fgetl(fid); end
    NCT = str2num(fgetl(fid)); fgetl(fid);
    YCT = zeros(NCT,2);
    % Get equivalent plastic strain:
    for i=1:1:NCT
        YCT(i,1)=str2num(fgetl(fid));
    end
    fgetl(fid); fgetl(fid);
    % Get stresses:
    for i=1:1:NCT
        YCT(i,2)=str2num(fgetl(fid));
    end
    %
    %
    fclose(fid);
end
%% Read material with PARABOLOID PLASTICITY LAW
function [Dens,Enn,Ess,Ett,...
         sig_nn,sig_ss,sig_tt,...
         Gc_nn,Gc_ss,Gc_tt,...
         d_stab] = read_cohesive(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end % Skip the first 13 lines
    Dens = str2num(fgetl(fid));
    fgetl(fid);
    Enn = str2num(fgetl(fid));
    Ess = str2num(fgetl(fid));
    Ett = str2num(fgetl(fid)); 
    fgetl(fid);
    sig_nn = str2num(fgetl(fid));
    sig_ss = str2num(fgetl(fid));
    sig_tt = str2num(fgetl(fid));
    fgetl(fid);
    Gc_nn = str2num(fgetl(fid));
    Gc_ss = str2num(fgetl(fid));
    Gc_tt = str2num(fgetl(fid));
    fgetl(fid);
    d_stab = str2num(fgetl(fid));
    %
    fclose(fid);
end
%% Read hyperelasticity (arruda-boyce) properties for this material
function [Dens,mu,lambda,...
          D,alpha,T0] = read_hyperelas_arruda(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end % Skip the first 13 lines
    Dens = str2num(fgetl(fid)); % Density - Dens
    fgetl(fid);
    mu = str2num(fgetl(fid));
    fgetl(fid);
    lambda = str2num(fgetl(fid));
    fgetl(fid); 
    D = str2num(fgetl(fid));
    fgetl(fid); 
    alpha = str2num(fgetl(fid)); % Coefficient of thermal expansion
    fgetl(fid);
    T0 = str2num(fgetl(fid)); % Initial temperature (also the stress free temperature, where there is no thermal expansion) - T0
    %
    fclose(fid);
end
%% Read hyperelasticity (neo-hookean) properties for this material
function [Dens,C10,D1,...
          alpha,T0] = read_hyperelas_neo(input_file,analyses_dir)
    %
    inputfile_name = strcat(input_file,'.txt');
    fid = fopen(inputfile_name,'r');
    if fid == -1
        disp(strcat('Could not open file ',inputfile_name));
        errorfile_name = strcat(analyses_dir,'/ERROR_materials');
        fid = fopen(errorfile_name,'a+');
        fprintf(fid,'\nCould not open file %s',inputfile_name);
        fclose(fid);
    end
    %
    for i=1:1:13, fgetl(fid); end % Skip the first 13 lines
    Dens = str2num(fgetl(fid)); % Density - Dens
    fgetl(fid);
    C10 = str2num(fgetl(fid));
    fgetl(fid);
    D1 = str2num(fgetl(fid));
    fgetl(fid); 
    alpha = str2num(fgetl(fid)); % Coefficient of thermal expansion
    fgetl(fid);
    T0 = str2num(fgetl(fid)); % Initial temperature (also the stress free temperature, where there is no thermal expansion) - T0
    %
    fclose(fid);
end
