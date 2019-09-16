%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Generate mesh of the RVEs                                              %
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
function generate_mesh(DoE_point,DiscInput_point,UserSettings,jDoE,iDiscInp,jDoE_dir)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

C1 = DoE_point(1);

C2 = DoE_point(2);

% DiscInput_points(1);

Lx = UserSettings.analyses.RVE_dimensions(1);

Ly = UserSettings.analyses.RVE_dimensions(2);

% LM = max(Lx,Ly)*0.5; % Margin of the RVE for construction purposes

RVEcenter = UserSettings.analyses.RVE_center;

Mesh_size = UserSettings.analyses.mesh_size;

refine_factor = UserSettings.analyses.mesh_refinement_factor;

niter = UserSettings.analyses.mesh_iter;

dummy1 = UserSettings.analyses.dummy1;

dummy2 = UserSettings.analyses.dummy2;

abaqus_path = UserSettings.analyses.abaqus_path;

nmat_phases = UserSettings.analyses.Number_Mat_Phases;

% Vector with the IDs of the "nmaterials" material phases
mat_phases_IDs = [DiscInput_point{1}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Python File for Abaqus CAE                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Started creation of Python file');
%
scriptfile_name = strcat(jDoE_dir,'/script_RVE_meshing.py');
fid = fopen(scriptfile_name,'wt');
%
% Heading
%
fprintf(fid,'# Abaqus/CAE script\n');
fprintf(fid,'# Created by M.A. Bessa and M. Fernandes on %s\n',datestr(now));
fprintf(fid,'#\n');
fprintf(fid,'from abaqus import *\n');
fprintf(fid,'from abaqusConstants import *\n');
fprintf(fid,'from caeModules import *\n');
fprintf(fid,'from driverUtils import executeOnCaeStartup\n');
% fprintf(fid,'from numpy import *\n');
fprintf(fid,'import numpy as np\n');
fprintf(fid,'from math import *\n');
fprintf(fid,'executeOnCaeStartup()\n');
fprintf(fid,'session.viewports[''Viewport: 1''].partDisplay.geometryOptions.setValues(\n');
fprintf(fid,'    referenceRepresentation=ON)\n');
fprintf(fid,'Mdb()\n');
% A new model database has been created.
% The model "Model-1" has been created.
fprintf(fid,'os.chdir(r''%s'')\n',jDoE_dir);
% Problem parameteres needs to be defined and double check before running the code
fprintf(fid,'\n');
fprintf(fid,'Lx=%7.5f #sample dimention in um VF9 and 28%%\n',Lx);
fprintf(fid,'Ly=%7.5f #sample dimention in um VF9 and 28%%\n',Ly);
% fprintf(fid,'LM=%7.5f # Margin width\n',LM);
fprintf(fid,'RVEcenter = [%7.5f,%7.5f] # Center position of RVE\n',RVEcenter(1),RVEcenter(2));
% Mesh parameter. Characteristic length needs to be small enough.
fprintf(fid,'Mesh_size=%7.5f # Mesh parameter\n',Mesh_size);
%
%
fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=None)\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', sheetSize=%7.5f,)\n',3.0*max(Lx,Ly));
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.rectangle(point1=(RVEcenter[0]-Lx/2,RVEcenter[1]-Ly/2), point2=(RVEcenter[0]+Lx/2, RVEcenter[1]+Ly/2))\n');

%% User code defining geometry of this RVE:
fprintf(fid,'### MODEL PARAMETERS ###\n');
fprintf(fid,'NUMBER_OF_POINTS=100 # NUMBER OF POINTS TO GENERATE THE CENTER HOLE THROUGH THE PARAMETRIC FUNCTION\n');
fprintf(fid,'RO=1.\n');
fprintf(fid,'C1=%1.8e\n',C1);
fprintf(fid,'C2=%1.8e\n',C2);
% fprintf(fid,'Lx=3.5\n');
% fprintf(fid,'Ly=3.5\n');
fprintf(fid,'TOL=1E-5\n');


fprintf(fid,'### GENERATING INNER SHAPE OF THE GEOMETRY ###\n');
fprintf(fid,'\n');
fprintf(fid,'THETAALL=np.linspace(0.,2.*pi,NUMBER_OF_POINTS)\n');
fprintf(fid,'POINTS=[]\n');
fprintf(fid,'for i in xrange(NUMBER_OF_POINTS):\n');
fprintf(fid,'    THETA=THETAALL[i]\n');
fprintf(fid,'    rr=RO*(1.+C1*cos(4.*THETA)+C2*cos(8.*THETA))\n');
fprintf(fid,'    POINTS.append((RVEcenter[0]+rr*cos(THETA),RVEcenter[1]+rr*sin(THETA)))\n');
fprintf(fid,'    if i==0: xFirst=RVEcenter[0]+rr*cos(THETA);yFirst=RVEcenter[1]+rr*sin(THETA)\n');
fprintf(fid,'    if i==NUMBER_OF_POINTS-1: POINTS.append((xFirst,yFirst))\n');
fprintf(fid,'\n');
fprintf(fid,'### GENERATE SPLINE SHAPE ###\n');
fprintf(fid,'s.Spline(points=POINTS)\n');

fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''FinalRVE'', dimensionality=TWO_D_PLANAR,\n');
fprintf(fid,'    type=DEFORMABLE_BODY)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''FinalRVE'']\n');
fprintf(fid,'p.BaseShell(sketch=s)\n');
fprintf(fid,'s.unsetPrimaryObject()\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''FinalRVE'']\n');
fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=p)\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
fprintf(fid,'\n');




% Create PART
% fprintf(fid,'mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', sheetSize=200.0)\n');
% fprintf(fid,'mdb.models[''Model-1''].Part(dimensionality=TWO_D_PLANAR, name=''Part-1'', type=\n');
% fprintf(fid,'    DEFORMABLE_BODY)\n');
% fprintf(fid,'\n');
% fprintf(fid,'mdb.models[''Model-1''].sketches[''__profile__''].rectangle(point1=(-Lx/2., -Ly/2.), \n');
% fprintf(fid,'    point2=(Lx/2., Ly/2.))\n');
% fprintf(fid,'mdb.models[''Model-1''].sketches[''__profile__''].Spline(points=POINTS)\n');
% fprintf(fid,'\n');
% fprintf(fid,'\n');
% fprintf(fid,'mdb.models[''Model-1''].parts[''Part-1''].BaseShell(sketch=\n');
% fprintf(fid,'    mdb.models[''Model-1''].sketches[''__profile__''])\n');
% fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
% fprintf(fid,'\n');
% fprintf(fid,'#PART DEFINITION\n');
% fprintf(fid,'Part_Full=mdb.models[''Model-1''].parts[''FinalRVE'']\n');
% fprintf(fid,'\n');
% fprintf(fid,'\n');

% CREATE ASSEMBLY:
fprintf(fid,'### CREATING ASSEMBLY ###\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.DatumCsysByDefault(CARTESIAN)\n');
fprintf(fid,'#INSTANCE DEFINITION\n');
fprintf(fid,'Instance_Full=mdb.models[''Model-1''].rootAssembly.Instance(dependent=ON, name=''FinalRVE-1'', part=p)\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

% CREATE SETS FOR EACH MATERIAL PHASE:
fprintf(fid,'p = mdb.models[''Model-1''].parts[''FinalRVE'']\n');
% Elements of Material Phase 1:
fprintf(fid,'# All faces:\n');
fprintf(fid,'f = p.faces\n');
fprintf(fid,'faces = f[:]\n');
fprintf(fid,'p.Set(faces=faces, name=''Phase_1'')\n');
fprintf(fid,'\n');
%
% Add here other material phases:

%
%
%% Common to every RVE:

fprintf(fid,'# Create sets useful for meshing\n');
fprintf(fid,'delta = min((min(Lx,Ly)/1000),Mesh_size/10)\n');
%
% Edges:
fprintf(fid,'# edges\n');
fprintf(fid,'s = p.edges\n');
fprintf(fid,'edgesLEFT = s.getByBoundingBox(RVEcenter[0]-Lx/2-delta,RVEcenter[1]-Ly/2-delta,0,RVEcenter[0]-Lx/2+delta,RVEcenter[1]+Ly/2+delta,0)\n');
fprintf(fid,'p.Set(edges=edgesLEFT, name=''LeftEdge'')\n');
fprintf(fid,'edgesRIGHT = s.getByBoundingBox(RVEcenter[0]+Lx/2-delta,RVEcenter[1]-Ly/2-delta,0,RVEcenter[0]+Lx/2+delta,RVEcenter[1]+Ly/2+delta,0)\n');
fprintf(fid,'p.Set(edges=edgesRIGHT, name=''RightEdge'')\n');
fprintf(fid,'edgesTOP = s.getByBoundingBox(RVEcenter[0]-Lx/2-delta,RVEcenter[1]+Ly/2-delta,0,RVEcenter[0]+Lx/2+delta,RVEcenter[1]+Ly/2+delta,0)\n');
fprintf(fid,'p.Set(edges=edgesTOP, name=''TopEdge'')\n');
fprintf(fid,'edgesBOT = s.getByBoundingBox(RVEcenter[0]-Lx/2-delta,RVEcenter[1]-Ly/2-delta,0,RVEcenter[0]+Lx/2+delta,RVEcenter[1]-Ly/2+delta,0)\n');
fprintf(fid,'p.Set(edges=edgesBOT, name=''BotEdge'')\n');
% Vertices:
fprintf(fid,'# vertices\n');
fprintf(fid,'v = p.vertices\n');
fprintf(fid,'vertexLB = v.getByBoundingBox(RVEcenter[0]-Lx/2-delta,RVEcenter[1]-Ly/2-delta,0,RVEcenter[0]-Lx/2+delta,RVEcenter[1]-Ly/2+delta,0)\n');
fprintf(fid,'p.Set(vertices=vertexLB, name=''VertexLB'')\n');
fprintf(fid,'vertexRB = v.getByBoundingBox(RVEcenter[0]+Lx/2-delta,RVEcenter[1]-Ly/2-delta,0,RVEcenter[0]+Lx/2+delta,RVEcenter[1]-Ly/2+delta,0)\n');
fprintf(fid,'p.Set(vertices=vertexRB, name=''VertexRB'')\n');
fprintf(fid,'vertexRT = v.getByBoundingBox(RVEcenter[0]+Lx/2-delta,RVEcenter[1]+Ly/2-delta,0,RVEcenter[0]+Lx/2+delta,RVEcenter[1]+Ly/2+delta,0)\n');
fprintf(fid,'p.Set(vertices=vertexRT, name=''VertexRT'')\n');
fprintf(fid,'vertexLT = v.getByBoundingBox(RVEcenter[0]-Lx/2-delta,RVEcenter[1]+Ly/2-delta,0,RVEcenter[0]-Lx/2+delta,RVEcenter[1]+Ly/2+delta,0)\n');
fprintf(fid,'p.Set(vertices=vertexLT, name=''VertexLT'')\n');
fprintf(fid,'\n');
%
% Make MESH:
fprintf(fid,'##########################################################\n');
fprintf(fid,'#Mesh the FinalRVE part\n');
fprintf(fid,'import os\n');
fprintf(fid,'import fileinput\n');
fprintf(fid,'import shutil\n');
fprintf(fid,'#\n');
fprintf(fid,'niter=1 # iteration number for meshing procedure\n');
fprintf(fid,'status_mesh = 0 # flag signaling if mesh was created\n');
fprintf(fid,'refine_factor = %7.5f # Parameter used to refine mesh (should be larger than 1)\n',refine_factor);
%
% START WHILE LOOP (until mesh is found or maximum number of iterations is reached)
fprintf(fid,'while status_mesh == 0:\n');
fprintf(fid,'    session.viewports[''Viewport: 1''].partDisplay.setValues(sectionAssignments=OFF,\n');
fprintf(fid,'        engineeringFeatures=OFF, mesh=ON)\n');
fprintf(fid,'    session.viewports[''Viewport: 1''].partDisplay.meshOptions.setValues(\n');
fprintf(fid,'        meshTechnique=ON)\n');
fprintf(fid,'    p = mdb.models[''Model-1''].parts[''FinalRVE'']\n');
fprintf(fid,'    p.seedPart(size=Mesh_size, deviationFactor=0.4, minSizeFactor=0.4)\n');
fprintf(fid,'    p.seedEdgeBySize(edges=edgesLEFT, size=Mesh_size, deviationFactor=0.4,\n');
% fprintf(fid,'        constraint=FINER)\n');
fprintf(fid,'        constraint=FIXED)\n');
fprintf(fid,'    p.seedEdgeBySize(edges=edgesRIGHT, size=Mesh_size, deviationFactor=0.4,\n');
% fprintf(fid,'        constraint=FINER)\n');
fprintf(fid,'        constraint=FIXED)\n');
fprintf(fid,'    p.seedEdgeBySize(edges=edgesTOP, size=Mesh_size, deviationFactor=0.4,\n');
% fprintf(fid,'        constraint=FINER)\n');
fprintf(fid,'        constraint=FIXED)\n');
fprintf(fid,'    p.seedEdgeBySize(edges=edgesBOT, size=Mesh_size, deviationFactor=0.4,\n');
% fprintf(fid,'    constraint=FINER)\n');
fprintf(fid,'    constraint=FIXED)\n');
fprintf(fid,'    p.generateMesh()\n');
%
% Export mesh input file that does not include Periodic Boundary Conditions
fprintf(fid,'    # Export mesh to .inp file\n');
fprintf(fid,'    mdb.Job(name=''include_mesh_FinalRVE_noPBCs'', model=''Model-1'', type=ANALYSIS, explicitPrecision=SINGLE,\n');
fprintf(fid,'        nodalOutputPrecision=SINGLE, description='''',\n');
fprintf(fid,'        parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT,\n');
fprintf(fid,'        numDomains=1, userSubroutine='''', numCpus=1, memory=90,\n');
fprintf(fid,'        memoryUnits=PERCENTAGE, scratch='''', echoPrint=OFF, modelPrint=OFF,\n');
fprintf(fid,'        contactPrint=OFF, historyPrint=OFF)\n');
fprintf(fid,'    mdb.jobs[''include_mesh_FinalRVE_noPBCs''].writeInput(consistencyChecking=OFF)\n');
fprintf(fid,'    #\n');
%
% Include basic information about the sections for each material:
% Copy mesh file to write a new file with Periodic Boundary Conditions
fprintf(fid,'    # Write the Section information for each material\n');
fprintf(fid,'    insert_line = False\n');
fprintf(fid,'    count_lines = 0\n');
fprintf(fid,'    vertex_count = 0\n');
fprintf(fid,'    ncount = len(p.nodes) # Count number of nodes in the Part\n');
fprintf(fid,'    finput=fileinput.input(''include_mesh_FinalRVE_noPBCs.inp'', inplace=1)\n');
fprintf(fid,'    for line in finput:\n');
fprintf(fid,'      if line.startswith(''*Nset, nset=Vertex''):\n');
fprintf(fid,'        vertex_count += 1\n');
fprintf(fid,'        if vertex_count == 4:\n');
fprintf(fid,'          insert_line = True\n');
fprintf(fid,'      else:\n');
fprintf(fid,'        if insert_line:\n');
fprintf(fid,'          if count_lines == 1:\n');
fprintf(fid,'            # START to insert text:\n');
fprintf(fid,'            print ''**''\n');
% Add sections for each material phase:
for jmat_phase=1:nmat_phases
    %
    phase_name = strcat('Phase_',num2str(jmat_phase));
    material_name = strcat('Mat_Phase_',num2str(jmat_phase));
    material_type = mat_phases_IDs(jmat_phase);
    %
        fprintf(fid,'            print ''** Add Section for %s''\n',phase_name);
    if material_type ~= 5 % Then it is not a cohesive section
        fprintf(fid,'            print ''*Solid Section, elset=%s, material=%s, controls=ctrls_%s''\n',phase_name,material_name,phase_name);
        % If you want to include orientation comment line above and uncomment the next 4 lines:
%         fprintf(fid,'            print ''*Solid Section, elset=%s, material=%s, orientation=orient_%s, controls=ctrls_%s''\n',phase_name,material_name,phase_name,phase_name);
%         fprintf(fid,'            print ''*Orientation, Name=orient_%s''\n',phase_name);
%         fprintf(fid,'            print ''1, 0, 0, 0, 1, 0''\n');
%         fprintf(fid,'            print ''3, 0.''\n');
%
    else % It is a cohesive section
        fprintf(fid,'            print ''COHESIVE SECTION, RESPONSE=TRACTION SEPARATION, MATERIAL=%s, ELSET=%s''\n',material_name,phase_name);
        fprintf(fid,'            print ''1.000''\n');
%
    end
end
%
%
% Stop inserting the sections.
fprintf(fid,'            print ''**''\n');
fprintf(fid,'            s = ''** NOTE 1: Global mesh size parameter used as seed is Mesh_size = ''+repr(Mesh_size)\n');
fprintf(fid,'            print s\n');
fprintf(fid,'            s = ''** NOTE 2: Number of attempts to get the mesh niter = ''+repr(niter)\n');
fprintf(fid,'            print s\n');
fprintf(fid,'            print ''**''\n');
fprintf(fid,'            # STOP to insert text.\n');
fprintf(fid,'            insert_line = False\n');
fprintf(fid,'          count_lines += 1\n');
fprintf(fid,'      print line,\n');
fprintf(fid,'    \n');
%
% Copy mesh file to write a new file with Periodic Boundary Conditions
fprintf(fid,'    shutil.copyfile(''include_mesh_FinalRVE_noPBCs.inp'',''include_mesh_FinalRVE.inp'')\n');
% fprintf(fid,'    shutil.copyfile(''%s'',''%s'')\n',...
%     srtcat(jDoE_dir,'include_mesh_FinalRVE_noPBCs.inp'),srtcat(jDoE_dir,'include_mesh_FinalRVE.inp'));
fprintf(fid,'    # Include in ''include_mesh_FinalRVE.inp'' commands for Periodic Boundary Conditions\n');
fprintf(fid,'    insert_line = False\n');
fprintf(fid,'    count_lines = 0\n');
fprintf(fid,'    vertex_count = 0\n');
fprintf(fid,'    ncount = len(p.nodes) # Count number of nodes in the Part\n');
fprintf(fid,'    finput=fileinput.input(''include_mesh_FinalRVE.inp'', inplace=1)\n');
fprintf(fid,'    for line in finput:\n');
fprintf(fid,'      if line.startswith(''*Nset, nset=Vertex''):\n');
fprintf(fid,'        vertex_count += 1\n');
fprintf(fid,'        if vertex_count == 4:\n');
fprintf(fid,'          insert_line = True\n');
fprintf(fid,'      else:\n');
fprintf(fid,'        if insert_line:\n');
fprintf(fid,'          if count_lines == 1:\n');
fprintf(fid,'            # START to insert text:\n');
fprintf(fid,'            print ''**''\n');
fprintf(fid,'            print ''** Added the dummy nodes for PBCs (any coordinates work)''\n');
fprintf(fid,'            print ''*Node, nset=%s'' # create dummy node to apply PBCs on Left & Right edges of RVE\n',dummy1);
fprintf(fid,'            dummy1 = ncount+1;\n');
fprintf(fid,'            s = repr(dummy1)+'', ''+repr(RVEcenter[0]+Lx/2)+'', ''+repr(0.0)\n');
fprintf(fid,'            print s\n');
fprintf(fid,'            print ''*Node, nset=%s'' # create dummy node to apply PBCs on Top & Bottom edges of RVE\n',dummy2);
fprintf(fid,'            dummy2 = ncount+2;\n');
fprintf(fid,'            s = repr(dummy2)+'', ''+repr(0.0)+'', ''+repr(RVEcenter[1]+Ly/2)\n');
fprintf(fid,'            print s\n');
fprintf(fid,'            # Include Constraint equations with the Periodic Boundary Conditions\n');
fprintf(fid,'            print ''**''\n');
fprintf(fid,'            print ''** Include Periodic Boundary Conditions (PBCs)''\n');
fprintf(fid,'            print ''*include, input=include_2D_pbcs.inp''\n');
fprintf(fid,'            print ''**''\n');
fprintf(fid,'            # STOP to insert text.\n');
fprintf(fid,'            insert_line = False\n');
fprintf(fid,'          count_lines += 1\n');
fprintf(fid,'      print line,\n');
fprintf(fid,'    \n');
%
% Periodic Boundary Conditions:
%
%            edge4
%       V4 _________ V3
%         |         |
%   edge1 |         |edge2
%         |         |
%         |_________|
%       V1   edge3  V2
%
fprintf(fid,'    # Periodic Boundary conditions\n');
fprintf(fid,'    def get_node_x(node):\n');
fprintf(fid,'        return node.coordinates[0]\n');
fprintf(fid,'    \n');
fprintf(fid,'    BotEdge_nodes = p.sets[''BotEdge''].nodes\n');
fprintf(fid,'    BotEdge_nodes_sorted = sorted(BotEdge_nodes, key=get_node_x)\n');
fprintf(fid,'    TopEdge_nodes = p.sets[''TopEdge''].nodes\n');
fprintf(fid,'    TopEdge_nodes_sorted = sorted(TopEdge_nodes, key=get_node_x)\n');
fprintf(fid,'    #\n');
fprintf(fid,'    def get_node_y(node):\n');
fprintf(fid,'        return node.coordinates[1]\n');
fprintf(fid,'    \n');
fprintf(fid,'    LeftEdge_nodes = p.sets[''LeftEdge''].nodes\n');
fprintf(fid,'    LeftEdge_nodes_sorted = sorted(LeftEdge_nodes, key=get_node_y)\n');
fprintf(fid,'    RightEdge_nodes = p.sets[''RightEdge''].nodes\n');
fprintf(fid,'    RightEdge_nodes_sorted = sorted(RightEdge_nodes, key=get_node_y)\n');
fprintf(fid,'    #\n');
fprintf(fid,'    # Check if the node count in the Bottom and Top edges is the same:\n');
fprintf(fid,'    if len(BotEdge_nodes_sorted) != len(TopEdge_nodes_sorted):\n');
fprintf(fid,'        status_mesh = 0\n');
fprintf(fid,'    elif len(LeftEdge_nodes_sorted) != len(RightEdge_nodes_sorted):\n');
fprintf(fid,'        status_mesh = 0\n');
fprintf(fid,'    else:\n');
fprintf(fid,'        status_mesh = 1 # Found valid mesh\n');
fprintf(fid,'        # Create include_2D_pbcs.inp file\n');
fprintf(fid,'        BotTop_nodecount = len(BotEdge_nodes_sorted)\n');
fprintf(fid,'        LeftRight_nodecount = len(LeftEdge_nodes_sorted)\n');
fprintf(fid,'        #\n');
fprintf(fid,'        #file to be created\n');
fprintf(fid,'        include_pbcs_file = open(''include_2D_pbcs.inp'', ''w'')\n');
fprintf(fid,'        #writing the entered content to the file we just created\n');
fprintf(fid,'        # Vertex 3 and 1\n');
fprintf(fid,'        for ndof in range(1,3):\n');
fprintf(fid,'          include_pbcs_file.write(''*EQUATION\\n'')\n');
fprintf(fid,'          include_pbcs_file.write('' 4\\n'')\n');
fprintf(fid,'          s = repr(TopEdge_nodes_sorted[len(TopEdge_nodes_sorted)-1].label)+'', ''+repr(ndof)+'', ''+repr(1.0) + '', '' + repr(BotEdge_nodes_sorted[0].label)+'', ''+repr(ndof)+'', ''+repr(-1.0)+'', ''+''%s''+'', ''+repr(ndof)+'', ''+repr(-Lx)+'', ''+''%s''+'', ''+repr(ndof)+'', ''+repr(-Ly)\n',dummy1,dummy2);
fprintf(fid,'          include_pbcs_file.write(s+''\\n'')\n');
fprintf(fid,'        \n');
fprintf(fid,'        # Vertex 2 and 4\n');
fprintf(fid,'        for ndof in range(1,3):\n');
fprintf(fid,'          include_pbcs_file.write(''*EQUATION\\n'')\n');
fprintf(fid,'          include_pbcs_file.write('' 4\\n'')\n');
fprintf(fid,'          s = repr(RightEdge_nodes_sorted[0].label)+'', ''+repr(ndof)+'', ''+repr(1.0) + '', '' + repr(LeftEdge_nodes_sorted[len(LeftEdge_nodes_sorted)-1].label)+'', ''+repr(ndof)+'', ''+repr(-1.0)+'', ''+''%s''+'', ''+repr(ndof)+'', ''+repr(-Lx)+'', ''+''%s''+'', ''+repr(ndof)+'', ''+repr(Ly)\n',dummy1,dummy2);
fprintf(fid,'          include_pbcs_file.write(s+''\\n'')\n');
fprintf(fid,'        \n');
fprintf(fid,'        # Edge 1-2 (Left-Right)\n');
fprintf(fid,'        for enode in range(1,LeftRight_nodecount-1):\n');
fprintf(fid,'          for ndof in range(1,3):\n');
fprintf(fid,'            include_pbcs_file.write(''*EQUATION\\n'')\n');
fprintf(fid,'            include_pbcs_file.write('' 3\\n'')\n');
fprintf(fid,'            s = repr(RightEdge_nodes_sorted[enode].label)+'', ''+repr(ndof)+'', ''+repr(1.0) + '', '' + repr(LeftEdge_nodes_sorted[enode].label)+'', ''+repr(ndof)+'', ''+repr(-1.0)+'', ''+''%s''+'', ''+repr(ndof)+'', ''+repr(-Lx)\n',dummy1);
fprintf(fid,'            include_pbcs_file.write(s+''\\n'')\n');
fprintf(fid,'        \n');
fprintf(fid,'        # Edge 3-4 (Bottom-Top)\n');
fprintf(fid,'        for enode in range(1,BotTop_nodecount-1):\n');
fprintf(fid,'          for ndof in range(1,3):\n');
fprintf(fid,'            include_pbcs_file.write(''*EQUATION\\n'')\n');
fprintf(fid,'            include_pbcs_file.write('' 3\\n'')\n');
fprintf(fid,'            s = repr(TopEdge_nodes_sorted[enode].label)+'', ''+repr(ndof)+'', ''+repr(1.0) + '', '' + repr(BotEdge_nodes_sorted[enode].label)+'', ''+repr(ndof)+'', ''+repr(-1.0)+'', ''+''%s''+'', ''+repr(ndof)+'', ''+repr(-Ly)\n',dummy2);
fprintf(fid,'            include_pbcs_file.write(s+''\\n'')\n');
fprintf(fid,'          \n');
fprintf(fid,'        \n');
fprintf(fid,'        #\n');
fprintf(fid,'        # Select two nodes to avoid rigid body motion\n');
fprintf(fid,'        allnodes = p.nodes\n');
fprintf(fid,'        xx_min = RVEcenter[0]-Lx/2\n');
fprintf(fid,'        xx_max = RVEcenter[0]+Lx/2\n');
fprintf(fid,'        yy_min = RVEcenter[1]-Ly/2\n');
fprintf(fid,'        yy_max = RVEcenter[1]+Ly/2\n');
fprintf(fid,'        ysupportnode = allnodes.getByBoundingBox((0.1*xx_max+1.9*xx_min)/2,(0.1*yy_max+1.9*yy_min)/2,0,(0.3*xx_max+1.7*xx_min)/2,(0.3*yy_max+1.7*yy_min)/2,0)\n');
fprintf(fid,'        xsupportnode = allnodes.getByBoundingBox((1.7*xx_max+0.3*xx_min)/2,(1.7*yy_max+0.3*yy_min)/2,0,(1.9*xx_max+0.1*xx_min)/2,(1.9*yy_max+0.1*yy_min)/2,0)\n');
fprintf(fid,'        #\n');
fprintf(fid,'        include_pbcs_file.write(''**\\n'')\n');
fprintf(fid,'        include_pbcs_file.write(''** Inlcude sets with two nodes within the RVE to restrict rigid body motion:\\n'')\n');
fprintf(fid,'        include_pbcs_file.write(''*Nset, nset=ySupport\\n'')\n');
fprintf(fid,'        include_pbcs_file.write(str(ysupportnode[0].label)+''\\n'')\n');
fprintf(fid,'        include_pbcs_file.write(''*Nset, nset=xSupport\\n'')\n');
fprintf(fid,'        include_pbcs_file.write(str(xsupportnode[0].label)+''\\n'')\n');
fprintf(fid,'        #\n');
fprintf(fid,'        #writing task is completed\n');
fprintf(fid,'        include_pbcs_file.close()\n');
fprintf(fid,'        \n');
fprintf(fid,'    # end if (comparing number of nodes on the edges)\n');
fprintf(fid,'    #\n');
fprintf(fid,'    # If we have too many iterations then print error to ERROR file\n');
fprintf(fid,'    if niter <= %i:\n',niter);
fprintf(fid,'        niter = niter + 1\n');
fprintf(fid,'        Mesh_size = Mesh_size/refine_factor # refine mesh\n');
fprintf(fid,'    else:\n');
fprintf(fid,'        status_mesh = 2 # Did not find valid mesh...\n');
fprintf(fid,'        ERROR_file = open(''ERROR_FILE'', ''a'')\n');
fprintf(fid,'        #writing the entered content to the end of the ERROR_FILE\n');
fprintf(fid,'        ERROR_file.write(''Failed to mesh the RVE for this sample\\n'')\n');
fprintf(fid,'        ERROR_file.close()\n');
fprintf(fid,'    # end if niter\n');
fprintf(fid,'\n');
fprintf(fid,'# end while status_mesh = 0\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'# End of python script\n');
%
%
fclose(fid);
%
% disp(' ');
% disp('Creation of Python file COMPLETED');
% disp('Elapsed Time [min]: ');
% disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Python Script in ABAQUS CAE                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Started generation of mesh files'); disp(' ');
unix(strcat(abaqus_path,' cae noGUI=',jDoE_dir,'/script_RVE_meshing.py'));
% disp(' ');
% disp('Generation of mesh COMPLETED');
% disp('Elapsed Time [min]: ');
% disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
end
