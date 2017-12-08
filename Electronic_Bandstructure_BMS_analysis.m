%% 3D_Electronic_Bandstructure_FE_Code
% ======================================================================= %
% Dimitri Krattiger
% 10-12-2012

close all; clear; clc; 
format compact
format shortg
warning off all
profile clear
profile on
tic

addpath('Global_Files')

%% Check for & Create Save directories
% ======================================================================= %

save_results = true;
load_results = true;

% if save_data folder does not exist, make it
if save_results
    if ~exist('save_data','dir')
        mkdir('save_data');
    end
    if ~exist('save_data/models','dir')
        mkdir('save_data/models');
    end
    if ~exist('save_data/solutions','dir')
        mkdir('save_data/solutions');
    end
end

%% User Selected Parameters
% ======================================================================= %

% pattern is a 3D array. The number of terms in the 1st, 2nd, and 3rd
% indices correspond to the number of elements in the x,y, and z directions


% pat_size = 1;
% pat_size = 2;
pat_size = 3;
% pat_size = 4;
% pat_size = 5;
% pat_size = 8;
% pat_size = 10;
% pat_size = 20;

n=6;

% frac = pat_size*2/3;
pattern = 4*ones(pat_size,pat_size,pat_size);

% Solution Options
n_curves     = 8;       % number of dispersion curves to plot
full_disp    = true;  	% 1 = compute full dispersion, 0 = don't
WFE          = true; 	% 1 = Bloch Boundary Conditions, 0 = Bloch Operator
RBME_disp    = false;  	% Compute RBME dispersion
BMS_disp     = true;    % compute BMS dispersion
plot_disp    = false;

% n_FIs = [10,15,20,30,50,100,200,500,1000,1500];
% n_FIs = [10,20,50,100,200,500,1000];
% n_FIs = [10,20,50];
% n_FIs = [0:50];
% n_FIs = [1:10,15,20,30,40,50,75,100,150,250,500]
% n_FIs = [1:5];
n_FIs = [10,20,50,100];
% n_FIs = [20];

plot_disp = true;
n_FIs = 10;

% size of k-point grid for monkhorst pack
use_ibz = false;
if ~plot_disp
    n_kap_1 = 10;
%     n_kap_1 = 12;  
%     n_kap_1 = 16; 
%     n_kap_1 = 20;  
    use_ibz = false;
end

% Mode Animations
% mode_plot = [1];     % select modes to plot. [] =  don't animate any modes
k_sel = 1;          % k-point to plot modes at


% units are in Ry to begin with (from pseudopotential)
eV_Per_Ha = 27.2114; 
eV_Per_Ry = eV_Per_Ha/2;
Ha_Per_Ry = 1/2;


cmap = viridis;
colormap(cmap)

% number of atoms in base primitive cell
n_atoms_base = 2;

%% Define Model dimensions
% ======================================================================= %
n_cellx = 1;        % number of repetitions of unit cell in x-direction
n_celly = 1;        % number of repetitions of unit cell in y-direction
n_cellz = 1;        % number of repetitions of unit cell in z-direction
pattern = repmat(pattern,[n_cellx,n_celly,n_cellz]);


% n = 2;              % number of element nodes per edge
n_bz = 1;           % number of brillouin Zones to cover
n_atoms = n_atoms_base*n_cellx*n_celly*n_cellz;

unit_select = 'atomic_units';
% length_select = 'meters';
switch unit_select
    case 'angstroms'
        a0 = 5.431;         % angstroms
    case 'meters'
        a0 = 5.431e-10;     % meters
    case 'nanometers'
        a0 = 0.543;         % nanometers
    case 'atomic_units'
        a0 = 10.261;        % atomic units     
        % THIS IS THE ONE THAT WORKS CURRENTLY
end

plot_units = 'nanometers';
switch plot_units
    case 'angstroms'
        conversionFactor = 5.431/a0;
        plot_unit_string = '\r{A}';
    case 'meters'
        conversionFactor = 5.431e-10/a0;
        plot_unit_string = 'm';
    case 'nanometers'
        conversionFactor = 0.543/a0;
        plot_unit_string = 'nm';
    case 'atomic_units'
        conversionFactor = 10.261/a0;
        plot_unit_string = 'a.u.';     
end



% lattice vectors
a1 = (a0/2)*[0;1;1]*n_cellx;
a2 = (a0/2)*[1;0;1]*n_celly;
a3 = (a0/2)*[1;1;0]*n_cellz;
R = [a1,a2,a3];

% reciprocal lattice vectors
b1 = 2*pi*cross(a2,a3)/dot(a1,cross(a2,a3));
b2 = 2*pi*cross(a3,a1)/dot(a2,cross(a3,a1));
b3 = 2*pi*cross(a1,a2)/dot(a3,cross(a1,a2));

%% Define Potential
% ======================================================================= %

% potential 1
V1 = 0*6.5;

% potential 2
V2 = 1*6.5; 

% potential 3
V3 = 2*6.5; 

% potential 4
V4 = 3*6.5; 

% silicon pseudo-potential function handle
Vhandle_Si_Pseudo = @(r)Si_Pseudo(r,a0);

% collect possible potentials into a cell array
Vs = {V1,V2,V3,V4,Vhandle_Si_Pseudo};

%% Discretize Brillouin zone
% ======================================================================= %

if plot_disp
    kp = 2*pi/(a0)*1;
    GammaPt     = [0;0;0];
    XPt         = [1;0;0]*kp;
    XprimePt    = [1;1;0]*kp;
    LPt         = [1/2;1/2;1/2]*kp;
    KPt         = [3/4;3/4;0]*kp;
    WPt         = [1;1/2;0]*kp;
    UPt         = [1;1/4;1/4]*kp;

    % KappaPts = [GammaPt,XPt,MPt,RPt,GammaPt];
    KappaPts = [LPt,GammaPt,XPt,nan(3,1),XprimePt,KPt ,GammaPt];
%     KappaPts = [LPt,GammaPt,XPt,LPt,KPt,GammaPt];

    % determine number of BZ segments
    i_real = ~isnan(KappaPts(1,:));
    n_segs = sum(i_real(1:end-1) & i_real(2:end));

    % n_segs = sum(~isnan(KappaPts(1,:)))-1;
    n_kap = (n_segs)*2^4+1;    % number of kappa points

    kappa = zeros(3,n_kap);
    kappa_plot = zeros(1,n_kap);
    count = 0;
    for i = 1:size(KappaPts,2)-1
        if i_real(i) && i_real(i+1)
            count = count+1;
            n_seg_pts = (n_kap-1)/n_segs;
            i_kap = ((count-1)*n_seg_pts+1):(count*n_seg_pts+1);
            kappa(:,i_kap) = ...
                KappaPts(:,i)*ones(1,n_seg_pts+1)+(KappaPts(:,i+1)-KappaPts(:,i))*...
                linspace(0,1,n_seg_pts+1);
            segment_length = sqrt(sum(KappaPts(:,i+1)-KappaPts(:,i)).^2);
            kappa_plot(i_kap) = kappa_plot(i_kap(1))+linspace(0,segment_length,n_seg_pts+1);
        end
    end
    
    kappa_use = kappa;    
    n_kap = size(kappa_use,2);
    w_use = ones(1,size(kappa_use,2))/size(kappa_use,2);
else
    
%     [kappa2,i_ibz,w_ibz] = MonkhorstPackGrid(b1,b2,b3,n_kap_1);

% profile clear
% profile on
    %n_kap_1 = 10;

    n_MP = n_kap_1*ceil([1/n_cellx,1/n_celly,1/n_cellz]);
%     n_MP = 10*ceil([1/n_cellx,1/n_celly,1/n_cellz]);
%     n_MP = n_kap_1*[1,1,1];
%     [kappa2,i_ibz,w_ibz] = MonkhorstPackGrid(b1,b2,b3,n_MP);
    [~,w_ibz,kappa2,i_ibz,~] = MonkhorstPackGrid(b1,b2,b3,n_MP);
    length(i_ibz)
    n_kap_1^3/length(i_ibz)
    
    if use_ibz
        kappa_use = kappa2(:,i_ibz);
        w_use = w_ibz;
    else
        kappa_use  = kappa2;
        w_use = ones(1,size(kappa_use,2))/size(kappa_use,2);
    end
    
    n_kap = size(kappa_use,2);
    
    %%
    figure(889);clf
    xpatchBZ1 = 2*pi/a0 *	[1,   1,    1,    1;...
                        -1,  -1,   -1,   -1;...
                        1/2,   0,      -1/2,   0;...    
                        1/2,   0,      -1/2,   0;...
                        0,     1/2,    0,      -1/2;...    
                        0,     1/2,    0,      -1/2];%;...
    
        
    ypatchBZ1 = 2*pi/a0 *	[1/2,   0,      -1/2,   0;...    
                        1/2,   0,      -1/2,   0;...
                        1,   1,    1,    1;...
                        -1,  -1,   -1,   -1;...
                        1/2,   0,      -1/2,   0;...    
                        1/2,   0,      -1/2,   0];%;...   
    
    
    zpatchBZ1 = 2*pi/a0 *	[0,     1/2,    0,      -1/2;...    
                        0,     1/2,    0,      -1/2;...  
                        0,     1/2,    0,      -1/2;...    
                        0,     1/2,    0,      -1/2;...
                        1,   1,    1,    1;...
                        -1,  -1,   -1,   -1];%;... 
    
    xpatchBZ2 = 2*pi/a0 *	[1,   1/2,    0,      0,      1/2,    1;...
                        -1,  -1/2,   0,      0,      -1/2,   -1;...
                        -1,  -1/2,   0,      0,      -1/2,   -1;...
                        1,   1/2,    0,      0,      1/2,    1;...
                        1,   1/2,    0,      0,      1/2,    1;...
                        -1,  -1/2,   0,      0,      -1/2,   -1;...
                        -1,  -1/2,   0,      0,      -1/2,   -1;...
                        1,   1/2,    0,      0,      1/2,    1];%;...
                 
    
    ypatchBZ2 = 2*pi/a0 *	[1/2,   1,   	1, 	1/2, 	0,      0;...
                    	1/2,   1,   	1, 	1/2, 	0,      0;...
                      	-1/2,  -1,   -1, 	-1/2, 	0,      0;...
                        -1/2,  -1,   -1, 	-1/2, 	0,      0;...
                        1/2,   1,   	1, 	1/2, 	0,      0;...
                        1/2,   1,   	1, 	1/2, 	0,      0;...
                        -1/2,  -1,   -1, 	-1/2, 	0,      0;...
                        -1/2,  -1,   -1, 	-1/2, 	0,      0];%;...
        
    
    zpatchBZ2 = 2*pi/a0 *	[0,     0,      1/2,      1, 	1,    1/2;...
                        0,     0,      1/2,      1, 	1,    1/2;...
                        0,     0,      1/2,      1, 	1,    1/2;...
                        0,     0,      1/2,      1, 	1,    1/2;...
                        0,     0,      -1/2,     -1, -1,   -1/2;...
                        0,     0,      -1/2,     -1, -1,   -1/2;...
                        0,     0,      -1/2,     -1, -1,   -1/2;...
                        0,     0,      -1/2,     -1, -1,   -1/2];%;...
                    
    xpatchBZ3 = 2*pi/a0 *	[1,     1      1;...
                             0,     1/2    3/4];
    ypatchBZ3 = 2*pi/a0 *	[0,     1/2    1/4;...
                             0,     1/2    3/4];
    zpatchBZ3 = 2*pi/a0 *	[0,     0      1/4;...
                             0,     1/2    0];
                         
     xpatchBZ4 = 2*pi/a0 *	[1,     0,     3/4,      1;...
                             1/2,   1,      1,      3/4;...
                             0,     1,      1,      1/2];
     ypatchBZ4 = 2*pi/a0 *	[0,     0,     3/4,      1/2;...
                             1/2,   1/4,    1/2,     3/4;...
                             0,     0,      1/4,     1/2];
     zpatchBZ4 = 2*pi/a0 *	[0,     0,     0,       0;...
                             1/2,   1/4,    0,      0;...
                             0,     0,      1/4,    1/2];
    view(3)
    axis equal
    %camproj('perspective')
    
%     plot3(kappa2(1,:),kappa2(2,:),kappa2(3,:),'k.');hold on
    plot3(kappa2(1,i_ibz),kappa2(2,i_ibz),kappa2(3,i_ibz),'k.');hold on
    h1 = patch(xpatchBZ1',ypatchBZ1',zpatchBZ1','r','facealpha',0.5);hold on%,'edgecolor','none')
    h2 = patch(xpatchBZ2',ypatchBZ2',zpatchBZ2','b','facealpha',0.5);hold on%,'edgecolor','none')
    h3 = patch(xpatchBZ3',ypatchBZ3',zpatchBZ3','g','facealpha',0.5);hold on%,'edgecolor','none')
    h4 = patch(xpatchBZ4',ypatchBZ4',zpatchBZ4','g','facealpha',0.5);hold on%,'edgecolor','none')
    B = [b1,b2,b3];
    plot3([[0,0,0];B(1,:)],[[0,0,0];B(2,:)],[[0,0,0];B(3,:)],'k-')
    axis off
    axis equal
    
    i_in_BZ = (abs(kappa_use(1,:)) +  abs(kappa_use(2,:)) + abs(kappa_use(3,:)) < 3*pi/a0)  & ...
              (abs(kappa_use(1,:)) < 2*pi/a0)& ...
              (abs(kappa_use(2,:)) < 2*pi/a0)& ...
              (abs(kappa_use(3,:)) < 2*pi/a0);
    xlabel('\it{k_x}')
    ylabel('\it{k_y}')
    zlabel('\it{k_z}')
    
%     camlight
%     set(h1,'facelighting','gouraud')
%     set(h2,'facelighting','gouraud')
    view(120,30)
    grid on
    
end

% define mesh
eleorder = n-1;
[coordinates,elenodes] = MeshParallelPiped(pattern,eleorder,R);

n_eles = size(elenodes,1);

%% Create patches for element faces
% ======================================================================= %

bottom = [1:n-1,...
          n:n:n^2-n,...
          n^2:-1:n^2-n+2,...
          n^2-n+1:-n:n+1];   
top = (n-1)*n^2+fliplr(bottom);

% front and back side patches
side_f = [1:n:n^2-2*n+1,...
          n^2-n+1:n^2:n^3-n^2-n+1,...
          n^3-n+1:-n:n^3-n^2+n+1,...
          n^3-n^2+1:-n^2:n^2+1];
side_b = n-1+fliplr(side_f);

% left and right side patches
side_l = [1:n^2:n^3-2*n^2+1,...
          n^3-n^2+1:n^3-n^2+n-1,...
          n^3-n^2+n:-n^2:n^2+n,...
          n:-1:2];
side_r = n^2-n+fliplr(side_l);

six_patches = [bottom;top;side_f;side_l;side_b;side_r];

% patchx = [];    patchy = [];    patchz = [];

faces = zeros(n_eles*6,4);
for i = 1:n_eles
    for j = 1:size(six_patches,1);
        patch_faces((i-1)*6+j,:) = elenodes(i,six_patches(j,:));
    end
end

%% plot isosurface contours for potential
% ======================================================================= %

if true
    
    % 3d grid in lattice vector coordinates
    ncells1 = 1;
    ncells2 = 1;
    ncells3 = 1;

    ndisc = 40;
    ndisc1 = ndisc*ncells1;
    ndisc2 = ndisc*ncells2;
    ndisc3 = ndisc*ncells3;

    a1vec2 = linspace(0,ncells1,ndisc1);
    a2vec2 = linspace(0,ncells2,ndisc2);
    a3vec2 = linspace(0,ncells3,ndisc3);
%     [a1grid2,a2grid2,a3grid2] = meshgrid(a1vec2,a2vec2,a3vec2);
    [a1grid2,a2grid2,a3grid2] = ndgrid(a1vec2,a2vec2,a3vec2);

    % 3d grid in cartesian coordinates
    xgrid2 = ((a1grid2*a1(1) + a2grid2*a2(1) + a3grid2*a3(1)));
    ygrid2 = ((a1grid2*a1(2) + a2grid2*a2(2) + a3grid2*a3(2)));
    zgrid2 = ((a1grid2*a1(3) + a2grid2*a2(3) + a3grid2*a3(3)));

    C2 = zeros(size(xgrid2));
    for i  = 1:numel(C2)
        C2(i) = Si_Pseudo([xgrid2(i);ygrid2(i);zgrid2(i)],a0);
    end
    % convert from Ry to Ha
    C2 = C2*(1/2);


    figure(1);clf;hold on
    h0 = patch('faces',patch_faces,'vertices',coordinates*conversionFactor);
    set(h0,'FaceColor','none',...
           'EdgeColor','none',...
           'edgealpha',0.5);

    maxC = max(C2(:));
    minC = min(C2(:));
    n_isosteps = 5;
    isosteps = linspace(0,1,n_isosteps+2);isosteps = isosteps(2:end-1); 
    isovals = minC + isosteps*(maxC-minC);

    legend_entries = cell(1,length(isosteps));
    h = zeros(1,length(isosteps));
    for i = 1:length(isosteps)
        isoSurfPatchStruct = isosurface(conversionFactor*xgrid2,...
                                        conversionFactor*ygrid2,...
                                        conversionFactor*zgrid2,...
                                        C2,isovals(i));
        h(i) = patch(isoSurfPatchStruct);
        i_color = round(1+(i-1)*(size(cmap,1)-1)/(length(isosteps)-1));
        set(h(i),'FaceColor',cmap(i_color,:),...
                 'EdgeColor','none',...
                 'facealpha',isosteps(i)^(0));
             
         roundedval = round(isovals(i)*1e2)/1e2;
         if roundedval == 0
             roundedval = 0;    % this seems pointless, but it is 
                                % necessary to avoid -0.00 in legend
         end
        legend_entries{i} = sprintf('%3.2f Ha',roundedval);
    end
    
    % plot unit cell edges
    xedge = a0*conversionFactor*[0, 1/2, 1/2, 0,   0, nan, 1/2, 1,   1, 1/2, 1/2, nan, 0, 1/2, nan, 1/2, 1,    nan, 1/2, 1, nan, 0,   1/2];
    yedge = a0*conversionFactor*[0, 0,   1/2, 1/2, 0, nan, 1/2, 1/2, 1, 1,   1/2, nan, 0, 1/2, nan, 0,   1/2,  nan, 1/2, 1, nan, 1/2, 1];
    zedge = a0*conversionFactor*[0, 1/2, 1,   1/2, 0, nan, 0,   1/2, 1, 1/2, 0,   nan, 0, 0,   nan, 1/2, 1/2,  nan, 1,   1, nan, 1/2, 1/2];
    h_boundary = plot3(xedge,yedge,zedge,'k-','linewidth',1);hold on
    
    %Add legend
    legend(h(end:-1:1),legend_entries(end:-1:1))
    
    daspect([1,1,1])
    view([5,20]); 
    view(3);     
    view([-75,15])
    
    axis tight 
    axis equal 
    
    grid on
    
    
    xlabel(['x (',plot_unit_string,')']);    
    ylabel(['y (',plot_unit_string,')']);
    zlabel(['z (',plot_unit_string,')']);
    set(gcf,'color','w')
    

    camlight        
    camlight left
    lighting gouraud

    drawnow
end


%% Form Master Mass and Stiffness
% ======================================================================= %
        
% define model save path
pathstring = 'save_data/models/';

% define element order string
orderstrings = {'1st','2nd','3rd','4th','5th',...
            '6th','7th','8th','9th','10th'};                    
orderstring = orderstrings{n-1};  
n_dof_per = prod(size(pattern)*(n-1));

% full model description string
modeldescription = [sprintf('Si_%ix%ix%iPrimCell_%iDOF_',...
                n_cellx,n_celly,n_cellz,n_dof_per),...
                orderstring,'OrderEles'];

model_savestring = [pathstring,modeldescription,'_FreeModel'];% form master mass and stiffness matrices
if exist([model_savestring,'.mat'],'file') && load_results
    load([model_savestring,'.mat'])

else
    % store mesh info in a structure
    Mesh_info.pattern_vec   = pattern(:);
    Mesh_info.Vs            = Vs;
    Mesh_info.n             = n;
    Mesh_info.elenodes      = elenodes;
    Mesh_info.coordinates   = coordinates;

    % compute system matrices
    WFE = true;
    if WFE
        [H0f,Sf] = ...
            master_matrices_3D_Electronic(Mesh_info);
    else
        [H0f,Hxf,Hyf,Hzf,Sf] = ...
            master_matrices_3D_Electronic_Bloch_Operator(kappa_use,Mesh_info);
    end

    H0f = sparse(H0f);
    Sf = sparse(Sf);
    
    if save_results
        save(model_savestring,'H0f','Sf','coordinates')
    end
end

%% Create Degree-of-freedom sets
% ======================================================================= %

% coordinates = [xlocs,ylocs,zlocs];
% R = [a1,a2,a3];
node_sets = find_node_sets(coordinates,R);
dof_sets = node_sets;
T_per = Periodic_Boundary_Conditions(dof_sets);

%% spy matrix
% ======================================================================= %
i_i = 1:length(dof_sets.i);
i_b = length(dof_sets.i):n_dof_per;
% dof_set_fields_per = {'l','f','d','lf','fd','ld','lfd'};

% for q = [];%1:n_kap
% wave vector at current k-point
kvec = kappa_use(:,1);

% phase modification at boundaries due to wavevector
lam  = exp(-1i*kvec'*R);
lam(end+1:3) = 0;

% Obtain full model result
T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
        + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
        + T_per.s123*lam(1)*lam(2)*lam(3);

Hhat = T_per_k'*H0f*T_per_k;

Hii = Hhat; Hii(i_b,:) = 0; Hii(:,i_b) = 0;
Hbb = Hhat; Hbb(i_i,:) = 0; Hbb(:,i_i) = 0;

[i_all,j_all,Hhatvec] = find(flipud(Hhat([i_i,i_b],[i_i,i_b])));
[i_all_ii,j_all_ii,Hiivec] = find(flipud(Hii([i_i,i_b],[i_i,i_b])));
[i_all_bb,j_all_bb,Hbbvec] = find(flipud(Hbb([i_i,i_b],[i_i,i_b])));

figure(222);clf
Hhatcolors = (abs(imag(Hhatvec))./abs(Hhatvec))*[1,1,1]*0.9;
Hhatcolors = [1,1,1]*0.75;
scatter(j_all,i_all,10,Hhatcolors,'.');hold on

Hiicolors = (abs(imag(Hiivec))./abs(Hiivec))*[1,1,1]*0.9; Hiicolors(:,3) = 1;
Hiicolors = [0,0,1];
scatter(j_all_ii,i_all_ii,10,Hiicolors,'.');hold on

Hbbcolors = (abs(imag(Hbbvec))./abs(Hbbvec))*[1,1,1]*0.9; Hbbcolors(:,1) = 1;
Hbbcolors = [1,0,0];
scatter(j_all_bb,i_all_bb,10,Hbbcolors,'.');hold on

xlim([0,size(Hhat,1)]);
ylim([0,size(Hhat,2)]);
axis equal
axis off
drawnow
rez = 150;
% saveas(gcf,sprintf('Hamiltonian_Sparsity_kpt_%i_of_%i.png',q,n_kap))


%%

ib_groups_color1 = 'l';
ib_groups_color2 = 'r';
color1 = [0,0,1];
color2 = [0,.25,1];

ib_groups_color1 = 'f';
ib_groups_color2 = 'b';
color1 = [0,1,0];
color2 = [0,1,.25];


ib_groups_color1 = 'd';
ib_groups_color2 = 't';
color1 = [1,0,0];
color2 = [1,0,.25];


i_b_color1 = [];
i_b_color2 = [];     
dof_set_names = fieldnames(dof_sets);
for i = 1:length(dof_set_names)
    if any(dof_set_names{i}==ib_groups_color1)
        i_b_color1 = [i_b_color1,dof_sets.(dof_set_names{i})];
    elseif any(dof_set_names{i}==ib_groups_color2)
        i_b_color2 = [i_b_color2,dof_sets.(dof_set_names{i})];
    end
end

if true
i_b_color = i_b;
i_b_black = true(1,size(coordinates,1));
i_b_black(i_b_color) = false;
else
% i_b_color1 = i_b

i_b_color = [i_b_color1,i_b_color2];
i_b_black = setdiff(i_b,i_b_color);
end

[~,i_keep_patch_faces] = unique(sort(patch_faces,2),'rows','stable');
patch_faces2 = patch_faces(i_keep_patch_faces,:);


i_patch_color = all(ismember(patch_faces2,i_b_color),2);

figure(11);clf;hold on
h0 = patch('faces',patch_faces2(~i_patch_color,:),'vertices',coordinates*conversionFactor);
set(h0,'FaceColor','b',...
       'EdgeColor','k',...
       'edgealpha',1,...
       'facealpha',0.5);
h1 = patch('faces',patch_faces2(i_patch_color,:),'vertices',coordinates*conversionFactor);
set(h1,'FaceColor',color1,...
       'EdgeColor','k',...
       'edgealpha',1,...
       'facealpha',0.5);

plot3(coordinates(i_b_black,1)*conversionFactor,coordinates(i_b_black,2)*conversionFactor,coordinates(i_b_black,3)*conversionFactor,...
    'marker','.','markerfacecolor','w','color','k','linestyle','none');hold on
plot3(coordinates(i_b_color,1)*conversionFactor,coordinates(i_b_color,2)*conversionFactor,coordinates(i_b_color,3)*conversionFactor,...
    'marker','.','markerfacecolor','w','color',color1,'linestyle','none')
% plot3(coordinates(i_b_color2,1)*conversionFactor,coordinates(i_b_color2,2)*conversionFactor,coordinates(i_b_color2,3)*conversionFactor,...
%     'marker','.','markerfacecolor','w','color',color2,'linestyle','none')

view([-75,15])
axis equal
axis off
zlim([0.15,0.543])
zlim([0,0.4])
zlim([0.4,0.543])
% pause

%% Dispersion Solution
% ======================================================================= %

if full_disp == 1

%     pathstring = 'save_data/solutions/';
%     modelstring = sprintf('%i_DOF_%i_NodesPerEle_%i_',...
%                     n_dof_per,n^3,n_kap);
% 
%     % BZ or IBZ in file names
%     if use_ibz;     
%         BZ_or_IBZ = 'IBZ';  
%     else
%         BZ_or_IBZ = 'BZ';  
%     end
%     
%     full_soln_savestring = [pathstring,modelstring,BZ_or_IBZ,'kpts_Silicon_Full_Soln'];
%     if exist([full_soln_savestring,'.mat'],'file') && load_results
%         load([full_soln_savestring,'.mat'],'E_disp','t_kloop','PHI')
% disp(['element order = ',num2str(n-1),...
%       ', n_ele_edge = ',num2str(n_ele_edge),...
%       ', n_kap_edge = ',num2str(n_kap_edge)])
              
    % BZ or IBZ in file names
    if use_ibz;     
        BZ_or_IBZ = 'IBZ';  
    else
        BZ_or_IBZ = 'BZ';  
    end

    % define model save path
    solutionpathstring = 'save_data/solutions/';

    solutiondescription = [sprintf('%i',n_kap),BZ_or_IBZ,sprintf('kpts_%iBands',n_curves)];

    solution_savestring = [solutionpathstring,modeldescription,'_',solutiondescription];     

    if exist([solution_savestring,'.mat'],'file') && load_results
%                     load([full_soln_savestring,'.mat'])
        load([solution_savestring,'.mat'])
    else
        t_loop = 0;
        
        solve_options.n_curves = n_curves;
        [omega,PHI,t_kloop] = dispersion_solver_w_k(kappa_use,H0f,Sf,dof_sets,R,solve_options);
        E_disp = omega.^2;
        
        if save_results
            save(solution_savestring,'E_disp','PHI','t_kloop')
        end
        
    end
    t_full = sum(t_kloop);
end

for i = 1:n_kap
    [E_disp(:,i),isort] = sort(E_disp(:,i));
    PHI(:,:,i) = PHI(:,isort,i);
end

%% Plot full model dispersion
% ======================================================================= %
if plot_disp
    figure(87);clf
    plot(kappa_plot,E_disp(1:n_curves,:),'.-k','linewidth',1);
    % Plot Dashed lines at high symmetry points
    hold on
    i_high_sym = linspace(1,n_kap,n_segs+1);

    % tick_labels = {'Gamma';'X';'M';'R';'Gamma'};
    % set('gca','xtick',
    tick_labels = {'L';'\Gamma';'X';'K';'\Gamma'};
    % plot([1;1]*kappa_plot(i_high_sym),ylimits'*ones(1,n_segs+1),'k:','linewidth',linewidths)
    set(gca,'xtick',kappa_plot(i_high_sym))
    set(gca,'xticklabels',tick_labels);
    set(gca,'xgrid','on')
    xlim(kappa_plot([1,end]))
    hold off

    % Create axis labels and figure title
    ylabel('Energy (Ha)');xlabel('wavenumber \kappa')
end

%% choose a model to use for error reference
% ======================================================================= %
use_kill_model_soln_as_reference = false;
if use_kill_model_soln_as_reference
%     ref_model = ['save_data/Bilal_Run/','15625_DOF_216NodesPerEle_91kpts_Silicon_Full_Soln'];
    ref_model = ['save_data/solutions/','Si_1x1x1PrimCell_15625DOF_5thOrderEles_47IBZkpts_8Bands'];
    ref_soln = load(ref_model);
    
    E_ref = ref_soln.E_disp;
    PHI_ref = ref_soln.PHI;
else
    E_ref = E_disp;
    PHI_ref = PHI;
end

% % convert from Ry to Ha
% ref_soln.E_disp = ref_soln.E_disp*2;


%% Do some plotting of unique solutions vs k-points
% ======================================================================= %

if false
    % determine unique energies
    comparison_type = 'Modal-esque';

    switch comparison_type
        case 'direct_rounded'
            [~,ia,~] = unique(round(sort(real(E_disp)),2)','rows'); 

        case 'Modal-esque'
            testmat = E_disp*diag(diag(E_disp'*E_disp).^(-1/2));
            testmat = 1-testmat'*testmat<1e-8;

            [~,ia,~] = unique(testmat,'rows');
    end



    rhotemp = zeros(n_dof_per,n_kap);
    % for i = 1:4
    for j = 1:n_kap
        PHItemp = squeeze(PHI(:,1:4,j));
        [PHItemp,Rmat] = qr(PHItemp,0);
        PHItemp = PHItemp*diag(diag(PHItemp'*PHItemp).^(-1/2));
        rhoj = sum(PHItemp.*conj(PHItemp),2);
        rhotemp(:,j) = rhotemp(:,j) + rhoj;
    end

    disp(['number of unique energy values (within some tolerance): ',num2str(length(ia))]);
    disp(['number of unique points predicted by IBZ code: ',num2str(length(i_ibz))]);


    plotit = false;
    if plotit
        for i = 1:length(ia)
            ibsets{i} = find(ib==i);
        end

        Rb = [b1,b2,b3];
        Rb = eye(3);
        gvals = Rb\kappa_use;e
        figure(88);clf
        % plot3(kappa(1,:),kappa(2,:),kappa(3,:),'r.');hold on
        plot3(gvals(1,:),gvals(2,:),gvals(3,:),'k.');hold on
        % i_unique = [87,    88,    51,    44,    45,    52,    58,     8,    46,    15,     9,     2,     1,    53,    16,    22,     3,...
        %      6,    10,    29,    47,    23,    12,    17,     5,     4,    18,    11];
        plot3(gvals(1,ia),gvals(2,ia),gvals(3,ia),'go','linewidth',1.5);
        plot3(gvals(1,i_ibz),gvals(2,i_ibz),gvals(3,i_ibz),'m*','linewidth',1.5);
        for i = 1:n_kap
            text(gvals(1,i),gvals(2,i),gvals(3,i),num2str(i))
        end
        axis equal

        % plot3([0,b1(1)],[0,b1(2)],[0,b1(3)],'m-o','linewidth',4)
        % plot3([0,b2(1)],[0,b2(2)],[0,b2(3)],'m-o','linewidth',4)
        % plot3([0,b3(1)],[0,b3(2)],[0,b3(3)],'m-o','linewidth',4)
        for i = 1:length(ia)
            h1 = plot3(gvals(1,ibsets{i}),gvals(2,ibsets{i}),gvals(3,ibsets{i}),'r*-');
            pause
            delete(h1)
        end
    end
end


%% Use BMS reduction to compute band structure
% ======================================================================= %

% w_cut = -1;
% w_cut = max(max(sqrt(E_disp)))*1.5;
% n_LI = 70;
n_n_FIs = length(n_FIs);

E_BMS = zeros(n_curves,n_kap,n_n_FIs);
% PHI_BMS = zeros(n_dof,n_curves,n_kap,n_n_FIs);
t_BMS = zeros(n_n_FIs,1);
n_dof_BMS = zeros(n_n_FIs,1);

if BMS_disp
    for i = 1:n_n_FIs

        clear options_BMS
        options_BMS.n_FI = n_FIs(i);
        options_BMS.InteriorMethod      = 'CB+';
        
%         options_BMS.BoundaryMethod      = 'weak';
%         options_BMS.preOrthoTypeLIRWeak = 'qr';
%         options_BMS.svdTolLIRWeak       = 1e-1;       
%         options_BMS.n_CC                = 100;

        options_BMS.BoundaryMethod      = 'exact';
%         options_BMS.BoundaryMethod      = 'none';
        %options_BMS.orthoTypeExact      = 'svd';
        options_BMS.n_CC                = 40;
        options_BMS.n_CC                = 35;

        t_start = tic;
        
        [H_BMS{i},S_BMS{i},dof_sets_BMS{i},t_up_front,T_BMS{i}] = ...
                BMS(H0f,Sf,coordinates,R,options_BMS);

%         [H_BMS,S_BMS,dof_sets_BMS,t_up_front,T_BMS{i}] = ...
%                 BMS(H0f,Sf,coordinates,w_cut,R,n_FI,n_LI);
            
        solve_options.n_curves = n_curves;
        [omega_BMS,phi_BMS,t_kloop_BMS] = ...
            dispersion_solver_w_k(kappa_use,H_BMS{i},S_BMS{i},dof_sets_BMS{i},R,solve_options);
        
        
        t_BMS(i) = toc(t_start);

        T_per_BMS{i} = Periodic_Boundary_Conditions(dof_sets_BMS{i});
        n_dof_BMS(i) = size(T_per_BMS{i}.s0,2);
        
%         E_BMS = omega_BMS.^2;
        PHI_BMS{i} = phi_BMS;
        E_BMS(:,:,i) = omega_BMS.^2; % convert from Ry to Ha
%         E_BMS(:,:,i) = omega_BMS.^2*eV_Per_Ry;  % convert from Ry to eV
        disp([num2str(n_FIs(i)),' fixed interface modes, ',...
              num2str(t_BMS(i)),' s'])
    end
end


%% Plot Sparsity diagram for BMS model
% wave vector at current k-point
kvec = kappa_use(:,1);

% phase modification at boundaries due to wavevector
lam  = exp(-1i*kvec'*R);
lam(end+1:3) = 0;

% Obtain full model result
T_per_BMS_k = T_per_BMS{1}.s0 + T_per_BMS{1}.s1*lam(1) + T_per_BMS{1}.s2*lam(2) + T_per_BMS{1}.s3*lam(3)...
        + T_per_BMS{1}.s12*lam(1)*lam(2) + T_per_BMS{1}.s23*lam(2)*lam(3) + T_per_BMS{1}.s13*lam(1)*lam(3)...
        + T_per_BMS{1}.s123*lam(1)*lam(2)*lam(3);

Hhat = T_per_BMS_k'*S_BMS{1}.w0*T_per_BMS_k;

i_i = 1:length(dof_sets_BMS{1}.i);
i_b = (length(dof_sets_BMS{1}.i)+1):size(T_per_BMS_k,2);

Hii = Hhat; Hii(i_b,:) = 0; Hii(:,i_b) = 0;
Hbb = Hhat; Hbb(i_i,:) = 0; Hbb(:,i_i) = 0;

[i_all,j_all,Hhatvec] = find(flipud(Hhat([i_i,i_b],[i_i,i_b])));
[i_all_ii,j_all_ii,Hiivec] = find(flipud(Hii([i_i,i_b],[i_i,i_b])));
[i_all_bb,j_all_bb,Hbbvec] = find(flipud(Hbb([i_i,i_b],[i_i,i_b])));

figure(222);clf
Hhatcolors = (abs(imag(Hhatvec))./abs(Hhatvec))*[1,1,1]*0.9;
Hhatcolors = [1,1,1]*0.75;
scatter(j_all,i_all,10,Hhatcolors,'.');hold on

Hiicolors = (abs(imag(Hiivec))./abs(Hiivec))*[1,1,1]*0.9; Hiicolors(:,3) = 1;
Hiicolors = [0,1,1];
scatter(j_all_ii,i_all_ii,10,Hiicolors,'.');hold on

Hbbcolors = (abs(imag(Hbbvec))./abs(Hbbvec))*[1,1,1]*0.9; Hbbcolors(:,1) = 1;
Hbbcolors = [1,0,0];
scatter(j_all_bb,i_all_bb,10,Hbbcolors,'.');hold on

xlim([0,size(Hhat,1)]);
ylim([0,size(Hhat,2)]);
axis equal
axis off
drawnow
rez = 150;
% saveas(gcf,sprintf('Hamiltonian_Sparsity_BMS_kpt_%i_of_%i.png',q,n_kap))

%% Evaluate charge density errors
% ======================================================================= %

% index showing which bands are occupied
f_occ = zeros(n_curves,1);
f_occ(1:4) = 1;

% if ~plot_disp
    
    % preallocate error arrays
    e_rho_BMS_ij    =  zeros(n_kap,n_n_FIs);
    e_E_B_BMS_ij 	=  zeros(n_kap,n_n_FIs);
    e_E_BMS         =  zeros(n_curves,n_kap,n_n_FIs);
    
    % preallocate charge densities
    rho_full = zeros(n_dof_per,n_n_FIs);
    rho_BMS = zeros(n_dof_per,n_n_FIs);
    
    % preallocate band energies
    E_B_BMS = zeros(1,n_n_FIs);
    E_B_ref = zeros(1,n_n_FIs);
    
    % preallocate error arrays
    e_E_B_BMS = zeros(1,n_n_FIs);
    e_rho_BMS = zeros(1,n_n_FIs);
    
    for j = 1:n_n_FIs
        rho_full_j = zeros(n_dof_per,1);
        rho_BMS_j = zeros(n_dof_per,1);
        
        for i = 1:n_kap
            i
            
            % BMS mode shapes
            if isstruct(T_BMS{j})
                PHI_test = BMS_Plus_Mode_Expansion(PHI_BMS{j}(:,:,i),...
                 dof_sets_BMS{j},kappa_use(:,i),R,T_BMS{j},H_BMS{j},S_BMS{j});
            else
                PHI_test = T_BMS{j}*PHI_BMS{j}(:,:,i);
            end
            
            % full model mode shapes
            PHI_full = PHI(:,:,i);
            
            % Reduce Modes to periodic representation
            i_per = ([dof_sets.i,...
                     dof_sets.l,...
                     dof_sets.f,...
                     dof_sets.d,...
                     dof_sets.lf,...
                     dof_sets.ld,...
                     dof_sets.fd,...
                     dof_sets.lfd]);
            PHI_test = PHI_test(i_per,:);
            
            % Orthogonalize Mode Shapes
            %...IS THIS APPROPRIATE? yes I think so
            %...IS THIS THE CORRECT IMPLEMENTATION? (I.E. SHOULD DEGENERATE
            %   SETS BE HANDLED SEPARATELY?) It seems to be fine 
            [PHI_test,~] = qr(PHI_test,0);
            [PHI_full,~] = qr(PHI_full,0);

            % normalize mode shapes
            PHI_test = PHI_test*diag(diag(PHI_test'*PHI_test).^(-1/2));
            PHI_full = PHI_full*diag(diag(PHI_full'*PHI_full).^(-1/2)); 
                        
%             % align phases and normalize again
%             PHI_test = PHI_test*diag(diag(PHI_full'*PHI_test).^(-1));
%             PHI_test = PHI_test*diag(diag(PHI_test'*PHI_test).^(-1/2));

            % charge density
            rho_BMS_i  = (conj(PHI_test).*PHI_test)*f_occ;
            rho_full_i  = (conj(PHI_full).*PHI_full)*f_occ;
            
            rho_BMS_j  = rho_BMS_j + rho_BMS_i*w_use(i);
            rho_full_j  = rho_full_j + rho_full_i*w_use(i);
                    
            % charge density error
%             e_rho_BMS(i,j) = 1-(rho_test'*rho_full)/(norm(rho_test)*norm(rho_full));
            e_rho_BMS_ij(i,j) = norm(rho_BMS_i-rho_full_i)/norm(rho_full_i);
            
            % band energy for ith k-point
            E_B_BMS_i   = sum(E_BMS(:,i,j).*f_occ)/n_atoms;
            E_B_ref_i   = sum(E_ref(:,i).*f_occ)/n_atoms;
            
            % summed band energy up to ith k-point
            E_B_BMS(j)  = E_B_BMS(j) +  E_B_BMS_i*w_use(i);
            E_B_ref(j)  = E_B_ref(j) +  E_B_ref_i*w_use(i);
            
            % energy error
            % e_E_BMS(:,i,j) = (E_BMS(:,i,j)-E_disp(:,i))./E_disp(:,i);
            
            % band energy error stored by k-point
            e_E_B_BMS_ij(i,j) = (E_B_BMS_i-E_B_ref_i);
            

        end
        
        
        
        
        % store jth density function
        rho_BMS(:,j) = rho_BMS_j;
        rho_full = rho_full_j; 
        
        e_rho_BMS(j) = norm(rho_BMS(:,j)-rho_full*ones)/norm(rho_full);
        e_E_B_BMS(j) = (E_B_BMS(j)-E_B_ref(j));
        
        disp(['Band energy full        = ',num2str(E_B_ref),' Ha'])
        disp(['Band energy BMS         = ',num2str(E_B_BMS),' Ha'])
        disp(['Band energy difference  = ',num2str(E_B_BMS-E_B_ref),' Ha'])
        disp(['Band energy % error     = ',num2str(100*(E_B_BMS-E_B_ref)./E_B_ref),'%'])        
        disp(['Charge Density full norm     = ',num2str(sqrt(sum(rho_full.^2))),' ?'])
        disp(['Charge Density BMS norm      = ',num2str(sqrt(sum(rho_BMS.^2))),' ?'])
        disp(['Charge Density difference    = ',num2str(sqrt(sum((rho_BMS-rho_full*ones(1,n_n_FIs)).^2))),' ?'])
        disp(['Charge Density error         = ',num2str(sqrt(sum((rho_BMS-rho_full*ones(1,n_n_FIs)).^2))./sqrt(sum((rho_full*ones(1,n_n_FIs)).^2))),' ?'])
    end 
% end

%% Plot Error
if ~plot_disp
    e_E_B_BMS_max = squeeze(max(abs(e_E_B_BMS_ij)));
    e_rho_BMS_max = squeeze(max(abs(e_rho_BMS_ij)));
    % kappa = 1;


    figure(121);clf
    h1 = semilogy(n_dof_BMS,e_rho_BMS_max,'ro-','markerfacecolor','r');hold on
    h2 = semilogy(n_dof_BMS,e_rho_BMS,'ro--','markerfacecolor','w');hold on
    h3 = semilogy(n_dof_BMS,e_E_B_BMS_max,'bv-','markerfacecolor','b');
    h4 = semilogy(n_dof_BMS,e_E_B_BMS,'bv--','markerfacecolor','w');
    xlabel('# of BMS Degrees of Freedom')
    ylabel('Error')
    legend([h1,h2,h3,h4],'max charge density error over all k-pts, max_k(e_{\rho}(k))',...
                         'total charge density error, e_{\rho}',...
                         'max energy error over all k-pts, max_k(e_{E}(k))',...
                         'total band energy error, e_{Eb}')

    figure(122);clf
    h1 = semilogy(t_BMS/t_full,e_rho_BMS_max,'ro-','markerfacecolor','r');hold on
    h2 = semilogy(t_BMS/t_full,e_rho_BMS,'ro--','markerfacecolor','w');hold on
    h3 = semilogy(t_BMS/t_full,e_E_B_BMS_max,'bv-','markerfacecolor','b');
    h4 = semilogy(t_BMS/t_full,e_E_B_BMS,'bv--','markerfacecolor','w');
    xlabel('comptation time ratio t_{BMS}/t_{full}')
    ylabel('Error')
    legend([h1,h2,h3,h4],'max charge density error over all k-pts, max_k(e_{\rho}(k))',...
                         'total charge density error, e_{\rho}',...
                         'max energy error over all k-pts, max_k(e_{E}(k))',...
                         'total band energy error, e_{Eb}')
end

%% Plot Dispersion
% ======================================================================= %

if plot_disp
    % Create Dispersion Plot
    linewidths = 1.5;
    figure(3);clf
    hold on 

    % plot full dispersion (black lines)
    h_legend = [];
    legendstrings = {};
    if full_disp == 1
        h1 = plot(kappa_plot,E_disp(1:n_curves,:)-E_disp(4,17),'-k','linewidth',linewidths);
        h_legend = [h_legend,h1(1)];
        legendstrings = [legendstrings,['Full Dispersion (',num2str(t_full,4),' sec.)']];
    end
    if BMS_disp
        h2 = plot(kappa_plot,E_BMS(1:n_curves,:,end)-E_BMS(4,17,end),'--r','linewidth',linewidths);
        h_legend = [h_legend,h2(1)];
        legendstrings = [legendstrings,['BMS Dispersion (',num2str(t_BMS(end),4),' sec.)']];
    end
    
    % Determine Axis Limits
    xlimits = [0,kappa_plot(end)];
    ylimits = ylim;
    % ylimits = [-3,17];
    xlim(xlimits)
    ylim(ylimits)

    % Plot Dashed lines at high symmetry points
    hold on
    i_high_sym = linspace(1,n_kap,n_segs+1);

    % tick_labels = {'Gamma';'X';'M';'R';'Gamma'};
    tick_labels = {'L';'\Gamma';'X';'K';'\Gamma'};
%     plot([1;1]*kappa_plot(i_high_sym),ylimits'*ones(1,n_segs+1),'k:','linewidth',linewidths)
    set(gca,'xtick',kappa_plot(i_high_sym))
    set(gca,'xticklabels',tick_labels);
    set(gca,'xgrid','on');
    grid on
    hold off

    % Create axis labels and figure title
    ylabel('Energy (Ha)');xlabel('\kappa')
    title('Dispersion Curves')

    % create legend
    legend(h_legend,legendstrings)

    if save_results
        clockvec = clock;
        savestring = sprintf('figures/3d_electronic_structure_dispersion_%iDOF_%i_%i_%i',...
        n_dof_per,clockvec(2),clockvec(3),clockvec(1));
        saveas(gcf,savestring);
    end
end

%% Prepare Mode Shapes for plotting
% ======================================================================= %
% 
% mode_plot = [1];
write_video = false;  

% % L pt
% i_kappa_plot = [1]; 

% % gamma pt
i_kappa_plot = [17];      

% X pt
% i_kappa_plot = [33];

% % K pt
% i_kappa_plot = [49]; 

kappa = kappa_use;

n_nodes_dir = size(pattern)*(n-1)+1;
n_node_a1 = n_nodes_dir(1);
n_node_a2 = n_nodes_dir(2);
n_node_a3 = n_nodes_dir(3);
n_dof = prod(n_nodes_dir);

xgrid = reshape(coordinates(:,1),n_nodes_dir);
ygrid = reshape(coordinates(:,2),n_nodes_dir);
zgrid = reshape(coordinates(:,3),n_nodes_dir);

% e_rho_val = zeros(n_n_FIs,1);
% rho_full_save = zeros(n_dof,n_n_FIs);
% rho_BMS_save = zeros(n_dof,n_n_FIs);
Phi_full_save = zeros(n_dof,n_curves,n_n_FIs);
Phi_BMS_save = zeros(n_dof,n_curves,n_n_FIs);
for j = 1:n_n_FIs

    T_per = Periodic_Boundary_Conditions(dof_sets); 

    % wave vector at current k-point
    kvec = kappa(:,i_kappa_plot);



    % phase modification at boundaries due to wavevector
    phase  = exp(1i*kvec'*[xgrid(:),ygrid(:),zgrid(:)]').';

    % phase modification at boundaries due to wavevector
    lam  = exp(-1i*kvec'*R);
    lam(end+1:3) = 0;

    % Obtain full model result
    T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
            + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
            + T_per.s123*lam(1)*lam(2)*lam(3);

    PHI_temp = PHI(:,:,i_kappa_plot);
    PHI_temp = T_per_k*PHI_temp;
    [PHI_temp,~] = qr(PHI_temp,0);

    % charge density
    Phi_full_save(:,:,j)  = PHI_temp.*repmat(phase,[1,n_curves]);

    % BMS+ result                
    PHI_temp2 = BMS_Plus_Mode_Expansion(PHI_BMS{j}(:,:,i_kappa_plot),...
             dof_sets_BMS{j},kappa_use(:,i_kappa_plot),R,T_BMS{j},H_BMS{j},S_BMS{j});
    [PHI_temp2,~] = qr(PHI_temp2,0);

    % BMS charge density
%     rho_BMS_save(:,j) = (conj(PHI_temp2).*PHI_temp2)*f_occ;
    Phi_BMS_save(:,:,j) = PHI_temp2.*repmat(phase,[1,n_curves]);

end

% Prepare Charge Densities for plotting

T_per_Gamma = T_per.s0 + T_per.s1 + T_per.s2 + T_per.s3...
        + T_per.s12 + T_per.s23 + T_per.s13 + T_per.s123;

rho_full_save  = T_per_Gamma*rho_full;
rho_BMS_save  = T_per_Gamma*rho_BMS;


%% pattern Mode Shapes (and charge densities) to show multiple unit cells
% ======================================================================= %


[a1grid,a2grid,a3grid] = meshgrid(linspace(0,1,n_nodes_dir(1)),...
                                linspace(0,1,n_nodes_dir(2)),...
                                linspace(0,1,n_nodes_dir(3)));

n_cellx = 1;
n_celly = 1;
n_cellz = 1;
n_nodes_dir2 = n_nodes_dir.*[n_cellx,n_celly,n_cellz];
[a1grid2,a2grid2,a3grid2] = meshgrid(linspace(0,n_cellx,n_nodes_dir2(1)),...
                                    linspace(0,n_celly,n_nodes_dir2(2)),...
                                    linspace(0,n_cellz,n_nodes_dir2(3)));
                                
% PHIrep = Phi_BMS_save(:,mode_plot,n_n_FIs);
mode_type = 'rhoBMS';
% mode_type = 'PHIBMS';
switch mode_type
    case 'PHIfull'
        % select mode and kpoint
        mode_plot = [5];
%         i_kappa_plot = [33];  
        
        % specify colormap
        cmap = viridis;
        colorbox = [0,0,0];
        colorboxlinestyle = '-';
        
        % extract mode
        PHIrep = Phi_full_save(:,mode_plot,n_n_FIs);
        PHIrep = PHIrep.*conj(PHIrep);
        
    case 'rhofull' 
        
        % specify colormap
        cmap = viridis;
        colorbox = [0,0,0];
        colorboxlinestyle = '-';
        
        % extract mode
        PHIrep = rho_full_save;
        
    case 'PHIBMS'
        % select mode and kpoint
        mode_plot = [5];
%         i_kappa_plot = [5];   
        
        % specify colormap
        cmap = viridis;
        colorbox = [0,0,0];
        colorboxlinestyle = '-';
        
        % extract mode
        PHIrep = Phi_BMS_save(:,mode_plot,n_n_FIs);
        PHIrep = PHIrep.*conj(PHIrep);
        
    case 'rhoBMS' 
        
        % specify colormap
        cmap = viridis;
        colorbox = [0,174,239]/256;
        colorboxlinestyle = '--';
        % extract mode
        PHIrep = rho_BMS_save(:,1);
        
    case 'TCB'
        
        % extract mode
        PHIrep = T_BMS{1}.w0(:,dof_sets_BMS{1}.d(12));
        
        % specify colormap
        cmap = colormap(autumn);
        cmap = colormap(winter);
        colorbox = cmap(1,:)*0.6;
end
        
        
        
PHIrep = reshape(PHIrep,n_nodes_dir);

Phi_exp = zeros(n_nodes_dir2);
for i = 1:n_cellx
    for j = 1:n_celly
        for k = 1:n_cellz
                % silly j,i,k ordering because this must coincide with
                % meshgrid
                Phi_exp((j-1)*n_nodes_dir(2) + (1:n_nodes_dir(2)),...
                        (i-1)*n_nodes_dir(1) + (1:n_nodes_dir(1)),...
                        (k-1)*n_nodes_dir(3) + (1:n_nodes_dir(3))) = ...
                       PHIrep;
        end
    end
end


% create a set of patches to outline the domain
z0 = [0;0;0];
Na = [nan;nan;nan];
boxpatch = [z0,     a1,     a1+a2,      a2,         Na,...
            a3,     a1+a3,  a1+a2+a3,   a2+a3,      Na,...
            z0,     a1,     a1+a3,      a3,         Na,...
            a2,     a1+a2,  a1+a3+a2,   a3+a2,      Na,...
            z0,     a3,     a3+a2,      a2,         Na,...
            a1,     a1+a3,  a1+a3+a2,   a1+a2,      Na];    

boxoutline = [z0,       a1,         a1+a2,      a2,         z0,         Na,...
              a3,       a1+a3,      a1+a2+a3,   a2+a3,      a3,         Na,...
              z0,       a3,         Na,...
              a1,       a1+a3,      Na,...
              a1+a2,    a1+a2+a3,   Na,...
              a2,       a2+a3,      Na];
                  
          
boxpatchx = reshape(boxpatch(1,:),[5,size(boxpatch,2)/5]);    boxpatchx(end,:) = [];
boxpatchy = reshape(boxpatch(2,:),[5,size(boxpatch,2)/5]);    boxpatchy(end,:) = [];
boxpatchz = reshape(boxpatch(3,:),[5,size(boxpatch,2)/5]);    boxpatchz(end,:) = [];

% Create a separate box outline for EVERY primitive cell
% boxoutline2 = [];
% for i = 1:n_cellx
%     for j = 1:n_celly
%         for k = 1:n_cellz
%                boxoutlinetemp = boxpatch;
%                boxoutlinetemp(1,:) = boxoutlinetemp(1,:)+a1(1)*(i-1)+a2(1)*(j-1)+a3(1)*(k-1);
%                boxoutlinetemp(2,:) = boxoutlinetemp(2,:)+a1(2)*(i-1)+a2(2)*(j-1)+a3(2)*(k-1);
%                boxoutlinetemp(3,:) = boxoutlinetemp(3,:)+a1(3)*(i-1)+a2(3)*(j-1)+a3(3)*(k-1);
%                boxoutline2 = [boxoutline2,boxoutlinetemp];
%         end
%     end
% end
boxoutline2 = [];
for i = 1:n_cellx
    for j = 1:n_celly
        for k = 1:n_cellz
               boxoutlinetemp = boxoutline;
               boxoutlinetemp(1,:) = boxoutlinetemp(1,:)+a1(1)*(i-1)+a2(1)*(j-1)+a3(1)*(k-1);
               boxoutlinetemp(2,:) = boxoutlinetemp(2,:)+a1(2)*(i-1)+a2(2)*(j-1)+a3(2)*(k-1);
               boxoutlinetemp(3,:) = boxoutlinetemp(3,:)+a1(3)*(i-1)+a2(3)*(j-1)+a3(3)*(k-1);
               boxoutline2 = [boxoutline2,boxoutlinetemp];
        end
    end
end


%%% Interpolate results to a finer grid
% ======================================================================= %

n_nodes_dir3 = 5*n_nodes_dir2;
[a1grid3,a2grid3,a3grid3] = meshgrid(linspace(0,n_cellx,n_nodes_dir3(1)),...
                                     linspace(0,n_celly,n_nodes_dir3(2)),...
                                     linspace(0,n_cellz,n_nodes_dir3(3)));

       
% interpolate BMS charge density grid;
Phi_exp3 = interp3(a1grid2,a2grid2,a3grid2,Phi_exp,a1grid3,a2grid3,a3grid3,'spline');


xgrid3 = a1grid3*a1(1) + a2grid3*a2(1) + a3grid3*a3(1);
ygrid3 = a1grid3*a1(2) + a2grid3*a2(2) + a3grid3*a3(2);
zgrid3 = a1grid3*a1(3) + a2grid3*a2(3) + a3grid3*a3(3);
%%% Plot Mode shapes
% ======================================================================= %

nt = 1;
tvec = linspace(0,2*pi,nt+1);
    
% wave vector at current k-point
kvec = kappa(:,i_kappa_plot);

% phase modification at boundaries due to wavevector
phase  = exp(-1i*kvec'*[xgrid3(:),ygrid3(:),zgrid3(:)]').';

% mode_type = 'TCB';
box_faces = false;
box_lines = true;
switch mode_type
    case 'PHIfull'
        Vec_phase = reshape((Phi_exp3(:).*phase),n_nodes_dir3);
        maxC = max(abs(Vec_phase(:)));
        minC = -maxC;
    case 'rhofull'
        Vec_phase = reshape(Phi_exp3(:),n_nodes_dir3);
        maxC = max(Vec_phase(:));
        minC = min(Vec_phase(:));
    case 'PHIBMS'
        Vec_phase = reshape((Phi_exp3(:).*phase),n_nodes_dir3);
        maxC = max(abs(Vec_phase(:)));
        minC = -maxC;
    case 'rhoBMS'
        Vec_phase = reshape(Phi_exp3(:),n_nodes_dir3);
        maxC = max(Vec_phase(:));
        minC = min(Vec_phase(:));
    case 'TCB'
        logarithmic = false;
        if logarithmic
            Vec_phase = reshape(log10(abs(Phi_exp3(:))),n_nodes_dir3);
            Vec_phase(Phi_exp3(:)==0) = -1000;
            maxC = 0;
            minC = -2;
        else
            Vec_phase = reshape(Phi_exp3(:),n_nodes_dir3);
            maxC = max(real(Vec_phase(:)));
            minC = min(real(Vec_phase(:)));
        end
        box_faces = true;
end

n_isosteps = 5;
isosteps = linspace(0,1,n_isosteps+2);isosteps = isosteps(2:end-1); 
isovals = minC + isosteps*(maxC-minC);
n_isosteps = length(isovals);
clear h legend_entries

 for p = 1:nt 
    
    figure(601);clf;  
    
    % harmonic term
    harmonic = exp(1i*tvec(p));
    Vec_plot = real(Vec_phase*harmonic);
     
    for i = 1:n_isosteps
        legend_entries{i} = [num2str(isovals(i),3)];
        isoSurfPatchStruct = isosurface(conversionFactor*xgrid3,...
                                        conversionFactor*ygrid3,...
                                        conversionFactor*zgrid3,...
                                        Vec_plot,isovals(i));
        h(i) = patch(isoSurfPatchStruct);
        i_color = round(1+(i-1)*(size(cmap,1)-1)/(length(isosteps)-1));
        set(h(i),'FaceColor',cmap(i_color,:),...
                 'EdgeColor','none',...
                 'facealpha',isosteps(i)^(0));
    end
    
    % plot outline
    if box_lines
        hold on
        hboxline = plot3(conversionFactor*boxoutline2(1,:),...
                         conversionFactor*boxoutline2(2,:),...
                         conversionFactor*boxoutline2(3,:));
        set(hboxline,'color',colorbox)
        set(hboxline,'linestyle',colorboxlinestyle)
        set(hboxline,'linewidth',3)
    end
        
    % plot box patches
    if box_faces
        hold on
        hbox = patch(conversionFactor*boxpatchx,...
              conversionFactor*boxpatchy,...
              conversionFactor*boxpatchz,colorbox);
        set(hbox,'edgecolor','none')
        set(hbox,'facealpha',0.25)
    end
    
%     plot_box = true;
%     if plot_box
%         hboxline = plot3(conversionFactor*boxoutlinex([1:end,1],:),...
%                          conversionFactor*boxoutliney([1:end,1],:),...
%                          conversionFactor*boxoutlinez([1:end,1],:));
%         set(hboxline,'color',colorbox)
%         set(hboxline,'linestyle',colorboxlinestyle)
%         set(hboxline,'linewidth',3)
%     end
    
    
    % Add legend
    if false
        legend(h(end:-1:1),legend_entries(end:-1:1),'location','northeast')
    end
    
    % view options
    view(3);     
    view([-15,15])
    axis tight
    axis equal
    axis off
    
    % axis labels
    xlabel(['x (',plot_unit_string,')']);    
    ylabel(['y (',plot_unit_string,')']);
    zlabel(['z (',plot_unit_string,')']);
    set(gcf,'color','w')
    
    rendering = 'vector';
    rendering = 'pixel';
    switch rendering
        case 'pixel'
            camlight        
            camlight left
            lighting gouraud
        case 'vector'
            set(gcf,'renderer','painters')
            
    end
    
    drawnow
    disp(['timestep ',num2str(p),' of ',num2str(nt)])

    if write_video
        frames(p) = getframe(gcf);
    end  

end

if save_results & false
    clockvec = clock;
    if q == 1
        savestring = sprintf('figures/3d_electronic_full_density_%iDOF_%i_%i_%i',...
        n_dof,clockvec(2),clockvec(3),clockvec(1));
    elseif q == 2
        savestring = sprintf('figures/3d_electronic_BMS_density_%iDOF_%iFIM_%i_%i_%i',...
        n_dof,n_FIs(i_FI_plot),clockvec(2),clockvec(3),clockvec(1));
    end
    saveas(gcf,savestring);
end

% if j < size(e_rho3save,2)
%     delete(h_dot)
% end

%%
if write_video
%     n_frames = round(logspace(log10(5),0,n_n_FIs));
    n_frames = ones(1,length(frames));
%     video_object = VideoWriter('Silicon_Charge_Error.avi','Uncompressed AVI');
%     video_object = VideoWriter('Silicon_Charge_Bloch_Wave.avi');
%     video_object = VideoWriter('Silicon_Plane_Wave.avi');
    video_object.FrameRate = 30;
    open(video_object);
    for j = 1:length(frames)
        for k = 1:n_frames(j)
            writeVideo(video_object,frames(j))
        end
    end
    close(video_object);
end


%% Save Results
% ======================================================================= %
if save_results
    clockvec = clock;
    savestring = sprintf('Silicon_BMS_Analysis_%iDOF_%i_%i_%i',...
        n_dof_per,clockvec(2),clockvec(3),clockvec(1));
    save(savestring);
end