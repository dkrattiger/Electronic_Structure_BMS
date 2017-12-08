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

addpath(genpath('Global_Files'))

%% User Selected Parameters
% ======================================================================= %

% Solution Options
save_results    = true;
save_figures    = true;
load_results    = true;
n_curves        = 8;       % number of dispersion curves to plot
% n_kap_1      = 3;
use_ibz         = true;
% end

% Mode Animations
mode_plot = [1];     % select modes to plot. [] =  don't animate any modes
k_sel = 1;          % k-point to plot modes at

%% Check for & Create Save directories
% ======================================================================= %
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
    
% if figures does not exist, make it
if save_figures
    if ~exist('figures','dir')
        mkdir('figures');
    end
end
%% Define Model dimensions
% ======================================================================= %
n_atoms_base = 2;
n_cellx = 1;        % number of repetitions of unit cell in x-direction
n_celly = 1;        % number of repetitions of unit cell in y-direction
n_cellz = 1;        % number of repetitions of unit cell in z-direction
n_atoms = n_atoms_base*n_cellx*n_celly*n_cellz;

% n = 5;              % number of element nodes per edge
% n_bz = 1;           % number of brillouin Zones to cover

length_select = 'atomic_units';
% length_select = 'meters';
switch length_select
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
% a0 = 1;
% a0 = 5.43;


% lattice vectors
a1 = (a0/2)*[0;1;1]*n_cellx;
a2 = (a0/2)*[1;0;1]*n_celly;
a3 = (a0/2)*[1;1;0]*n_cellz;
R = [a1,a2,a3];

% reciprocal lattice vectors
b1 = 2*pi*cross(a2,a3)/dot(a1,cross(a2,a3));
b2 = 2*pi*cross(a3,a1)/dot(a2,cross(a3,a1));
b3 = 2*pi*cross(a1,a2)/dot(a3,cross(a1,a2));

% units are in Ry to begin with (from pseudopotential)
eV_Per_Ha = 27.2114; 
eV_Per_Ry = eV_Per_Ha/2;
Ha_Per_Ry = 1/2;

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


%% Loop through different model sizes
% ======================================================================= %

runselect = 'big';

switch runselect
    
    
    case 'test'
        stretch = 3;
        n_dof_ideal_list = ( 10.^(linspace(2^(stretch),...
                            log10(18000)^(stretch),10).^(1/stretch)) );
        elem_order_list = [5];
        MPrule_list = [12];
    
    case 'big'
        stretch = 3;
        n_dof_ideal_list = ( 10.^(linspace(2^(stretch),...
                            log10((4*6)^3)^(stretch),10).^(1/stretch)) );
%         elem_order_list = [1,2,3,4,5];
        elem_order_list = [3,4,5,6,7];
        MPrule_list = [7:10];

    case 'small'
        stretch = 3;
        n_dof_ideal_list = ( 10.^(linspace(2^(stretch),...
                            log10(5000)^(stretch),5).^(1/stretch)) );
        elem_order_list = [1,2];
        MPrule_list = [12];
        
    case 'MPsmalltest'
        n_dof_ideal_list = 12^3;%=1728
        elem_order_list = [1,2,3];
        MPrule_list = [1:15];

    case 'MP1'
        n_dof_ideal_list = 12^3;%=1728
        elem_order_list = [1,2,3,4,6];
        MPrule_list = [1:15];
        
    case 'MP2'
        n_dof_ideal_list = 15^3;%=8000
        elem_order_list = [1,3,5];
        MPrule_list = [1:15];
        
    case 'MP3'
        n_dof_ideal_list = 24^3;%=8000
        elem_order_list = [1,2,3,4];
        MPrule_list = [1:15];
        
    case 'MPCollect'
        n_dof_ideal_list = [(12^3)*ones(1,5),(15^3)*ones(1,3),(20^3)*ones(1,4)];
        elem_order_list = [[1,2,3,4,6],[1,3,5],[1,2,4,5]];
        MPrule_list = [1:15];
        
    case 'kill'
        n_dof_ideal_list = (6*5)^3; %=27000
        elem_order_list = [6];
        MPrule_list = [20];
end

for q1 = 1:length(elem_order_list)
    
    % nodes per element edge
    n = elem_order_list(q1)+1;
    eleorder = elem_order_list(q1);

    % number of dof to aim for in each model
    if length(n_dof_ideal_list) == 1
        n_dof_ideal_list = n_dof_ideal_list*ones(1,length(elem_order_list));        
    end
    n_dof_ideal = n_dof_ideal_list;
    
%     n_dof_ideal = n_dof_ideal_list(q1);
    
    % number of elements per edge
    n_ele_edge_list{q1} = unique(round((n_dof_ideal.^(1/3))./(eleorder)));
    n_ele_edge_list{q1}(n_ele_edge_list{q1}==0) = [];
    
    % preallocate results storage arrays
    E_b_store{q1} = zeros(length(n_ele_edge_list{q1}),length(MPrule_list));
    n_dof_store{q1} = zeros(length(n_ele_edge_list{q1}),length(MPrule_list));
    t_full_store{q1} = zeros(length(n_ele_edge_list{q1}),length(MPrule_list));
    E_b_max_store{q1} = zeros(length(n_ele_edge_list{q1}),length(MPrule_list));

            
            
    for q2 = 1:length(n_ele_edge_list{q1})

        % number of elements per model edge
        n_ele_edge = n_ele_edge_list{q1}(q2);

        % frac = pat_size*2/3;
        pattern = 4*ones(n_ele_edge,n_ele_edge,n_ele_edge);
        pattern = repmat(pattern,[n_cellx,n_celly,n_cellz]);

        
        %% Form Model (Mesh, and FE matrices)
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
            
            % define mesh
%             eleorder = n-1;
            [coordinates,elenodes] = MeshParallelPiped(pattern,eleorder,R);
            
            % store mesh info in a structure
            Mesh_info.pattern_vec   = pattern(:);
            Mesh_info.Vs            = Vs;
            Mesh_info.n             = n;
            Mesh_info.elenodes      = elenodes;
            Mesh_info.coordinates   = coordinates;

            % compute system matrices
            [H0f,Sf] = master_matrices_3D_Electronic(Mesh_info);

            H0f = sparse(H0f);
            Sf = sparse(Sf);

            if save_results
                save(model_savestring,'H0f','Sf','coordinates')
            end
        end

        %% Create Degree-of-freedom sets
        % ======================================================================= %
        
        node_sets = find_node_sets(coordinates,R);
        dof_sets = node_sets;
        T_per = Periodic_Boundary_Conditions(dof_sets);

        %% Discretize Brillouin zone
        % ======================================================================= %
        for q3 = 1:length(MPrule_list)
            n_kap_edge = MPrule_list(q3);
            [kappa_ibz,w_ibz,kappa_full] = MonkhorstPackGrid(b1,b2,b3,n_kap_edge);

            if use_ibz
                kappa_use = kappa_ibz;
                w_use = w_ibz;
            else
                kappa_use = kappa_full;
                w_use = ones(1,size(kappa_use,2))/size(kappa_use,2);
            end
            n_kap = size(kappa_use,2);

            %% Dispersion Solution
            % ======================================================================= %

            % E_disp = zeros(n_curves,n_kap);
            % PHI = zeros(n_dof,n_curves,n_kap);
            % pause
            disp(['element order = ',num2str(n-1),...
                  ', n_ele_edge = ',num2str(n_ele_edge),...
                  ', n_kap_edge = ',num2str(n_kap_edge)])
              
            % BZ or IBZ in file names
            if use_ibz     
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
                load([solution_savestring,'.mat'],'E_disp','t_kloop')
            else
                t_loop = 0;

                [omega,PHI,t_kloop] = dispersion_solver_w_k(kappa_use,H0f,Sf,dof_sets,R,n_curves);
                E_disp = omega.^2;

                if save_results
                    save([solution_savestring,'.mat'],'E_disp','PHI','t_kloop')
                end

            end
            t_full = sum(t_kloop);
            t_kap = sum(t_kloop)/length(t_kloop);


            % total band energy
            f_occ = [ones(4,1);zeros(4,1)];          
            E_b = bandenergy(E_disp,f_occ,w_use)/n_atoms;

            % compute band energy for each k-point
            for i = 1:n_kap
                E_b_k_store{q1}{q3}(q2,i) = f_occ'*E_disp(:,i)/n_atoms;
            end
            E_b_store{q1}(q2,q3) = E_b;
            n_dof_store{q1}(q2,q3) = n_dof_per;
            t_full_store{q1}(q2,q3) = t_full;
            

        end
    end
end

%% choose a model to use for error reference
% ======================================================================= %
use_kill_model_soln_as_reference = true;
if use_kill_model_soln_as_reference
    
    ref_type = 'pw';
    switch ref_type
        case 'pw'
            ref_model = ['save_data/solutions/','Si_1x1x1PrimCell_PW_6859DOF_256IBZkpts_8Bands'];
            ref_soln = load(ref_model,'E');
            n_atoms_ref = n_atoms;
            E_ref = ref_soln.E;

            options = [];
            [~,w_ref] = MonkhorstPackGrid(b1,b2,b3,20);
%     [~,~,w_ref] = MonkhorstPackGrid(b1,b2,b3,8); length(w_ref)
        case 'fe'
            ref_model = ['save_data/solutions/','Si_1x1x1PrimCell_27000DOF_6thOrderEles_256IBZkpts_8Bands'];
            [~,w_ref] = MonkhorstPackGrid(b1,b2,b3,20);
    end
else
    
    % use last result (should be most converged one)
    E_ref = E_disp;
    w_ref = w_use;
end

% total band energy
f_occ = [ones(4,1);zeros(n_curves-4,1)];          
E_b_ref = bandenergy(E_ref,f_occ,w_ref)/n_atoms_ref;
n_kap = size(E_ref,2);
for i = 1:n_kap
    E_b_ref_k(i) = f_occ'*E_ref(:,i)/n_atoms_ref;
end


%% Plotting Parameters
% ======================================================================= %
figure(1)
% colorbank  = get(gca,'ColorOrder');
colorbank =[...
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840]; % this is the newest default color set

linestyles = {'-','--',':','-.'};
markers  = {'.','s','v','o'};
orders = {'1st','2nd','3rd','4th','5th','6th','7th'};
close

%% Plot E_b convergence vs model size
% ======================================================================= %
figure(1);clf
for q1 = 1:length(elem_order_list)
    for q3 = 1:length(MPrule_list)
    
        h1 = plot(n_dof_store{q1}(:,q3),E_b_store{q1}(:,q3));hold on
        set(h1,'color',colorbank(rem(q1-1,size(colorbank,1))+1,:))
        set(h1,'linestyle',linestyles{rem(q3-1,length(linestyles))+1})
        set(h1,'linewidth',2)
        set(h1,'marker','.')
        set(h1,'markersize',12)
        
        legendstrings{q1} = ['Ele. order = ',num2str(elem_order_list(q1)),', ',...
            num2str(MPrule_list(q3)),'\times',...
            num2str(MPrule_list(q3)),'\times',...
            num2str(MPrule_list(q3)),'k-pts'];
    end
end

xlabel('n_{DOF}')
ylabel('Band energy, E_b (Ry)')
legend(legendstrings)

%% Compute and plot E_b error vs model size
% ======================================================================= %
e_B_type = 'avg'; % 'max' or 'avg'

if ~strcmp(runselect([1,2]),'MP')
%     n_atoms = 2;

    % compute error
    for q1 = 1:length(elem_order_list)
        for q3 = 1:length(MPrule_list)
            switch e_B_type
                case 'max'
                    n_kap = size(E_ref,2);size(E_ref)
                    for i = 1:n_kap
                        e_E_b_k{q1}{q3}(:,i) = E_b_k_store{q1}{q3}(:,i) - E_b_ref_k(:,i);
                    end
                    e_E_b{q1}(:,q3) = max(abs(e_E_b_k{q1}{q3}(:,:)),[],2);
                case 'avg'
                    e_E_b{q1}(:,q3) = (E_b_store{q1}(:,q3)-E_b_ref);
            end
        end
    end

    % plot error vs model size for each MP grid and element type
    for i = 1:4    
        figure(200+i);clf
        % linewidth = 2;
        for q1 = 1:length(elem_order_list)
            for q3 = 1:length(MPrule_list)
                switch i
                    case 1
                        xlab = 'n_{DOF}';
                        h1 = loglog(n_dof_store{q1}(:,q3),e_E_b{q1}(:,q3)*(Ha_Per_Ry),'.-');hold on        
                        save_string_end = 'NDOF';
                    case 2
                        h1 = semilogy(n_dof_store{q1}(:,q3).^(1/3),e_E_b{q1}(:,q3)*(Ha_Per_Ry),'.-');hold on
                        xlab = 'nodes per edge';                
                        save_string_end = 'NodesPerEdge';
                    case 3
                        h1 = loglog(n_ele_edge_list{q1},e_E_b{q1}(:,q3)*(Ha_Per_Ry),'.-');hold on
                        xlab = 'elements per edge';
                        set(gca,'xtick',[1,2,3,4,5,8,12,18,30])                
                        save_string_end = 'ElesPerEdge';
                    case 4
                        h1 = loglog(t_full_store{q1}(:,q3),e_E_b{q1}(:,q3)*(Ha_Per_Ry),'.-');hold on
                        xlab = 'computation time';
                        save_string_end = 'Time';
                end

                set(h1,'color',colorbank(rem(q1-1,size(colorbank,1))+1,:))
                set(h1,'linestyle',linestyles{rem(q3-1,length(linestyles))+1})
                set(h1,'linewidth',2)
                set(h1,'marker','.')
                set(h1,'markersize',12)

            %     legendstrings{q1} = [num2str(n_node_edge_list(q1)^3),' nodes/element'];
                legendstrings{(q1-1)*length(MPrule_list)+q3} = [orders{elem_order_list(q1)},'-order brick elements, ',...
                                     num2str(MPrule_list(q3)),'\times',...
                                     num2str(MPrule_list(q3)),'\times',...
                                     num2str(MPrule_list(q3)),'k-pts'];
            end
        end

        grid on
        grid minor,grid minor
        xlabel(xlab)
        switch e_B_type
            case 'max'
                ylabel('Max band energy error, e_{E_b} (Ha/atom)')
            case 'avg'
                ylabel('Avg band energy error, e_{E_b} (Ha/atom)')
        end
        legend(legendstrings)

        if save_figures
            name=['Silicon_Electronic_BandEnergyError_vs_',save_string_end]; %what you want the file to be called
            path= ['figures/']; %the folder where you want to put the file

            % save .FIG file
            saveas(gcf,[path,name])
        end

    end
end

if true
    save('save_data/Silicon_Electronic_FE_Error_vs_Time','t_full_store','e_E_b')
end

%% plot error vs MP grid size for each element type
% ======================================================================= %

if strcmp(runselect([1,2]),'MP')

    % compute error
    for q1 = 1:length(elem_order_list)
        for q2 = 1:length(n_ele_edge_list{q1})

            error_type = 'total';
            % error_type = 'MP';
            switch error_type
                case 'total'
%             e_E_b2{q1}(q2,:) = (E_b_store{q1}(q2,:)-E_b_store{end}(end,end))/(Ha2eV*n_atom);            
                    e_E_b2{q1}(q2,:) = (E_b_store{q1}(q2,:)-E_b_ref)*Ha_Per_Ry;
                case 'MP'
                    e_E_b2{q1}(q2,:) = (E_b_store{q1}(q2,:)-E_b_store{q1}(q2,end))*Ha_Per_Ry;
            end
        end
    end
    
    
    figure(3);clf
    for q1 = 1:length(elem_order_list)
        for q2 = 1:length(n_ele_edge_list{q1})
            
            
            n = elem_order_list(q1);

            x_axis_val = 'Time';
            x_axis_val = 'MPgridsize';
            switch x_axis_val
                case 'Time'
                    x = t_full_store{q1}(q2,:);
                    xlab = 'Computation Time (s)';
                case 'MPgridsize'
                    x = MPrule_list;
                    xlab = '(n \times n \times n) MP-rule';
            end            
            
            h1 = semilogy(x,abs(e_E_b2{q1}(q2,:)),'.-');hold on        
%             h1 = loglog(n_kap_edge_list,abs(e_E_b2{q1}(q2,:)),'.-');hold on        
    %         save_string_end = 'NDOF';            

            set(h1,'color',colorbank(rem(n-elem_order_list(1),size(colorbank,1))+1,:))
%             set(h1,'linestyle',linestyles{rem(q2-1,length(linestyles))+1})
            set(h1,'linestyle',linestyles{find(unique(n_dof_ideal_list) == n_dof_ideal_list(q1))})
            set(h1,'linewidth',2)
            set(h1,'marker','.')
            set(h1,'markersize',12)

        %     legendstrings{q1} = [num2str(n_node_edge_list(q1)^3),' nodes/element'];
            legendstrings{(q1-1)*length(n_ele_edge_list{q1})+q2} = [orders{elem_order_list(q1)},'-order brick elements, ',...
                                 num2str(n_ele_edge_list{q1}(q2)),'\times',...
                                 num2str(n_ele_edge_list{q1}(q2)),'\times',...
                                 num2str(n_ele_edge_list{q1}(q2)),'eles/prim. cell'];
        end
    end

    grid on
    grid minor,grid minor
    xlabel(xlab)
    ylabel('Absolute Band energy error, e_{E_b} (Ha/atom)')
    legend(legendstrings)
    
%     xlim([0,15])
%     ylim([

    if save_figures
        name=['Silicon_Electronic_',error_type,'_BandEnergyError_vs_',x_axis_val]; %what you want the file to be called
        path='figures/'; %the folder where you want to put the file

        % save .FIG file
        saveas(gcf,[path,name])
    end
end
