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


% Solution Options
save_results = true;
save_figures = true;
load_results = true;
n_curves     = 8;       % number of dispersion curves to plot
full_disp    = true;  	% 1 = compute full dispersion, 0 = don't
WFE          = true; 	% 1 = Bloch Boundary Conditions, 0 = Bloch Operator
RBME_disp    = false;  	% Compute RBME dispersion
BMS_disp     = true;    % compute BMS dispersion
plot_disp    = false;

% unit conversions
eV_Per_Ry = 13.6;
eV_Per_Ha = 2*eV_Per_Ry;
Ha_Per_Ry = 1/2;


runtype = 'big_run';
%runtype = 'compare_LIR';
switch runtype
    case 'small_run'

        % these vary along each curve
        n_FIs           = [25,50,100];
        
        % these vary between curves
        n_LIs           = [40,40,40,40];
        pat_size_list   = [3,3,3,3];
        n_list          = [4,6,6,6];
        n_kap_edge_list = [3,7,8,9];
        Boundary_Methods = {'exact','exact','exact','exact'};
       
        
    case 'big_run'
        
        % these vary along each curve
        n_FIs           = [5,10,25,50,100,200];        
        
        % these vary between curves
        n_LIs           = [40,40,40,40];
        pat_size_list   = [3,3,3,3];
        n_list          = [6,6,6,6];
        n_kap_edge_list = [7,8,9,10];
        Boundary_Methods = {'exact','exact','exact','exact'};
        
    case 'compare_LIR'
        
        % these vary along each curve
        n_FIs           = [5,10,25,50,100];
        
        % these vary between curves
        n_LIs           = [40,70];
        pat_size_list   = [3,3];
        n_list          = [6,6];
        n_kap_edge_list = [9,9];
        Boundary_Methods = {'exact','hybrid'};
end
        

n_n_FIs = length(n_FIs);
% n_n_LIs = length(n_LIs);
list_size = length(pat_size_list);

e_E_B_BMS = zeros(n_n_FIs,list_size);
e_E_B_full= zeros(1,list_size);
t_BMS = zeros(n_n_FIs,list_size);
n_dof_BMS = zeros(n_n_FIs,list_size);

for q1 = 1:list_size
%     for q2 = 1:n_n_LIs

        %% User Selected Parameters
        % ======================================================================= %

        % pattern is a 3D array. The number of terms in the 1st, 2nd, and 3rd
        % indices correspond to the number of elements in the x,y, and z directions


        % pat_size = 1;
        % pat_size = 2;
        % pat_size = 4;
        pat_size = pat_size_list(q1);
        % pat_size = 8;
        % pat_size = 10;
        % pat_size = 20;


        % frac = pat_size*2/3;
        pattern = 4*ones(pat_size,pat_size,pat_size);

        % size of k-point grid for monkhorst pack
        if ~plot_disp
            n_kap_1 = n_kap_edge_list(q1);

            use_ibz = true;
        end

        % Mode Animations
        mode_plot = [1];     % select modes to plot. [] =  don't animate any modes
        k_sel = 1;          % k-point to plot modes at

        %% Define Model dimensions
        % ======================================================================= %
        n_cellx = 1;        % number of repetitions of unit cell in x-direction
        n_celly = 1;        % number of repetitions of unit cell in y-direction
        n_cellz = 1;        % number of repetitions of unit cell in z-direction
        n = n_list(q1);              % number of element nodes per edge
        n_bz = 1;           % number of brillouin Zones to cover

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
            KappaPts = [LPt,GammaPt,XPt,nan(3,1),XprimePt,KPt,GammaPt];
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

            % Try All high sym points and their permutations
            try_high_sym_pts_and_perms = true;
            if try_high_sym_pts_and_perms 
                kappa_use = [(XprimePt),(XPt)];
                n_kap = size(kappa_use,2);
                kappa_plot = 1:n_kap;
            end

            n_kap = size(kappa_use,2);
            w_use = ones(1,size(kappa_use,2))/size(kappa_use,2);
        else

            %[kappa2,i_ibz,w_ibz] = MonkhorstPackGrid(b1,b2,b3,n_kap_1);
            [~,w_ibz,kappa2,i_ibz,~] = MonkhorstPackGrid(b1,b2,b3,n_kap_1);


            if use_ibz
                kappa_use = kappa2(:,i_ibz);
                w_use = w_ibz;
            else
                kappa_use  = kappa2;
                w_use = ones(1,size(kappa_use,2))/size(kappa_use,2);
            end

            n_kap = size(kappa_use,2);

        end


%         %% Create Mesh
%         % ======================================================================= %
%         R = [a1,a2,a3];
%         pattern = repmat(pattern,[n_cellx,n_celly,n_cellz]);
%         eleorder = n-1;
%         [coordinates,elenodes] = MeshParallelPiped(pattern,eleorder);
% 
%         %% Create Degree-of-freedom sets
%         % ======================================================================= %
% 
%         coordinates = [xlocs,ylocs,zlocs];
%         R = [a1,a2,a3];
%         node_sets = find_node_sets(coordinates,R);
%         dof_sets = node_sets;
%         T_per = Periodic_Boundary_Conditions(dof_sets);

        % pause
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
            
            % define mesh
            pattern = repmat(pattern,[n_cellx,n_celly,n_cellz]);
            eleorder = n-1;
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
                save([model_savestring],'H0f','Sf','coordinates')
            end
        end
        

        %% Create Degree-of-freedom sets
        % ======================================================================= %

        node_sets = find_node_sets(coordinates,R);
        dof_sets = node_sets;
        T_per = Periodic_Boundary_Conditions(dof_sets);
        [n_dof_per_store(q1),n_dof_store(q1)] = size(T_per.s0);
        
        %% Dispersion Solution
        % ======================================================================= %

        if full_disp == 1

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
                load([solution_savestring,'.mat'],'E_disp','t_kloop')
            else
                t_loop = 0;

                [omega,PHI,t_kloop] = dispersion_solver_w_k(kappa_use,H0f,Sf,dof_sets,R,n_curves);
                E_disp = omega.^2;

                if save_results
                    save(solution_savestring,'E_disp','PHI','t_kloop')
                end

            end
            t_full(q1) = sum(t_kloop);
        end

%         E_disp = E_disp*eV_Per_Ry; % Ry to eV 

        %% choose a model to use for error reference
        % ======================================================================= %
        use_kill_model_soln_as_reference = true;
        if use_kill_model_soln_as_reference
            
            ref_type = 'pw';
            switch ref_type
                
                case 'pw'
                    
                    ref_model = 'save_data/solutions/Si_1x1x1PrimCell_PW_3375DOF_120IBZkpts_8Bands';
                    ref_soln = load(ref_model);
                    E_ref = ref_soln.E;
                    [~,w_ref] = MonkhorstPackGrid(b1,b2,b3,15);
                    
                case 'fe'

                    ref_model = ['save_data/solutions/','Si_1x1x1PrimCell_13824DOF_6thOrderEles_256IBZkpts_8Bands'];
                    ref_soln = load(ref_model);
                    E_ref = ref_soln.E_disp;
                    [~,w_ref] = MonkhorstPackGrid(b1,b2,b3,20);
                    
            end

        else
            E_ref = E_disp;
            w_ref = w_use;
        end



        %% Evaluate charge density errors
        % ======================================================================= %

        % index showing which bands are occupied
        f_occ = zeros(n_curves,1);
        f_occ(1:4) = 1;
        n_atom = 2;

        % reference model band energy
        E_B_ref = bandenergy(E_ref,f_occ,w_ref)/(n_atom); 

        % full model band energy and band energy error
        E_B_full = bandenergy(E_disp,f_occ,w_use)/n_atom;   
        e_E_B_full(q1) = (E_B_full-E_B_ref);
        
        
        %% Use BMS reduction to compute band structure
        % ======================================================================= %

        w_cut = -1;
        w_cut = max(max(sqrt(E_disp/eV_Per_Ry)))*1.5;
        % n_LI = 70;
    %     n_n_FIs = length(n_FIs);

        E_BMS = zeros(n_curves,n_kap,n_n_FIs);
        % PHI_BMS = zeros(n_dof,n_curves,n_kap,n_n_FIs);
    %     t_BMS = zeros(n_n_FIs,1);
    %     n_dof_BMS = zeros(n_n_FIs,1);
        E_B_BMS = zeros(1,n_n_FIs);
        if BMS_disp
            for i = 1:n_n_FIs

                clear options_BMS
                options_BMS.n_FI = n_FIs(i);
                options_BMS.InteriorMethod      = 'CB+';

        %         options_BMS.BoundaryMethod      = 'weak';
        %         options_BMS.preOrthoTypeLIRWeak = 'qr';
        %         options_BMS.svdTolLIRWeak       = 1e-1;       
        %         options_BMS.n_CC                = 100;

                options_BMS.BoundaryMethod      = Boundary_Methods{q1};
%                 options_BMS.orthoTypeExact      = 'svd'; 
%                 options_BMS.svdTolLIRExact      = 1e-5;         
                options_BMS.n_CC                = n_LIs(q1);

                t_start = tic;

                [H_BMS{i},S_BMS{i},dof_sets_BMS{i},t_up_front,T_BMS{i}] = ...
                        BMS(H0f,Sf,coordinates,R,options_BMS);

                [omega_BMS,phi_BMS,t_kloop_BMS{i,q1}] = ...
                    dispersion_solver_w_k(kappa_use,H_BMS{i},S_BMS{i},dof_sets_BMS{i},R,n_curves);


                t_BMS(i,q1) = toc(t_start);

                T_per_BMS{i} = Periodic_Boundary_Conditions(dof_sets_BMS{i});
                n_dof_BMS(i,q1) = size(T_per_BMS{i}.s0,2);

                PHI_BMS{i} = phi_BMS;
                E_BMS(:,:,i) = omega_BMS.^2;
                
                % band energy error                
                E_B_BMS(i) = bandenergy(E_BMS(:,:,i),f_occ,w_use)/n_atom;
                e_E_B_BMS(i,q1) = (E_B_BMS(i)-E_B_ref);
                
                disp([num2str(n_FIs(i)),' fixed interface modes, ',...
                      num2str(t_BMS(i)),' s'])
            end
        end
       

        disp(['Band energy ref        = ',num2str(E_B_ref*Ha_Per_Ry),' Ha/atom'])
        disp(['Band energy full         = ',num2str(E_B_BMS*Ha_Per_Ry),' Ha/atom'])
        disp(['Band energy BMS         = ',num2str(E_B_BMS*Ha_Per_Ry),' Ha/atom'])
        disp(['full Band energy error  = ',num2str(e_E_B_full*Ha_Per_Ry),' Ha/atom'])     
        disp(['BMS Band energy error  = ',num2str(e_E_B_BMS(:,q1)'*Ha_Per_Ry),' Ha/atom']) 

%     end
end

%% Plot Error vs # Degrees of Freedom
% ======================================================================= %

figure(121);clf
colorbank = get(gca,'colororder');
linestyles = {'-','--',':'};
markers = '.';
count = 1;
for q1 = 1:list_size
    h1 = semilogy(n_dof_BMS(:,q1),e_E_B_BMS(:,q1)*Ha_Per_Ry);hold on
    set(h1,'color',colorbank(q1,:));
    set(h1,'marker','.')
    set(h1,'linestyle',linestyles{1})
    set(h1,'linewidth',1.5)

    legendstrings{count} = ['BMS Band Energy Error Model ',num2str(q1),', ',num2str(n_LIs(q1)), 'CC modes'];
    count=count+1;

    h2 = semilogy(n_dof_per_store(q1),e_E_B_full(q1)*Ha_Per_Ry,'x');
    set(h2,'linewidth',1.5)
    
    legendstrings{count} = ['FE Band Energy Error',num2str(q1)];
    count=count+1;
    set(h2,'color',colorbank(q1,:));
end    
legend(legendstrings)
xlabel('# of Degrees of Freedom')
ylabel('Error')


%% Plot Error vs Computation Time
% ======================================================================= %

figure(122);clf
colorbank = get(gca,'colororder');
lightcolorbank1 = 1-(1-colorbank)*0.2;
lightcolorbank2 = 1-(1-colorbank)*0.15;
linestyles = {'-','--',':','-.'};
markers = '.';


% load and plot PW results
if true
    PW_results = load('save_data/Silicon_Electronic_PW_Error_vs_Time');
    
    for i = 1:size(PW_results.t_soln_save,2)
        h1 = loglog(PW_results.t_soln_save(:,i),PW_results.e_E_b(:,i)*Ha_Per_Ry,'.-');hold on
        set(h1,'color',lightcolorbank1(i,:))
        %set(h1,'color',[0.9,0.9,0.9])
        set(h1,'linestyle','-')
        set(h1,'linewidth',3)
        set(h1,'marker','.')
    end
end

% load and plot FE results
if true
    FE_results = load('save_data/Silicon_Electronic_FE_Error_vs_Time');
    
    % FE results are computed for 3rd - 7th order elements. Keep just
    % 4th-6th.
    FE_results.t_full_store = FE_results.t_full_store(2:end-1);
    FE_results.e_E_b = FE_results.e_E_b(2:end-1);
    
%     FE_results.t_full_store = FE_results.t_full_store(1);
%     FE_results.e_E_b = FE_results.e_E_b(1);
    
    for i = 1:size(FE_results.t_full_store,2)
        for j = 1:size(FE_results.t_full_store{i},2) 
            h1 = loglog(FE_results.t_full_store{i}(:,j),FE_results.e_E_b{i}(:,j)*(Ha_Per_Ry),'.-');hold on
            set(h1,'color',lightcolorbank2(j,:))
            %set(h1,'color',[0.95,0.95,0.95])
            set(h1,'linestyle',linestyles{i+1})
            set(h1,'linewidth',4)
            set(h1,'marker','.')
        end
    end
end

count = 1;
for q1 = 1:list_size
    h1 = loglog(t_BMS(:,q1),e_E_B_BMS(:,q1)*Ha_Per_Ry);hold on
    set(h1,'color',colorbank(q1,:));
    set(h1,'marker','.')
    set(h1,'linestyle',linestyles{1})
    set(h1,'linewidth',1.5)
    legendstrings{count} = ['BMS Band Energy Error Model ',num2str(q1),', ',num2str(n_LIs(q1)), 'CC modes'];
    count=count+1;
    
    h2 = semilogy(t_full(q1),e_E_B_full(q1)*Ha_Per_Ry,'x','markersize',15);
    set(h2,'linewidth',1.5)
    legendstrings{count} = ['FE Band Energy Error',num2str(q1)];
    count=count+1;
    set(h2,'color',colorbank(q1,:));
end    
%legend(legendstrings)
xlabel('computation time (s)')
ylabel('Band Energy Error')
grid on;grid minor;grid minor
xlim([10^0,10^3]);ylim([10^-6,10^-2])

if save_figures
    name='Silicon_Electronic_BMS_BandEnergyError_vs_Time'; %what you want the file to be called
    path='figures/'; %the folder where you want to put the file

    % save .FIG file
    saveas(gcf,[path,name])
end