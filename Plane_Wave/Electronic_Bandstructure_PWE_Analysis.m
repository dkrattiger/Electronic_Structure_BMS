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


% Solution Options
n_bands      = 8;       % number of dispersion curves to plot
full_disp    = true;  	% 1 = compute full dispersion, 0 = don't
RBME_disp    = false;  	% Compute RBME dispersion
plot_disp    = false;


% size of k-point grid for monkhorst pack
use_ibz = false;
if ~plot_disp
    n_kap_1 = 8;
%     n_kap_1 = 12;  
%     n_kap_1 = 16; 
%     n_kap_1 = 20;  
    use_ibz = true;
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

% a1 = (a0)*[1;0;0]*n_cellx;
% a2 = (a0)*[0;1;0]*n_celly;
% a3 = (a0)*[0;0;1]*n_cellz;


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

%% Real-Space Potential Grid
% ======================================================================= %

% real space grid dimensions
n = 41;
nx = n*n_cellx;
ny = n*n_celly;
nz = n*n_cellz;

% define grid
a1vec = linspace(0,1,nx+1);a1vec = a1vec(1:end-1);
a2vec = linspace(0,1,ny+1);a2vec = a2vec(1:end-1);
a3vec = linspace(0,1,nz+1);a3vec = a3vec(1:end-1);

[a1grid,a2grid,a3grid] = ndgrid(a1vec,a2vec,a3vec);

coordinates = [a1grid(:),a2grid(:),a3grid(:)]*R';
xgrid = a1(1)*a1grid + a2(1)*a2grid + a3(1)*a3grid;
ygrid = a1(2)*a1grid + a2(2)*a2grid + a3(2)*a3grid;
zgrid = a1(3)*a1grid + a2(3)*a2grid + a3(3)*a3grid;
Vgrid = reshape(Si_Pseudo(coordinates',a0),[nx,ny,nz]);


%% Form PW discretization of Hamiltonian

runselect = 'big';
switch runselect
    case 'small'
        nPWs = [7:2:9];
        n_nPWs = length(nPWs);
        n_kaps = [6,7,8];
        n_MPs = length(n_kaps);
    case 'big'
        nPWs = [5:1:11];
        n_nPWs = length(nPWs);
        n_kaps = [7,8,9,10];
        n_MPs = length(n_kaps);
    case 'ref'
        nPWs = [19];
        n_nPWs = length(nPWs);
        n_kaps = [20];
        n_MPs = length(n_kaps);
end

% preallocate timing results;
t_model_save = zeros(n_nPWs,n_MPs);
t_soln_save = zeros(n_nPWs,n_MPs);

for qq2 = 1:n_MPs
    
    % Discretize Brillouin zone
    % =================================================================== %
    n_kap_1 = n_kaps(qq2);

    n_MP = n_kap_1*ceil([1/n_cellx,1/n_celly,1/n_cellz]);
    [~,w_ibz,kappa_full,i_ibz,~] = MonkhorstPackGrid(b1,b2,b3,n_MP);
    %length(i_ibz)
    %n_kap_1^3/length(i_ibz)

    if use_ibz
        kappa = kappa_full(:,i_ibz);
        w_use = w_ibz;
    else
        kappa  = kappa_full;
        w_use = ones(1,size(kappa,2))/size(kappa,2);
    end

    n_kap = size(kappa,2);

    % plotting parameters
    kappa_plot = 1:n_kap;
    n_segs = 1;
    
    % preallocate Energy solutions
    E_save{qq2} = zeros(n_bands,n_kap,n_nPWs);

    for qq1 = 1:n_nPWs

        % Number of Plane-Wave Basis functions along each dimension
        P = nPWs(qq1);
        Q = nPWs(qq1);
        R = nPWs(qq1);
        n_dof = P*Q*R;
        
        % full model description string
        modeldescription = [sprintf('Si_%ix%ix%iPrimCell_PW_%iDOF_',...
                        n_cellx,n_celly,n_cellz,n_dof)];
                    
        % BZ or IBZ in file names
        if use_ibz     
            BZ_or_IBZ = 'IBZ';  
        else
            BZ_or_IBZ = 'BZ';  
        end

        % define model save path
        solutionpathstring = 'save_data/solutions/';

        solutiondescription = [sprintf('%i',n_kap),BZ_or_IBZ,sprintf('kpts_%iBands',n_bands)];

        solution_savestring = [solutionpathstring,modeldescription,solutiondescription];
        
        % Note, Model creation and solution are not separated because
        % models are not stored (due to memory) after they are evaluated
        if exist([solution_savestring,'.mat'],'file') && load_results
%             load([solution_savestring,'.mat'],'E_disp','t_kloop')
            load([solution_savestring,'.mat'],'E','PHI','n_dof','t_soln_k','t_model_k')
%             'E','PHI','n_dof','t_soln_k','t_model_k'
        else
            Vf = (fftn(Vgrid)) / (nx*ny*nz);

            % center of fft array
            p = [-floor(P/2):floor(P/2)];
            q = [-floor(Q/2):floor(Q/2)];
            r = [-floor(R/2):floor(R/2)];

            % centered location
            p0 = 1+floor(nx/2);
            q0 = 1+floor(ny/2);
            r0 = 1+floor(nz/2);

            H = zeros(n_dof,n_dof);
            E = zeros(n_bands,n_kap);
            PHI = zeros(n_dof,n_bands,n_kap);
            t_soln_k = zeros(n_kap,1);

            for i = 1:n_kap
                tstartmodel = tic;
                for rrow = 1:R
                for qrow = 1:Q
                for prow = 1:P
                    row = (rrow-1)*Q*P + (qrow-1)*P + prow;

                    if false

                        % Looped assembly
                        for rcol = 1:R
                        for qcol = 1:Q
                        for pcol = 1:P
                            col = (rcol-1)*Q*P + (qcol-1)*P + pcol;
                            pfft = p(prow)-p(pcol);
                            qfft = q(qrow)-q(qcol);
                            rfft = r(rrow)-r(rcol);
                            H(row,col) = Vf(p0+pfft,q0+qfft,r0+rfft);
                            if row==col
                                qvec = kappa(:,i)+ p(prow)*b1 + q(qrow)*b2 + r(rrow)*b3;
                                H(row,col) = H(row,col) + qvec'*qvec;
                            end        
                        end
                        end
                        end
                    else

                        % Vectorized Indexing to form entire row of Hamiltonian matrix
                        % at once
                        Vkeep = Vf(p0+p(prow)-p,q0+q(qrow)-q,r0+r(rrow)-r);
                        H(row,:) = Vkeep(:);
                        qvec = kappa(:,i)+ p(prow)*b1 + q(qrow)*b2 + r(rrow)*b3;
                        H(row,row) = H(row,row) + qvec'*qvec;
                    end
                end
                end
                end        
                t_model_k(i) = toc(tstartmodel);


                % solve for eigenvalues (energies)
                % ======================================================= %
                tstartsoln = tic;
                [PHI(:,:,i),L] = eigs(H,n_bands,'sm');
        %         [PHI,L] = eig(H(:,:,i));
                [L,isort] = sort(diag(L));
                E(:,i) = L(1:n_bands);
                PHI(:,:,i) = PHI(:,isort,i);
                t_soln_k(i) = toc(tstartsoln);
                
                fprintf('k point %i of %i, model time: %4.2fs, solution time: %4.2s\n',...
                    i,n_kap,t_model_k(i),t_soln_k(i))
            end
            
            if save_results
                save(solution_savestring,'E','PHI','n_dof','t_soln_k','t_model_k')
            end
        end

        t_model_save(qq1,qq2) = sum(t_model_k(:));
        t_soln_save(qq1,qq2) = sum(t_soln_k(:));


        % store results
        f_occ = [ones(4,1);zeros(n_bands-4,1)];     
        E_save{qq2}(:,:,qq1) = E;
        E_b_save(qq1,qq2) = bandenergy(E,f_occ,w_use)/n_atoms;    
        PHI_save{qq1,qq2} = PHI;
        n_dof_save(qq1,qq2) = n_dof;

        % Display loop info
        fprintf('nPWs = %i, modeling time = %4.2f, solution time = %4.2f\n\n',nPWs(qq1),t_model_save(qq1),t_soln_save(qq1));
    end
end

%% choose a model to use for error reference
% ======================================================================= %

ref_model = 'save_data/solutions/Si_1x1x1PrimCell_PW_6859DOF_256IBZkpts_8Bands';
if exist([ref_model,'.mat'],'file')    
    n_atoms_ref = 2;
    n_kap_ref_1 = 20;
    n_MP = n_kap_ref_1*ceil([1/n_cellx,1/n_celly,1/n_cellz]);
    [~,w_use,~,~,~] = MonkhorstPackGrid(b1,b2,b3,n_MP);
    load([ref_model,'.mat'],'E','PHI','n_dof','t_soln_k','t_model_k')
end

% use last result (should be most converged one)
n_atoms_ref = 2;
E_ref = E;
w_ref = w_use;

% total band energy
f_occ = [ones(4,1);zeros(n_bands-4,1)];          
E_b_ref = bandenergy(E_ref,f_occ,w_ref)/n_atoms_ref;
n_kap_ref = size(E_ref,2);
E_b_ref_k = zeros(1,n_kap);
for i = 1:n_kap_ref
    E_b_ref_k(i) = f_occ'*E_ref(:,i)/n_atoms_ref;
end

%% Plot Last dispersion
% ======================================================================= %

figure(2);clf
plot(kappa_plot,sort(E_save{end}(:,:,end)))
xlim(kappa_plot([1,end]))
drawnow

%% Compute and plot E_b error vs model size
% ======================================================================= %
e_E_b = zeros(qq1,qq2);
for qq1 = 1:n_nPWs
    % compute error
    e_E_b(qq1,:) = E_b_save(qq1,:)-E_b_ref;
end

colorbank = get(groot,'defaultaxescolororder');

% plot error vs model size for each MP grid and element type
for i = 1:3
    figure(200+i);clf
    switch i
        case 1
            xlab = 'n_{DOF}';
            h1 = loglog(n_dof_save,e_E_b*Ha_Per_Ry,'.-');hold on        
            save_string_end = 'NDOF';
        case 2
            h1 = semilogy(n_dof_save.^(1/3),e_E_b*Ha_Per_Ry,'.-');hold on
            xlab = 'nodes per edge';                
            save_string_end = 'NodesPerEdge';
        case 3
            h1 = loglog(t_soln_save,e_E_b*Ha_Per_Ry,'.-');hold on
            xlab = 'computation time';
            save_string_end = 'Time';
    end
    
    grid on

    %set(h1,'color',colorbank(rem(qq2-1,size(colorbank,1))+1,:))
    set(h1,'linestyle','--')
    set(h1,'linewidth',2)
    set(h1,'marker','.')
    set(h1,'markersize',8)

%     legendstrings{(q1-1)*length(MPrule_list)+q3} = [orders{elem_order_list(q1)},'-order brick elements, ',...
%                          num2str(MPrule_list(q3)),'\times',...
%                          num2str(MPrule_list(q3)),'\times',...
%                          num2str(MPrule_list(q3)),'k-pts'];

    grid on
    grid minor,grid minor
    xlabel(xlab)
%     switch e_B_type
%         case 'max'
%             ylabel('Max band energy error, e_{E_b} (Ha/atom)')
%         case 'avg'
            ylabel('Avg band energy error, e_{E_b} (Ha/atom)')
%     end
%     legend(legendstrings)

end
% end

if false
    save('save_data/Silicon_Electronic_PW_Error_vs_Time','t_soln_save','e_E_b')
end


%% plot isosurface contours for potential
% ======================================================================= %

if true
    C2 = Si_Pseudo(coordinates',a0);
    C2 = C2*(1/2);


    figure(1);clf;hold on

    maxC = max(C2(:));
    minC = min(C2(:));
    n_isosteps = 5;
    isosteps = linspace(0,1,n_isosteps+2);isosteps = isosteps(2:end-1); 
    isovals = minC + isosteps*(maxC-minC);

    legend_entries = cell(1,length(isosteps));
    h = zeros(1,length(isosteps));
    for i = 1:length(isosteps)
        isoSurfPatchStruct = isosurface(conversionFactor*reshape(coordinates(:,1),[nx,ny,nz]),...
                                        conversionFactor*reshape(coordinates(:,2),[nx,ny,nz]),...
                                        conversionFactor*reshape(coordinates(:,3),[nx,ny,nz]),...
                                        reshape(C2,[nx,ny,nz]),isovals(i));
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
    coordsedge = conversionFactor*[...
        [zeros(3,1),a1,a1+a2,a2,zeros(3,1)],[nan;nan;nan]...
        a3*ones(1,5)+[zeros(3,1),a1,a1+a2,a2,zeros(3,1)],[nan;nan;nan],...
        zeros(3,1),a3,[nan;nan;nan],...
        a1,a3+a1,[nan;nan;nan],...
        a1+a2,a3+a1+a2,[nan;nan;nan],...
        a2,a3+a2,[nan;nan;nan]];
                        
%                         
%                         a1,a1+a2,a2,zeros(3,1)],[nan;nan;nan],...
                        
    
    h_boundary = plot3(coordsedge(1,:),coordsedge(2,:),coordsedge(3,:),'k-','linewidth',1);hold on
    
    
    
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

