clear;
%Example demonstration of the NNDSVD initialization for disjoint community detection

%====================
%Load the network datasets
%file = load('data\washington.mat');
%edges = file.washington.adj;
%adj = Compute_adjmatrix(edges); %Adjacent matrix
%att = file.washington.content; %Network attribute matrix
%gnd = file.washington.groundtruth(:,2); %Ground-truth of community membership
%==========
file = load('data\texas.mat');
edges = file.texas.adj;
adj = Compute_adjmatrix(edges); %Adjacent matrix
att = file.texas.content; %Adjacent matrix
gnd = file.texas.groundtruth(:,2); %Adjacent matrix
%==========
%file = load('data\cora.mat');
%adj = file.G; %Adjacent matrix
%att = file.content; %Network attribute matrix
%gnd = file.labels; %Ground-truth of community membership
%==========
%file = load('data\cite.mat');
%adj = file.G; %Adjacent matrix
%att = file.content; %Network attribute matrix
%gnd = file.labels; %Ground-truth of community membership
%==========
%Get the network parameters
num_nodes = size(adj, 1); %Number of nodes
num_atts = size(att, 2); %Number of node attributes
num_topo_clus = max(gnd); %Number of topology clusters
num_att_clus = num_topo_clus; %Number of attribtue clusters

%====================
max_iter = 1e4; %Maximum number of optimization iteratons
min_error = 1e-5; %Minimum relative error to determine the convergence of optimization
%==========
num_runs = 10; %Number of independent runs of DHCD
num_ch_runs = 10; %Number of independent runs of channel selection
sample_rate = 0.10; %Sampling rate of Ch-NMI
%==========
lambd = 0;

%====================
%Remove self-connectd edges
for i=1:num_nodes
   adj(i, i) = 0; 
end
%==========
adj = sparse(adj);
att = sparse(att);

%p = parpool;
%====================
params = [0.1:0.1:0.9, 1:1:10];
[~, num_params] = size(params); %Number of parameter settings
for l=1:num_params
    %==============================
    %T-A Channel, i.e., DHCD T-A
    alpha = params(l); 
    fprintf('Alpha: %f\n', alpha);
    %====================
    [topo_mem_init, ~] = NNDSVD(adj, num_topo_clus, 0); %Initialize the topology cluster membership matrix, i.e., X
    [A, att_desc_init] = NNDSVD(att, num_att_clus, 0); %Initialize the attribute description matrix, i.e., Z
    att_desc_init = att_desc_init';
    [~, trans_TA_init] = NNDSVD(A, num_topo_clus, 0); %Initialize the T-A transition matrix, i.e., U
    [topo_mem_TA,att_desc_TA,trans_TA,obj] = DHCD_TA(adj,att,topo_mem_init,att_desc_init,trans_TA_init,alpha,lambd,max_iter,min_error);
    fprintf('T-A Obj. %8.4f\n', obj);
    %==========
    [~, labels] = max(topo_mem_TA, [], 2); %Extract the community detection result from X
    chNMI = get_Ch_NMI(labels, gnd, sample_rate, num_ch_runs); %Derive the corresponding Ch-NMI
    fprintf('Ch-NMI %f\n', chNMI);
    %====================
    %Evalute the performance of the best community detection result
    [~, labels] = max(topo_mem_TA, [], 2);
    NMI_TA = compute_NMI(gnd, labels);
    res = bestMap(gnd, labels);
    AC_TA = length(find(gnd == res))/length(gnd);
    fprintf('(Obj. Min.)T-A Ch. NMI: %8.4f; AC: %8.4f\n', NMI_TA, AC_TA);
    fprintf('====================\n');
    %===========
	fid = fopen('DHCD_T-A(dis_SVD).txt', 'at');
	fprintf(fid, 'Obj. Min. T-A Ch. Alpha: %8.4f; Obj: %8.4f; Ch-NMI: %8.4f; NMI: %8.4f; AC: %8.4f\n', [alpha, obj, chNMI, NMI_TA, AC_TA]);
    fprintf(fid, '====================\n');
	fclose(fid);

    %==============================
    %A-T Channel, i.e., DHCD A-T
    beta = alpha;
    fprintf('Beta: %f\n', beta);
    %====================
    %Initialize the attribute cluster membership matrix & attribute description matrix, i.e., Y & Z
    [att_mem_init, att_desc_init] = NNDSVD(att, num_att_clus, 0);
    att_desc_init = att_desc_init';
    [A, ~] = NNDSVD(adj, num_topo_clus, 0);
    [~, trans_AT_init] = NNDSVD(A, num_att_clus, 0); %Initialize the A-T transition matrix, i.e., V
    [att_desc_AT,att_mem_AT,trans_AT,obj] = DHCD_AT(adj,att,att_mem_init,att_desc_init,trans_AT_init,beta,lambd,max_iter,min_error);
    fprintf('A-T Obj. %8.4f\n', obj);
    %==========
    [~, labels] = max(att_mem_AT, [], 2); %Extract the community detection result from Y
    chNMI = get_Ch_NMI(labels, gnd, sample_rate, num_ch_runs); %Derive the corresponding Ch-NMI
    fprintf('Ch-NMI %f\n', chNMI);
    %======================
    %Evalute the performance of the best community detection result
    [~, labels] = max(att_mem_AT, [], 2);
    NMI_AT = compute_NMI(gnd, labels);
    res = bestMap(gnd, labels);
    AC_AT = length(find(gnd == res))/length(gnd);
    fprintf('(Obj. Min.)A-T Ch. NMI: %8.4f; AC: %8.4f\n', NMI_AT, AC_AT);
    fprintf('====================\n');
    %==========
    fid = fopen('DHCD_A-T(dis_SVD).txt', 'at');
	fprintf(fid, 'Obj. Min. A-T Ch. Beta: %8.4f; Obj: %8.4f; Ch. NMI: %8.4f; NMI: %8.4f; AC: %8.4f\n', [beta, obj, chNMI, NMI_AT, AC_AT]);
    fprintf(fid, '=====================\n');
	fclose(fid);
end
%delete(p);
