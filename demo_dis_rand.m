clear;
%Example demonstration of the ranodm initialization for disjoint communoty detection

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
    [topo_mem_init, ~] = NMF_rand_init(adj, num_topo_clus); %Initialize the topology cluster membership matrix, i.e., X
	[~, att_desc_init] = NMF_rand_init(att, num_att_clus); %Initialize the attribute description matrix, i.e., Z
	trans_TA_init = rand(num_topo_clus, num_att_clus); %Initialize the T-A transition matrix, i.e., U
    [topo_mem_TA,att_desc_TA,trans_TA,obj] = DHCD_TA(adj,att,topo_mem_init,att_desc_init,trans_TA_init,alpha,lambd,max_iter,min_error);
    fprintf('T-A Obj. %8.4f\n', obj);
    %==========
    [~, labels] = max(topo_mem_TA, [], 2); %Extract the community detection result from X
    chNMI = get_Ch_NMI(labels, gnd, sample_rate, num_ch_runs); %Derive the corresponding Ch-NMI
    fprintf('Ch-NMI %f\n', chNMI);
    %==========
    %Evaluate the performance of current community detection result
    NMI_TA = compute_NMI(gnd, labels);
    res = bestMap(gnd, labels);
    AC_TA = length(find(gnd == res))/length(gnd);
    fprintf('T-A Ch. NMI: %8.4f; AC: %8.4f\n', NMI_TA, AC_TA);
    fid = fopen('DHCD_T-A(dis_rand).txt', 'at');
    fprintf(fid, 'T-A Ch. Alpha: %8.4f; Obj: %8.4f; Ch-NMI: %8.4f; NMI: %8.4f; AC: %8.4f\n', [alpha, obj, chNMI, NMI_TA, AC_TA]);
    fclose(fid);
    %====================
    %Independently run the DHCD T-A algorihtm multiple times
    for t=2:num_runs
        [topo_mem_init, ~] = NMF_rand_init(adj, num_topo_clus); %Initialize the topology cluster membership matrix, i.e., X
        [~, att_desc_init] = NMF_rand_init(att, num_att_clus); %Initialize the attribute cluster membership matrix, i.e., Z
        trans_TA_init = rand(num_topo_clus, num_att_clus); %Initialize the T-A transition matrix, i.e., U
        [cur_topo_mem_TA,cur_att_desc_TA,cur_trans_TA,cur_obj] = DHCD_TA(adj,att,topo_mem_init,att_desc_init,trans_TA_init,alpha,lambd,max_iter,min_error);
        fprintf('T-A Obj. %8.4f\n', cur_obj);
        %==========
        [~, labels] = max(cur_topo_mem_TA, [], 2); %Extract the community detection result from X
        chNMI = get_Ch_NMI(labels, gnd, sample_rate, num_ch_runs); %Derive the corresponding Ch-NMI
        %==========
        %Evaluate the performance of current community detection result
        fprintf('Ch-NMI %f\n', chNMI);
        NMI_TA = compute_NMI(gnd, labels);
        res = bestMap(gnd, labels);
        AC_TA = length(find(gnd == res))/length(gnd);
        fprintf('T-A Ch. NMI: %8.4f; AC: %8.4f\n', NMI_TA, AC_TA);
        %==========
        fid = fopen('DHCD_T-A(dis_rand).txt', 'at');
        fprintf(fid, 'T-A Ch. Alpha: %8.4f; Obj: %8.4f; Ch-NMI: %8.4f; NMI: %8.4f; AC: %8.4f\n', [alpha, cur_obj, chNMI, NMI_TA, AC_TA]);
        fclose(fid);
        %====================
        %Update the best community detection result based on the value of objective function
        if cur_obj<obj
            topo_mem_TA = cur_topo_mem_TA;
            att_desc_TA = cur_att_desc_TA;
            trans_TA = cur_trans_TA;
            obj = cur_obj;
        end
    end
    %====================
    %Evalute the performance of the best community detection result
    [~, labels] = max(topo_mem_TA, [], 2);
    NMI_TA = compute_NMI(gnd, labels);
    res = bestMap(gnd, labels);
    AC_TA = length(find(gnd == res))/length(gnd);
    fprintf('(Obj. Min.)T-A Ch. NMI: %8.4f; AC: %8.4f\n', NMI_TA, AC_TA);
    fprintf('====================\n');
    %===========
	fid = fopen('DHCD_T-A(dis_rand).txt', 'at');
	fprintf(fid, 'Obj. Min. T-A Ch. Alpha: %8.4f; Obj: %8.4f; NMI: %8.4f; AC: %8.4f\n', [alpha, obj, NMI_TA, AC_TA]);
    fprintf(fid, '====================\n');
	fclose(fid);

    %==============================
    %A-T Channel, i.e., DHCD A-T
    beta = alpha;
    fprintf('Beta: %f\n', beta);
    %====================
    %Initialize the attribute cluster membership matrix & attribute description matrix, i.e., Y & Z
    [att_mem_init, att_desc_init] = NMF_rand_init(att, num_att_clus);
    trans_AT_init = rand(num_att_clus, num_topo_clus); %Initialize the A-T transition matrix, i.e., V
    [att_desc_AT,att_mem_AT,trans_AT,obj] = DHCD_AT(adj,att,att_mem_init,att_desc_init,trans_AT_init,beta,lambd,max_iter,min_error);
    fprintf('A-T Obj. %8.4f\n', obj);
    %==========
    [~, labels] = max(att_mem_AT, [], 2); %Extract the community detection result from Y
    chNMI = get_Ch_NMI(labels, gnd, sample_rate, num_ch_runs); %Derive the corresponding Ch-NMI
    fprintf('Ch-NMI %f\n', chNMI);
    %==========
    %Evaluate the performance of current community detection result
    NMI_AT = compute_NMI(gnd, labels);
    res = bestMap(gnd, labels);
    AC_AT = length(find(gnd == res))/length(gnd);
    fprintf('A-T Ch. NMI: %8.4f; AC: %8.4f\n', NMI_AT, AC_AT);
    fid = fopen('DHCD_A-T(dis_rand).txt', 'at');
    fprintf(fid, 'A-T Ch. Beta: %8.4f; Obj: %8.4f; Ch-NMI: %8.4f; NMI: %8.4f; AC: %8.4f\n', [beta, obj, chNMI, NMI_AT, AC_AT]);
    fclose(fid);
    %====================
    %Independently run the DHCD A-T algorihtm multiple times
    for t=2:num_runs
        %Initialize the attribute cluster membership matrix & attribute description matrix, i.e., Y & Z
        [att_mem_init, att_desc_init] = NMF_rand_init(att, num_att_clus);
        trans_AT_init = rand(num_att_clus, num_topo_clus); %Initialize the A-T transition matrix, i.e., V
        [cur_att_desc_AT,cur_att_mem_AT,cur_trans_AT,cur_obj] = DHCD_AT(adj,att,att_mem_init,att_desc_init,trans_AT_init,beta,lambd,max_iter,min_error);
        fprintf('A-T Obj. %8.4f\n', cur_obj);
        %==========
        [~, labels] = max(cur_att_mem_AT, [], 2); %Extract the community detection result from Y
        chNMI = get_Ch_NMI(labels, gnd, sample_rate, num_ch_runs); %Derive the corresponding Ch-NMI
        fprintf('Ch-NMI %f\n', chNMI);
        %==========
        %Evaluate the performance of current community detection result
        NMI_AT = compute_NMI(gnd, labels);
        res = bestMap(gnd, labels);
        AC_AT = length(find(gnd == res))/length(gnd);
        fprintf('A-T Ch. NMI: %8.4f; AC: %8.4f\n', NMI_AT, AC_AT);
        fid = fopen('DHCD_A-T(dis_rand).txt', 'at');
        fprintf(fid, 'A-T Ch. Beta: %8.4f; Obj: %8.4f; Ch-NMI: %8.4f; NMI: %8.4f; AC: %8.4f\n', [beta, cur_obj, chNMI, NMI_AT, AC_AT]);
        fclose(fid);
        %===================
        %Update the best community detection result based on the value of objective function
        if cur_obj<obj
            att_mem_AT = cur_att_mem_AT;
            att_desc_AT = cur_att_desc_AT;
            trans_AT = cur_trans_AT;
            obj = cur_obj;
        end
    end
    %======================
    %Evalute the performance of the best community detection result
    [~, labels] = max(att_mem_AT, [], 2);
    NMI_AT = compute_NMI(gnd, labels);
    res = bestMap(gnd, labels);
    AC_AT = length(find(gnd == res))/length(gnd);
    fprintf('(Obj. Min.)A-T Ch. NMI: %8.4f; AC: %8.4f\n', NMI_AT, AC_AT);
    fprintf('====================\n');
    %==========
    fid = fopen('DHCD_A-T(dis_rand).txt', 'at');
	fprintf(fid, 'Obj. Min. A-T Ch. Beta: %8.4f; Obj: %8.4f; NMI: %8.4f; AC: %8.4f\n', [beta, obj, NMI_AT, AC_AT]);
    fprintf(fid, '=====================\n');
	fclose(fid);
end
%delete(p);
