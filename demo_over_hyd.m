clear;
%Example demonstration of the hybrid initialization (random & NNDSVD initialization) for overlapping communoty detection

%====================
%Load the network datasets
adj = load('data/Reddit26_adj.txt'); %Adjacency matrix
att = load('data/Reddit26_node_attri.txt'); %Network attribtue matrix
gnd_path = 'data/Reddit26_gnd.txt'; %Path of ground-truth
res_path_TA = 'Reddit26_res_TA(over_hyd).txt'; %Path to save the (overlapping community detection) result of DHCD T-A
res_path_AT = 'Reddit26_res_AT(over_hyd).txt'; %Path to save the (overlapping community detection) reuslt of DHCD A-T
%==========
%adj = load('data/Reddit25_adj.txt');
%att = load('data/Reddit25_node_attri.txt');
%gnd_path = 'data/Reddit25_gnd.txt';
%res_path_TA = 'Reddit25_res_TA(over_hyd).txt';
%res_path_AT = 'Reddie25_res_AT(over_hyd).txt';

%====================
%Get the network parameters
num_nodes = size(adj,1); %Number of nodes
num_atts = size(att, 2); %Number of node attributes
num_topo_clus = 3; %Number of topology clusters
num_att_clus = num_topo_clus; %Number of attribute clusters

%====================
max_iter = 1e4; %Maximum number of optimization iteratons
min_error = 1e-5; %Minimum relative error to determine the convergence of optimization
%==========
num_runs = 10; %Number of independent runs of DHCD
lambd = 1e3;

%====================
%Remove self-connectd edges
for i=1:num_nodes
    adj(i, i) = 0; 
end
%==========
adj = sparse(adj);
att = sparse(att);

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
    [~, att_desc_init] = NNDSVD(att, num_att_clus, 0); %Initialize the attribute description matrix, i.e., Z
    att_desc_init = att_desc_init';
    trans_TA_init = rand(num_topo_clus, num_att_clus); %Initialize the T-A transition matrix, i.e., U
    [topo_mem_TA,att_desc_TA,trans_TA,obj] = DHCD_TA(adj,att,topo_mem_init,att_desc_init,trans_TA_init,alpha,lambd,max_iter,min_error);
    fprintf('T-A Obj. %8.4f\n', obj);
    %==========
    %Extract the disjoint community detection result
    [~, dis_labels] = max(topo_mem_TA*trans_TA, [], 2);
    %Extract and save the overlapping community detection result
    save_over_res(dis_labels, trans_TA, res_path_TA);
    %Evaluate the performance of current community detection result
    [Fsc, Jac] = over_eva(gnd_path, res_path_TA);
    fprintf('T-A Ch. F-score: %8.4f; Jaccard: %8.4f\n', [Fsc, Jac]);
	fid = fopen('DHCD_T-A(over_hyd).txt', 'at');
	fprintf(fid, 'T-A Ch. Alpha: %8.4f; Obj: %8.4f; F-Score: %8.4f; Jaccard: %8.4f\n', [alpha, obj, Fsc, Jac]);
	fclose(fid);
	%====================
    %Independently run the DHCD T-A algorihtm multiple times
    for t=2:num_runs
        trans_TA_init = rand(num_topo_clus, num_att_clus); %Initialize the T-A transition matrix, i.e., U
        [cur_topo_mem_TA,cur_att_desc_TA,cur_trans_TA,cur_obj] = DHCD_TA(adj,att,topo_mem_init,att_desc_init,trans_TA_init,alpha,lambd,max_iter,min_error);
        fprintf('T-A Obj. %8.4f\n', cur_obj);
        %==========
        %Extract the disjoint community detection result
        [~, dis_labels] = max(cur_topo_mem_TA*cur_trans_TA, [], 2);
        %Extract and save the overlapping community detection result
        save_over_res(dis_labels, cur_trans_TA, res_path_TA);
        %Evaluate the performance of current community detection result
        [Fsc, Jac] = over_eva(gnd_path, res_path_TA);
        fprintf('T-A Ch. F-score: %8.4f; Jaccard: %8.4f\n', [Fsc, Jac]);
        fid = fopen('DHCD_T-A(over_hyd).txt', 'at');
        fprintf(fid, 'T-A Ch. Alpha: %8.4f; Obj: %8.4f; F-Score: %8.4f; Jaccard: %8.4f\n', [alpha, cur_obj, Fsc, Jac]);
        fclose(fid);
        %==========
        %Update the best community detection result based on the value of objective function
        if cur_obj<obj
            topo_mem_TA = cur_topo_mem_TA;
            att_desc_TA = cur_att_desc_TA;
            trans_TA = cur_trans_TA;
            obj = cur_obj;
        end
    end
    %======================
    %Evalute the performance of the best community detection result
    [~, dis_labels] = max(topo_mem_TA*trans_TA, [], 2);
    save_over_res(dis_labels, trans_TA, res_path_TA);
    [Fsc, Jac] = over_eva(gnd_path, res_path_TA);
    fprintf('Min. Obj. T-A Ch. F-score: %8.4f; Jaccard: %8.4f\n', [Fsc, Jac]);
    fprintf('====================\n');
	fid = fopen('DHCD_T-A(over_hyd).txt', 'at');
	fprintf(fid, 'Obj. Min. T-A Ch. Alpha: %8.4f; Obj: %8.4f; F-Score: %8.4f; Jaccard: %8.4f\n', [alpha, obj, Fsc, Jac]);
    fprintf(fid, '====================\n');
	fclose(fid);

    %==============================
    %T-A Channel, i.e., DHCD T-A
    beta = alpha;
    fprintf('Beta: %f\n', beta);
    %====================
    %Initialize the attribute cluster membership matrix & attribute description matrix, i.e., Y & Z
    [att_mem_init, att_desc_init] = NNDSVD(att, num_att_clus, 0);
    att_desc_init = att_desc_init';
    trans_AT_init = rand(num_att_clus, num_topo_clus); %Initialize the A-T transition matrix, i.e., V
    [att_desc_AT,att_mem_AT,trans_AT,obj] = DHCD_AT(adj,att,att_mem_init,att_desc_init,trans_AT_init,beta,lambd,max_iter,min_error);
    fprintf('A-T Obj. %8.4f\n', obj);
    %==========
    %Extract the disjoint community detection result
    [~, dis_labels] = max(att_mem_AT*trans_AT, [], 2);
    %Extract and save the overlapping community detection result
    save_over_res(dis_labels, trans_AT, res_path_AT);
    %Evaluate the performance of current community detection result
    [Fsc, Jac] = over_eva(gnd_path, res_path_AT);
    fprintf('A-T Ch. F-score: %8.4f; Jaccard: %8.4f\n', [Fsc, Jac]);
	fid = fopen('DHCD_A-T(over_hyd).txt', 'at');
	fprintf(fid, 'A-T Ch. Beta: %8.4f; Obj: %8.4f; F-Score: %8.4f; Jaccard: %8.4f\n', [beta, obj, Fsc, Jac]);
	fclose(fid);
    %====================
	%Independently run the DHCD A-T algorihtm multiple times
    for t=2:num_runs
        trans_AT_init = rand(num_att_clus, num_topo_clus); %Initialize the A-T transition matrix, i.e., V
        [cur_att_desc_AT,cur_att_mem_AT,cur_trans_AT,cur_obj] = DHCD_AT(adj,att,att_mem_init,att_desc_init,trans_AT_init,beta,lambd,max_iter,min_error);
        fprintf('A-T Obj. %8.4f\n', cur_obj);
        %==========
        %Extract the disjoint community detection result
        [~, dis_labels] = max(cur_att_mem_AT*cur_trans_AT, [], 2);
        %Extract and save the overlapping community detection result
        save_over_res(dis_labels, cur_trans_AT, res_path_AT);
        %Evaluate the performance of current community detection result
        [Fsc, Jac] = over_eva(gnd_path, res_path_AT);
        fprintf('A-T Ch. F-score: %8.4f; Jaccard: %8.4f\n', [Fsc, Jac]);
        fid = fopen('DHCD_A-T(over_hyd).txt', 'at');
        fprintf(fid, 'A-T Ch. Beta: %8.4f; Obj: %8.4f; F-Score: %8.4f; Jaccard: %8.4f\n', [beta, cur_obj, Fsc, Jac]);
        fclose(fid);
        %==========
        %Update the best community detection result based on the value of objective function
        if cur_obj<obj
            att_mem_AT = cur_att_mem_AT;
            att_desc_AT = cur_att_desc_AT;
            trans_AT = cur_trans_AT;
            obj = cur_obj;
        end
    end
    %====================
    [~, dis_labels] = max(att_mem_AT*trans_AT, [], 2);
    save_over_res(dis_labels, trans_AT, res_path_AT);
    [Fsc, Jac] = over_eva(gnd_path, res_path_AT);
    fprintf('Min. Obj. A-T Ch. F-score: %8.4f; Jaccard: %8.4f\n', [Fsc, Jac]);
    fprintf('====================\n');
	fid = fopen('DHCD_A-T(over_hyd).txt', 'at');
	fprintf(fid, 'Obj. Min. A-T Ch. Beta: %8.4f; Obj: %8.4f; F-Score: %8.4f; Jaccard: %8.4f\n', [beta, obj, Fsc, Jac]);
    fprintf(fid, '====================\n');
	fclose(fid);
end
