function save_over_res(dis_labels, trans, res_path)
%Function to extract & save the overlapping community detection result
%dis_labels: label sequence of disjoint community detection
%trans: transition matrix
%res_path: path to save the overlapping community detection result

    %====================
    [num_nodes, ~] = size(dis_labels); %Number of nodes
    num_clus = max(dis_labels); %Number of communities
    %==========
    %Normalize each column of the transition matrix
    [N, M] = size(trans);
    for j=1:M
        trans(:, j) = trans(:, j)/max(sum(trans(:, j)), realmin);
        if norm(trans(:, j))==0
            for i=1:N
                trans(i, j) = 1.0/N;
            end
        end
    end
    
    %====================
    %Label sequences of overlapping community detection
    over_labels = cell(num_clus); 
    %Initialize the lebale sequences
    for c=1:num_clus
        over_labels{c} = []; 
    end
    mem_thres = 1.0/num_clus; %Threshold to determine the overlapping community membership 
    %====================
    %Derive the label sequence w.r.t. each community
    for i=1:num_nodes
        dis_label = dis_labels(i);
        for c=1:num_clus
            if trans(c, dis_label)>=mem_thres
                over_labels{c} = [over_labels{c}, i];
            end
        end
    end
    %==========
    for c=1:num_clus
        if length(over_labels{c})==0
            for k=1:num_nodes
                if dis_labels(k)==c
                    over_labels{c} = [over_labels{c}, k];
                end
            end
        end
    end
    
    %====================
    %Save the extracted overlapping community detection result
    fid = fopen(res_path, 'w+');
    for c=1:num_clus
       [~, len] = size(over_labels{c});
       for j=1:len
          fprintf(fid, '%d ', over_labels{c}(j)); 
       end
       fprintf(fid, '\n'); 
    end
    fclose(fid);
end

