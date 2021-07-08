function [topo_mem,att_desc,trans_TA,obj] = DHCD_TA(adj,att,topo_mem,att_desc,trans_TA,alpha,lambd,max_iter,min_error)
%Funtion to implement DHCD T-A
%adj: adjacency matrix, i.e., A
%att: network attribute matrix, i.e., C
%topo_mem: topology cluster membership matrix, i.e., X 
%att_desc: attribute description matrix, i.e., Z
%trans_TA: T-A transition matrix, i.e., U
%alpha, lambd: hyper-parameters of the model
%max_iter: maximum number of iterations
%min_error: minimum relative error to determine convergence
%obj: objective function value

    %====================
    %num_topo_clus: number of topology clusters
    %num_att_clus: number of attribute clusters
    [num_topo_clus, num_att_clus] = size(trans_TA);
    iter_cnt = 0; %Counter of iterations
    error = 0; %Relative error of the objective fucntion
    %====================
    %Calculate the objective function value
    obj = norm(adj - topo_mem*topo_mem','fro')^2;
    obj = obj + alpha*norm(att - topo_mem*trans_TA*att_desc', 'fro')^2;
    obj = obj + lambd*norm(ones(1, num_topo_clus)*trans_TA-ones(1, num_att_clus), 'fro')^2;
    %====================
    while iter_cnt==0 || (error >= min_error && iter_cnt <= max_iter)
        pre_obj = obj; %Objective function value in previous iteraton
        %==========
        %X-Process: update the topology cluster membership matrix X
        aux1 = trans_TA*att_desc'; %Auxiliary matrix, i.e., Q
        numer = 2*adj*topo_mem + alpha*att*aux1'; %Numerator
        denom = 2*topo_mem*(topo_mem'*topo_mem) + alpha*topo_mem*(aux1*aux1'); %Denominator
        topo_mem = topo_mem.*(numer./max(denom, realmin));
        %======
        topo_mem = topo_mem.*(0.5 + (adj*topo_mem)./max(2*(topo_mem*topo_mem')*topo_mem, realmin));
        
        %==========
        %Z-Process: update the attribute description matrix Z
        aux2 = topo_mem*trans_TA; %Auxiliary matrix, i.e., P
        numer = att'*aux2; %Numerator
        denom = att_desc*(aux2'*aux2); %Denominator
        att_desc = att_desc.*(numer./max(denom, realmin));
        
        %===========
        %U-Process: update the T-A transition matrix U
        numer = alpha*topo_mem'*att*att_desc + lambd*ones(1, num_topo_clus)'*ones(1, num_att_clus); %Numberator
        denom = alpha*(topo_mem'*topo_mem)*trans_TA*(att_desc'*att_desc) + lambd*ones(1, num_topo_clus)'*ones(1, num_att_clus)*trans_TA; %Denominator
        trans_TA = trans_TA.*(numer./max(denom, realmin));
        
        %=====================
        iter_cnt = iter_cnt + 1;
        %Update the objective function value
        obj = norm(adj - topo_mem*topo_mem','fro')^2;
        obj = obj + alpha*norm(att - topo_mem*trans_TA*att_desc', 'fro')^2;
        obj = obj + lambd*norm(ones(1, num_topo_clus)*trans_TA-ones(1, num_att_clus), 'fro')^2;
        %==========
        error = abs(obj - pre_obj)/pre_obj; % Relative error of the objective function
        %fprintf('Obj. Value %8.4f; Error: %8.8f\n', [obj, error]);
    end
    %fprintf('Obj. Value: %8.4f\n', obj);
    %fprintf('Total Iterations %d\n', iter_cnt);
end

