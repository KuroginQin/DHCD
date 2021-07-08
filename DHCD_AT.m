function [att_desc,att_mem,trans_AT,obj] = DHCD_AT(adj,att,att_mem,att_desc,trans_AT,beta,lambd,max_iter,min_error)
%Funtion to implement DHCD A-T
%adj: adjacency matrix, i.e., A
%att: network attribute matrix, i.e., C
%att_desc: attribute description matrix, i.e., Z
%att_mem: attribute cluster membership matrix, i.e., X
%trans_AT: A-T transition matrix, i.e., U
%beta, lambd: hyper-parameters of the model
%max_iter: maximum number of iterations
%min_error: minimum relative error to determine convergence
%obj: objective function value

    %====================
    %num_att_clus: number of attribute clusters
    %num_topo_clus: number of topology clusters
    [num_att_clus, num_topo_clus] = size(trans_AT);
    iter_cnt = 0; %Counter of iterations
    error = 0; %Relative error of the objective fucntion
    %====================
    %Calculate the objective function value
    obj = beta*norm(adj - (att_mem*trans_AT)*(trans_AT'*att_mem'), 'fro')^2;
    obj = obj + norm(att - att_mem*att_desc', 'fro')^2;
    obj = obj + lambd*norm(ones(1, num_topo_clus)*trans_AT - ones(1, num_att_clus), 'fro')^2;
    %====================
    while iter_cnt==0 || (error >= min_error && iter_cnt <= max_iter)
        pre_obj = obj; %Objective function value in previous iteraton
        %==========
        %Y-Process: update the attribute cluster membership matrix Y
        aux1 = trans_AT*trans_AT'; %Auxiliary matrix, i.e., S
        aux2 = att*att_desc; % Auxiliary matrix, i.e., CZ
        numer = 2*beta*adj*att_mem*aux1 + aux2; %Numerator
        denom = 2*beta*att_mem*aux1*(att_mem'*att_mem)*aux1 + att_mem*(att_desc'*att_desc); %Denominator
        att_mem = att_mem.*(numer./max(denom, realmin));
        %att_mem = att_mem.*(numer./max(denom, 1e-10)); % 1e-50
        
        %==========
        %Z-Process: update the attribute description matrix Z
        numer = att'*att_mem; %Numerator
        denom = att_desc*(att_mem'*att_mem); %Denominator
        att_desc = att_desc.*(numer./max(denom, realmin));
        %att_desc = att_desc.*(numer./max(denom, 1e-10)); % 1e-50
        
        %==========
        %V-Process: update the A-T transition matrix V
        aux1 = adj*att_mem; %Auxiliary matrix, i.e., AY
        aux2 = att_mem'*att_mem; %Auxiliary matri, i.e., R
        numer = beta*att_mem'*aux1*trans_AT + beta*aux1'*att_mem*trans_AT + lambd*ones(1, num_att_clus)'*ones(1, num_topo_clus); %Numerator
        denom = 2*beta*(aux2*trans_AT)*trans_AT'*(aux2*trans_AT) + lambd*ones(1, num_att_clus)'*ones(1, num_att_clus)*trans_AT; %Denominator
        trans_AT = trans_AT.*(numer./max(denom, realmin));
        %trans_AT = trans_AT.*(numer./max(denom, 1e-10)); % 1e-50
        
        %====================
        iter_cnt = iter_cnt + 1;
        %Update the objective function value
        obj = beta*norm(adj - (att_mem*trans_AT)*(trans_AT'*att_mem'), 'fro')^2;
        obj = obj + norm(att - att_mem*att_desc', 'fro')^2;
        obj = obj + lambd*norm(ones(1, num_topo_clus)*trans_AT - ones(1, num_att_clus), 'fro')^2;
        %==========
        error = abs(obj - pre_obj)/pre_obj; % Relative error of the objective function
        %fprintf('Obj. Value %8.4f; Error: %8.8f\n', [obj, error]);
    end
    %fprintf('Obj. Value: %8.4f\n', obj);
    %fprintf('Total Iterations %d\n', iter_cnt);
end

