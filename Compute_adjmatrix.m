function [ adjmatrix,nodenum ] = Compute_adjmatrix( txtmatrix )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% 从txt中读入网络后，计算邻接矩阵
% txtmatrix: 读取的txt中的节点网络，一般为n*2矩阵，n是边数
% adjmatrix: 返回的邻接矩阵
% nodenum: 返回网络节点的个数

nodemin=min(min(txtmatrix));
if nodemin==0
    txtmatrix=txtmatrix+1;
end
nodenum=max(max(txtmatrix));
b=sparse(txtmatrix(:,1),txtmatrix(:,2),ones(size(txtmatrix,1),1),nodenum,nodenum);
Aij=full(b);
Bij=Aij';
adjmatrix=Aij+Bij;%社区的邻接矩阵
end

