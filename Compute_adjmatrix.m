function [ adjmatrix,nodenum ] = Compute_adjmatrix( txtmatrix )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% ��txt�ж�������󣬼����ڽӾ���
% txtmatrix: ��ȡ��txt�еĽڵ����磬һ��Ϊn*2����n�Ǳ���
% adjmatrix: ���ص��ڽӾ���
% nodenum: ��������ڵ�ĸ���

nodemin=min(min(txtmatrix));
if nodemin==0
    txtmatrix=txtmatrix+1;
end
nodenum=max(max(txtmatrix));
b=sparse(txtmatrix(:,1),txtmatrix(:,2),ones(size(txtmatrix,1),1),nodenum,nodenum);
Aij=full(b);
Bij=Aij';
adjmatrix=Aij+Bij;%�������ڽӾ���
end

