function [Fsc, Jac] = over_eva(gnd_path, res_path)
%Function to compute the F-score and Jaccard of overlapping community detection result
%gnd_path: path of the ground-truth
%res_path: path of the community detection result

for alli = 0:10
    gndfile = textread(gnd_path,'%s','delimiter','\n','whitespace','','bufsize',80000);
    gmlfile = textread(res_path,'%s','delimiter','\n','whitespace','','bufsize',80000);
    
    [grow, gcol] = size(gmlfile);
    [gndrow, gndcol] = size(gndfile);
    gndnum = 0;
    for i=1: gndrow
        str = deblank(gndfile{i,1});
        S1 = regexp(str, '\s+', 'split');
        S2 = str2num(char(S1));
        S3 = sort(S2);
        gndnum = gndnum + length(S3);
        gmlnum = 0;
        for j = 1:grow
            strgnd = deblank(gmlfile{j,1});
            S1gnd = regexp(strgnd, '\s+', 'split');
            S2gnd = str2num(char(S1gnd));
            gmlnum = gmlnum + length(S2gnd);
            score(i,j) = Fscore(S3,S2gnd);
            jindex(i,j) = Jaccard(S3,S2gnd);
        end
    end
    Fmeasure = max(score');
    Fmeasure1 = max(score);
    Jindex = max(jindex');
    Jindex1 = max(jindex);
    
    first_f = sum(Fmeasure)/(2*gndrow);
    second_f = sum(Fmeasure1)/(2*grow);
    first_j = sum(Jindex)/(2*gndrow);
    second_j = sum(Jindex1)/(2*grow);
    
    Fmeasures(alli+1) = first_f + second_f;
    J(alli+1) = first_j + second_j;
    Fsc = max(Fmeasures);
    Jac = max(J);
end
i = 1;

function Fscore_harmonic = Fscore(detected, gnd)

% temp = gnd;
% gnd = detected;
% detected = temp;

inter = intersect(detected, gnd);
precision = length(inter)/length(detected);
recall = length(inter)/length(gnd);
Fscore_harmonic = harmmean([precision recall]);

function Jindex = Jaccard(detected, gnd)
inter = intersect(detected, gnd);
all = union(detected, gnd);
Jindex = length(inter)/length(all);

function Pindex = Purity(detected, gnd)
inter = intersect(detected, gnd);
Pindex = length(inter);
