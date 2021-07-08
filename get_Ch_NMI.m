function [ch_NMI] = get_Ch_NMI(labels, gnd, rate, num_ch_runs)
%Get the channel NMI (Ch-NMI) in a semi-supervised way
%labels: label sequence generated by a certain channel
%gnd: ground-truth sequence of the dataset
%rate: sampling rate
%testNum: number of independent runs
%ch_NMI: the Ch-NMI result

    %===================
    ch_NMI = 0.0;
    cnt = 0;
    [num_nodes, ~] = size(labels); %Get the number of nodes
    num_samp = int8(num_nodes*rate); %Number of sampled nodes
    %===========
    for t=1:num_ch_runs
       [~, indexs] = sort(rand(1, num_nodes));
       labels_sample = zeros(1, num_samp);
       gnd_sample = zeros(1, num_samp);
       for n=1:num_samp
           labels_sample(n) = labels(indexs(n));
           gnd_sample(n) = gnd(indexs(n));
       end
       cur_NMI = compute_NMI(labels_sample, gnd_sample);
       if isnan(cur_NMI)==1
          t = t-1;
          %fprintf('NaN skip\n');
          continue;
       end
       ch_NMI = ch_NMI + cur_NMI;
       cnt = cnt + 1;
    end
    %===========
    ch_NMI = ch_NMI/cnt;
end
