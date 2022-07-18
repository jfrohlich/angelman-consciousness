function[P_val] = permclusttest(A,B)

% PERMUTATION TEST

% Do permutation test for each point in channel-frequency or channel-scale space 

assert(sum(size(A)==size(B)) == ndims(A),'Sizes of A and B do not match')
C = cat(3,A,B); % concatenate both groups/conditions
n = size(A,3); % number of subjects in each group/condition
N = 10^4; % number of permutations
p_thresh = 0.01;

load 1020_neighbors.mat

% unshuffled
[h,p,ci,stat] = ttest2(C(:,:,n+1:end), C(:,:,1:n),'dim',3,'vartype','unequal','tail','both');
% separate positive and negative differences so they don't cluster together
pos = (stat.tstat>0 & p<p_thresh);
[lbl_pos,n_pos,idx_n_pos] = ro_cluster_spacefreq(pos',neighbours2);
neg = (stat.tstat<0 & p<p_thresh);
[lbl_neg,n_neg,idx_n_neg] = ro_cluster_spacefreq(neg',neighbours2);
% pool together clusters for both tails
n2tail = [n_pos,n_neg]; 
L = length(n2tail);
empSize = zeros(N,1);

rng('default')
for k = 1:N   
    if mod(k,100) == 0
        fprintf('%1.0f percent done\n',k/N*100)
    end
    rng(k+342398) % keep the same seed each time we run the script
    DX = randperm(n*2); % permute conditions
    [h,p,ci,stat] = ttest2(C(:,:,DX(n+1:end)),C(:,:,DX(1:n)),'dim',3);
    pos_perm = (stat.tstat>0 & p<p_thresh); % positive t-statistics 
    [~,n_pos_perm,~] = ro_cluster_spacefreq(pos',neighbours2);
    neg_perm = (stat.tstat<0 & p<p_thresh); % negative t-statistics 
    [~,n_neg_perm,~] = ro_cluster_spacefreq(neg',neighbours2);
    % pool together clusters for both tails
    n_perm = [n_pos_perm,n_neg_perm]; 
    biggest = max(n_perm); % largest cluster size 
    empSize(k) = biggest; % 
end

X = repmat(n2tail,[N,1]);
Y = repmat( empSize, [1,L]);
howMany = sum( Y > X ); % how many permuted clusters are bigger than the 'real' clusters?
P_val = howMany./N; % convert to p-vals 

keyboard

end


%% Older code (archived) 
% 
% %% Do permutation test for each point in channel-frequency or channel-scale space 
% 
% assert(sum(size(A)==size(B)) == ndims(A),'Sizes of A and B do not match')
% C = cat(3,A,B); % concatenate both groups/conditions
% n = size(A,3); % number of subjects in each group/condition
% N = 10^4; % number of permutations 
% 
% D = nan(size(A,1),size(A,2),N); % mean difference allocation
% T = nan(size(A,1),size(A,2),N); % t-stats allocation
% P = nan(size(A,1),size(A,2),N); % p-values allocation
% 
% % unshuffled
% [h,p,ci,stat] = ttest2(C(:,:,n+1:end), C(:,:,1:n),'dim',3);
% true_D = mean(ci,3);
% true_T = stat.tstat; % t-statistic 
% true_P = p; % true p-value 
% rng('default')
% for k = 1:N   
%     if mod(k,100) == 0
%         fprintf('     %1.0f percent done permuting labels\n',k/N*100)
%     end
%     rng(k+342398)
%     DX = randperm(size(C,3)); % permute labels
%     [h,p,ci,stat] = ttest2(C(:,:,DX(n+1:end)), C(:,:,DX(1:n)),'dim',3);
%     D(:,:,k) = mean(ci,3);
%     T(:,:,k) = stat.tstat;
%     P(:,:,k) = p;
% end
% 
% % number of comparsions greater than true value
% howMany = sum(abs(repmat(true_T,[1,1,N]))>repmat((max(max(T),[],2)),[size(A,1),size(A,2),1]),3); 
% % divide by number of permutations to get the p-value
% pemp = howMany./N; 
% 
% % multiple comarison across space and freq
% pemp = (sum(abs(repmat(true_T,[1,1,N]))<repmat((max(max(T),[],2)),[size(A,1),size(A,2),1]),3)/N);
% temp=pemp<0.05;
% 
% if plotting
%     figure, imagesc(-log10(pemp),[-4,4]), colorbar, colormap jet
%     figure, imagesc(temp,[0,1])
% end


                

