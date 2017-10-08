addpath( './external-lib/bnt/graph' )
addpath( './external-lib/bnt/KPMtools')
addpath( '../')

filename = 'graphlab_test.mat';

% synthetic_data
% 
% %% This is a bit ugly, but you can easily plug your data in
% 
% 
% % initial parameters for the algorithm
% % Width and height of the input matrix, which is to be decomposed 
% [W,K] = size(X);
% I = 5; % number of factors
% %
% % Shapes and scales of gammas
% a_ve = ones(I,nr_params_to_initialize)*2;
% b_ve = 10./a_ve;
% %
% a_utility = repmat(1000,[I,1]);
% %
% a_tm = ones(W,I);
% b_tm = ones(W,I); % rightfactor will be mean of exponential distros
% 
% % Dirichlet - related parameters. Wherever we have number of parents i > 1,
% % we'll need to put a vector of size i inside the corresponding cell 
% % Note: other Dirichlet parameters, if defined, will be ignored.
% u_dirichlet = cell(K,1);
% relevant_indices = 1:K;
% relevant_indices = relevant_indices(numparents>1);
% for(i = relevant_indices)
%     u_dirichlet{i} = ones(I,numparents(i)).*0.5;
% end
% %

% Populate adjacency matrix with utility gamma nodes, initialize moments
% and stuff - all of these need to be provided to graphlab, because
% currently that implementation does not know how to generate initial data
x=X;
% generate_initial_data
load('graphlab_test_matlab.mat');
%
u_dirichlet = [u_dirichlet; cell(size(gamma_chain_adjacency,1)-K,1)]; % append to reflect the added aux gammas
%
order = topological_sort(gamma_chain_adjacency);
[ditch, backward_mapping] = sort(order);
gamma_chain_adjacency = gamma_chain_adjacency(order, order);
% to revert topological ordering, use gamma_chain_adjacency(backward_mapping, backward_mapping);

expectations_markov = expectations_markov(order,:,:);
u_dirichlet = u_dirichlet(order);
numparents = full(sum(gamma_chain_adjacency,1));
numchildren = full(sum(gamma_chain_adjacency,2));
%
% Translating some of the variables to what ghaphlab implementation accepts
[rows_adjacency_colwise, cols_adjacency_colwise] = find(gamma_chain_adjacency);
rows_adjacency_colwise = rows_adjacency_colwise-1;  % zero-based indexing
cols_adjacency_colwise = cols_adjacency_colwise-1;  % zero-based indexing
% Also reordering the input matrix and mask kM
orderX = order(order<=size(X,2));
[ditch backward_mapping] = sort(orderX);
X = X(:,orderX);
M = M(:,orderX);

%
colwise_num = sum(M);
[colwise_mask_rows, colwise_mask_cols] = find(M);
colwise_mask_rows = colwise_mask_rows - 1;
colwise_mask_cols = colwise_mask_cols - 1;
x_raw = X(M);

% Parameter initialization
%
expectations_discrete = cell(network_size,1);
expectations_dirichlet = cell(network_size,1);
u_dirichlet = u_dirichlet(~cellfun('isempty',u_dirichlet));
indeces = find(numparents>1);
cnt = 1;
for(i = indeces)
    parent_indices = find( gamma_chain_adjacency(:,i) );
    expectations_discrete{i} = ones(numparents(i),I)./numparents(i);
    expectations_dirichlet{i} = psi(u_dirichlet{cnt}) - ...
        ones(size(u_dirichlet{cnt}))*diag( psi(sum(u_dirichlet{cnt},1)) );
    cnt = cnt + 1;
end
clear initial;

% [u_dirichlet(~cellfun('isempty',u_dirichlet)), expectations_dirichlet(~cellfun('isempty',expectations_dirichlet)), expectations_discrete(~cellfun('isempty',expectations_discrete))]

markov_statistic_1 = expectations_markov(:,:,1)';
markov_statistic_2 = expectations_markov(:,:,2)';
E_t_transposed = E_t';
L_t_transposed = L_t';
a_tm_transposed = a_tm';
b_tm_transposed = 1./b_tm'; % Note that Graphlab will accept inverses. I guess I should use shapes and scales, not this legacy thing...
numparents = uint32(numparents);
numchildren = uint32(numchildren);
order = uint32(order);
a_ve_auxillary = a_utility;
rows_adjacency_colwise = uint32(rows_adjacency_colwise);
cols_adjacency_colwise = uint32(cols_adjacency_colwise);
colwise_num = uint32(colwise_num);
colwise_mask_cols = uint32(colwise_mask_cols);
colwise_mask_rows = uint32(colwise_mask_rows);
% All of these need to be provided to graphla
a_ve = a_ve(order,:)';
b_ve = 1./b_ve(order,:)'; % Again, I should probably use shapes and scales, not inverse means... Let's have this as a todo.
b_ve = b_ve(:,numparents == 0);
%
order = order-1; % zero-based indexing
%
%
%
save(filename...
    , 'W', 'K', 'I'...
    , 'rows_adjacency_colwise', 'cols_adjacency_colwise'...
    , 'numparents', 'numchildren', 'order'...
    , 'a_ve', 'b_ve', 'a_ve_auxillary'...
    , 'a_tm_transposed', 'b_tm_transposed'...
    , 'u_dirichlet'...
    ... % Initialization data - those are mandatory at the moment
    , 'expectations_discrete', 'expectations_dirichlet'...
    , 'markov_statistic_1', 'markov_statistic_2'...
    , 'E_t_transposed', 'L_t_transposed'...
    , 'colwise_num', 'colwise_mask_rows', 'colwise_mask_cols', 'x_raw'...
    , 'backward_mapping'...
);