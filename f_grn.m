function [ network ] = f_grn( Y, name )
%F_GRN Summary of this function goes here
% Y - m by n matrix. m is # of genes and n time point's expression profile  
% name - n_genes by 1 cell, name for each gene
%% Box1
figure();hold on
plot(Y);legend(name);
[T, n_genes] = size(Y);

%% Discretize genes
n_levels = 2; % discretized onto 2 levels
m0 = mean(mean(Y));
genes = Y;
for i = 1 : n_genes
    for j = 1 : T
        if Y(j,i) > m0(1)
            genes(j,i) = 1;
        else
            genes(j,i) = 0;
        end
    end
end
genes = genes';

%% Box2
ESS = 15; % Equivalent sample size (hyperparameters of the Dirichlet
 % distribution)
 BDei_matrix =  -Inf * ones(n_genes,n_genes);
for i = 1 : n_genes%TAR
    for j = 1 : n_genes%REG
        if i ~= j
            REG = genes(j,1:end-1)'; % Subset of regulators (all time points except
             % for the last)
            TAR = genes(i,2:end)'; % Target gene (all time points except for the
             % first)
            %Calculate scores:
            BDei = 0;
            for r = 1:n_levels % # states of each target or regulator
                Nj = REG == r-1;
                 BDei = BDei + log(gamma(ESS/n_levels)/gamma(sum(Nj)+ESS/n_levels));
                for d = 1:n_levels
                     Njk = REG == r-1 & TAR == d-1;
                     BDei = BDei + log(gamma(sum(Njk)+ESS/(n_levels*n_levels))/...
                     gamma(ESS/(n_levels*n_levels)));
                end
            end
            BDei_matrix(j,i) = BDei;
        end
    end
end
%score matrix;

%% for i find hight j from BDei_matrix
network = BDei_matrix == repmat(max(BDei_matrix')',1,n_genes);

 %% Box3
 bg1 = biograph(network,name,'ShowWeights', 'off','EdgeFontSize',12);
for i = 1:size(bg1.edges)
 if get(bg1.edges(i),'Weight') < 0
 set(bg1.edges(i),'LineColor',[1 0 0]);
 end
 set(bg1.edges(i),'LineWidth',2);
end
for i = 1:size(bg1.nodes)
 set(bg1.nodes(i),'FontSize',12)
 set(bg1.nodes(i),'size',[30,30])
 bg1.Nodes(i).Shape = 'circle';
 set(bg1.nodes(i),'Color',[0.71,0.8,1])
 set(bg1.nodes(i),'LineColor',[0,0,0])
end
view(bg1)

end

