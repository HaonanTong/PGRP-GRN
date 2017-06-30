%%
% Read Data
[ fig, ngene, expr, plotData, agis, agis_new ] = ...
    f_plotTable( './Process1/Profiles-DEGs-nodes-networks.csv', ...
    './OriginData/NodesForSmallNetworks.txt');

[ network1 ] = f_grn( expr([1 2 4 5 ],1:end)', agis_new([1 2 4 5]) );
[ network2 ] = f_grn( plotData([1 2 4 5 ],1:end)', agis_new([1 2 4 5]) );
[ network3 ] = f_grn( plotData([1 2 4 5 ],2:end)', agis_new([1 2 4 5]) );


%%
[ fig, ngene, expr, plotData, agis, agis_new ] = ...
    f_plotTable( './Process1/Profiles-nodes-networks.csv', ...
    './OriginData/NodesForSmallNetworks.txt');

[ network1 ] = f_grn( expr([1 2 3 4 6 7 ],1:end)', agis_new([1 2 3 4 6 7 ]) );
[ network2 ] = f_grn( plotData([1 2 3 4 6 7 ],1:end)', agis_new([1 2 3 4 6 7 ]) );
[ network3 ] = f_grn( plotData([1 2 3 4 6 7 ],2:end)', agis_new([1 2 3 4 6 7 ]) );
