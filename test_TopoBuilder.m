clc
disp(' ');
disp('#############################################');
%                ne ke ni ki
tb = TopoBuilder(10, 5,10, 5);
tb = tb.connect_e2i_rand(3);
tb.postidxMat
tb = tb.connect_i2e_rand(3);
tb.postidxMat

%                ne ke ni ki
tb = TopoBuilder(10, 5,10, 5);
tb = tb.connect_i2e_rand([]);
tb.postidxMat

%                ne ke ni ki
tb = TopoBuilder(10, 5,10, 5);
tb = tb.connect_i2i_rand([]);
i2i_entries = tb.postidxMat';
assert(all(i2i_entries(1:50)   == 0),'i2i links are wrong A')
assert(all(i2i_entries(51:100) > 10),'i2i links are wrong B')

%                ne ke ni ki
tb = TopoBuilder(10, 5,10, 5);
tb = tb.connect_e2e_rand([]);
e2e_entries = tb.postidxMat';
assert(all(e2e_entries(1:50)  <= 10),'e2e links are wrong A')
assert(all(e2e_entries(51:100) == 0),'e2e links are wrong B')
%%
%                ne ke ni ki
tb = TopoBuilder(3, 5,3, 5);
tb = tb.connect_e2i_fromTo([1 2 3 1 2 3],[2 3 1 3 1 2])
tb.postidxMat
%%
%                ne ke ni ki
tb = TopoBuilder(3, 5,3, 5);
tb = tb.connect_i2e_fromTo([1 2 3 1 2 3],[2 3 1 3 1 2])
tb.postidxMat
%%
%                ne ke ni ki
tb = TopoBuilder(3, 5,3, 5);
tb = tb.connect_e2e_fromTo([1 2 3 1 2 3],[2 3 1 3 1 2])
tb.postidxMat
%%
%                ne ke ni ki
tb = TopoBuilder(3, 5,3, 5);
tb = tb.connect_i2i_fromTo([1 2 3 1 2 3],[2 3 1 3 1 2])
tb.postidxMat
%% Test neighbor connetions
%                ne ke ni ki
matSizeE = [5 5];
matSizeI = [5 5];

tb = TopoBuilder(prod(matSizeE), 8,prod(matSizeI),8,matSizeE,matSizeI );
fromN = 13;
toN   = 13;
tb = tb.connect_i2i_fromToNeighbor(fromN,toN,1);

tb.plotMat_2D();
%% Test big neighbor connetions
tic
matSizeE = [100 100];
matSizeI = [100 100];

tb = TopoBuilder(prod(matSizeE), 8,prod(matSizeI),8,matSizeE,matSizeI );
nE = prod(matSizeE);
nI = prod(matSizeI);
fromN = randperm(nI,5);
toN   = randperm(nE,5);
tb = tb.connect_i2e_fromToNeighbor(fromN,toN,[6 8]);
tb.plotMat_2D();
toc
