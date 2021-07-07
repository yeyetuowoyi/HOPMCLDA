
lncSim = load ('.\Dataset\lncRNAsimilarity.txt');
load('.\Dataset\interaction.mat');
disSim = load('.\Dataset\diseasesimilarity.txt');



%%parameters
maxiter = 500;
beta=10;
alpha=0.5;
omega=0.9;
tol = 0.001;
k1=88;
k2=58;
order = 2;
weights = 0.1;


%%



Aori=interaction;

DS= HOP(disSim,k1,order,weights);
LS= HOP(lncSim,k2,order,weights);

Wll = LS;
Wdd = DS;
Wdl = Aori;
Wld = Aori';
[dn,dr] = size(Wdl);
T = [Wll, Wld; Wdl, Wdd];


[t1, t2] = size(T);
trIndex = double(T ~= 0);
[WW,iter] = MC(omega,alpha, beta, T, trIndex, tol, maxiter, 0,1);


M_recovery = WW((t1-dn+1) : t1, 1 : dr);
%% Loocv
index_1 = find(1 == Aori);
pp = length(index_1);
for t = 1 : pp
        DL=Aori;
        DL(index_1(t)) = 0;

        Wll = LS;
        Wdd =  DS;   
        Wdl = DL;
        Wld = DL';

        [dn,dr] = size(Wdl);
        T = [Wll, Wld; Wdl, Wdd];

        
        [t1, t2] = size(T);
        trIndex = double(T ~= 0);

       [WW,iter] = MC(omega,alpha, beta, T, trIndex, tol, maxiter, 0,1);

        result = WW((t1-dn+1) : t1, 1 : dr);
        M_recovery(index_1(t)) = result(index_1(t));        
        t
        iter
        
end       
        pre_label_score = M_recovery(:);
        label_y = Aori(:);
        auc = roc(pre_label_score,label_y);

