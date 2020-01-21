function pred = knntest(Xtest,class_t,X,class,K,dist_type,pret_type)

% prediction of new samples with calculated model
%
% pred = knnpred(Xtest,X,class,K,dist_type,pret_type)
%
% ------------ INPUT ---------------------------------------------------
% Xtest:        dataset to be predicted [n_test x p] n objects, p variables
% X:            training data matrix (n x p)
% class:        training class vector (n x 1)
% K:            number of neighbors
% dist_type:    'euclidean' Euclidean distance
%               'mahalanobis' Mahalanobis distance
%               'cityblock' City Block metric
%               'minkowski' Minkowski metric
%               'sm' Sokal-Michener 
%               'jt' Jaccard Tanimoto
%               'gle' Gleason-Dice
%               'ct4' Consonni-Todeschini
%               'ac' Austin-Colwell
% pret_type:    'cent' cenering
%               'scal' variance scaling
%               'auto' for autoscaling (centering + variance scaling)
%               'rang' range scaling (0-1)
%               'fp'   fingerprints
%
% ------------ OUTPUT --------------------------------------------------
% pred is a structure conyaining
% class_pred    predicted class vector [n_test x 1]
% neighbors     list of k neighbors for each predicted sample [n_test x k]
% 
% version 1.0 - september 2009
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% www.disat.unimib.it/chm

% version 2.0 - February 2012
% Kamel Mansouri
% Milano Chemometrics and QSAR Research Group
% www.disat.unimib.it/chm

% data check
if length(class)~=size(X,1)
    disp('the class input should be for the training set')
    %class_tr=input('class tr');
    %class=evalin(WS,);
    %keyboard
end

[n,p] = size(Xtest);
[X_scal_train,param] = data_pretreatment(X,pret_type);
X_scal = test_pretreatment(Xtest,param);
Xd = [X_scal;X_scal_train];
% D = pdist(Xd,model.set.dist_type);
% D = squareform(D);
D = knn_calc_dist(X_scal_train,X_scal,dist_type,pret_type);
neighbors = zeros(n,K);
for i=1:n
    D_in = D(i,:);
    [d_tmp,n_tmp] = sort(D_in);
    neighbors(i,:) = n_tmp(1:K);
    d_neighbors = d_tmp(1:K);
    class_calc(i) = knnclass(class(neighbors(i,:)),d_neighbors,max(class),K);
    dc(i,:)=d_neighbors;
end

pred.neighbors  = neighbors;
pred.class_pred = class_calc';
pred.class_param = calc_class_param(pred.class_pred ,class_t);
pred.D=D;
pred.dc=dc;