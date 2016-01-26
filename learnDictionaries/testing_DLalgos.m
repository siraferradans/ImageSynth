%Testing which learning algorithm suits better 

load('./Database3000.mat');

%we choose the dataset with highest std (for instance) 
% linearindx = @(lambda)(lambda(1)-1).*L+lambda(2);
% 
% for j2=2:J
%     for l2=1:L
%         lambda2=[j2 l2];
%         
%         covarX=X{linearindx(lambda2)}*X{linearindx(lambda2)}';
%         
%         var(linearindx(lambda2)) = trace(covarX);
%     end 
% end 
% 
% [~,varI]=max(var);

varI = size(X,2);

%contiguous patches are very similar, thus we need to randomize their
%position

XX = X{varI};
IndexP = randperm(size(XX,2));

trainX = double(XX(:,IndexP(1:round(end*4/5))));
testX  = double(XX(:,IndexP(round(end*4/5)+1:end)));

%testing on the size of the dataset
params.modeD=0; %atoms with norm 1
params.mode=1; %l1 norm on the coefs alpha
params.lambda=0.00;
params.K=size(trainX,1); %squared dictionary
params.iter = 1000;
params.batchsize=512;
params.numThreads=-1; % number of threads
params.verbose=false;
%param.modeParam=0;
Energy=@(X,D,alpha)mean(sum((X-D*alpha).^2,1)./sum(X.^2,1)); %
R=@(X,D,alpha)mean( 0.5*sum((X-D*alpha).^2,1)+params.lambda*sum(abs(alpha),1));

%%
rast = 1;
d = size(XX,1);
i = [d^2 d^2*5 d^2*10 d^2*100 min(d^2*1000,size(trainX,2))];%[0.001 0.05 0.1 0.2 0.4 0.6 0.8 1];
t=[];
for rast=1:length(i) %percentage

    %% training data
    data = trainX(:,1:i(rast));
    
    %get dictionary 
 %   [D,alpha] = nmf(data,params);
    
    tic
    D = mexTrainDL(data,params);
    t(rast)=toc
%     %associated coefs
     alpha=mexLasso(data,D,params);
     %and compute the error
    Etraining(rast) = Energy(data,double(D),double(alpha));
    reminder_training(rast) = R(data,double(D),double(alpha));
    
    %% testing data set (does it generalize?)
    alpha2=mexLasso(testX,D,params);   
    Etesting(rast) = Energy(testX,double(D),double(alpha2));
    reminder_testing(rast) = R(testX,double(D),double(alpha2));
  
%    rast=rast+1;
end 
 
subplot(121);plot(i,Etraining); hold on; plot(i,Etesting,'r');
subplot(122);plot(i,reminder_training); hold on; plot(i,reminder_testing,'r');


