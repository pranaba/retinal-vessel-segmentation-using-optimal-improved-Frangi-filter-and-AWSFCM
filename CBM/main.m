% Reference: "A novel framework for retinal vessel segmentation using optimal improved
% Frangi filter and adaptive weighted spatial FCM" To be published in
% Computers in biology and medicine, Elsevier.
% Authors:S. Mahapatra, S. Agrawal, P.K. Mishro, R.B. Pachori
% Please cite the paper on using
% 
% In this paper, we suggest a novel framework for retinal vessel segmentation using optimal improved
% Frangi filter and adaptive weighted spatial FCM. The optimal improved Frangi-based multi-scale filter 
% is developed for vessel enhancement. The parameters of the Frangi filter are optimized using a modified 
% enhanced leaderparticle swarm optimization (MELPSO). The enhanced image is segmented using a novel adaptive 
% weighted spatial fuzzyc-means (AWSFCM) clustering technique.


clc; clear all; close all;


%%%%%%%%%%%%% READ IMAGE %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% HRF %%%%%%%%%%%%%%%

 I = imread('13_g.jpg'); % Original Color Image
 ref_im=imread('13_gt.tif'); % Groundtruth
 bw_mask =imread('13_g_mask.tif'); % Mask
 bw_mask=bw_mask(:,:,2);
%%%%%%%%%%%%%% Display Image
 
%  figure;imshow(I);
%  figure;imshow(ref_im);
%  figure;imshow(bw_mask);

bw_mask=logical(bw_mask);
ref_bw=im2bw(ref_im,0.5);

I=(I(:,:,2)); %Extract Green Channel
[H L]=size(I); % Find size of Image
I=imresize(I,[H/4 L/4]); % reduce image size For fast computation
ref_bw=imresize(ref_bw,size(I));
bw_mask=imresize(bw_mask,size(I));


 %%%%%%%%%%%%% MELPSO START %%%%%%%%%%%
  nVar=2;      %number of  variables
  VarSize=[1,nVar]; 
     
% % % %% Parameters of MELPSO
  MaxIt=10;             % % %maximum number of iterations
  nPop=30;              % % %population size (swarm size)
  w=2;                  % % %Inertia coefficient
  wdamp=0.5;            % % %damping factor
  c1=1.5;               % % personal acceleration coefficient
  c2=1.5;               % % social acceleration coefficient

% % % Initialization
% % %creating a template for the particles
 
empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
% % % % % % % %creating a population array
 particle=repmat(empty_particle,nPop,1);
% % % % % % % %creating a global best for refernce
GlobalBest.Cost=0;
GlobalBest.Position=[];
c1L=1;c1H=3;c2L=0.4;c2H=1.0;
constraints=[1 3;0.3 1.0];	
a=0.2;
b=1.5;

  Gsigma=2;
  SigCauchy=2;
  c4=0.2;%%%% new mutation constat etta

for i=1:nPop
     %Generate random solution
     
      particle(i).Position(1)=round(c1L+(c1H-c1L)*rand);
      particle(i).Position(2)=c2L+(c2H-c2L)*rand;   
      sigmas=a+(b-a)*(rand(1,5));
% %     %initialize velocity
      particle(i).Velocity=zeros(VarSize);
% %     
% %     %Evaluation
       particle(i).Cost=evaluateFitnessFRANGI(particle(i).Position,I,sigmas);
    
% %     %update the personal cost
     particle(i).Best.Position=particle(i).Position;
     particle(i).Best.Cost=particle(i).Cost;
%     %update the global best
    if particle(i).Best.Cost>GlobalBest.Cost
        GlobalBest.Position=particle(i).Best.Position;
        GlobalBest.Cost=particle(i).Best.Cost;
       
    end
    
end

% % % % array to hold best costs at each step
 BestCosts=zeros(MaxIt,1);
% % %% Main Loop of PSO
  for it=1:MaxIt
      for i=1:nPop
 %         %update velocity
          particle(i).Velocity=w*particle(i).Velocity+c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position)+c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
% %         %update position
          particle(i).Position=particle(i).Position+particle(i).Velocity;
% %         %%%%% Check within constrains
 
          particle(i).Position(1)=round(particle(i).Position(1));
         if particle(i).Position(1)<c1L     
            particle(i).Position(1)=c1L;
         else if particle(i).Position(1)>c1H
                particle(i).Position(1)=c1H; 
             end
         end
             
             if particle(i).Position(2)<c2L
            particle(i).Position(2)=c2L;
         else if particle(i).Position(2)>c2H
                particle(i).Position(2)=c2H; 
             end
              end 
           
 		end % for k
% %         
% %         %evaluation
          particle(i).Cost=evaluateFitnessFRANGI(particle(i).Position,I,sigmas);
% %         %update personal best
        if particle(i).Cost>particle(i).Best.Cost
            particle(i).Best.Position=particle(i).Position; %Update particle best
            particle(i).Best.Cost=particle(i).Cost;
%             %update global best
% 
            if particle(i).Best.Cost>GlobalBest.Cost
                GlobalBest.Position=particle(i).Best.Position; %Update Global Best
        GlobalBest.Cost=particle(i).Best.Cost;
                 
            end
        end
  %%%%%%%%%%%%%% APPLY MUTATION %%%%%%
    
    
    
%     %%%%%%%%%%%%%%%% 1. GAUSSIAN MUTATION %%%%%%%%%%%
%     
      [GlobalBest.Position GlobalBest.Cost]=gausianMutation(I,constraints,GlobalBest.Position,GlobalBest.Cost,Gsigma,sigmas);
   Gsigma=Gsigma-(1/MaxIt);
     
% %     %%%%%%%%%%% 2.CAUCHY Mutation %%%%%%%%%%%%%%%%
   
    [GlobalBest.Position GlobalBest.Cost ]=CauchyMutation(I,constraints,GlobalBest.Position,GlobalBest.Cost,SigCauchy,sigmas);
      SigCauchy=SigCauchy-(1/MaxIt);

% %      %%%%%%%%%%% DE based Mutation %%%%%%%%%%%%%%%%
%    
           R=randi(nPop,1);S=randi(nPop,1);
        
      [GlobalBest.Position GlobalBest.Cost]=DEMutation(I,constraints,GlobalBest.Position,GlobalBest.Cost,particle(R).Position,particle(S).Position,sigmas);
%     
      [GlobalBest.Position GlobalBest.Cost]=  mutation_new(I,constraints,GlobalBest.Position,GlobalBest.Cost,sigmas,it,MaxIt,c4);  

c4=c4-(1/MaxIt);
% 
     BestCosts(it)=GlobalBest.Cost;
    disp(['Iterations' num2str(it) 'BestCost--' num2str(BestCosts(it))])
    w=2-((2/MaxIt)*it);
 
  end

% figure;
% plot(BestCosts,'LineWidth',2)
xlabel('Iterations');
ylabel('Best Cost');
grid on;

% % % % % % %    final output Image  % % % % % % %
    II=finaloutputImageFRANGI(GlobalBest.Position,I,sigmas); %%% compute using GlobalBest (Leader) value
  figure;imshow((II))

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION %%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%% AWSFCM %%%%%%%%%%%%%%%

  I3=(II);
  imn=mat2gray((I3));
 data=imn(:);
 [CENTER, U, OBJ_FCN] = awsfcm(data);%%%%%%%%%%%%% 

  [A I]= sort (CENTER);
 img=I3;
member=U';
%%%%%%%%%%%
DA1=U(I(1),:);
 DA1=reshape(DA1,[size(img)]).*255;
  DA1=uint8(DA1);

%%%%%%%%
DA3=U(I(3),:);

DA3=reshape(DA3,[size(img)]).*255;
DA3=uint8(DA3);
 figure;imshow(DA3);
DA=(DA3);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outImg=DA;
 Img=outImg;
 figure;imshow(DA);


% %%%%%%%%%%% Post Processing %%%%%%%%%%%%%%%%%

 outImg=im2bw(outImg);
    bw_result=bwareaopen(outImg,30);
  bw_result=bw_result.*bw_mask;
figure;imshow(bw_result);

%%%%%% Compute Performance index %%%%%%%%%%
 r1=evaluate_metrics(bw_result,ref_bw,bw_mask)

 
 
 
 function z= evaluateFitnessFRANGI(x,I,sigmas)

I=uint8(I);
 IA=I/(max(max(I)));

c=x(1);
d=x(2);


V1 = fun_2D_vessel(I, sigmas, [c;c], d, false);


   z=entropy(V1);
 end


 function [GlobalBest,BestCost] = gausianMutation(I,constraints,GlobalBest,BestCost,sigma,sigmas)


pd=makedist('normal','mu',0,'sigma',2);
x=(-20:1:20);
y=pdf(pd,x);
    
for i=1:size(constraints, 1)
    c1L=constraints(i,1);
    c1H=constraints(i,2);
    GaussiGlobalBest(i)=round(c1L+(c1H-c1L)*y(21));
   
    end
   
   
    bestcost2=evaluateFitnessFRANGI(GaussiGlobalBest,I,sigmas);
    if bestcost2 > BestCost
       GlobalBest =GaussiGlobalBest;
     
        BestCost=bestcost2;
    end

 end

 function [GlobalBest,BestCost] =CauchyMutation(I,constraints,GlobalBest,BestCost,SigCauchy,sigmas)



pd=makedist('tLocationScale','mu',0,'sigma',2,'nu',1);
x=(-20:1:20);
y=pdf(pd,x);
    
for i=1:size(constraints, 1)
    c1L=constraints(i,1);
    c1H=constraints(i,2);
    cauchyGlobalBest(i)=round(c1L+(c1H-c1L)*y(21));
   
    end
   
   
    bestcost2=evaluateFitnessFRANGI(cauchyGlobalBest,I,sigmas);
    if bestcost2 > BestCost
       GlobalBest =cauchyGlobalBest;
     
        BestCost=bestcost2;
    end
 end

 function [G C] = DEMutation(I,constraint,Gp,GpCost,pop1,pop2,sigmas)

F=0.8;
       
       
    MuP=Gp+F*(pop1-pop2);


[D r]=size(constraint);
for d=1:D
    
     if   MuP(d)<constraint(d,1)
         MuP(d)=constraint(d,1);
     else if MuP(d)>constraint(d,2)
           MuP(d)=constraint(d,2);  
         end
     end
end


MuP(1)=round(MuP(1));

MuPCost=evaluateFitnessFRANGI(MuP,I,sigmas);
if MuPCost>GpCost
    G=MuP;
    C=MuPCost
else 
    G=Gp;
    C=GpCost;
end


 end

 function [G C] = mutation_new(I,constraint,Gp,GpCost,sigmas,t,T,c4)
%  c4=0.2;
[D r]=size(constraint);
rand1=rand(1,1);
rand2=rand(1,1);
for d=1:D
MuP(d)=Gp(d)+c4*(rand1-rand2)*(1-(t/T));
end


MuP(1)=round(MuP(1));
for d=1:D
   
     if   MuP(d)<constraint(d,1)
         MuP(d)=constraint(d,1);
     else if MuP(d)>constraint(d,2)
           MuP(d)=constraint(d,2);  
         end
     end
end
MuPCost=evaluateFitnessFRANGI(MuP,I,sigmas);
if MuPCost>GpCost
    G=MuP;
    C=MuPCost
else 
    G=Gp;
    C=GpCost;
end

 end

 function [V1] = finaloutputImageFRANGI(x,I,sigmas)

c=x(1);d=x(2);

V1 = fun_2D_vessel(I,sigmas, [c;c], d, false);

% imwrite(g,'ELPSOEnhance7.jpg','jpg');
 end
 
 
 
function [CENTER, U, OBJ_FCN] = awsfcm(img)

imn=mat2gray(img);
m=2; % degree of fuzziness
n=3; % no of cluster centers
max_iter=100;
min_impro=0.0001;
obj_fcn=zeros(max_iter,1);


[CENTER, U, OBJ_FCN] = fcm(img,n);
uij=U;
data=imn(:);
img1=wiener2(img,5);
[rn,cn]=size(img);
imgsiz=rn*cn;
imgv=reshape(img1,imgsiz,1);
imgv=double(imgv);


for ii=1:max_iter

mf = U.^m; % MF matrix after exponential modification
center = mf*data./((ones(size(data, 2),1)*sum(mf'))'); % new center
dist = distfcm(center, data); % fill the distance matrix
obj_fcn(ii) = sum(sum((dist.^2).*mf)); % objective function
tmp = dist.^(-2/(m-1)); % calculate new U, suppose expo != 1
U_new = tmp./(ones(n, 1)*sum(tmp));
nwin=5;
dims=[rn cn]; % dimention of image
tempwin=ones(nwin);
mfwin=zeros(size(U_new));
for i=1:size(U_new,1)
tempmf=reshape(U_new(i,:), dims);
tempmf=imfilter(tempmf,tempwin,'conv');
mfwin(i,:)=reshape(tempmf,1,size(U_new,2));
end
spw=2; % p
mfw=2; % q
mfwin=mfwin.^spw;
U_new=U_new.^mfw;
tmp=mfwin.*U_new;
U_new=tmp./(ones(n, 1)*sum(tmp));

if ii > 1,
		if abs(obj_fcn(ii) - obj_fcn(ii-1)) < min_impro, break; end,
	end
end
uij=U_new;

U=uij;
obj_fcn;
center;
end
 
 function r=evaluate_metrics(bw,ref_bw,bw_mask)
% INPUTS
%     bw: segmented image-(logical)
%     ref_bw: ground-truth-(logical)
%     bw_mask: FOV mask-(logical)
% OUTPUTS
%     r=[TPR;FPR;accuracy;precision];

TP_image=ref_bw&bw;
TP=sum(TP_image(:))% # of hits (True Positive)
FN_image=ref_bw&~bw;
FN=sum(FN_image(:))% # of misses (False Negative\Type 2 Error)

FP_image=~ref_bw&bw;
FP=sum(FP_image(:))% # of false alarms (False Positive/Type 1 Error)
TN_image=~ref_bw&~bw;
TN=sum(TN_image(:))% # of correct rejections (True Negative)

accuracy=(TP+TN)/(TP+FN+FP+TN);
TPR=TP/(TP+FN);% True Positive Rate (sensitivity/recall/hit rate)
FPR=FP/(FP+TN);% False Positive Rate (specificity=1-FPR)
PPV=TP/(TP+FP);%positive predictive value (precision)
specificity=1-FPR;
% NU=(TP*TN)-(FP*TN);
% DE=(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN);
% MCC=NU/(sqrt(DE));
% AUC=0.5*((TN/(TN+FP))+TPR);
%  bw=(bw.*bw_mask);
 ref_bw=double(ref_bw);
j=jaccard(bw,ref_bw);
D=dice(bw,ref_bw);
accuracy
sensitivity=TPR
specificity
jaccard_coef=j
Dice_coef=D
 r=[accuracy;TPR;specificity;j;D;PPV];
%r1=[TPR;accuracy;specificity];
 end
