clear all;
clc;
close all;

%************************(1)***************************
%read the fMRI data
spm_defaults;
ans='D:\renzhi\TSCICA\simulated data';%data path
cd(ans); 
img_file=dir('*.img');
img_file_cell=struct2cell(img_file);
filename=img_file_cell(1,:)';
Timepoints=size(filename,1);
V=spm_vol(filename);
for ii=1:135
data(:,:,ii)=spm_read_vols(V{ii});
end
cd('D:\renzhi\TSCICA');%The path of the source code
level=2;%Threshold

%data stores the fMRI data
%***********************************************************************
a=size(data,1);
b=size(data,1);
t=size(data,3);
Data=zeros(a*b,t);
for i=1:t
    Data(:,i)=reshape(data(:,:,i),a*b,1);
end
voxel=a*b;
result=zeros(a,b); 
total=zeros(voxel,1); 

% Centering
  [q,n]=size(Data);
  meanx = mean(Data);
  e = ones(q,1);
  Data = Data-(e*meanx);
  
% Covariance matrix
  CC = cov(Data);
% Whitening by the classical PCA
  [EE,DD]= eig(CC);
  ID=inv(sqrtm(DD));
  
  x=(EE*ID*EE')*Data';
  Data=x;
  
w0 = rand(size(Data,1),1);w0=w0/norm(w0);
mu1 = 1;mu2=1;
lambda = 1;
OverValue=0.001;  maxIter = 200;%the number of iteration steps
templatedata00=ICA_rtfmri_selectroi_readtemplate(5); %load the spatial reference 
ref= templatedata00';  %spatial reference
    RT=2;
    p=[8 10 2 2 6 0 32];
    hrf=spm_hrf(RT,p);
    block_design1=zeros(135,1);
    block_design1(16:30,1)=ones(15,1);
    block_design1(46:60,1)=ones(15,1);
    block_design1(76:90,1)=ones(15,1);
    block_design1(106:120,1)=ones(15,1);
    ref_t=conv(block_design1,hrf);
    ref_t=ref_t(1:size(block_design1,1),1);
    reft=ref_t-ones(size(ref_t))*mean(ref_t);
    reft=reft/norm(reft);%temporal reference
threshold = 0.9; 
[IC, w] = TSCICA(Data, ref, reft,threshold, w0, mu1,mu2, lambda, maxIter, OverValue);

 %show the spatial template
        result0=ICA_rtfmri_two2three(templatedata00,a,b,1,1);
        figure;
        templatedata0=ICA_rtfmri_show_active(result0,a,b,1) ; 
        montage(templatedata0);
        
%********************the first **************************
        result1=ICA_rtfmri_two2three(IC,a,b,1,level);
        figure;
        IC1=ICA_rtfmri_show_active(result1,a,b,1) ; 
        montage(IC1);

