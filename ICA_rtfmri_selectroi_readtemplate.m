function Data2=ICA_rtfmri_selectroi_readtemplate(n)
filename='template.img';
Timepoints=size(filename,1);
V=spm_vol(filename);
Data=spm_read_vols(V);
[a,b]=size(Data);
data=Data;
%%%%%%%%%%%%%%%%%%%%%%%%%%% overlap rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n==1
data(:,57:b)=zeros(a,b-56);data(1:161,58)=zeros(161,1);% 1%
elseif n==2
data(:,59:b)=zeros(a,b-58);data(1:161,58)=zeros(161,1);% 2%
elseif n==3
data(:,60:b)=zeros(a,b-59);data(1:162,59)=zeros(162,1);% 3%
elseif n==4
data(:,60:b)=zeros(a,b-59);% 4%
elseif n==5
 data(:,61:b)=zeros(a,b-60);data(1:142,60)=zeros(142,1);% 5%
elseif n==6
data(:,62:b)=zeros(a,b-61);data(1:149,61)=zeros(149,1);% 6%
elseif n==7
data(:,63:b)=zeros(a,b-62);data(1:158,62)=zeros(158,1);% 7%
elseif n==8
data(:,63:b)=zeros(a,b-62);data(1:134,62)=zeros(134,1);% 8%
elseif n==9
data(:,64:b)=zeros(a,b-63);data(1:145,63)=zeros(145,1);% 9%
elseif n==10
 data(:,65:b)=zeros(a,b-64);data(1:157,64)=zeros(157,1);% 10%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 错误率 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if m==1
% data(163:167,52:56)=ones(5,5);  % 1%
% elseif m==2
%  data(163:169,50:56)=ones(7,7);  % 2%
% elseif m==3
% data(163:169,46:56)=ones(7,11);  % 3%
% elseif m==4
% data(163:172,46:55)=ones(10,10);  % 4%
% elseif m==5
% data(163:171,43:56)=ones(9,14);  % 5%
% elseif m==6
% data(163:177,46:55)=ones(15,10);  % 6%
% elseif m==7
%  data(163:177,47:58)=ones(15,12);  % 7%
% elseif m==8
% data(163:176,42:55)=ones(14,14);  % 8%
% elseif m==9
%  data(163:177,42:56)=ones(15,15);  % 9%
% elseif m==10
% data(163:178,43:58)=ones(16,16);  % 10%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%end

Data1=[];Data2=[];
n1=0;n2=0;

    Data1(1:n1,:,:)=data((a-n1+1):a,:,:);
    Data1((n1+1):a,:,:)=data(1:(a-n1),:,:);

    Data2(:,1:n2,:)=Data1(:,(b-n2+1):b,:);
    Data2(:,(n2+1):b,:)=Data1(:,1:(b-n2),:);


  Data=reshape(Data,a*b,1);
  Data2=reshape(Data2,a*b,1);
  result=ICA_rtfmri_two2three(Data,a,b,1,1);
  result2=ICA_rtfmri_two2three(Data2,a,b,1,1);
  
  
  m=find(result==1);
  m2=find(result2==1);
  mm=intersect(m,m2);%取交集
  fugailv=size(mm,1)/size(m,1);
  
  n=find(result==0);
  nn=intersect(n,m2);
  cuowujihuolv=size(nn,1)/size(m,1);
  
%   figure(1);
%   templatedata1=ICA_rtfmri_show_active(result,a,b,1) ; 
%   montage(templatedata1);
%   figure(2);
%   templatedata2=ICA_rtfmri_show_active(result2,a,b,1) ; 
%   montage(templatedata2)

