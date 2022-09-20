function template= ICA_rtfmri_show_structimg(filename,x,y,z)
%filename is the structure image's name
%y is the image's row number
%x is the image's collum number
%first is the beginning slice
%last is the endding slice

fid=fopen(filename,'r');
n=z;
%template=zeros(x,y,1,n);
template=zeros(y,x,3,n);
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = ICA_rtfmri_spm_hread(filename);
%dtype=spm_type(TYPE);
dtype='int16';
for i=1:z
  voxel=y*x;
  switch (dtype)
      case 'uint8'
          off=(i-1)*voxel;
      case 'int16'
          off=(i-1)*voxel*2;
      case 'int32'
          off=(i-1)*voxel*4;
      case 'float'
          off=(i-1)*voxel*4;
  end
 fseek(fid,off,-1);
 out=fread(fid,[x,y],dtype); 
%  off=(i-1)*voxel*2;
% % off=(i-1)*voxel*1;
%   fseek(fid,off,-1);
%  out=fread(fid,[x,y],'int16');
%   out=fread(fid,[x,y],'uint8');
  %out=invert(out',row);
  %out=out';
  out=flipud(out');
  out=mat2gray(out);
  
  filename0='one_ROI.img';
  Timepoints=size(filename0,1);
  V=spm_vol(filename);
  Data=spm_read_vols(V);
  [a,b]=size(Data);
  data=zeros(a,b);
for ii=1:a
    for jj=1:b
        if Data(ii,jj)<=100
          data(ii,jj)=0;
        elseif Data(ii,jj)<820 && Data(ii,jj)>100
          data(ii,jj)=0.7;
        else
          data(ii,jj)=1;
        end 
    end
end
DATA=zeros(a,b);
  for ij=1:a
      DATA(:,ij)=data(ij,:)';
  end
  out=DATA;
  
  [X,map]=gray2ind(out);
  RGB=ind2rgb(X,map);
  template(1:y,1:x,1:3,i)=RGB;
end
fclose(fid);

