function result= ICA_rtfmri_two2three(total,x,y,z,level);
voxel=x*y;

%将数据转换为0/1序列
pos1=find(isnan(total));
total(pos1)=0;
pos2=find(total<level);
total(pos2)=0; 
pos3=find(total>=level);
total(pos3)=1;


%将数据还原为x*y*z
for j=1:z
    k=(j-1)*voxel+1;
    off=j*voxel;
    z_slice=total(k:off,:);
    out=reshape(z_slice,x,y);
    result(:,:,j)=out;
end