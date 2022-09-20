function result= ICA_rtfmri_two2three(total,x,y,z,level);
voxel=x*y;

%������ת��Ϊ0/1����
pos1=find(isnan(total));
total(pos1)=0;
pos2=find(total<level);
total(pos2)=0; 
pos3=find(total>=level);
total(pos3)=1;


%�����ݻ�ԭΪx*y*z
for j=1:z
    k=(j-1)*voxel+1;
    off=j*voxel;
    z_slice=total(k:off,:);
    out=reshape(z_slice,x,y);
    result(:,:,j)=out;
end