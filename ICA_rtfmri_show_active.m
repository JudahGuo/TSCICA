function template=ICA_rtfmri_show_ICs(total,x,y,z);
filename='one_ROI.img';
%filename='srf20071107_01_YL_LiYongnian-0002-00129-000129.img';
%filename='wsingle_subj_T1.img';%smp内的模板
%filename='fuck.img';%smp内的模板
%filename='template20.img';%smp内的模板
template=ICA_rtfmri_show_structimg(filename,x,y,z);
%经ws画激活图

 N=z;
voxel=y*x;
for slice=1:N
    out=total(:,:,slice);
    out=flipud(out');
    for i=1:y
        for j=1:x
            if (out(i,j)==1)
                template(i,j,1,slice)=1;
                template(i,j,2,slice)=0;
                template(i,j,3,slice)=0;
            end  
            if(out(i,j)==2)
                template(i,j,1,slice)=0;
                template(i,j,2,slice)=1;
                template(i,j,3,slice)=0;
            end
            if(out(i,j)==3)
                template(i,j,1,slice)=0;
                template(i,j,2,slice)=0;
                template(i,j,3,slice)=1;
            end
            if(out(i,j)==4)
                template(i,j,1,slice)=1;
                template(i,j,2,slice)=0;
                template(i,j,3,slice)=1;
            end
            if(out(i,j)==5)
                template(i,j,1,slice)=1;
                template(i,j,2,slice)=1;
                template(i,j,3,slice)=0;
            end
            if(out(i,j)==6)
                template(i,j,1,slice)=0;
                template(i,j,2,slice)=1;
                template(i,j,3,slice)=1;
            end
            if(out(i,j)==7)
                template(i,j,1,slice)=1;
                template(i,j,2,slice)=1;
                template(i,j,3,slice)=1;
            end

        end
    end  
end


