function transformation_array = imgregistration(img_stack,position)

[m,n,num_t]=size(img_stack);

r=20;

h=fspecial('gaussian',[3,3],1);


ymin=max(round(position(1,2)-r),1);
ymax=min(round(position(1,2)+r),m);
xmin=max(round(position(1,1)-r),1);
xmax=min(round(position(1,1)+r),n);

img_pre=img_stack(ymin:ymax,xmin:xmax,1); 


transformation_array=cell(num_t,2);

transformation_array{1,1}=maketform('affine',eye(3));
transformation_array{1,2}=[xmin ymin xmax ymax];

   

 for j=2:num_t
     
     %pre-alignment using cross correlation
     
     ymin=max(round(position(j,2)-r),1);
     ymax=min(round(position(j,2)+r),m);
     xmin=max(round(position(j,1)-r),1);
     xmax=min(round(position(j,1)+r),n);
     
     img=img_stack(ymin:ymax,xmin:xmax,j); 
     
     img_filtered=imfilter(img,h);
     img_pre_filtered=imfilter(img_pre,h);

   
     
     [optimizer, metric] = imregconfig('monomodal');
     optimizer.MaximumIterations=300;
     optimizer.MinimumStepLength = 5e-4;
     optimizer.MinimumStepLength = 0.02;
     tform = imregtform(double(img_pre_filtered),double(img_filtered),'affine',optimizer,metric);
     
     T1=maketform('affine',double(tform.T));
     
     [ysize,xsize]=size(img_pre_filtered);
     
     img_reg=imtransform(img_pre_filtered,T1,'XData',[1 xsize],'YData',[1 ysize]); 
     
     image_overlay(:,:,1)=imadjust(uint16(img_filtered));
     image_overlay(:,:,2)=imadjust(uint16(img_reg));
     image_overlay(:,:,3)=0;
     figure (3); imshow(image_overlay);
     
     
     %figure (4); imshowpair(img_filtered,img_reg,'Scaling','joint');
     
     transformation_array{j,1}=T1;
     transformation_array{j,2}=[xmin ymin xmax ymax];
     
     img_pre=img;
     
 end
 
 
end
 

 
 
         
          
     
     
     
     
     
     
 
 
