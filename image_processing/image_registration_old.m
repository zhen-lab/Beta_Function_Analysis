function transformation_array = image_registration(img_stack,template,corner)

[m,n,num_t]=size(img_stack);

h=fspecial('gaussian',[3,3],1);

img_pre=template;

[ysize,xsize]=size(template);


ymin=corner(1);
ymax=corner(2);
xmin=corner(3);
xmax=corner(4);

transformation_array=cell(num_t,2);

transformation_array{1,1}=maketform('affine',eye(3));
transformation_array{1,2}=[xmin ymin xmax ymax];

   

 for j=2:num_t
     
     %pre-alignment using cross correlation
     
     img=img_stack(ymin:ymax,xmin:xmax,j); 
     
     img_filtered=imfilter(img,h);
     %thres=mean(img_filtered(:))+std(double(img_filtered(:)));
     %img_filtered_bw=(img_filtered>thres);
     %img_filtered_bw=bwareaopen(img_filtered_bw,4);
     
     img_pre_filtered=imfilter(img_pre,h);
     %thres=mean(img_pre_filtered(:))+std(double(img_pre_filtered(:)));
     %img_pre_filtered_bw=(img_pre_filtered>thres);
     %img_pre_filtered_bw=bwareaopen(img_pre_filtered_bw,4);
     
     %c=xcorr2(double(img_filtered_bw),double(img_pre_filtered_bw));
     %c=xcorr2(img_filtered,img_pre_filtered);
     %[~,imax]=max(abs(c(:)));

     %[ypeak,xpeak]=ind2sub(size(c),imax(1));
     

     %ymin=max(ypeak-ysize+ymin,1);
     %xmin=max(xpeak-xsize+xmin,1);
     
     
     %ymax=min(ymin+ysize-1,m);
     %xmax=min(xmin+xsize-1,n);
     
     %fine alignement that allows rotation of the images
     
     %img=img_stack(ymin:ymax,xmin:xmax,j);
     %img_filtered=imfilter(img,h);
     %thres=mean(img_filtered(:))+std(double(img_filtered(:)));
     %img_filtered_bw=(img_filtered>thres);
     %img_filtered_bw=bwareaopen(img_filtered_bw,4);
     
     
     [optimizer, metric] = imregconfig('monomodal');
     optimizer.MaximumIterations=400;
     tform = imregtform(double(img_pre_filtered),double(img_filtered),'rigid',optimizer,metric);
     
     T1=maketform('affine',double(tform.T));
     
     img_reg=imtransform(img_pre_filtered,T1,'XData',[1 xsize],'YData',[1 ysize]); 
     
     %image_overlay(:,:,1)=imadjust(uint16(img_filtered));
     %image_overlay(:,:,2)=imadjust(uint16(img_reg));
     %image_overlay(:,:,3)=0;
     %figure (3); imshow(image_overlay);
     
     
     figure (3); cla; imshowpair(img_filtered,img_reg,'Scaling','joint');
     
     transformation_array{j,1}=T1;
     transformation_array{j,2}=[xmin ymin xmax ymax];
     
     img_pre=img;
     
 end
 
 
end
 

 
 
         
          
     
     
     
     
     
     
 
 