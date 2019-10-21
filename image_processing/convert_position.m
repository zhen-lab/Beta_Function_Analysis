%This function uses the timing information to calculate the real position
%of the object on the stage.



function position_on_plate = convert_position( position_data,t_nd2,t_stage,metadata_stage,istart,iend)
scale_x=0.1; %one motor step in x is 0.1 micron   
scale_y=0.1; %one motor step in y is 0.1 micron
pixel_scale=3.33; % one pixel is 3.33 micron at 4x

 
position_on_plate=zeros(length(position_data),2);

N=iend-istart+1;


for i=1:N
    
    
    
   
    j=i+istart-1;
    
    [~,idx]=min(abs(t_nd2(j)-t_stage));
    
    if i==1

        shift_x=metadata_stage(idx,2);
        shift_y=metadata_stage(idx,3);
       
    end
    

    
    stage_x=(metadata_stage(idx,2)-shift_x)*scale_x;
    stage_y=(metadata_stage(idx,3)-shift_y)*scale_y;
    
    position_on_plate(i,1)=-position_data(i,1)*pixel_scale-stage_x; %positive when worm moves right
    position_on_plate(i,2)=position_data(i,2)*pixel_scale+stage_y; %positive when worm moves up
    

end

end

