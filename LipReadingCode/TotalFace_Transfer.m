clc;clear;
load('TotalFace');

    for user = 1 : 540
        for image = 1 : 20
            total_image{user,image}=rgb2gray(total_image{user,image});
        end 
    end
   
  total_image=reshape(total_image,1,10800);
  total_image=cell2mat(total_image); %123456
 total_image=reshape(total_image,96,96,20,540);
 
 for i = 1 : 6
     for j = 1 : 90
         label((i-1)*90+j) = i;
         if j>81
            set((i-1)*90+j) = 2;
         else
            set((i-1)*90+j) = 1;
         end
     end
 end

 mitdb.data=total_image;
 mitdb.label=label;
 mitdb.set=set;