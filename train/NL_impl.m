function [output,anomaly]=NL_impl(ref,src,search_windiow,patch_size,eps)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  ref : reference image (image without anomaly)
 %  src : source image (image to check)
 %  input: image to be filtered
 %  search_windiow: radio of search window
 %  patch_size: radio of similarity window
 %  eps: degree of filtering
 %
 %  Author: Jose Vicente Manjon Herrera & Antoni Buades
 %  Date: 09-03-2006
 %
 %  Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
 %  "A non-local algorithm for image denoising"
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Size of the image
 [m , n]=size(ref);
 
 
 
 % Memory for the output
 output=zeros(m,n);
 anomaly = double(zeros(m,n));

 % Replicate the boundaries
 ref_pad = padarray(ref,[patch_size patch_size],'symmetric');
 src_pad = padarray(src,[patch_size patch_size],'symmetric');
 
 % Used kernel
 kernel = make_kernel(patch_size);
 kernel = kernel / sum(sum(kernel));
 
 eps=eps*eps;
 
 %For every pixel of the src image
 for i=1:m
 for j=1:n
                 
         i1 = i+ patch_size;
         j1 = j+ patch_size;
         
         %Obtain the patch around the pixel
         W1= src_pad(i1-patch_size:i1+patch_size , j1-patch_size:j1+patch_size);
         
         wmax=0; 
         average_pixel=0;
         sweight=0;
         
         %Find the bound of the search radius around the pixel
         rmin = max(i1-search_windiow,patch_size+1);
         rmax = min(i1+search_windiow,m+patch_size);
         smin = max(j1-search_windiow,patch_size+1);
         smax = min(j1+search_windiow,n+patch_size);
         
         %For every pixel in the search radius
         for r=rmin:1:rmax
         for s=smin:1:smax
                 
                %avoid overweight of pixel (i,j) with itself
                if(r==i1 && s==j1) 
                    continue; 
                end;
                
                %patch around the pixel in the serach radius
                W2= ref_pad(r-patch_size:r+patch_size , s-patch_size:s+patch_size);                
                
                %compute distance between the patches
                d = sum(sum(kernel.*(W1-W2).*(W1-W2)));
                 
                %compute weight
                w=exp(-d/eps);                 
                                 
                if w>wmax                
                    wmax=w;
                end
                
                %sum all the weight
                sweight = sweight + w;
                %compute the value of the output without normalization
                average_pixel = average_pixel + w*ref_pad(r,s); 
         end 
         end
             
	    %wmax is a sort of estimation of the weight of a pixel with itself
        %pixel estimation
        average_pixel = average_pixel + wmax*ref_pad(i1,j1);
        %patch estimation

        sweight = sweight + wmax;

        %normalize the output value
        output(i,j) = average_pixel / sweight;           
 end
 end
 

output_pad = padarray(output,[patch_size patch_size],'symmetric');
for x1=1:m
    for x2=1:n
         x1_ = x1+ patch_size;
         x2_ = x2+ patch_size;
         P_out= output_pad(x1_-patch_size:x1_+patch_size , x2_-patch_size:x2_+patch_size);
         P = src_pad(x1_-patch_size:x1_+patch_size , x2_-patch_size:x2_+patch_size);
         anomaly(x1,x2) = norm(P_out - P);
    end
end
 
function [kernel] = make_kernel(f)              
 
kernel=zeros(2*f+1,2*f+1);   
for d=1:f    
  value= 1 / (2*d+1)^2 ;    
  for i=-d:d
  for j=-d:d
    kernel(f+1-i,f+1-j)= kernel(f+1-i,f+1-j) + value ;
  end
  end
end
kernel = kernel ./ f;

        