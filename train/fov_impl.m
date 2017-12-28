function [y_hat,patch_hat]=fov_impl(ref,src,search_radius,U_radius,sigma,disableFov)
%
%  Foveated Nonlocal Means denoising algorithm (ver 1.12, Feb. 16, 2012)
%  Author:  Alessandro Foi, Tampere University of Technology, Finland
% -------------------------------------------------------------------------
%  SYNTAX
%
%     [y_hat] = FovNLM(z,search_radius,U_radius,sigma,disableFov)
%
%
%  z :              input noisy image to be filtered
%  search_radius :  radius of search window     (e.g., 5 <= search_radius <= 13 )
%  U_radius :       radius of patch             (e.g., 1 <= U_radius <= 10)
%  sigma :          standard deviation of the additive white Gaussian noise
%  disableFov :     disables foveation when set to 1  (optional)
%
%  y_hat :          output denoised image
%
%
%  This software implements the Foveated Nonlocal Means algorithm
%  proposed in the paper A. Foi and G. Boracchi, "Foveated self-similarity
%  in nonlocal image filtering", in Proc. IS&T/SPIE EI2012 - Human Vision
%  and Electronic Imaging XVII, 8291-32, Burlingame (CA), USA, Jan. 2012.
%
%  By setting disableFov=1, this software can replicate the results of the
%  NLmeansfilter.m  code by J.V. Manjon Herrera & A. Buades (31.07.2008),
%  which implements the NL-means filter proposed by A. Buades, B. Coll and
%  J.M. Morel in the paper "A non-local algorithm for image denoising".
%  http://www.mathworks.com/matlabcentral/fileexchange/13176-non-local-means-filter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% windowing k_kernel  
% (withing the Foveated NL-means this kernel is not used for windowing,
% but to determine the standard-deviation of the blurring kernels used by the foveation operator)
[u2,u1]=meshgrid(-U_radius:U_radius,-U_radius:U_radius);
linfdistance=max(abs(u1),abs(u2));
K_tilde=cumsum((2*(U_radius:-1:0)+1).^-2);
K_tilde=K_tilde([U_radius,U_radius:-1:1]);
k_kernel=interp1(0:U_radius,K_tilde,linfdistance,'linear','extrap');

k_kernel=k_kernel/sum(k_kernel(:));


%% construction of the Gaussian blurring kernels used by the foveation operator

% the set of values assumed by the windowing kernel k_kernel, sorted from smallest to largest
K=unique(k_kernel);  
% sorts from largest to smallest (the sorting is not essential, but placed to maintain consistent 
% indexing with the examples in the paper)
K=K(end:-1:1);    

%   1-exp(-2*pi) = 0.9981  % PSF mass over disk of unit radius
p=1-exp(-2*pi);
% standard-deviations of the Gaussian blurring kernels 
%(here max(K) corresponds to \kappa_0 in the paper)
varsigma=sqrt(1/(-2*log(1-p))*max(K)./K);  
% support widths of the gaussian blurring kernels
v_kernel_widths=2*ceil(3*varsigma)+1;   

v_kernels=cell(size(K));
jj=0;
for kappa=K',
    jj=jj+1;
    % generate k_kernel as scaled bivariate Gaussian p.d.U_radius.
    v_kernels{jj}=(sqrt(pi*kappa)*2*varsigma(jj))*...
        fspecial('gaussian',v_kernel_widths(jj),varsigma(jj));   
end


%% option: disable foveation and implement windowing by setting all blurring kernels as suitably scaled Dirac impulses. In this way the rest of the code will operate equivalent to the standard NL-means (with windowed patch distance).

if exist('disableFov','var')
    if disableFov==1
        disp('  * * *  foveation is disabled   * * *  ');
        v_kernels=cell(size(K));
        jj=0;
        for kappa=K',
            jj=jj+1;
            v_kernels{jj}=sqrt(kappa);
            v_kernel_widths(jj)=1;
        end
    end
end


%% foveation operator

refBlurred = cell(size(v_kernels));
srcBlurred = cell(size(v_kernels));
for jj=1:numel(v_kernels)
    if 0
        zBlurred{jj}=conv2(...
            padarray(...
            z,[U_radius U_radius]+(v_kernel_widths(jj)-1)/2,'symmetric'),v_kernels{jj},'valid');
    else   % realize 2D convolution as cascaded 1D convolutions in separable manner
        % 1D section (2D k_kernel is separable)
        v_kernel1D=v_kernels{jj}(:,(size(v_kernels{jj},2)+1)/2);  
        v_kernel1D=v_kernel1D/sqrt(max(v_kernel1D(:)));
        refBlurred{jj}=conv2(...
            padarray(...
            ref,[U_radius U_radius]+(v_kernel_widths(jj)-1)/2,'symmetric'),v_kernel1D,'valid');
        refBlurred{jj}=conv2(refBlurred{jj},v_kernel1D','valid');
        
        srcBlurred{jj}=conv2(...
            padarray(...
            src,[U_radius U_radius]+(v_kernel_widths(jj)-1)/2,'symmetric'),v_kernel1D,'valid');
        srcBlurred{jj}=conv2(srcBlurred{jj},v_kernel1D','valid');
    end
end

[size_z_1, size_z_2]=size(ref);
numel_k_kernel=numel(k_kernel);

% big array with all foveated patches;
BIGfov_ref=zeros(numel_k_kernel,size_z_1,size_z_2); 
BIGfov_src = zeros(numel_k_kernel,size_z_1,size_z_2); 
for jj=1:numel(v_kernels)
    W=k_kernel==K(jj);
    for x1=1:size_z_1
        for x2=1:size_z_2
            W_ref= refBlurred{jj}(x1:x1+2*U_radius,x2:x2+2*U_radius);
            %  BIGfov(:,x1,x2) contains the foveated patch centered at 
            % (x1,x2), reshaped to a column vector.
            BIGfov_ref(W,x1,x2)=W_ref(W);   
            
            W_src= srcBlurred{jj}(x1:x1+2*U_radius,x2:x2+2*U_radius);
            %  BIGfov(:,x1,x2) contains the foveated patch centered at 
            % (x1,x2), reshaped to a column vector.
            BIGfov_src(W,x1,x2)=W_src(W); 
        end
    end
end


%% foveated NL-means

y_hat = zeros(size_z_1,size_z_2);
sigma2 = sigma^2;
for x1=1:size_z_1
    for x2=1:size_z_2
        % set boundaries of search window
        search1min=max(x1-search_radius,1);  search1max=min(x1+search_radius,size_z_1);  % vertical (left; right)
        search2min=max(x2-search_radius,1);  search2max=min(x2+search_radius,size_z_2);  % horizontal (top; bottom)
        
        % compute foveated differences and foveated distance for all foveated patches within the search window
        dfov=0;
        for u=1:numel_k_kernel
            %compute foveated distance        
            dfov = dfov+(BIGfov_ref(u,search1min:search1max,search2min:search2max)-BIGfov_src(u,x1,x2)).^2; 
        end
        
        % weights
        % permutation is used to 'squeeze' away the first singleton dimension
        w=exp(-permute(dfov,[2 3 1])/sigma2);  
        w(x1-search1min+1,x2-search2min+1)=0; 
        wmax=max(w(:)); 
        % weight for z(x1,x2) is set equal to the maximum weight 
        %(unless all weights are zero, in which case it is set to 1)
        w(x1-search1min+1,x2-search2min+1)=wmax+~wmax; 
        w=w/sum(w(:));
        
        % weighted averaging
        y_hat(x1,x2)=sum(sum(w.*ref(search1min:search1max,search2min:search2max)));
        

    end
end

patch_hat = zeros(size_z_1,size_z_2);
y_hat_pad = padarray(y_hat,[U_radius U_radius],'symmetric');
src_pad= padarray(src,[U_radius U_radius],'symmetric');
for x1=1:size_z_1
    for x2=1:size_z_2
         x1_ = x1+ U_radius;
         x2_ = x2+ U_radius;
         P_hat= y_hat_pad(x1_-U_radius:x1_+U_radius , x2_-U_radius:x2_+U_radius);
         P = src_pad(x1_-U_radius:x1_+U_radius , x2_-U_radius:x2_+U_radius);
         patch_hat(x1,x2) = norm(P_hat - P);
    end
end