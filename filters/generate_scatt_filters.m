function [filters,lpal] = generate_scatt_filters(options)

N  = getoptions(options,'N',256);
J{1} = getoptions(options,'J1',5);
L{1} = getoptions(options,'L1',8);
Q{1} = getoptions(options,'Q1',1);
J{2} = getoptions(options,'J2',5);
L{2} = getoptions(options,'L2',8);
Q{2} = getoptions(options,'Q2',1);
J{3} = getoptions(options,'J3',5);
L{3} = getoptions(options,'L3',8);
Q{3} = getoptions(options,'Q3',1);
spiral = getoptions(options,'spiral',0);
stretchspiral = getoptions(options,'stretchspiral',0);
isGPU= getoptions(options,'isGPU',0);
isComplex=getoptions(options,'isComplex',0);
mirror=0;
onedim = getoptions(options,'onedim',0);
spline = getoptions(options,'spline',0);
l1renorm = getoptions(options,'l1renorm',0);

if onedim
    Q{1} = L{1};
    Q{2} = L{2};
    Q{3} = L{3};
%     L{1} = 1;
%     L{2} = 1;
%     L{3} = 1;
end
gpu = getoptions(options,'gpu',0);

for m=1:size(J,2)
    
    filtopts.L=L{m};
    filtopts.Q=Q{m};
    filtopts.J=J{m}*filtopts.Q;
    if onedim
        filtopts.filter_format='fourier_multires';
        filtopts.boundary = 'per';
        filtopts.shrink = 1;
%         if m==2
%             filtopts.shrink = L{m};
%         end
    end
    if mirror
        filtopts.min_margin = 0.25 * 2^(J{m}) ;
    else
        filtopts.min_margin = 0 * 2^(J{m}) ;
    end
    if onedim
        if spline
            filnew=spline_filter_bank_1d(N, filtopts);
        else
            filnew=morlet_filter_bank_1d(N, filtopts);
        end
    else
        if spiral
       %     filnew=morlet_filter_bank_2d_spiral([N N], filtopts,isGPU,isComplex);
            filnew=morlet_filter_bank_2d_spiral(N, filtopts,isGPU,isComplex);
        else if stretchspiral
               filtopts.coef_stretching = getoptions(options,'coef_stretching',1);
               filnew=morlet_filter_bank_2d_stretchspiral(N, filtopts);
            else
                filnew=morlet_filter_bank_2d(N, filtopts);
            end 
        end
    end
    for r=1:size(filnew.psi.filter{1}.coefft,2)%J{m}
        for j=1:J{m}
            for l=1:L{m}
                %if onedim
                %    filters{m}.psi{r}{j}{l}=filnew.psi.filter{j}.coefft{r};
                %else
                    filters{m}.psi{r}{j}{l}=filnew.psi.filter{(j-1)*L{m}+l}.coefft{r};
                    if l1renorm
                        if onedim
                            rien = ifft(filters{m}.psi{r}{j}{l});
                            filters{m}.psi{r}{j}{l} = filters{m}.psi{r}{j}{l}/(sum(abs(rien(:)))+eps);
                        else
                            rien = ifft2(filters{m}.psi{r}{j}{l});
                            filters{m}.psi{r}{j}{l} = filters{m}.psi{r}{j}{l}/(sum(abs(rien(:)))+eps);
                        end
                    end
                %end
            end
        end
        %if onedim
        %filters{m}.phi{r} = filnew.phi{r};
        %else
        filters{m}.phi{r} = filnew.phi.filter.coefft{r};
      %  filters{m}.Q = 1;%Q{m};
        filters{m}.Q = Q{m};
        
        if onedim && ~spline
             filters{m}.meta.wo = max(filnew.psi.meta.center);
        end 
      
        clear tmp;
        tmp.psi = filters{m}.psi{r};
        tmp.phi = filters{m}.phi{r};
        dtmp = generate_dualfilters(tmp);
        filters{m}.dphi{r} = dtmp.dphi;
        filters{m}.dpsi{r} = dtmp.dpsi;
       
        
    end
    if m==1
        lpal = littlewood_paley(filnew);
    end
    
    clear filtopts;
    
end

if gpu
    filts0 = filters;
    %clear filters;
    for m=1:size(J,2)
        for r=1:J{m}
            for j=1:J{m}*filts0{m}.Q
                for l=1:L{m}
                    filters{m}.psi{r}{j}{l} = gpuArray(filts0{m}.psi{r}{j}{l});
                    filters{m}.dpsi{r}{j}{l} = gpuArray(filts0{m}.dpsi{r}{j}{l});
                end
            end
            filters{m}.phi{r} = gpuArray(filts0{m}.phi{r});
            filters{m}.dphi{r} = gpuArray(filts0{m}.dphi{r});
        end
    end
end
