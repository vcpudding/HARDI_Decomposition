
function [weights fiber_directions final_func]=deconvolve_two_fibers(S, kernel, g2, tri_areas, GD, initialization, options)

% initialization

order=options.order; % order=12 works very well
delta=options.delta; % estimated from the data
lambda=options.lambda; % initial damping parameter value
tol=options.tol; % set to 1e-6
r=options.rank; % set to 2
t=options.modsel; % threshold for model selection, to avoid model selection set very low
max_iter=options.maxiter; % maximal number of iterations, default is 200
v_coef=initialization;

% if options.init 
% % Initialize with the signal minimum
%     [SS idx] = sort(S);
%     v_coef(:,1)=GD(idx(1),:); % initialize with the minimum of S
%     v_coef(:,2)=GD(idx(2),:);
% else
%     % random initialization
%     v_coef=randn(3,r);
%     %v_coef(:,1)=v_coef(:,1)/norm(v_coef(:,1));
%     %v_coef(:,2)=v_coef(:,2)/norm(v_coef(:,2));
% end

%n=length(S);

%lambda1=lambda;
%lambda2=lambda;

%v_coef=Stabilized_Gram_Schmidt(v_coef); % Orthonormalize initilization for better results

%v_coef_update=v_coef; %zeros(3,r);

dt=1; %1/S0;
iter=1;
c=10; % update factor for the LM damping parameter
%dFunc=1;
%dNorm=1; 

%func2=1e4;

lambda1=lambda;
lambda2=lambda;

[dF M] = compute_gradient_ALS_fast(S, kernel, v_coef, g2, tri_areas, GD, order, delta, 1); % area, cen

func=0.5*sum(dF.^2);

%while dNorm>tol

%v_coef_old=v_coef;

for iter=1:max_iter

v_coef_new(:,2)=v_coef(:,2);
    
v_coef_new(:,1)=v_coef(:,1)-dt*(M'*M+lambda1*diag(diag(M'*M)))\(M'*dF);
[new_dF new_M] = compute_gradient_ALS_fast(S, kernel, v_coef_new, g2, tri_areas, GD, order, delta, 2); % area, cen
new_func=0.5*sum(new_dF.^2);
if new_func<func
        v_coef(:,1)=v_coef_new(:,1);
        func=new_func;
        lambda1=lambda1/c;
        M=new_M;
        dF=new_dF;
        %final_lambda1(iter)=lambda(i);
else
      
      while new_func>func
            lambda1=lambda1*c;
            v_coef_new(:,1)=v_coef(:,1)-dt*(M'*M+lambda1*diag(diag(M'*M)))\(M'*dF);
            [new_dF new_M] = compute_gradient_ALS_fast(S, kernel, v_coef_new, g2, tri_areas, GD, order, delta, 2); % area, cen
            new_func=0.5*sum(new_dF.^2);
      end
      v_coef(:,1)=v_coef_new(:,1);
      M=new_M;
      dF=new_dF;
      func=new_func;
end

v_coef_new(:,1)=v_coef(:,1);

v_coef_new(:,2)=v_coef(:,2)-dt*(M'*M+lambda2*diag(diag(M'*M)))\(M'*dF);
[new_dF new_M] = compute_gradient_ALS_fast(S, kernel, v_coef_new, g2, tri_areas, GD, order, delta, 1); % area, cen
new_func=0.5*sum(new_dF.^2);
    if new_func<func
        v_coef(:,2)=v_coef_new(:,2);
        func=new_func;
        lambda2=lambda2/c;
        M=new_M;
        dF=new_dF;
        %final_lambda2(iter)=lambda(i);
    else
        
        while new_func>func
            lambda2=lambda2*c;
            v_coef_new(:,2)=v_coef(:,2)-dt*(M'*M+lambda2*diag(diag(M'*M)))\(M'*dF);
            [new_dF new_M] = compute_gradient_ALS_fast(S, kernel, v_coef_new, g2, tri_areas, GD, order, delta, 1); % area, cen
            new_func=0.5*sum(new_dF.^2);
        end
        v_coef(:,2)=v_coef_new(:,2);
        M=new_M;
        dF=new_dF;
        func=new_func;
    end


final_func(iter)=func;

%BIC11(iter)=n*log(norm(S'-kernel*(((v_coef(1,1)*g2(:,1)+v_coef(2,1)*g2(:,2)+v_coef(3,1)*g2(:,3)).^order).*tri_areas')).^2/n)+3*log(n);
%BIC12(iter)=n*log(norm(S'-kernel*(((v_coef(1,2)*g2(:,1)+v_coef(2,2)*g2(:,2)+v_coef(3,2)*g2(:,3)).^order).*tri_areas')).^2/n)+3*log(n);
%BIC2(iter)=n*log(norm(S'-kernel*(((v_coef(1,1)*g2(:,1)+v_coef(2,1)*g2(:,2)+v_coef(3,1)*g2(:,3)).^order).*tri_areas')-...
%    kernel*(((v_coef(1,2)*g2(:,1)+v_coef(2,2)*g2(:,2)+v_coef(3,2)*g2(:,3)).^order).*tri_areas')).^2/n)+6*log(n);

%BIC11(iter)=1-norm(S'-kernel*(((v_coef(1,1)*g2(:,1)+v_coef(2,1)*g2(:,2)+v_coef(3,1)*g2(:,3)).^order).*tri_areas'))^2/norm(S)^2;
%BIC12(iter)=1-norm(S'-kernel*(((v_coef(1,2)*g2(:,1)+v_coef(2,2)*g2(:,2)+v_coef(3,2)*g2(:,3)).^order).*tri_areas'))^2/norm(S)^2;
% % BIC2(iter)=norm(S'-kernel*(((v_coef(1,1)*g2(:,1)+v_coef(2,1)*g2(:,2)+v_coef(3,1)*g2(:,3)).^order).*tri_areas')-...
% %         kernel*(((v_coef(1,2)*g2(:,1)+v_coef(2,2)*g2(:,2)+v_coef(3,2)*g2(:,3)).^order).*tri_areas'))^2;
% 
% 
% if iter>6 & BIC11(iter)<t
%     weights=[0 1];
%     fiber_directions(:,1)=zeros(3,1);
%     fiber_directions(:,2)=v_coef(:,2)/norm(v_coef(:,2));
%     return
% elseif iter>6 & BIC12(iter)<t
%     weights=[1 0];
%     fiber_directions(:,1)=v_coef(:,1)/norm(v_coef(:,1));
%     fiber_directions(:,2)=zeros(3,1);
%     return
% end
    
    
if iter>1 & abs(final_func(iter)-final_func(iter-1))/final_func(iter)<tol
    break
end


%fnorm1(iter)=norm(v_coef(:,1));
%fnorm2(iter)=norm(v_coef(:,2));


end
 
if iter==max_iter
    display('slow convergence');
end
 
figure(1),plot(final_func)
%figure(2),plot(BIC11,'r')
%hold on,plot(BIC12,'g')
%hold on,plot(BIC2,'b')
%hold off

v1=v_coef(:,1);
v2=v_coef(:,2);

fnorm1=norm(v1);
fnorm2=norm(v2);

fiber_directions(:,1) = v1/fnorm1;
fiber_directions(:,2) = v2/fnorm2;

weights(1) = fnorm1^order/(fnorm1^order+fnorm2^order);
weights(2) = fnorm2^order/(fnorm1^order+fnorm2^order);




