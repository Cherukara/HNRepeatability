function x = ak_ConjGrad(A,b,tol,disp_norm)
%DESCRIPTION: x = ak_ConjGrad(A,b,tol)
%             Solves Ax = b for x using Conjugate Gradient Least Squares
%
%INPUTS:
%   A(function) - System 'matrix'
%   b(double vector) - RHS
%   tol(double) - Tolerance level for norm(residuals)
%   disp_norm(string) - Display norm: 'Yes' or 'No'
%
%OUTPUTS:
%   x(double vector) - Least Squares Solution
%
%DEPENDENCIES:
%   None
%
%AUTHOR:
%   Anita Karsa, University College London, 2016

r = b;
x = zeros(size(r));
p = r;

% res = 2;
% 
% K = 1;

while and(sqrt((r'*r)/(b'*b))>tol,0==0)%sqrt((r'*r)/(b'*b))<res) %sqrt((r'*r)/(b'*b)) %norm(r)
        
    res = sqrt((r'*r)/(b'*b));
    
    % Calculate Ap
    Ap = A(p);
    
    % Calculate alpha
    denom = p'*Ap;
    alpha = p'*r/denom;
    
    % Calculate x
    x = x + alpha*p;
    
    % Calculate r
    r = r - alpha*Ap;
    
    % Calculate Ar
    Ar = A(r);
    
    %Calculate beta
    beta = p'*Ar/denom;
    
    %Calculate p
    p = r - beta*p;
    
    % Stages
%     X = reshape(x,[164 205 205])/0.4309;
%     imshow(rot90(squeeze(X(:,96,:))),[-0.2 0.2]);
%     pause
% %    x_out(:,:,K) = rot90(squeeze(X(:,:,71)));
%     Params.Resolution = [1 1 1];
%     Freq = ak_FrequencyMap(X,Params);
%     B = reshape(b,[240 240 144]);
%     subplot(1,2,2); imshow(rot90(squeeze(Freq(120,:,:)-B(120,:,:))),[-0.1 0.1]);
%    res = reshape(r,[240 240 144]);
%    figure(1); imshow(rot90(squeeze(res(120,:,:))),[-0.01 0.01]);
%    f_out(:,:,K) = rot90(squeeze(Freq(:,:,71)));
%    K = K + 1;
    
    if strcmp(disp_norm,'Yes')
        disp(['"Residuals" = ',num2str(norm(sqrt((r'*r)/(b'*b))))]);
    end
end

end