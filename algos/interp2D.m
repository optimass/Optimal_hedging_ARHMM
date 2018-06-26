%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a wrapper for the function interp2 bc there is no extrapolate
% feature.
%
%  input:
%      [X,Y] = meshgrid(Xgrid,Ygrid);
%       Z : function to interpolate
%      X0 : where to interpolate
%      Y0 : idem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = interp2D(X, Y, Z, X0, Y0)


out = interp2(X, Y, Z, X0, Y0,'linear');

if isnan(out)
    
    Xgrid = X(1,:)';
    Ygrid = Y(:,1);
    
    X_too_small = X0 < min(min(X));
    X_too_big   = X0 > max(max(X));
    
    Y_too_small = Y0 < min(min(Y));
    Y_too_big   = Y0 > max(max(Y));
    
    
    if X_too_small && Y_too_small
        out = Z(1,1);
    elseif X_too_small && Y_too_big
        out = Z(end,1);
    elseif X_too_big && Y_too_small
        out = Z(1,end);
    elseif X_too_big && Y_too_big
        out = Z(end,end);
    elseif X_too_small
        out = interp1(Ygrid,Z(:,1),Y0,'linear','extrap');
    elseif X_too_big
        out = interp1(Ygrid,Z(:,end),Y0,'linear','extrap');
    elseif Y_too_small
        out = interp1(Xgrid,Z(1,:),X0,'linear','extrap');
    elseif Y_too_big
        out = interp1(Xgrid,Z(end,:),X0,'linear','extrap');
    end
    
end



% [X,Y] = ndgrid(1:10,1:10);

 
    