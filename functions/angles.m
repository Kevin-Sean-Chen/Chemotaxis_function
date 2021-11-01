function aa = angles(v,u)%u=target
    
%     CosTheta = dot(u,v)/sqrt(sum(u.^2)*sum(v.^2));%(norm(u)*norm(v));
%     ThetaInDegrees = acosd(CosTheta);
%     aa = ThetaInDegrees;
    
    %%%with angle sign
    v_3d = [v, 0];
    u_3d = [u, 0];
    c = cross(v_3d, u_3d);
    
    % calculate degrees
    if c(3) < 0
        aa = -atan2(norm(c),dot(v,u))*(180/pi);
    else
        aa = atan2(norm(c),dot(v,u))*(180/pi);
    end
end