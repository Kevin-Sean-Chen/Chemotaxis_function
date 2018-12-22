function aa = angles(v,u)%u=target
    
    CosTheta = dot(u,v)/sqrt(sum(u.^2)*sum(v.^2));%(norm(u)*norm(v));
    ThetaInDegrees = acosd(CosTheta);
    aa = ThetaInDegrees;

end