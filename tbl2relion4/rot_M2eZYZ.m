function eu = rot_M2eZYZ(R)

eu = zeros(1,3);
tol = 5e-5;

if( R(3,3) < (1-tol))
   
    if( R(3,3) > (tol-1) )
       
        % GENERAL CASE
        eu(2) = acos ( R(3,3) )*180/pi;
        eu(1) = atan2( R(2,3), R(1,3) )*180/pi;
        eu(3) = atan2( R(3,2),-R(3,1) )*180/pi;
       
    else
       
        % r22 <= -1
        eu(2) = 180;
        eu(1) = -atan2( R(2,1), R(2,2) )*180/pi;
        eu(3) = 0;
       
    end
else
   
    % r22 <= -1
    eu(2) = 0;
    eu(1) = atan2( R(2,1), R(2,2) )*180/pi;
    eu(3) = 0;
   
end

end