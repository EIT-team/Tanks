function Q = quat(vect,theta)
	Q = zeros(1,4);
	xyz_and_angle=[vect theta] ;   %displays inputed data and theta angle before deg to rad conversion
	%thetaconvert = (theta*pi/180)./norm(theta*pi/180) %convert deg to rad
	thetaconvert=theta*pi/180;
	px=vect(1,1);py=vect(1,2);pz=vect(1,3); 
	Q(1,1) = cos(thetaconvert/2); %w value
	Q(1,2)= px*sin(thetaconvert/2); % x value
	Q(1,3)= py*sin(thetaconvert/2); % y value
	Q(1,4)= pz*sin(thetaconvert/2); % z value
	
	q1w=Q(1,1);%qw value
	q1x=Q(1,2);%qx value
	q1y=Q(1,3);%qy value
	q1z=Q(1,4);%qz value
	
	%Normalize Quaternion due to floating point errors
     	qnorm=sqrt(q1w^2+q1x^2+q1y^2+q1z^2) ;
     	Q=Q./qnorm ;%normalize Q array
end