function [x,y,z, H]=pitch(in, a, L)
% vector; [in] are provide in earth reference frame
  % x,y,z are produced in body frame
  % a is the angle to move from inertial to body frame
  
%check to ensure that in is a 1x3 or vector
  M=size(in);
  if(M(1)~=3)
    printf('Error > Input vector dimensions are not nx3!\n');
  else
    out=zeros(size(in));
    for(i=1:L)
      H = [cos(a(i)), 0, -sin(a(i));
           0, 1, 0;
           sin(a(i)), 0, cos(a(i))];
      out(:,i) = H*in(:,i);
    end
    x=out(1, :);
    y=out(2, :);
    z=out(3, :);
  end
end