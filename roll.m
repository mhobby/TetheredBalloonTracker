function [x,y,z,H]=roll(in, B, L)
  % vector; [in] are provide in body reference frame
  % x,y,z are produced in inertial/earth reference frame
  % B is the angle to move from body to earth frame
  
  %check to ensure that in is a 1x3 or vector
  M=size(in);
  if(M(1)~=3)
    printf('Error > Input vector dimensions are not nx3!\n');
  else
    out=zeros(size(in));
    for i=1:L
      H = [1, 0, 0;
           0, cos(B(i)), sin(B(i));
           0, -sin(B(i)), cos(B(i));];
      out(:,i) = H*in(:,i);
    end
    x=out(1, :);
    y=out(2, :);
    z=out(3, :);
  end
end