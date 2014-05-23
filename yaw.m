function [x,y,z,H]=yaw(in, g, L)
  %check to ensure that in is a 1x3 or vector
  % vector; [in] are provide in earth reference frame
  % x,y,z are produced in body reference frame
  % g is the angle to move from earth to body frame
  
  M=size(in);
  if(M(1)~=3)
     error('Error > Input vector dimensions are not nx3!');
  else
    out=zeros(size(in));
    for i=1:L
      H = [cos(g(i)), sin(g(i)), 0;
           -sin(g(i)), cos(g(i)), 0
           0, 0, 1];
      out(:,i) = H*in(:,i);
    end
    x=out(1, :);
    y=out(2, :);
    z=out(3, :);
  end
end