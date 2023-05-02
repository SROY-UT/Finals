function [ B  ] = bidiag_francis_step(B)
% Compute the Givens rotation matrix Q to introduce a bulge at (i+1,i)
( m, n ) = size(B);
i = m
x = B(i);
y = B(i+1);
if y == 0
    c = sign(x);
    s = 0;
else
    if abs(y) > abs(x)
        t = -x / y;
        s = 1 / sqrt(1 + t^2);
        c = s * t;
    else
        t = -y / x;
        c = 1 / sqrt(1 + t^2);
        s = c * t;
    end
end
Q = eye(length(B));
Q(i:i+1,i:i+1) = [c,s;-s,c];

% Apply the rotation to the current and next rows of B
B(i:i+1,:) = Q' * B(i:i+1,:);

% Compute the Givens rotation matrix P to chase out the bulge at (i+2,i)
x = B(i+1);
y = B(i+2);
if y == 0
    c = sign(x);
    s = 0;
else
    if abs(y) > abs(x)
        t = -x / y;
        s = 1 / sqrt(1 + t^2);
        c = s * t;
    else
        t = -y / x;
        c = 1 / sqrt(1 + t^2);
        s = c * t;
    end
end
P = eye(length(B));
P(i+1:i+2,i+1:i+2) = [c,s;-s,c];

% Apply the rotation to the next and next-next columns of B
B(:,i+1:i+2) = B(:,i+1:i+2) * P;

% Zero out any small off-diagonal entries of B
tol = eps(class(B)) * max(size(B));
B(abs(B) < tol) = 0;

end
