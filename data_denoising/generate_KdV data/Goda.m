function v_new = Goda(v_old, dt, dx)
% run ths Goda scheme as described in
% question 3 of AMATH 741 question 3
% with input of previous time step v_old (column vector)
% and time step dt and spatial step dx
% output is column vector v_new with the values for the next
% time step
   N = size(v_old,1);
   e = ones(N, 1);
   P = dt/dx;
   Q = dt/(dx.^3);
   A = spdiags([-0.5 * Q * e, 0 * e, e, 0 * e, 0.5*Q*e], [-2, -1, 0, 1, 2], N, N);
   for i = 1:N
       if i == 1
           j = N;
       else
           j = i - 1;
       end
       if i == N
           k = 1;
       else
           k = i + 1;
       end
       %v^n_i + v^n_i+1
       beta = v_old(i) + v_old(k);
       %v^n_i + v^n_i-1
       alpha = v_old(i) + v_old(j);
       if i == 1
           A(1, N-1) = -0.5 * Q;
           A(1, N) = P*alpha + Q;
           A(1, 2) = -P*beta - Q;
       elseif i == 2
           A(2, N) = -0.5*Q;
           A(2, 1) = P*alpha + Q;
           A(2, 3) = -P * beta - Q;
       elseif i == N-1
           A(N-1, N-2) = P * alpha + Q;
           A(N-1, N) = -P*beta - Q;
           A(N-1, 1) = 0.5 * Q;
       elseif i == N
           A(N, N-1) = P*alpha + Q;
           A(N, 1) = -P*beta - Q;
           A(N, 2) = 0.5 * Q;
       else
           A(i, i - 1) = P*alpha + Q;
           A(i, i + 1) = -P * beta - Q;
       end
   end
   v_new = A\v_old;
end
